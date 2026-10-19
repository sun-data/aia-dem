"""
Build the rolling 48-hour AIA DEM site that solar-movie-viewer plays.

Each run fills in the 3-minute slots of the last 48 hours that the site does
not have yet, newest first, until a time budget runs out. For each slot it
downloads the six coronal channels of the near-real-time AIA synoptic images
from JSOC, inverts them for a differential emission measure (DEM) with
:func:`utu.dem.plowman`, colors the DEM with :mod:`colorsynth`, and saves the
frame as a JPEG at each size. Frames older than the window are deleted, and
the manifest the viewer reads is rewritten from the frames on hand.

The synoptic images are 1024 x 1024 at 2.4 arcsec per pixel, level 1.5, with
the Sun centered. Each pixel is the mean, in DN, of a 4 x 4 block of the full
image, so the rate is the value over the exposure time, and the noise is that
of a mean of 16 pixels.

Usage: python build.py SITE [--hours HOURS] [--budget MINUTES] [--end TIME]
"""

import argparse
import concurrent.futures
import datetime
import io
import json
import pathlib
import time
import typing
import urllib.error
import urllib.request

import astropy.units as u
import colorsynth
import matplotlib
import named_arrays as na
import numpy as np
import sdo
import utu
from astropy.io import fits
from PIL import Image

matplotlib.use("agg")
import matplotlib.pyplot as plt  # noqa: E402

base_url = "https://jsoc1.stanford.edu/data/aia/synoptic/nrt"
"""Where JSOC keeps the near-real-time synoptic images."""

wavelengths = (94, 131, 171, 193, 211, 335)
"""The channels inverted, in angstroms."""

cadence = datetime.timedelta(minutes=3)
"""The spacing of the slots, which is that of the JSOC images."""

grid = slice(20, 71)
"""Log T 5.0 to 7.5 on the 0.05 grid of the response, which starts at 4.0."""

shown = slice(12, 45)
"""Log T 5.6 to 7.2 within :data:`grid`, the range the channels constrain."""

n_sample = 16
"""The number of full-resolution pixels each synoptic pixel is the mean of."""

saturation = 15000
"""A block mean this high, in DN, is at or near the 16383 of the camera."""

systematic_94 = 1
"""
An uncertainty in DN/s added in quadrature to that of 94. The response of 94
to plasma at 1-2 MK is too small, so the fit explains the faint 94 of coronal
holes with plasma near 6 MK and paints them red. This lets the fit miss 94 by
about 1 DN/s there, which leaves the quiet Sun as it was and adds about a
third to the 2.5-5.6 MK emission of active regions.
"""

vmax = 195
"""
The DEM, in 10^26 cm^-5 per unit log T, shown at full brightness. It is fixed
for every frame and every temperature, so that the colors and brightness of
any two pixels compare. It is the median, over a dozen frames from
2026-10-15 to 2026-10-17, of the 99.9th percentile of the DEM over every
pixel and temperature of the frame.
"""

stretch_center = 6.3
"""The log T at which the colors change fastest, near the peak of most coronal DEMs."""

stretch_width = 0.15
"""
The width, in log T, of the stretch of the colors about :data:`stretch_center`.
The colors are spread evenly over arctan((log T - 6.3) / 0.15) rather than over
log T, which gives most of the spectrum to 1.4-3 MK, where most of the corona
is, and squeezes cooler plasma into violet and hotter into red. Spread evenly
over log T, 1.4-2.5 MK was all green.
"""

version = 2
"""
The version of the frames, to raise whenever a change alters how they look.
A run that finds frames of another version makes all of them again, so that
the window never mixes two. Version 1 spread the colors evenly over log T.
"""

plate_scale = 2.4
"""The arcsec per pixel of the synoptic images."""

size = 1024
"""The number of pixels across a synoptic image."""

levels = (1024, 512)
"""The sizes each frame is saved at, in pixels across."""

disk = 1.2 * 976
"""
The radius, in arcsec, within which a pixel the inversion failed on is
painted white. The solar radius changes by under 2% over the year, which
the margin covers.
"""


def _slots(end: datetime.datetime, hours: float) -> list[datetime.datetime]:
    """The slots of the window ending at the last slot at or before `end`, newest first."""
    step = cadence // datetime.timedelta(minutes=1)
    last = end.replace(minute=end.minute - end.minute % step, second=0, microsecond=0)
    return [last - i * cadence for i in range(round(datetime.timedelta(hours=hours) / cadence))]


def _url(slot: datetime.datetime, wavelength: int) -> str:
    """Where JSOC keeps the image of one channel at one slot."""
    return f"{base_url}/{slot:%Y/%m/%d}/H{slot:%H}00/AIA{slot:%Y%m%d_%H%M%S}_{wavelength:04d}.fits"


def _fetch(url: str) -> bytes | None:
    """Download a file, or return :obj:`None` if JSOC does not have it (yet)."""
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=60) as response:
                return response.read()
        except urllib.error.HTTPError as error:
            if error.code == 404:
                return None
            if attempt == 3:
                raise
        except (urllib.error.URLError, TimeoutError, ConnectionError):
            if attempt == 3:
                raise
        time.sleep(2**attempt)
    raise AssertionError("unreachable")


def _download(slot: datetime.datetime) -> list[bytes] | None:
    """Every channel of a slot, or :obj:`None` if any is missing."""
    with concurrent.futures.ThreadPoolExecutor(len(wavelengths)) as pool:
        files = list(pool.map(_fetch, [_url(slot, w) for w in wavelengths]))
    present = [f for f in files if f is not None]
    if len(present) < len(files):
        return None
    return present


def _read(files: list[bytes]) -> tuple[np.ndarray, np.ndarray]:
    """The counts and exposure times of the channels of a slot."""
    counts, exptime = [], []
    for data in files:
        with fits.open(io.BytesIO(data)) as hdul:
            header = hdul[1].header
            # the viewer places the frames by these, so check that they hold
            centered = header["CRPIX1"] == header["CRPIX2"] == (size + 1) / 2
            centered &= header["CRVAL1"] == header["CRVAL2"] == 0
            if not centered or abs(header["CDELT1"] - plate_scale) > 1e-4 or header["CROTA2"] != 0:
                raise ValueError(f"{header['T_OBS']} {header['WAVELNTH']} is not a centered synoptic image")
            counts.append(hdul[1].data.astype(float))
            exptime.append(header["EXPTIME"])
    return np.array(counts), np.array(exptime)


def _response(when: datetime.datetime) -> na.FunctionArray:
    """
    The SolarSoft V10 response, corrected for degradation and normalized to
    EVE at `when`, with the empirical correction of 94.
    """
    channels = na.ScalarArray(np.array(wavelengths) * u.AA, axes="wavelength")
    response = sdo.aia.temperature_response(
        wavelength=channels,
        time=when.strftime("%Y-%m-%dT%H:%M"),
        eve=True,
        chiantifix=True,
    )
    return typing.cast(na.FunctionArray, response[dict(temperature=grid)])


def _invert(counts: np.ndarray, exptime: np.ndarray, response: na.FunctionArray) -> tuple[np.ndarray, np.ndarray]:
    """
    The DEM of every pixel in 10^26 cm^-5 per unit log T, with axes
    ``(y, x, temperature)`` over :data:`shown`, and its reduced chi^2.
    """
    axes = ("wavelength", "y", "x")
    channels = na.ScalarArray(np.array(wavelengths) * u.AA, axes="wavelength")
    sigma = sdo.aia.uncertainty(
        counts=na.ScalarArray(counts, axes=axes) * u.DN,
        wavelength=channels,
        n_sample=n_sample,
    )
    sigma = u.Quantity(sigma.ndarray_aligned(axes)).to_value(u.DN)
    rate = counts / exptime[:, np.newaxis, np.newaxis]
    # a saturated channel leaves the pixel unfit, and painted white
    rate[counts > saturation] = np.nan
    sigma = sigma / exptime[:, np.newaxis, np.newaxis]
    i94 = wavelengths.index(94)
    sigma[i94] = np.hypot(sigma[i94], systematic_94)
    dem, chi2 = utu.dem.plowman(
        intensity=na.ScalarArray(rate, axes=axes) * u.DN / u.s / u.pix,
        uncertainty=na.ScalarArray(sigma, axes=axes) * u.DN / u.s / u.pix,
        response=response,
        axis_channel="wavelength",
        axis_temperature="temperature",
    )
    values = dem.outputs.to(1e26 / u.cm**5).ndarray_aligned(("y", "x", "temperature"))
    return u.Quantity(values).value[..., shown], chi2.ndarray_aligned(("y", "x"))


def _stretch(logt: u.Quantity | np.ndarray) -> np.ndarray:
    """The scale the colors are spread evenly over; see :data:`stretch_width`."""
    return np.arctan((u.Quantity(logt).value - stretch_center) / stretch_width)


def _rgb(dem: np.ndarray, logt: np.ndarray) -> np.ndarray:
    """The colorsynth RGB image of a DEM cube, on the fixed scale."""
    rgb = colorsynth.rgb(
        np.nan_to_num(dem),
        wavelength=logt * u.dimensionless_unscaled,
        axis=-1,
        spd_min=0,
        spd_max=vmax,
        wavelength_norm=_stretch,
    )
    return np.clip(np.asarray(rgb), 0, 1)


def _save(rgb: np.ndarray, site: pathlib.Path, name: str) -> None:
    """Save a frame at every size, with the color at full resolution."""
    # row 0 of the array is the bottom of the Sun, the top of an image
    image = Image.fromarray(np.round(255 * rgb[::-1]).astype(np.uint8))
    for level in levels:
        scaled = image if level == size else image.resize((level, level), Image.Resampling.LANCZOS)
        path = site / str(level) / f"{name}.jpg"
        part = path.with_suffix(".part")
        scaled.save(part, format="JPEG", quality=90, subsampling=0)
        part.replace(path)


def _colorbar(logt: np.ndarray, path: pathlib.Path, wide: bool) -> None:
    """
    Save the color key on a transparent background: DEM along x and
    temperature along y, or the other way around if `wide`.
    """
    fraction = np.linspace(0, 1, 101)[:, np.newaxis]
    x, y, rgb = colorsynth.colorbar(
        fraction * (vmax * np.ones(logt.size)),
        wavelength=logt * u.dimensionless_unscaled,
        axis=-1,
        spd_min=0,
        spd_max=vmax,
        wavelength_norm=_stretch,
    )
    # spaced as the colors are, so that each gets its share of the key
    y = _stretch(y)
    label_dem = "DEM ($10^{26}\\,\\mathrm{cm^{-5}}$)"
    label_t = "temperature (MK)"
    # labeled in MK, which reads more easily than log T
    ticks = np.array([1, 1.5, 2, 3, 5, 10])
    positions = _stretch(np.log10(ticks * 1e6))
    if wide:
        fig, ax = plt.subplots(figsize=(4, 1.4), dpi=100, constrained_layout=True)
        ax.pcolormesh(y, x, rgb)
        ax.set_xticks(positions, [f"{t:g}" for t in ticks])
        ax.set_xlabel(label_t, color="white", fontsize=9)
        ax.set_ylabel("DEM\n($10^{26}\\,\\mathrm{cm^{-5}}$)", color="white", fontsize=8)
    else:
        fig, ax = plt.subplots(figsize=(1.6, 4), dpi=100, constrained_layout=True)
        ax.pcolormesh(x, y, rgb)
        ax.set_yticks(positions, [f"{t:g}" for t in ticks])
        ax.set_xlabel(label_dem.replace(" ", "\n", 1), color="white", fontsize=8)
        ax.set_ylabel(label_t, color="white", fontsize=9)
    ax.tick_params(colors="white", labelsize=8)
    fig.savefig(path, transparent=True)
    plt.close(fig)


def _logt() -> np.ndarray:
    """The log T of each temperature of the DEM shown, from the response table alone."""
    response = sdo.aia.temperature_response(wavelength=[94] * u.AA)
    return np.log10(response.inputs.ndarray.to_value(u.K))[grid][shown]


def _fill(
    site: pathlib.Path,
    slots: dict[str, datetime.datetime],
    have: dict[str, dict],
    logt: np.ndarray,
    end: datetime.datetime,
    deadline: float,
) -> int:
    """
    Make the frames of the slots not in `have`, newest first, until
    `deadline` (of :func:`time.monotonic`), adding each to `have`. Return
    the number of slots that failed, not counting those JSOC lacks.
    """
    missing = [name for name in slots if name not in have]
    print(f"{len(have)} of {len(slots)} frames on hand, {len(missing)} to try", flush=True)
    response = _response(end)
    coordinate = (np.arange(size) - (size - 1) / 2) * plate_scale
    on_disk = np.hypot(*np.meshgrid(coordinate, coordinate)) < disk
    failures = 0
    for name in missing:
        if time.monotonic() > deadline:
            print("out of time; the rest are left for the next run", flush=True)
            break
        slot = slots[name]
        t0 = time.monotonic()
        try:
            files = _download(slot)
            if files is None:
                print(f"{name} not on JSOC (yet)", flush=True)
                continue
            t1 = time.monotonic()
            counts, exptime = _read(files)
            dem, chi2 = _invert(counts, exptime, response)
            t2 = time.monotonic()
            good = chi2 >= 0
            if not good.any():
                raise ValueError("the inversion failed on every pixel")
            rgb = _rgb(dem, logt)
            rgb[~good & on_disk] = 1
            _save(rgb, site, name)
        except Exception as error:
            failures += 1
            print(f"{name} failed: {error!r}", flush=True)
            continue
        have[name] = dict(
            time=f"{slot:%Y-%m-%dT%H:%M}:00Z",
            name=name,
            exptime_min=round(float(exptime.min()), 3),
            chi2_median=round(float(np.median(chi2[good])), 3),
            failed=int((~good).sum()),
        )
        print(
            f"{name} download {t1 - t0:.1f} s, dem {t2 - t1:.1f} s, render {time.monotonic() - t2:.1f} s; "
            f"median chi2 {have[name]['chi2_median']:.2f}, {have[name]['failed']} pixels failed",
            flush=True,
        )
    return failures


def _manifest(site: pathlib.Path, have: dict[str, dict], logt: np.ndarray, now: datetime.datetime) -> None:
    """Write the manifest the viewer reads, listing the frames in `have`."""
    manifest = dict(
        title="AIA DEM, last 48 hours",
        source=(
            "JSOC AIA near-real-time synoptic images; utu.dem.plowman with the SolarSoft V10 "
            f"response normalized to EVE, with chiantifix, and {systematic_94} DN/s added to the "
            f"uncertainty of 94; colorsynth, spread over arctan((log T - {stretch_center}) / {stretch_width})"
        ),
        version=version,
        updated=f"{now:%Y-%m-%dT%H:%M:%S}Z",
        center=[0, 0],
        scale={str(level): plate_scale * size / level for level in levels},
        size={str(level): level for level in levels},
        logt=[round(float(logt[0]), 2), round(float(logt[-1]), 2)],
        vmax=vmax,
        presets=[dict(
            id="dem",
            label="DEM (10^26 cm^-5)",
            colorbar="colorbar.png",
            colorbar_wide="colorbar_wide.png",
            path="{level}/{name}.jpg",
        )],
        frames=[have[name] for name in sorted(have)],
    )
    path = site / "manifest.json"
    part = path.with_suffix(".part")
    part.write_text(json.dumps(manifest, indent=1))
    part.replace(path)


def main() -> None:
    """Fill in the missing frames of the window, prune old ones, and write the manifest."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("site", type=pathlib.Path, help="the directory the site is built in")
    parser.add_argument("--hours", type=float, default=48, help="the length of the window")
    parser.add_argument("--budget", type=float, default=40, help="minutes after which no new frame is started")
    parser.add_argument("--end", help="the end of the window, UTC, if not now")
    args = parser.parse_args()
    start = time.monotonic()
    site: pathlib.Path = args.site
    for level in levels:
        (site / str(level)).mkdir(parents=True, exist_ok=True)

    now = datetime.datetime.now(datetime.UTC).replace(tzinfo=None)
    end = datetime.datetime.fromisoformat(args.end) if args.end else now
    slots = {f"{s:%Y%m%d_%H%M}": s for s in _slots(end, args.hours)}
    previous = {}
    current = True
    if (site / "manifest.json").exists():
        manifest = json.loads((site / "manifest.json").read_text())
        previous = {f["name"]: f for f in manifest["frames"]}
        current = manifest.get("version", 1) == version
        if not current:
            print(f"the frames are version {manifest.get('version', 1)}; making all of them again", flush=True)
    have = {
        name: previous.get(name, dict(time=f"{slot:%Y-%m-%dT%H:%M}:00Z", name=name))
        for name, slot in slots.items()
        if current and all((site / str(level) / f"{name}.jpg").exists() for level in levels)
    }
    logt = _logt()
    for wide in (False, True):
        _colorbar(logt, site / ("colorbar_wide.png" if wide else "colorbar.png"), wide)

    # however the run ends, the site it leaves is consistent
    failures = 0
    try:
        failures = _fill(site, slots, have, logt, end, start + 60 * args.budget)
    finally:
        for path in site.glob("*/*.*"):
            if path.suffix == ".part" or path.stem not in have:
                path.unlink()
        _manifest(site, have, logt, now)
        print(f"{len(have)} frames in the manifest, {failures} failed, {time.monotonic() - start:.0f} s", flush=True)
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
