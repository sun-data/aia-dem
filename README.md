# aia-dem

A rolling 48-hour movie of the differential emission measure (DEM) of the
solar corona, from the near-real-time images of
[SDO](https://sdo.gsfc.nasa.gov)/AIA, a frame every 3 minutes, updated
every 15 minutes.

**Watch it: <https://sun-data.github.io/solar-movie-viewer/>, the DEM chip.**

The DEM says how much plasma there is along each line of sight at each
temperature. Each frame shows it as color: hue is temperature, and brightness
is the amount of plasma. Most of the spectrum goes to 1.4-3 MK, where most
of the corona is: violet below 1.4 MK, blue near 1.6 MK, green near 2 MK,
yellow and orange near 2.5 MK, and red above 3 MK. One fixed scale covers
every temperature and every frame,
so the brightness and color of any two pixels compare, across the disk and
across the 48 hours.

## How it works

About every 15 minutes a GitHub Actions workflow runs [`build.py`](build.py),
which fills in every 3-minute slot of the last 48 hours the site does not
have yet, newest first:

1. It downloads the six coronal channels (94, 131, 171, 193, 211 and 335 Å)
   of the AIA synoptic images from the JSOC near-real-time archive. These are
   1024 x 1024 at 2.4 arcsec per pixel, each pixel the mean of a 4 x 4 block
   of the full-resolution image, and appear about 11 minutes after they are
   taken.
2. It inverts them for a DEM with
   [`utu.dem.plowman`](https://github.com/sun-data/utu), the method of
   Plowman & Caspi (2020), from log T 5.0 to 7.5, using the SolarSoft V10
   temperature response with the EVE normalization and the empirical
   correction of 94 Å
   ([`sdo.aia.temperature_response`](https://github.com/sun-data/solar-dynamics-observatory)),
   and the noise SolarSoft estimates for each channel (shot noise, read
   noise, dark subtraction, quantization and compression;
   `sdo.aia.uncertainty`), reduced for the 16 pixels each synoptic pixel
   averages.
3. It colors the DEM from log T 5.6 to 7.2 with
   [`colorsynth`](https://github.com/sun-data/colorsynth), and saves it as
   JPEG at 1024 and 512 pixels. The colors are spread evenly over
   arctan((log T - 6.3) / 0.15) rather than over log T, which would leave
   1.4-2.5 MK all green.

Frames older than 48 hours are deleted, and `manifest.json` lists the rest,
which is what the viewer reads. The frames are carried from run to run in the
Actions cache, and the site is published to GitHub Pages at
<https://sun-data.github.io/aia-dem/>. A frame takes about 30 seconds on
the runner, and each run stops starting new frames after 40 minutes, so
after a long gap the window fills in over several runs: all 48 hours, 960
frames, take about eight hours. The manifest gives the `version` of the
frames, which goes up in `build.py` whenever a change alters how they look;
a run that finds frames of another version deletes them and makes the whole
window again, newest first, rather than mix the two.

GitHub drops most runs of a frequent schedule, so each run starts the next
itself. A run that took under 10 minutes, as one does once the window is
filled in, has the next wait 10 minutes first, in an environment (`pacing`)
whose wait timer holds the run without holding a runner; a longer one has it
start at once. A schedule every 15 minutes restarts the chain if it breaks,
and a late run makes up whatever was missed. Disabling the workflow stops it.

## Choices

* **94 Å gets 1 DN/s of extra uncertainty.** The response of 94 Å to plasma
  at 1-2 MK is too small, so a fit that matches 94 Å exactly explains the
  faint 94 Å of coronal holes with plasma near 6 MK, and paints them red.
  Adding 1 DN/s in quadrature to its uncertainty lets the fit miss 94 Å by
  about that much. Coronal holes then come out dark, the quiet Sun is
  unchanged, and active regions gain about a third more emission at 2.5-5.6
  MK. Leaving 94 Å out altogether does as well for the coronal holes but
  puts seven times as much plasma at those temperatures in active regions.
* **The scale is fixed**, at 195 x 10^26 cm^-5 per unit log T: the median,
  over a dozen frames from 2026-10-15 to 2026-10-17, of the 99.9th
  percentile of the DEM over every pixel and temperature.
* **Unfit pixels are white.** Where a channel saturates, as in the core of a
  flare, the fit has nothing to go on and the pixel is painted white.

## Limitations

* The resolution is that of the synoptic images, 2.4 arcsec per pixel, a
  quarter of AIA's.
* The near-real-time images are preliminary. Their pointing and calibration
  are not final.
* A DEM from six channels is a smooth estimate, regularized toward a DEM
  whose logarithm changes slowly with temperature. It cannot resolve
  structure finer than a few tenths in log T.

## Running it locally

```sh
pip install -r requirements.txt
python build.py site --hours 6
```

`--end` builds the window ending at another time; the near-real-time
archive goes back to at least 2024. To watch the result, serve `site` with
CORS headers and open a local copy of the viewer with `?dem=` pointing at
it, for example `index.html?dem=http://localhost:8001/`.
