# Readout-adaptive event-based velocimetry (RAEBV)

Python code for reconstructing velocity fields from event-camera recordings under continuous illumination.

Code accompanying **Readout-adaptive event-based velocimetry under continuous illumination**, by **Dongmin Sun and Hyungmin Park**, Multiphase Flow & Flow Visualization Laboratory, Department of Mechanical Engineering, Seoul National University.

![Real-event jet reconstruction](assets/sample.gif)

*Round jet, Re = 3,500. (a) Events accumulated over 4 readout segments (px). (b) Selected layer: cyan = Layer 0, magenta = Layer 1. (c) Velocity vectors colored by speed (m/s; positions in m). Equal-length arrows indicate direction. Horizontally mirrored; slow-motion playback.*

## Installation

Download or clone this repository and run the commands below from its root folder. Tested with Python 3.9.13 on Windows 10; no GPU is required.

```bash
python -m pip install -r requirements.txt
```

For recordings in .raw or .hdf5 format recorded with Metavision SDK, also install the [Metavision SDK](https://docs.prophesee.ai/stable/installation/index.html) with compatible Python bindings (tested with SDK 4.6.2).

Sample recordings are not currently included.

## Run on your data

Set the ROI, spatial calibration, layers, cell size and stride in `CONFIG` in [configuration.py](configuration.py). Flow settings from the paper are provided as commented examples.

Set `polarity` to `"ON"`, `"OFF"` or `"BOTH"`. In `"BOTH"` mode, triplets are matched within each polarity and combined into cell velocities.

Successive matching windows are advanced by a prescribed stride, expressed as an integer number of readout segments. The minimum possible stride is one segment; larger strides may be selected based on the time interval required for a suitable particle displacement between successive matching windows.

First, calculate and save a specified number of instantaneous velocity fields:

```bash
python main.py --input /path/to/recording.raw --output results/my_run --windows 10000 --workers 4
```

Alternatively, use `--duration` to process a specified recording duration in seconds:

```bash
python main.py --input /path/to/recording.hdf5 --output results/my_run --duration 10 --workers 4
```

Each run saves instantaneous fields, timestamps and settings. See the [usage guide](docs/usage.md) for speed-gate calibration and shorter recordings.

Then compute flow statistics from the saved fields:

```bash
python postprocess.py --input results/my_run --workers 4
```

Mean velocities, RMS fluctuations and Reynolds shear stress are saved in `results/my_run/statistics/`. Postprocessing uses the saved cell velocities directly, without additional spatial averaging. `--workers` selects the number of CPU worker processes.

See the [usage guide](docs/usage.md) for input formats, additional options and output units.

## License

[MIT License](LICENSE).
