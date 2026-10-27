# Readout-adaptive event-based velocimetry (EBV)

Python code for reconstructing velocity fields from event-camera recordings under continuous illumination.

Code accompanying **Readout-adaptive event-based velocimetry under continuous illumination**, by **Dongmin Sun and Hyungmin Park**, Multiphase Flow & Flow Visualization Laboratory, Seoul National University.

![Real-event jet reconstruction](assets/sample.gif)

*Round jet, Re = 3,500. (a) Events accumulated over 4 readout segments (axes in px). (b) Selected layer: cyan = Layer 0, magenta = Layer 1. (c) Velocity field colored by speed (m/s; axes in m). Equal-length arrows show direction. All valid cells are displayed; panels are mirrored horizontally and played in slow motion.*

## Installation

Download or clone this repository and run the commands below from its root folder. Tested with Python 3.9.13 on Windows 10; no GPU is required.

```bash
python -m pip install -r requirements.txt
```

For recordings in .raw or .hdf5 format recorded with Metavision SDK, also install the [Metavision SDK](https://docs.prophesee.ai/stable/installation/index.html) with compatible Python bindings (tested with SDK 4.6.2).

Sample recordings are not currently included.

## Run on your data

All settings are in [configuration.py](configuration.py), with presets for pipe, jet and cylinder-wake flows. Leave one `ACTIVE_FLOW` line uncommented. Adjust the ROI, spatial calibration, layer settings, cell size and stride for your recording.

Successive matching windows are advanced by a prescribed stride, expressed as an integer number of readout segments. The minimum possible stride is one segment; larger strides may be selected based on the time interval required for a suitable particle displacement between successive matching windows.

First, calculate and save a specified number of instantaneous velocity fields:

```bash
python main.py --input /path/to/recording.raw --output results/my_run --windows 10000 --workers 4
```

Alternatively, use `--duration` to process a specified recording duration in seconds:

```bash
python main.py --input /path/to/recording.hdf5 --output results/my_run --duration 10 --workers 4
```

Choose one command. MAIN always saves the instantaneous fields, timestamps and settings. RAW and SDK-converted HDF5 files use the same calculation. Speed-gate calibration uses the first `calibration_windows` windows (5,000 in the presets), independently of the number or duration of fields saved.

Then compute flow statistics from the saved fields:

```bash
python postprocess.py --input results/my_run --workers 4
```

Mean velocities, RMS fluctuations and Reynolds shear stress are saved in `results/my_run/statistics/`. Postprocessing uses the saved cell velocities directly, without additional spatial averaging. `--workers` selects the number of CPU worker processes.

See the [usage guide](docs/usage.md) for input formats, additional options and output units.

## License

[MIT License](LICENSE). The example GIF uses this project's own experimental data.
