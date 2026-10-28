# Usage guide

See the [README](../README.md) for installation. Run commands from the repository root using your own recording.

## Input

Use a Metavision `.raw` recording, an SDK-converted `.hdf5` recording, or an `.npz` file containing a structured NumPy array named `events`:

| Field | Format |
|---|---|
| `x`, `y` | Sensor pixel coordinates |
| `t` | Timestamp in microseconds, in nondecreasing order |
| `p` | Polarity: 0 = OFF, 1 = ON |

Recordings must contain a detectable row-wise readout pattern. RAW and HDF5 files require the Metavision SDK; NPZ files do not. HDF5 input must use the SDK event format and the lowercase `.hdf5` extension.

To convert a RAW recording with the SDK and process it:

```bash
metavision_file_to_hdf5 -i recording.raw -o recording.hdf5
python main.py --input recording.hdf5 --output results/hdf5_run --windows 10000 --workers 4
```

## Instantaneous velocity fields

Enter your processing parameters in `CONFIG` in [configuration.py](../configuration.py). Spatial calibration is in m/px; flow settings from the paper are provided as commented examples.

The implementation uses two temporal layers. Set `polarity` to `"ON"` for ON events (`p = 1`), `"OFF"` for OFF events (`p = 0`), or `"BOTH"` to use both. In `"BOTH"` mode, each triplet contains events of the same polarity; accepted ON and OFF triplets contribute to the cell velocity. Readout segmentation uses all events within the ROI in every mode. Parameter units and ranges are described beside each setting in `configuration.py`.

Recalculate the speed gate after changing polarity. A saved `--gate-file` must match the selected polarity and other configuration values. Automatic calibration requires events of the selected polarity and valid cell velocities in both layers.

```bash
python main.py --input recording.raw --output results/my_run --windows 10000 --workers 4
```

| Option | Purpose |
|---|---|
| `--config path/to/configuration.py` | Use another configuration file; saved `config.json` files are also accepted |
| `--windows N` | Save exactly N instantaneous fields |
| `--duration S` | Save complete matching windows within the first S recorded seconds |
| `--workers N` | Set CPU worker count; default is 1 |
| `--gate-file path/to/gate.json` | Reuse calibration for a matching configuration and the same recording |

Choose either `--windows` or `--duration`; if neither is given, `CONFIG["max_windows"]` is used. Duration is measured from timestamp zero in the decoded recording, not wall-clock processing time. Only complete matching windows ending before the time limit are saved. The actual count and covered duration are recorded in `complete.json`. With `--duration`, a shorter recording ends at its last complete window; `--windows N` requires at least N complete windows.

`main.py` saves one velocity field per matching window. Run `postprocess.py` afterward to compute statistics.

Speed-gate calibration uses the first `calibration_windows` windows (default: 5,000). For a shorter output selection, calibration can read beyond the requested number or duration; only the requested instantaneous fields are saved. If the recording is too short for calibration, provide a compatible `--gate-file` from the same recording.

Rerun the same command to resume completed blocks. Use a new output folder when changing input, settings, selection or code. Use `python main.py --help` for all options.

Instantaneous-field outputs:

| File | Contents |
|---|---|
| `config.json` | Settings used |
| `frame_times.npy` | Window start and end timestamps in µs |
| `vectors/000000.npy`, ... | One `[x, y, u, v, layer_id]` array per window; positions in px, velocities in px/s |
| `calibration/gate.json` | Calibration for reuse |
| `complete.json` | Completion summary and analyzed duration |

Additional files in `calibration/` and `blocks/` store calibration data and processing progress. Each vector row represents one valid cell; empty cells remain missing. Coordinates are local to the configured ROI and follow the original image-axis directions.

## Separate postprocessing

```bash
python postprocess.py --input results/my_run --workers 4
```

This computes cell-wise statistics from the saved vectors. The original recording and Metavision SDK are not needed.

Results are written to `results/my_run/statistics/` by default. Use `--output path/to/statistics` to choose another folder. Repeating the same command resumes completed statistics blocks; use a new output folder if the saved vectors or settings change.

| File | Contents |
|---|---|
| `statistics_px.csv` | Mean velocity, RMS, shear-stress term, energy and valid count, in pixel-based units |
| `statistics_si.csv` | Same statistics in m, m/s and squared velocity units |
| `layer_fractions.npz` | Layer fractions per window and pooled vector counts |
| `moments.npz` | Accumulated counts, means and fluctuation moments |
| `config.json`, `complete.json` | Settings and completion summary |

Statistics use the saved cell velocities directly, with equal weight for each valid window and no additional spatial averaging. CSV positions are cell centers. SI conversion uses the saved spatial scale. `uv_stress` is negative covariance, −〈u′v′〉, and `ke` is ½(u_rms² + v_rms²). `run.json` and `blocks/` support resuming postprocessing.
