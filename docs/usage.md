# Usage guide

See the [README](../README.md) for installation. Run commands from the repository root using your own recording.

## Input

Use a Metavision `.raw` recording or an `.npz` file containing a structured NumPy array named `events`:

| Field | Format |
|---|---|
| `x`, `y` | Sensor pixel coordinates |
| `t` | Timestamp in microseconds, in nondecreasing order |
| `p` | Polarity: 0 = OFF, 1 = ON |

Recordings must contain a detectable row-wise readout pattern. RAW files require the Metavision SDK; NPZ files do not.

## Configuration and execution

In [configuration.py](../configuration.py), select one `ACTIVE_FLOW` and adjust `COMMON` and the corresponding `FLOWS` entry. Comments explain the ROI, layer settings, cell size, stride, window count and spatial scale in m/px.

```bash
python main.py --input recording.raw --output results/my_run --workers 4 --save-vectors
```

| Option | Purpose |
|---|---|
| `--flow pipe5300` | Override the selected flow |
| `--config path/to/configuration.py` | Use another configuration file; saved `config.json` files are also accepted |
| `--windows N` | Override the window count; use at least 5,000 for a new calibration with the supplied presets |
| `--workers N` | Set CPU worker count; default is 1 |
| `--save-vectors` | Save vectors for each window; omit for statistics only |
| `--gate-file path/to/gate.json` | Reuse calibration for a matching configuration and the same recording |

Rerun the same command to resume completed blocks. Use a new output folder when changing the input, settings, vector-output option or code. Use `python main.py --help` for all options.

## Outputs

| File | Contents |
|---|---|
| `config.json` | Settings used |
| `frame_times.npy` | Window start and end timestamps in µs |
| `vectors/*.npy` | Optional `[x, y, u, v, layer_id]`; positions in px, velocities in px/s |
| `statistics_px.csv` | Statistics in px and px/s |
| `statistics_si.csv` | Statistics in m and m/s, with squared velocity units for covariance and energy |
| `calibration/gate.json` | Calibration for reuse |
| `complete.json` | Completion summary and analyzed duration |

Additional files support calibration and resuming runs. Coordinates follow the original image axes; SI conversion uses the configured spatial scale. The vectors are cell estimates, and empty cells remain missing. In the statistics tables, `uv_stress` is negative covariance, −〈u′v′〉, and `ke` is ½(u_rms² + v_rms²).
