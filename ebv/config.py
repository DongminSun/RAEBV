"""Explicit, portable settings for the four final paper cases."""
import json
import runpy
from pathlib import Path


def load_config(path, flow=None):
    path = Path(path)
    if path.suffix.lower() == ".py":
        module = runpy.run_path(str(path))
        get_config = module.get("get_config")
        if not callable(get_config):
            raise ValueError("Python configuration must define get_config(flow=None)")
        c = get_config(flow)
    else:
        # Effective JSON snapshots from earlier runs remain readable.
        with path.open(encoding="utf-8") as fh:
            c = json.load(fh)
        if flow is not None and c.get("name") != flow:
            raise ValueError("Selected flow does not match the JSON configuration")
    required = {
        "name", "flow", "reynolds_number", "roi", "layers", "stride_segments",
        "cell_size_px", "trajectory_tolerance_px", "polarity",
        "row_jump_threshold_px", "calibration_windows", "histogram_max_px_per_ms",
        "histogram_bins", "histogram_smoothing_bins", "max_windows",
        "min_valid_samples", "spatial_scale_m_per_px",
    }
    if set(c) != required:
        raise ValueError(f"Configuration keys differ: missing={required-set(c)}, extra={set(c)-required}")
    if c["polarity"] != "ON":
        raise ValueError("The released paper pipeline uses ON-event matching.")
    for name in ("stride_segments", "cell_size_px", "row_jump_threshold_px",
                 "calibration_windows", "histogram_bins", "histogram_smoothing_bins",
                 "max_windows", "min_valid_samples"):
        if type(c[name]) is not int or c[name] <= 0:
            raise ValueError(f"{name} must be a positive integer")
    roi = c["roi"]
    if set(roi) != {"x", "y", "width", "height"}:
        raise ValueError("ROI requires x, y, width and height")
    if any(type(v) is not int for v in roi.values()) or min(roi["x"], roi["y"]) < 0 or min(roi["width"], roi["height"]) <= 0:
        raise ValueError("Invalid ROI")
    if len(c["layers"]) != 2:
        raise ValueError("Paper presets require exactly two layers")
    for layer in c["layers"]:
        if len(layer) != 6 or any(type(v) is not int for v in layer):
            raise ValueError("Each layer must be [s, g, xmin, xmax, ymin, ymax]")
        s, g, xmin, xmax, ymin, ymax = layer
        if s <= 0 or g < 0 or xmin > xmax or ymin > ymax:
            raise ValueError(f"Invalid layer: {layer}")
    if c["stride_segments"] > max(3*s + 2*g for s, g, *_ in c["layers"]):
        raise ValueError("Stride must not exceed the longest matching window")
    for key in ("trajectory_tolerance_px", "histogram_max_px_per_ms", "spatial_scale_m_per_px"):
        if not 0 < float(c[key]) < float("inf"):
            raise ValueError(f"Invalid {key}")
    return c


def write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temp.replace(path)
