"""Load and validate user-provided processing parameters."""
from copy import deepcopy
import json
import runpy
from pathlib import Path


def load_config(path):
    path = Path(path)
    if path.suffix.lower() == ".py":
        module = runpy.run_path(str(path))
        if "CONFIG" not in module:
            raise ValueError("Python configuration must define a CONFIG dictionary")
        c = module["CONFIG"]
    else:
        # Load settings saved by a processing run.
        with path.open(encoding="utf-8") as fh:
            c = json.load(fh)
    return validate_config(c)


def validate_config(c):
    if not isinstance(c, dict):
        raise ValueError("Configuration must be a dictionary")
    required = {
        "roi", "layers", "stride_segments",
        "cell_size_px", "trajectory_tolerance_px", "polarity",
        "row_jump_threshold_px", "calibration_windows", "histogram_max_px_per_ms",
        "histogram_bins", "histogram_smoothing_bins", "max_windows",
        "min_valid_samples", "spatial_scale_m_per_px",
    }
    # Ignore descriptive flow labels in saved configurations.
    legacy_labels = {"name", "flow", "reynolds_number"}
    missing, extra = required-set(c), set(c)-required-legacy_labels
    if missing or extra:
        raise ValueError(f"Configuration keys differ: missing={missing}, extra={extra}")
    c = deepcopy({k: v for k, v in c.items() if k not in legacy_labels})
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
        raise ValueError("Matching requires exactly two temporal layers")
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
