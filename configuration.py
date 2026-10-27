"""Shared settings and commented flow selection for the paper examples."""
from copy import deepcopy

# Select ONE flow: leave exactly one ACTIVE_FLOW line uncommented.
# ACTIVE_FLOW = "pipe5300"   # Pipe flow, Re = 5,300
# ACTIVE_FLOW = "pipe19000"  # Pipe flow, Re = 19,000
ACTIVE_FLOW = "jet3500"     # Round jet, Re = 3,500
# ACTIVE_FLOW = "wake3900"   # Cylinder wake, Re = 3,900

# Shared numerical settings. c is the cell side length in pixels.
# The statistics grid step equals c; no additional spatial averaging is applied.
COMMON = {
    "cell_size_px": 16,
    "trajectory_tolerance_px": 0.5,
    "polarity": "ON",
    "row_jump_threshold_px": 200,
    "calibration_windows": 5000,
    "histogram_max_px_per_ms": 30.0,
    "histogram_bins": 300,
    "histogram_smoothing_bins": 5,
}

# Each layer: [s, g, dx_min, dx_max, dy_min, dy_max].
# s = segments per event group; g = skipped segments between groups.
# Search offsets use the original image axes. A window has 3*s + 2*g segments.
# spatial_scale_m_per_px converts pixel distances to metres.
# max_windows is used when neither --windows nor --duration is supplied.
FLOWS = {
    # Pipe flow, Re = 5,300
    "pipe5300": {
        "flow": "pipe",
        "reynolds_number": 5300,
        "roi": {"x": 0, "y": 0, "width": 1280, "height": 360},
        "layers": [[22, 11, -1, 1, 2, 8], [22, 66, 0, 0, 1, 4]],
        "stride_segments": 20,
        "max_windows": 200000,
        "min_valid_samples": 1,
        "spatial_scale_m_per_px": 2.2222222222222223e-05,
    },
    # Pipe flow, Re = 19,000
    "pipe19000": {
        "flow": "pipe",
        "reynolds_number": 19000,
        "roi": {"x": 0, "y": 0, "width": 1280, "height": 360},
        "layers": [[6, 3, -1, 1, 2, 8], [6, 18, 0, 0, 1, 4]],
        "stride_segments": 6,
        "max_windows": 200000,
        "min_valid_samples": 1,
        "spatial_scale_m_per_px": 2.2222222222222223e-05,
    },
    # Round jet, Re = 3,500
    "jet3500": {
        "flow": "jet",
        "reynolds_number": 3500,
        "roi": {"x": 0, "y": 0, "width": 1280, "height": 720},
        "layers": [[4, 4, 0, 7, -2, 2], [8, 8, -1, 1, -1, 1]],
        "stride_segments": 8,
        "max_windows": 100000,
        "min_valid_samples": 10,
        "spatial_scale_m_per_px": 4.8e-05,
    },
    # Cylinder wake, Re = 3,900
    "wake3900": {
        "flow": "wake",
        "reynolds_number": 3900,
        "roi": {"x": 0, "y": 0, "width": 1280, "height": 720},
        "layers": [[4, 8, -4, 4, -4, 4], [12, 24, -2, 2, -2, 2]],
        "stride_segments": 12,
        "max_windows": 100000,
        "min_valid_samples": 2,
        "spatial_scale_m_per_px": 5.4945054945054945e-05,
    },
}


def get_config(flow=None):
    """Return an independent effective configuration for the selected flow."""
    selected = ACTIVE_FLOW if flow is None else flow
    if selected not in FLOWS:
        raise ValueError(f"Unknown flow {selected!r}; choose one of {', '.join(FLOWS)}")
    return deepcopy({**COMMON, "name": selected, **FLOWS[selected]})
