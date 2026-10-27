"""Edit CONFIG directly for your recording. Examples below are comments only."""

CONFIG = {
    "roi": {"x": 0, "y": 0, "width": 1280, "height": 720},
    # Cell side length c in pixels; the statistics grid step also equals c.
    "cell_size_px": 16,
    "trajectory_tolerance_px": 0.5,
    # Each layer: [s, g, dx_min, dx_max, dy_min, dy_max].
    # s = segments per event group; g = skipped segments between groups.
    # A matching window spans 3*s + 2*g segments; offsets use image axes.
    "layers": [[4, 4, 0, 7, -2, 2], [8, 8, -1, 1, -1, 1]],
    "stride_segments": 8,
    "polarity": "ON",
    "row_jump_threshold_px": 200,
    "spatial_scale_m_per_px": 4.8e-05,
    # Speed-gate calibration is independent of the number of saved fields.
    "calibration_windows": 5000,
    "histogram_max_px_per_ms": 30.0,
    "histogram_bins": 300,
    "histogram_smoothing_bins": 5,
    # Default field count when neither --windows nor --duration is supplied.
    "max_windows": 100000,
    # Minimum valid windows per cell for the statistics tables.
    "min_valid_samples": 10,
}

# Paper examples for reference. Copy the desired values into CONFIG above.
# These commented values are not loaded or selected by the program.
# The current CONFIG values correspond to the round-jet example.

# Pipe flow, Re = 5,300:
# {
#     "roi": {"x": 0, "y": 0, "width": 1280, "height": 360},
#     "layers": [[22, 11, -1, 1, 2, 8], [22, 66, 0, 0, 1, 4]],
#     "stride_segments": 20,
#     "max_windows": 200000,
#     "min_valid_samples": 1,
#     "spatial_scale_m_per_px": 2.2222222222222223e-05,
# }

# Pipe flow, Re = 19,000:
# {
#     "roi": {"x": 0, "y": 0, "width": 1280, "height": 360},
#     "layers": [[6, 3, -1, 1, 2, 8], [6, 18, 0, 0, 1, 4]],
#     "stride_segments": 6,
#     "max_windows": 200000,
#     "min_valid_samples": 1,
#     "spatial_scale_m_per_px": 2.2222222222222223e-05,
# }

# Round jet, Re = 3,500:
# {
#     "roi": {"x": 0, "y": 0, "width": 1280, "height": 720},
#     "layers": [[4, 4, 0, 7, -2, 2], [8, 8, -1, 1, -1, 1]],
#     "stride_segments": 8,
#     "max_windows": 100000,
#     "min_valid_samples": 10,
#     "spatial_scale_m_per_px": 4.8e-05,
# }

# Cylinder wake, Re = 3,900:
# {
#     "roi": {"x": 0, "y": 0, "width": 1280, "height": 720},
#     "layers": [[4, 8, -4, 4, -4, 4], [12, 24, -2, 2, -2, 2]],
#     "stride_segments": 12,
#     "max_windows": 100000,
#     "min_valid_samples": 2,
#     "spatial_scale_m_per_px": 5.4945054945054945e-05,
# }
