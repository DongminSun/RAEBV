"""Processing parameters for RAEBV."""

# Default settings: round jet, Re = 3,500.

CONFIG = {
    # Sensor ROI in pixels: top-left (x, y), width and height.
    # Integers; x/y >= 0, width/height > 0. Keep the ROI inside the sensor.
    "roi": {"x": 0, "y": 0, "width": 1280, "height": 720},

    # Cell side length c in pixels (positive integer), measured from the ROI origin.
    # Matching outputs and statistics use the same c-by-c grid.
    "cell_size_px": 16,

    # Radial distance tolerance in pixels (> 0) for accepting the oldest event
    # against the trajectory through the middle and current events.
    "trajectory_tolerance_px": 0.5,

    # Exactly two layers, ordered Layer 0 (shorter interval), Layer 1 (longer interval).
    # Each entry: [s, g, dx_min, dx_max, dy_min, dy_max], all integers.
    # s >= 1: segments per event group; g >= 0: skipped segments between groups.
    # Each layer spans 3*s + 2*g segments. dx/dy are inclusive search offsets in px
    # from a current event to a middle-group candidate, not velocity components.
    # Image axes: +x right, +y down; each minimum must be <= its maximum.
    "layers": [[4, 4, 0, 7, -2, 2], [8, 8, -1, 1, -1, 1]],

    # Matching-window stride in readout segments (integer >= 1).
    # Choose based on the desired sampling interval and particle displacement between successive windows.
    # Must not exceed max(3*s + 2*g) over the two layers.
    "stride_segments": 8,

    # "ON": p = 1 only; "OFF": p = 0 only; "BOTH": use both event polarities.
    # In BOTH mode, each triplet contains events of the same polarity.
    # Cell velocities combine the accepted ON and OFF triplets.
    "polarity": "ON",

    # Readout reset threshold in pixels (positive integer): split segments when
    # previous_y - next_y > this value. Adapt to the row-reset jump in your recording.
    # Reset detection uses all input events within the ROI, independent of matching polarity.
    "row_jump_threshold_px": 200,

    # Spatial calibration in metres per pixel (> 0), used for SI statistics export.
    # Saved instantaneous vectors remain in px and px/s.
    "spatial_scale_m_per_px": 4.8e-05,

    # First N matching windows used to estimate the layer-selection speed gate U_th.
    # Positive integer; can exceed the requested number/duration of saved fields.
    # Skipped when --gate-file supplies an existing calibration.
    "calibration_windows": 5000,

    # Upper speed limit (> 0) of the U_th histogram, in px/ms (1 px/ms = 1000 px/s).
    # Cover both layer distributions; the gate search starts at 0.5 px/ms.
    # Speeds above this value are excluded from the histogram, not from saved vectors.
    "histogram_max_px_per_ms": 30.0,

    # Number of equal-width histogram bins from 0 to histogram_max_px_per_ms (integer >= 2).
    "histogram_bins": 300,

    # Moving-average width in histogram bins: integer 1..histogram_bins; 1 disables it.
    # Smooths calibration histograms only, not velocity fields.
    "histogram_smoothing_bins": 5,

    # Default number of instantaneous fields to save (positive integer).
    # Used only when neither --windows nor --duration is supplied.
    "max_windows": 100000,

    # Minimum valid windows per cell included in the statistics CSVs (positive integer).
    # This output filter leaves saved vectors and accumulated moments unchanged.
    "min_valid_samples": 10,
}

# Flow settings from the paper. Copy values into CONFIG as needed.

# Pipe flow, Re = 5,300:
# {
#     "roi": {"x": 0, "y": 0, "width": 1280, "height": 360},
#     "layers": [[22, 11, -1, 1, 2, 8], [22, 66, 0, 0, 1, 4]],
#     "stride_segments": 20,
#     "max_windows": 200000,
#     "min_valid_samples": 1,
#     "spatial_scale_m_per_px": 2.22e-05,
# }

# Pipe flow, Re = 19,000:
# {
#     "roi": {"x": 0, "y": 0, "width": 1280, "height": 360},
#     "layers": [[6, 3, -1, 1, 2, 8], [6, 18, 0, 0, 1, 4]],
#     "stride_segments": 6,
#     "max_windows": 200000,
#     "min_valid_samples": 1,
#     "spatial_scale_m_per_px": 2.22e-05,
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
#     "spatial_scale_m_per_px": 5.49e-05,
# }
