"""Estimate layer velocities and select cell vectors."""
import numpy as np
from . import core
from .config import POLARITY_MODES


def estimate_layers(snap, c):
    w, h = c["roi"]["width"], c["roi"]["height"]
    polarity_mode = POLARITY_MODES[c["polarity"]]
    layers = []
    for s, g, xmin, xmax, ymin, ymax in c["layers"]:
        past, middle, current = core.build_triplet(snap, s, g)
        mx, my, mt, mp, mg = core.build_spatial_lookup(
            middle["x"], middle["y"], middle["t"], middle["p"], w, h)
        px, py, pt, pp, pg = core.build_spatial_lookup(
            past["x"], past["y"], past["t"], past["p"], w, h)
        tmin = np.min(past["t"]) if len(past["t"]) else 0
        tmax = np.max(past["t"]) if len(past["t"]) else 0
        layers.append(core.process_triplet_v8(
            current["x"], current["y"], current["t"], current["p"],
            mx, my, mt, mp, mg, px, py, pt, pp, pg, w, h,
            xmin, xmax, ymin, ymax, polarity_mode, 0.0, 7000e8,
            tmin, tmax, c["trajectory_tolerance_px"]))
    return layers


def select_vectors(layers, c, gate_pxs):
    v = core.fuse_layers_speed_gated_bin(
        layers, c["roi"]["width"], c["roi"]["height"], c["cell_size_px"], gate_pxs)
    return np.empty((0, 5), dtype=np.float64) if v is None else v


def histogram(layers, c):
    edges = np.linspace(0.0, c["histogram_max_px_per_ms"], c["histogram_bins"]+1)
    h = np.zeros((2, c["histogram_bins"]), dtype=np.int64)
    for j, layer in enumerate(layers):
        b = core._bin_one_layer_for_gated_fusion(
            layer, c["roi"]["width"], c["roi"]["height"], c["cell_size_px"])
        h[j] = np.histogram(np.hypot(b[:, 2], b[:, 3])/1000.0, bins=edges)[0]
    return h


def gate_from_histogram(hist, c):
    totals = hist.sum(axis=1)
    if np.any(totals == 0):
        raise ValueError("Both layers require valid calibration cell speeds")
    probability = hist/totals[:, None]
    smoothed = np.stack([core._smooth_hist_prob(p, c["histogram_smoothing_bins"]) for p in probability])
    edges = np.linspace(0.0, c["histogram_max_px_per_ms"], c["histogram_bins"]+1)
    gate, method = core._find_hist_intersection((edges[:-1]+edges[1:])/2, *smoothed)
    return dict(speed_gate_pxs=float(gate*1000.0), speed_gate_pxms=float(gate), method=method), edges, probability, smoothed
