"""Reproducible two-pass calibration and block-wise full-record processing."""
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import hashlib
import json
import time
import numpy as np
from .config import write_json
from .events import build_index, windows_range
from .estimator import estimate_layers, select_vectors, histogram, gate_from_histogram
from .statistics import Moments, grid_shape, direct_grid

BLOCK_SIZE = 5000  # Fixed to preserve the final reference's reduction order.


def implementation_hash():
    h = hashlib.sha256()
    for path in sorted(Path(__file__).parent.glob("*.py")):
        h.update(path.name.encode())
        h.update(path.read_bytes())
    return h.hexdigest()


def signature(source, c, windows, gate_file, save_vectors):
    path = Path(source).resolve()
    stat = path.stat()
    data = dict(input=str(path), size=stat.st_size, mtime_ns=stat.st_mtime_ns,
                config=c, windows=windows, save_vectors=save_vectors,
                implementation_sha256=implementation_hash())
    if gate_file:
        data["gate_file_sha256"] = hashlib.sha256(Path(gate_file).read_bytes()).hexdigest()
    digest = hashlib.sha256(json.dumps(data, sort_keys=True).encode()).hexdigest()
    return digest, data


def _calibration_block(args):
    source, c, index_path, folder, first, stop = args
    folder = Path(folder)
    target = folder/f"hist_{first:06d}_{stop:06d}.npy"
    if target.exists():
        return str(target)
    index = np.load(index_path, mmap_mode="r")
    hist = np.zeros((2, c["histogram_bins"]), dtype=np.int64)
    for _, snap, _ in windows_range(source, c, index, first, stop):
        hist += histogram(estimate_layers(snap, c), c)
    temporary = target.with_suffix(".tmp")
    with temporary.open("wb") as fh:
        np.save(fh, hist)
    temporary.replace(target)
    return str(target)


def _velocity_block(args):
    source, c, index_path, folder, first, stop, gate, save_vectors = args
    folder = Path(folder)
    block = folder/"blocks"/f"{first:06d}_{stop:06d}"
    block.mkdir(parents=True, exist_ok=True)
    if (block/"complete.json").exists():
        return first, str(block)
    index = np.load(index_path, mmap_mode="r")
    acc = Moments(grid_shape(c))
    totals = np.zeros(2, dtype=np.int64)
    fractions = np.zeros((stop-first, 2), dtype=np.float64)
    nonempty = np.zeros(stop-first, dtype=bool)
    t0 = time.perf_counter()
    for fi, snap, _ in windows_range(source, c, index, first, stop):
        layers = estimate_layers(snap, c)
        vectors = select_vectors(layers, c, gate)
        acc.add(direct_grid(vectors, c))
        counts = np.bincount(vectors[:, 4].astype(int), minlength=2)
        totals += counts
        if len(vectors):
            fractions[fi-first] = counts/len(vectors)
            nonempty[fi-first] = True
        if save_vectors:
            np.save(folder/"vectors"/f"{fi:06d}.npy", vectors)
        if (fi-first+1) % 1000 == 0:
            write_json(block/"progress.json", dict(first=first, processed=fi-first+1,
                        target=stop-first, elapsed_seconds=time.perf_counter()-t0))
    acc.save(block/"moments.npz")
    np.savez(block/"layers.npz", counts=totals, fractions=fractions, nonempty=nonempty)
    write_json(block/"complete.json", dict(first=first, stop=stop, windows=stop-first,
                elapsed_seconds=time.perf_counter()-t0))
    return first, str(block)


def _jobs(fn, jobs, workers, label):
    if workers == 1:
        for i, job in enumerate(jobs):
            yield fn(job)
            print(f"{label}: {i+1}/{len(jobs)} blocks", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(fn, job) for job in jobs]
            for i, future in enumerate(as_completed(futures)):
                yield future.result()
                print(f"{label}: {i+1}/{len(jobs)} blocks", flush=True)


def run(source, c, output, windows=None, workers=1, gate_file=None, save_vectors=False):
    source = str(Path(source).resolve())
    output = Path(output).resolve()
    windows = c["max_windows"] if windows is None else windows
    if windows <= 0 or workers <= 0:
        raise ValueError("windows and workers must be positive")
    if not gate_file and windows < c["calibration_windows"]:
        raise ValueError("Fewer windows than calibration_windows; supply a compatible --gate-file for short examples")
    output.mkdir(parents=True, exist_ok=True)
    sig, identity = signature(source, c, windows, gate_file, save_vectors)
    manifest = output/"run.json"
    if manifest.exists():
        if json.loads(manifest.read_text())["signature"] != sig:
            raise ValueError("Output belongs to a different input/configuration/version. Use a new output directory.")
    elif any(output.iterdir()):
        raise ValueError("Output directory is not empty and has no run manifest")
    write_json(manifest, dict(signature=sig, **identity))
    write_json(output/"config.json", c)
    if save_vectors:
        (output/"vectors").mkdir(exist_ok=True)
    index_path = output/"frame_times.npy"
    if not index_path.exists():
        print("Indexing readout windows from input events", flush=True)
        index = build_index(source, c, windows, lambda n: print(f"Index: {n}/{windows}", flush=True))
        temp = index_path.with_suffix(".tmp")
        with temp.open("wb") as fh:
            np.save(fh, index)
        temp.replace(index_path)
    index = np.load(index_path, mmap_mode="r")
    if index.shape != (windows, 2):
        raise ValueError("Invalid cached window index")
    calibration = output/"calibration"
    calibration.mkdir(exist_ok=True)
    if gate_file:
        loaded = json.loads(Path(gate_file).read_text(encoding="utf-8"))
        if loaded.get("config") != c:
            raise ValueError("Saved calibration configuration does not match this run")
        gate = float(loaded["speed_gate_pxs"])
        if not np.isfinite(gate) or gate < 0:
            raise ValueError("Invalid saved speed gate")
        gate_info = {**loaded, "source": "provided_calibration"}
    else:
        jobs = [(source, c, str(index_path), str(calibration), first,
                 min(first+1000, c["calibration_windows"]))
                for first in range(0, c["calibration_windows"], 1000)]
        hist = sum((np.load(path) for path in _jobs(_calibration_block, jobs, workers, "Calibration")),
                   np.zeros((2, c["histogram_bins"]), dtype=np.int64))
        gate_info, edges, prob, smooth = gate_from_histogram(hist, c)
        gate_info.update(config=c, windows=c["calibration_windows"], source="first_windows_of_input")
        gate = gate_info["speed_gate_pxs"]
        np.savez(calibration/"histogram.npz", histogram=hist, edges=edges,
                 probability=prob, smoothed_probability=smooth)
    write_json(calibration/"gate.json", gate_info)
    print(f"Speed gate: {gate:.12g} px/s ({gate_info['method']})", flush=True)
    jobs = [(source, c, str(index_path), str(output), first,
             min(first+BLOCK_SIZE, windows), gate, save_vectors)
            for first in range(0, windows, BLOCK_SIZE)]
    blocks = dict(_jobs(_velocity_block, jobs, workers, "Velocity"))
    acc = Moments(grid_shape(c))
    layer_counts = np.zeros(2, dtype=np.int64)
    fractions, masks = [], []
    for first in sorted(blocks):
        folder = Path(blocks[first])
        part = Moments(grid_shape(c))
        part.restore(folder/"moments.npz")
        acc.merge(part)
        with np.load(folder/"layers.npz") as data:
            layer_counts += data["counts"]
            fractions.append(data["fractions"])
            masks.append(data["nonempty"])
    acc.save(output/"moments.npz")
    acc.table(c).to_csv(output/"statistics_px.csv", index=False)
    acc.table(c, si=True).to_csv(output/"statistics_si.csv", index=False)
    fractions, masks = np.concatenate(fractions), np.concatenate(masks)
    np.savez(output/"layer_fractions.npz", fractions=fractions, nonempty=masks,
             pooled_vector_counts=layer_counts)
    mean_fraction = fractions[masks].mean(axis=0).tolist() if masks.any() else [None, None]
    result = dict(status="COMPLETE", windows=windows, signature=sig,
                  duration_s=float((index[-1, 1]-index[0, 0])/1e6),
                  speed_gate_pxs=gate, valid_cell_windows=int(acc.count.sum()),
                  layer_fraction_per_nonempty_window=mean_fraction,
                  pooled_layer_vector_counts=layer_counts.tolist(),
                  spatial_cell_px=c["cell_size_px"], grid_step_px=c["cell_size_px"],
                  temporal_batch_size=1, additional_spatial_average=False)
    write_json(output/"complete.json", result)
    return result
