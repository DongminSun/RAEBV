"""Calibrate the speed gate and save instantaneous cell velocity fields."""
from pathlib import Path
import hashlib
import json
import time
import numpy as np
from .config import write_json
from .events import build_index, windows_range
from .estimator import estimate_layers, select_vectors, histogram, gate_from_histogram
from .execution import BLOCK_SIZE, implementation_hash, jobs as _jobs


def signature(source, c, windows, duration_us, gate_file):
    path = Path(source).resolve()
    stat = path.stat()
    data = dict(input=str(path), size=stat.st_size, mtime_ns=stat.st_mtime_ns,
                config=c, windows=windows, duration_us=duration_us, save_vectors=True,
                stage="instantaneous_vectors",
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
    source, c, index_path, folder, first, stop, gate = args
    folder = Path(folder)
    block = folder/"blocks"/f"{first:06d}_{stop:06d}"
    block.mkdir(parents=True, exist_ok=True)
    if (block/"complete.json").exists() and all(
            (folder/"vectors"/f"{fi:06d}.npy").is_file() for fi in range(first, stop)):
        return first, json.loads((block/"complete.json").read_text())["valid_cell_windows"]
    index = np.load(index_path, mmap_mode="r")
    valid_cell_windows = 0
    t0 = time.perf_counter()
    for fi, snap, _ in windows_range(source, c, index, first, stop):
        layers = estimate_layers(snap, c)
        vectors = select_vectors(layers, c, gate)
        valid_cell_windows += len(vectors)
        target = folder/"vectors"/f"{fi:06d}.npy"
        temporary = target.with_suffix(".tmp")
        with temporary.open("wb") as fh:
            np.save(fh, vectors)
        temporary.replace(target)
        if (fi-first+1) % 1000 == 0:
            write_json(block/"progress.json", dict(first=first, processed=fi-first+1,
                        target=stop-first, elapsed_seconds=time.perf_counter()-t0))
    write_json(block/"complete.json", dict(first=first, stop=stop, windows=stop-first,
                valid_cell_windows=valid_cell_windows,
                elapsed_seconds=time.perf_counter()-t0))
    return first, valid_cell_windows


def run(source, c, output, windows=None, workers=1, gate_file=None, duration=None):
    source = str(Path(source).resolve())
    output = Path(output).resolve()
    if windows is not None and duration is not None:
        raise ValueError("Choose --windows or --duration, not both")
    if type(workers) is not int or workers <= 0:
        raise ValueError("workers must be a positive integer")
    duration_us = None
    if duration is not None:
        if not np.isfinite(duration) or duration <= 0:
            raise ValueError("duration must be a positive finite number of seconds")
        duration_us = int(round(duration*1e6))
        if duration_us < 1:
            raise ValueError("duration must be at least one microsecond")
    else:
        windows = c["max_windows"] if windows is None else windows
        if type(windows) is not int or windows <= 0:
            raise ValueError("windows must be a positive integer")
    output.mkdir(parents=True, exist_ok=True)
    sig, identity = signature(source, c, windows, duration_us, gate_file)
    manifest = output/"run.json"
    if manifest.exists():
        if json.loads(manifest.read_text())["signature"] != sig:
            raise ValueError("Output belongs to a different input/configuration/version. Use a new output directory.")
    elif any(output.iterdir()):
        raise ValueError("Output directory is not empty and has no run manifest")
    write_json(manifest, dict(signature=sig, **identity))
    write_json(output/"config.json", c)
    (output/"vectors").mkdir(exist_ok=True)
    index_path = output/"frame_times.npy"
    if not index_path.exists():
        print("Indexing readout windows from input events", flush=True)
        index = build_index(source, c, windows,
                            lambda n: print(f"Index: {n} windows", flush=True), stop_us=duration_us)
        temp = index_path.with_suffix(".tmp")
        with temp.open("wb") as fh:
            np.save(fh, index)
        temp.replace(index_path)
    index = np.load(index_path, mmap_mode="r")
    if (index.ndim != 2 or index.shape[1] != 2 or not len(index) or
            (windows is not None and len(index) != windows) or
            (duration_us is not None and np.any(index[:, 1] >= duration_us))):
        raise ValueError("Invalid cached window index")
    windows = len(index)
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
        calibration_index = index_path
        if windows < c["calibration_windows"]:
            calibration_index = calibration/"frame_times.npy"
            if not calibration_index.exists():
                print(f"Indexing {c['calibration_windows']} windows for speed-gate calibration", flush=True)
                try:
                    calibration_times = build_index(source, c, c["calibration_windows"])
                except RuntimeError as exc:
                    raise RuntimeError(
                        f"Speed-gate calibration needs {c['calibration_windows']} complete windows. "
                        "Use a longer recording or a compatible --gate-file from the same recording."
                    ) from exc
                temporary = calibration_index.with_suffix(".tmp")
                with temporary.open("wb") as fh:
                    np.save(fh, calibration_times)
                temporary.replace(calibration_index)
            if np.load(calibration_index, mmap_mode="r").shape != (c["calibration_windows"], 2):
                raise ValueError("Invalid cached calibration index")
        jobs = [(source, c, str(calibration_index), str(calibration), first,
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
             min(first+BLOCK_SIZE, windows), gate)
            for first in range(0, windows, BLOCK_SIZE)]
    blocks = dict(_jobs(_velocity_block, jobs, workers, "Velocity"))
    result = dict(status="COMPLETE", stage="instantaneous_vectors", windows=windows, signature=sig,
                  duration_s=float((index[-1, 1]-index[0, 0])/1e6),
                  requested_duration_s=None if duration_us is None else duration_us/1e6,
                  speed_gate_pxs=gate, valid_cell_windows=sum(blocks.values()),
                  spatial_cell_px=c["cell_size_px"], grid_step_px=c["cell_size_px"],
                  temporal_batch_size=1, additional_spatial_average=False)
    write_json(output/"complete.json", result)
    return result
