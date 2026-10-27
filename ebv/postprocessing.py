"""Compute flow statistics from saved cell velocity fields."""
from pathlib import Path
import hashlib
import json
import numpy as np
from .config import load_config, write_json
from .execution import BLOCK_SIZE, implementation_hash, jobs
from .statistics import Moments, grid_shape, direct_grid


def _statistics_block(args):
    source, c, output, first, stop = args
    source, output = Path(source), Path(output)
    block = output/"blocks"/f"{first:06d}_{stop:06d}"
    block.mkdir(parents=True, exist_ok=True)
    if (block/"complete.json").exists() and all(
            (block/name).is_file() for name in ("moments.npz", "layers.npz")):
        return first, str(block)
    acc = Moments(grid_shape(c))
    totals = np.zeros(2, dtype=np.int64)
    fractions = np.zeros((stop-first, 2), dtype=np.float64)
    nonempty = np.zeros(stop-first, dtype=bool)
    for fi in range(first, stop):
        path = source/"vectors"/f"{fi:06d}.npy"
        vectors = np.load(path, allow_pickle=False)
        if (vectors.ndim != 2 or vectors.shape[1] != 5 or
                not np.isfinite(vectors).all() or
                not np.isin(vectors[:, 4], (0, 1)).all()):
            raise ValueError(f"Invalid [x, y, u, v, layer_id] array: {path}")
        acc.add(direct_grid(vectors, c))
        counts = np.bincount(vectors[:, 4].astype(int), minlength=2)
        totals += counts
        if len(vectors):
            fractions[fi-first] = counts/len(vectors)
            nonempty[fi-first] = True
    acc.save(block/"moments.npz")
    np.savez(block/"layers.npz", counts=totals, fractions=fractions, nonempty=nonempty)
    write_json(block/"complete.json", dict(first=first, stop=stop, windows=stop-first))
    return first, str(block)


def run(source, output=None, workers=1):
    source = Path(source).resolve()
    output = (source/"statistics" if output is None else Path(output)).resolve()
    if type(workers) is not int or workers <= 0:
        raise ValueError("workers must be a positive integer")
    complete = json.loads((source/"complete.json").read_text(encoding="utf-8"))
    if complete.get("status") != "COMPLETE":
        raise ValueError("Finish MAIN processing before computing statistics")
    c = load_config(source/"config.json")
    index_path = source/"frame_times.npy"
    index = np.load(index_path, allow_pickle=False)
    windows = complete["windows"]
    if (type(windows) is not int or windows <= 0 or index.shape != (windows, 2) or
            not np.issubdtype(index.dtype, np.integer) or np.any(index[:, 1] < index[:, 0])):
        raise ValueError("Invalid source window index")
    vector_identity = hashlib.sha256()
    for fi in range(windows):
        path = source/"vectors"/f"{fi:06d}.npy"
        if not path.is_file():
            raise FileNotFoundError(f"Missing instantaneous field: {path}. MAIN must save every window.")
        stat = path.stat()
        vector_identity.update(f"{fi}:{stat.st_size}:{stat.st_mtime_ns}\n".encode())
    identity = dict(input=str(source), stage="flow_statistics", windows=windows, config=c,
                    frame_times_sha256=hashlib.sha256(index_path.read_bytes()).hexdigest(),
                    vector_files_sha256=vector_identity.hexdigest(),
                    implementation_sha256=implementation_hash())
    sig = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    output.mkdir(parents=True, exist_ok=True)
    manifest = output/"run.json"
    if manifest.exists():
        if json.loads(manifest.read_text())["signature"] != sig:
            raise ValueError("Statistics output belongs to different input/version. Use a new output directory.")
    elif any(output.iterdir()):
        raise ValueError("Statistics output directory is not empty and has no run manifest")
    write_json(manifest, dict(signature=sig, **identity))
    write_json(output/"config.json", c)
    items = [(str(source), c, str(output), first, min(first+BLOCK_SIZE, windows))
             for first in range(0, windows, BLOCK_SIZE)]
    blocks = dict(jobs(_statistics_block, items, workers, "Statistics"))
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
    result = dict(status="COMPLETE", stage="flow_statistics", windows=windows, signature=sig,
                  duration_s=float((index[-1, 1]-index[0, 0])/1e6),
                  valid_cell_windows=int(acc.count.sum()),
                  layer_fraction_per_nonempty_window=mean_fraction,
                  pooled_layer_vector_counts=layer_counts.tolist(),
                  spatial_cell_px=c["cell_size_px"], grid_step_px=c["cell_size_px"],
                  temporal_batch_size=1, additional_spatial_average=False)
    write_json(output/"complete.json", result)
    return result
