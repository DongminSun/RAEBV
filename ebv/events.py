"""Read RAW/HDF5/NPZ events and form readout segments without fixed-time matching."""
from collections import deque
from pathlib import Path
import numpy as np

CHUNK_US = 10_000


def event_chunks(path, roi, start_us=0, stop_us=None):
    """Yield ROI-local structured events without resetting chunk timestamps."""
    path = Path(path)
    if path.suffix.lower() == ".raw" or path.suffix == ".hdf5":
        try:
            from metavision_core.event_io import EventsIterator
        except ImportError as exc:
            raise RuntimeError(
                "RAW/HDF5 input requires the Prophesee Metavision SDK in this Python environment. "
                "Structured-event NPZ input works without the SDK."
            ) from exc
        args = dict(start_ts=int(start_us), delta_t=CHUNK_US, relative_timestamps=False)
        if stop_us is not None:
            args["max_duration"] = int(stop_us - start_us)
        source = EventsIterator(str(path), **args)
    elif path.suffix.lower() == ".npz":
        with np.load(path, allow_pickle=False) as data:
            events = data["events"]
        if events.dtype.names is None or not {"x", "y", "t", "p"} <= set(events.dtype.names):
            raise ValueError("NPZ must contain structured 'events' with x, y, t, p fields")
        if np.any(events["t"][1:] < events["t"][:-1]):
            raise ValueError("Events must be ordered by timestamp")
        lo = np.searchsorted(events["t"], start_us)
        hi = len(events) if stop_us is None else np.searchsorted(events["t"], stop_us)

        def batches():
            pos = lo
            while pos < hi:
                boundary = ((int(events["t"][pos]) // CHUNK_US) + 1) * CHUNK_US
                end = min(hi, int(np.searchsorted(events["t"], boundary)))
                yield events[pos:end]
                pos = end
        source = batches()
    else:
        raise ValueError("Supported inputs are Metavision .raw/.hdf5 and structured-event .npz files")
    x0, y0, w, h = (roi[k] for k in ("x", "y", "width", "height"))
    for ev in source:
        keep = ((ev["x"] >= x0) & (ev["x"] < x0+w) &
                (ev["y"] >= y0) & (ev["y"] < y0+h))
        cropped = ev[keep].copy()
        if x0:
            cropped["x"] -= x0
        if y0:
            cropped["y"] -= y0
        yield cropped


def segments(chunks, jump=200):
    """Preserve the final-run boundary convention, including a shared event."""
    previous_y = None
    pending = []
    started = False
    for ev in chunks:
        if not len(ev):
            continue
        ys = ev["y"].astype(np.int32)
        jumps = np.flatnonzero(-np.diff(ys) > jump) + 1
        if previous_y is not None and previous_y - int(ys[0]) > jump:
            jumps = np.r_[0, jumps]
        pos = 0
        for j in jumps:
            if started:
                seg = np.concatenate(pending + [ev[pos:j + 1]])
                yield {k: seg[k].astype(np.int64 if k == "t" else np.int32)
                       for k in ("x", "y", "t", "p")}
            started = True
            pending = []
            pos = j
        pending.append(ev[pos:].copy())
        previous_y = int(ys[-1])


def window_length(c):
    return max(3*s + 2*g for s, g, *_ in c["layers"])


def build_index(path, c, count=None, progress=None, stop_us=None):
    """Index a window count or all complete windows ending before stop_us."""
    if (count is None) == (stop_us is None):
        raise ValueError("Specify either a window count or a stop timestamp")
    if count is not None and (type(count) is not int or count <= 0):
        raise ValueError("Window count must be a positive integer")
    if stop_us is not None and (type(stop_us) is not int or stop_us <= 0):
        raise ValueError("Stop timestamp must be a positive integer in microseconds")
    buf, times = deque(), []
    need = window_length(c)
    for seg in segments(event_chunks(path, c["roi"], stop_us=stop_us), c["row_jump_threshold_px"]):
        buf.append(seg)
        if len(buf) < need:
            continue
        times.append((int(buf[0]["t"][0]), int(buf[-1]["t"][-1])))
        if progress and len(times) % 10000 == 0:
            progress(len(times))
        if count is not None and len(times) == count:
            return np.asarray(times, dtype=np.int64)
        for _ in range(c["stride_segments"]):
            buf.popleft()
    if stop_us is not None:
        if not times:
            raise RuntimeError("No complete matching windows within the requested duration")
        return np.asarray(times, dtype=np.int64)
    raise RuntimeError(f"Input contains only {len(times)} complete windows; requested {count}")


def windows_range(path, c, index, first, stop):
    """Replay exactly indexed windows after a bounded seek."""
    if not 0 <= first < stop <= len(index):
        raise ValueError((first, stop, len(index)))
    seek_us = max(0, int(index[first, 0])-CHUNK_US) // CHUNK_US * CHUNK_US
    end_us = int(index[stop-1, 1]) + 2*CHUNK_US
    need = window_length(c)
    buf, aligned, fi = deque(), False, first
    stream = event_chunks(path, c["roi"], seek_us, end_us)
    for seg in segments(stream, c["row_jump_threshold_px"]):
        buf.append(seg)
        if len(buf) < need:
            continue
        snap = list(buf)
        times = (int(snap[0]["t"][0]), int(snap[-1]["t"][-1]))
        expected = tuple(index[fi])
        if not aligned and times != expected:
            if times[0] > expected[0] or times[1] > expected[1]:
                raise RuntimeError(f"Cannot align window {fi}: {times} != {expected}")
            buf.popleft()
            continue
        if times != expected:
            raise RuntimeError(f"Window timestamp mismatch at {fi}: {times} != {expected}")
        aligned = True
        yield fi, snap, times
        fi += 1
        if fi == stop:
            return
        for _ in range(c["stride_segments"]):
            buf.popleft()
    raise RuntimeError(f"Input exhausted at window {fi}, requested stop={stop}")
