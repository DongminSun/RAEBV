"""Shared deterministic block execution and implementation identity."""
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import hashlib

BLOCK_SIZE = 5000  # Preserve the final reference's statistics reduction order.


def implementation_hash():
    h = hashlib.sha256()
    for path in sorted(Path(__file__).parent.glob("*.py")):
        h.update(path.name.encode())
        h.update(path.read_bytes())
    return h.hexdigest()


def jobs(fn, items, workers, label):
    if workers == 1:
        for i, item in enumerate(items):
            yield fn(item)
            print(f"{label}: {i+1}/{len(items)} blocks", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(fn, item) for item in items]
            for i, future in enumerate(as_completed(futures)):
                yield future.result()
                print(f"{label}: {i+1}/{len(items)} blocks", flush=True)
