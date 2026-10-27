"""Compute statistics from saved instantaneous fields; no event recording is needed."""
import os
for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_name, "1")

import argparse
import json
from ebv.postprocessing import run


def main():
    parser = argparse.ArgumentParser(description="Flow statistics from saved RAEBV cell vectors")
    parser.add_argument("--input", required=True, help="Completed main.py output folder containing vectors")
    parser.add_argument("--output", help="Statistics folder (default: INPUT/statistics)")
    parser.add_argument("--workers", type=int, default=1, help="Parallel block workers (default: 1)")
    args = parser.parse_args()
    print(json.dumps(run(args.input, args.output, args.workers), indent=2))


if __name__ == "__main__":
    main()
