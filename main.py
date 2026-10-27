"""Command-line entry point; see README.md for usage."""
import os
for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_name, "1")

import argparse
import json
from pathlib import Path
from ebv.config import load_config
from ebv.pipeline import run


def main():
    parser = argparse.ArgumentParser(description="Readout-adaptive triplet velocimetry")
    parser.add_argument("--input", required=True, help="Metavision .raw/.hdf5 or structured-event .npz")
    parser.add_argument("--config", default=str(Path(__file__).with_name("configuration.py")),
                        help="Configuration file (default: configuration.py beside main.py)")
    parser.add_argument("--flow", help="Override ACTIVE_FLOW, e.g. pipe5300 or jet3500")
    parser.add_argument("--output", required=True, help="New output folder, or matching run to resume")
    parser.add_argument("--windows", type=int, help="Override the preset analysis window count")
    parser.add_argument("--workers", type=int, default=1, help="Parallel block workers (default: 1)")
    parser.add_argument("--gate-file", help="Reuse a configuration-matched calibration JSON")
    parser.add_argument("--save-vectors", action="store_true", help="Write [x,y,u,v,layer_id] per window")
    args = parser.parse_args()
    result = run(args.input, load_config(args.config, args.flow), args.output,
                 args.windows, args.workers, args.gate_file, args.save_vectors)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
