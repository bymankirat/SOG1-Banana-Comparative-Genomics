"""Run phase 3: SOG1 under stress across crops, and CRISPR guides for banana SOG1.

    python run_phase3.py           # all steps (needs internet)
    python run_phase3.py --from 3  # start at step 3 from saved downloads

The GitHub workflow in .github/workflows/phase3.yml runs it and saves the results here.
"""
import argparse
import importlib
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "pipeline"))

STEPS = [
    ("p3_01_atlas", "find and download stress experiments (Expression Atlas)"),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--from", dest="start", type=int, default=1)
    args = ap.parse_args()
    t0 = time.time()
    for i, (mod, what) in enumerate(STEPS, 1):
        if i < args.start:
            continue
        print(f"\n=== phase 3 step {i}: {what} ===", flush=True)
        importlib.import_module(mod).main()
    print(f"\nfinished in {(time.time() - t0) / 60:.1f} min", flush=True)


if __name__ == "__main__":
    main()
