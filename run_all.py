"""Run the whole SOG1 analysis.

    python run_all.py            # everything (needs internet, HMMER via pyhmmer, MAFFT, trimAl, IQ-TREE)
    python run_all.py --offline  # steps 4-8 only, from the public data already in data/derived/

The GitHub workflow in .github/workflows/pipeline.yml runs the full version.
"""
import argparse
import importlib
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "pipeline"))

STEPS = [
    ("s01_fetch", "download proteomes, CDS and the NAC domain model", True),
    ("s02_nac_family", "NAC family and SOG1 clade in every genome", True),
    ("s03_annotate", "gene structure, Compara, promoters, AlphaFold", True),
    ("s04_phylogeny", "alignment and maximum-likelihood tree", False),
    ("s05_proteins", "SQ/TQ motifs, disorder, AlphaFold confidence", False),
    ("s06_selection", "dN/dS of NAC domain vs tail", False),
    ("s07_promoters", "SOG1 binding motif in target promoters", False),
    ("s08_export", "write docs/results.js", False),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--offline", action="store_true", help="skip the steps that need internet")
    ap.add_argument("--from", dest="start", type=int, default=1, help="first step to run (1-8)")
    args = ap.parse_args()
    t0 = time.time()
    for i, (mod, what, online) in enumerate(STEPS, 1):
        if i < args.start or (online and args.offline):
            continue
        print(f"\n=== step {i}: {what} ===", flush=True)
        importlib.import_module(mod).main()
    print(f"\nfinished in {(time.time() - t0) / 60:.1f} min", flush=True)


if __name__ == "__main__":
    main()
