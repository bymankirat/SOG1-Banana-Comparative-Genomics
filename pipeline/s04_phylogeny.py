"""Step 4: alignment and maximum-likelihood tree of the SOG1 clade.

MAFFT L-INS-i alignment -> trimAl (automated1) -> IQ-TREE 2 with ModelFinder
(best model by BIC) and 1000 ultrafast bootstrap replicates, rooted on three
Arabidopsis NACs from other subfamilies. The 2023 report used a neighbour-joining
tree of five proteins; this is the model-based version with the whole clade.

If the tools are missing (e.g. a laptop without IQ-TREE), the step reuses the
results already in results/tree/ as long as the input sequences have not changed.
"""
import re
import shutil
import subprocess

from util import DERIVED, RESULTS, config, log, md5, read_fasta, save_json, write_fasta

OUT = RESULTS / "tree"
OUT.mkdir(exist_ok=True)


def safe(name):
    return re.sub(r"[^A-Za-z0-9_.-]", "_", name.replace("|", "__"))


def run(cmd, **kw):
    log("$", " ".join(cmd))
    return subprocess.run(cmd, check=True, **kw)


def main():
    src = DERIVED / "sog1_clade_proteins.faa"
    stamp = OUT / "input.md5"
    if stamp.exists() and stamp.read_text().strip() == md5(src) and (OUT / "sog1_clade.treefile").exists():
        log("tree is up to date")
        return
    iqtree = shutil.which("iqtree2") or shutil.which("iqtree") or shutil.which("iqtree3")
    if not (shutil.which("mafft") and shutil.which("trimal") and iqtree):
        raise SystemExit("mafft, trimal and iqtree2 are needed for this step (they are installed in the GitHub workflow)")

    recs = read_fasta(src)
    names = {safe(h): h for h, _ in recs}
    write_fasta(OUT / "input.faa", [(safe(h), s) for h, s in recs])
    with open(OUT / "sog1_clade.aln.faa", "w") as fh:
        run(["mafft", "--localpair", "--maxiterate", "1000", "--thread", "-1", "--quiet", str(OUT / "input.faa")], stdout=fh)
    run(["trimal", "-in", str(OUT / "sog1_clade.aln.faa"), "-out", str(OUT / "sog1_clade.trim.faa"), "-automated1"])
    outgroup = [safe(f"At|{g}") for g in config.AT_OUTGROUP if safe(f"At|{g}") in names]
    for f in OUT.glob("sog1_clade.*"):
        if f.suffix not in (".faa",):
            f.unlink()
    run([iqtree, "-s", str(OUT / "sog1_clade.trim.faa"), "-m", "MFP", "-B", "1000", "-T", "AUTO", "--threads-max", "4",
         "--seed", str(config.SEED), "-o", ",".join(outgroup), "-pre", str(OUT / "sog1_clade"), "-quiet", "-redo"])

    report = (OUT / "sog1_clade.iqtree").read_text()
    model = re.search(r"Best-fit model according to BIC:\s*(\S+)", report)
    trimmed = read_fasta(OUT / "sog1_clade.trim.faa")
    full = read_fasta(OUT / "sog1_clade.aln.faa")
    save_json(OUT / "summary.json", {
        "model": model.group(1) if model else None, "sequences": len(recs),
        "aligned_columns": len(full[0][1]), "trimmed_columns": len(trimmed[0][1]),
        "outgroup": [names[o] for o in outgroup], "names": names, "iqtree": iqtree.rsplit("/", 1)[-1]})
    for f in ("input.faa",):
        (OUT / f).unlink(missing_ok=True)
    stamp.write_text(md5(src))
    log("best model:", model.group(1) if model else "?")


if __name__ == "__main__":
    main()
