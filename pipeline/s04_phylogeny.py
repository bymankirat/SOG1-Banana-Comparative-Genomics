"""Step 4: alignment, maximum-likelihood tree, and the final call on who is in the SOG1 clade.

MAFFT L-INS-i -> trimAl (automated1) -> IQ-TREE with ModelFinder (BIC) and 1000
ultrafast bootstraps. The tree holds every candidate from step 2, the three
Arabidopsis clade genes, and the Arabidopsis NACs that compete with them as best
match, and is rooted on NST1 (a distant subfamily).

Clade rule: walk up the tree from a candidate until the first group that contains
any Arabidopsis gene. If those Arabidopsis genes are only SOG1/NAC044/NAC085, the
candidate is in the SOG1 clade; it is then labelled by the closest of the three
(patristic distance). This works for moss and liverwort genes whose best
sequence match is a different NAC but which group with SOG1 in the tree.

If the tools are missing, the existing tree in results/tree/ is reused as long as
the input sequences are unchanged; the clade call is always recomputed.
"""
import re
import shutil
import subprocess

import pandas as pd

from util import DERIVED, RESULTS, config, log, md5, read_fasta, save_json, write_fasta, load_json

OUT = RESULTS / "tree"
OUT.mkdir(exist_ok=True)


def safe(name):
    return re.sub(r"[^A-Za-z0-9_.-]", "_", name.replace("|", "__"))


def run(cmd, **kw):
    log("$", " ".join(cmd))
    return subprocess.run(cmd, check=True, **kw)


# ------------------------------------------------------------------ tree utilities
def parse_newick(s):
    s = s.strip().rstrip(";")
    pos = 0

    def node(parent):
        nonlocal pos
        n = {"children": [], "parent": parent, "label": "", "len": 0.0}
        if s[pos] == "(":
            pos += 1
            while True:
                n["children"].append(node(n))
                if s[pos] == ",":
                    pos += 1
                    continue
                break
            pos += 1  # ')'
        start = pos
        while pos < len(s) and s[pos] not in ",():":
            pos += 1
        n["label"] = s[start:pos]
        if pos < len(s) and s[pos] == ":":
            pos += 1
            start = pos
            while pos < len(s) and s[pos] not in ",()":
                pos += 1
            n["len"] = float(s[start:pos])
        return n

    root = node(None)

    def depth(n, d):
        n["depth"] = d + n["len"]
        for c in n["children"]:
            depth(c, n["depth"])
    depth(root, -root["len"])
    return root


def tips(n):
    return [n] if not n["children"] else [t for c in n["children"] for t in tips(c)]


def classify(root, names):
    """names: tree label -> 'Sp|gene'. Returns {key: (in_clade, label, support)}."""
    clade_keys = {f"At|{g}": s for g, s in config.AT_CLADE.items()}
    leaves = {names.get(t["label"], t["label"]): t for t in tips(root)}
    at_leaves = {k for k in leaves if k.startswith("At|")}

    def ancestors(n):
        while n["parent"] is not None:
            n = n["parent"]
            yield n

    def lca_depth(a, b):
        anc = {id(x) for x in ancestors(a)}
        for x in ancestors(b):
            if id(x) in anc:
                return x["depth"]
        return 0.0

    out = {}
    for key, t in leaves.items():
        if key.startswith("At|"):
            out[key] = (key in clade_keys, clade_keys.get(key), None)
            continue
        verdict = (False, None, None)
        for a in ancestors(t):
            at_here = {names.get(x["label"], x["label"]) for x in tips(a)} & at_leaves
            if at_here:
                ok = at_here <= set(clade_keys)
                sup = float(a["label"]) if re.fullmatch(r"[\d.]+", a["label"] or "") else None
                if ok:
                    near = min(clade_keys, key=lambda c: t["depth"] + leaves[c]["depth"] - 2 * lca_depth(t, leaves[c]))
                    verdict = (True, clade_keys[near], sup)
                else:
                    verdict = (False, None, sup)
                break
        out[key] = verdict
    return out


def build(src):
    iqtree = shutil.which("iqtree2") or shutil.which("iqtree") or shutil.which("iqtree3")
    if not (shutil.which("mafft") and shutil.which("trimal") and iqtree):
        raise SystemExit("mafft, trimal and IQ-TREE are needed to rebuild the tree (installed in the GitHub workflow)")
    recs = read_fasta(src)
    names = {safe(h): h for h, _ in recs}
    write_fasta(OUT / "input.faa", [(safe(h), s) for h, s in recs])
    for f in OUT.glob("sog1_clade.*"):
        f.unlink()
    with open(OUT / "sog1_clade.aln.faa", "w") as fh:
        run(["mafft", "--localpair", "--maxiterate", "1000", "--thread", "-1", "--quiet", str(OUT / "input.faa")], stdout=fh)
    run(["trimal", "-in", str(OUT / "sog1_clade.aln.faa"), "-out", str(OUT / "sog1_clade.trim.faa"), "-automated1"])
    run([iqtree, "-s", str(OUT / "sog1_clade.trim.faa"), "-m", "MFP", "-B", "1000", "-T", "AUTO", "--threads-max", "4",
         "--seed", str(config.SEED), "-o", safe(f"At|{config.TREE_ROOT}"), "-pre", str(OUT / "sog1_clade"), "-quiet", "-redo"])
    report = (OUT / "sog1_clade.iqtree").read_text()
    model = re.search(r"Best-fit model according to BIC:\s*(\S+)", report)
    trimmed = read_fasta(OUT / "sog1_clade.trim.faa")
    full = read_fasta(OUT / "sog1_clade.aln.faa")
    (OUT / "input.faa").unlink(missing_ok=True)
    for f in OUT.glob("sog1_clade.*"):
        if f.suffix not in (".faa", ".treefile", ".iqtree", ".contree"):
            f.unlink()
    return {"model": model.group(1) if model else None, "sequences": len(recs),
            "aligned_columns": len(full[0][1]), "trimmed_columns": len(trimmed[0][1]),
            "root": f"At|{config.TREE_ROOT}", "names": names, "iqtree": iqtree.rsplit("/", 1)[-1]}


def main():
    src = DERIVED / "sog1_candidates.faa"
    stamp = OUT / "input.md5"
    if stamp.exists() and stamp.read_text().strip() == md5(src) and (OUT / "sog1_clade.treefile").exists():
        log("tree is up to date")
        summ = load_json(OUT / "summary.json")
    else:
        summ = build(src)
        stamp.write_text(md5(src))
        log("best model:", summ["model"])

    root = parse_newick((OUT / "sog1_clade.treefile").read_text())
    calls = classify(root, summ["names"])
    nac = pd.read_csv(DERIVED / "nac_genes.csv")
    nac["key"] = nac.short + "|" + nac.gene
    nac["sog1_clade"] = nac.key.map(lambda k: calls.get(k, (False,))[0]).fillna(False).astype(bool)
    nac["clade_label"] = nac.key.map(lambda k: calls.get(k, (False, None))[1])
    nac["clade_support"] = nac.key.map(lambda k: calls.get(k, (False, None, None))[2])
    nac.drop(columns=["key"]).to_csv(DERIVED / "nac_genes.csv", index=False)

    seqs = dict(read_fasta(src))
    members = nac[nac.sog1_clade]
    write_fasta(DERIVED / "sog1_clade_proteins.faa", [(f"{r.short}|{r.gene}", seqs[f"{r.short}|{r.gene}"]) for r in members.itertuples()])
    agree = nac[nac.candidate | nac.sog1_clade]
    summ["rule_agreement"] = {"both": int((agree.best_hit_rule & agree.sog1_clade).sum()),
                              "tree_only": sorted((agree[~agree.best_hit_rule & agree.sog1_clade].short + "|" + agree.gene).tolist()),
                              "best_hit_only": sorted((agree[agree.best_hit_rule & ~agree.sog1_clade].short + "|" + agree.gene).tolist())}
    summ["other_at"] = {f"At|{r.gene}": (r.symbol if isinstance(r.symbol, str) and r.symbol else r.gene)
                        for r in nac[(nac.species == config.REFERENCE_SPECIES) & nac.in_tree & ~nac.sog1_clade].itertuples()}
    save_json(OUT / "summary.json", summ)

    counts = pd.read_csv(RESULTS / "nac_counts.csv")
    per = members.groupby("species").size()
    counts["sog1_clade"] = counts.species.map(per).fillna(0).astype(int)
    counts.to_csv(RESULTS / "nac_counts.csv", index=False)
    log("SOG1 clade per species:", {k: int(v) for k, v in per.items()})
    log("tree vs best-hit rule:", summ["rule_agreement"])


if __name__ == "__main__":
    main()
