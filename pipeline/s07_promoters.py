"""Step 7: do SOG1 target genes carry the SOG1 binding motif CTT(N)7AAG, and do their banana orthologs?

Positive control first: in Arabidopsis, promoters of reported SOG1 targets should
hold more motifs than 1,000 random gene promoters. If the control works, the same
test on the banana orthologs (found by reciprocal best hit) asks whether the
regulatory wiring is conserved, which is what would make MaSOG1 a true
functional equivalent.

A composition control (CTT(N)7CTT: same letters, not a palindrome) checks that
any enrichment is about the motif and not just AT-rich promoters.
"""
import re

import numpy as np
import pandas as pd
from scipy import stats

from util import DERIVED, RESULTS, config, log, read_fasta, save_json

MOTIFS = {"SOG1 motif CTT(N)7AAG": config.SOG1_MOTIF, "control CTT(N)7CTT": "CTT.......CTT"}
PROXIMAL = 500
BIN = 100


def hits(seq, pattern):
    """Distances (bp upstream of the gene start) of every motif match, overlapping allowed, both strands."""
    rx = re.compile(f"(?=({pattern}))")
    L = len(seq)
    pos = {m.start() for m in rx.finditer(seq)}
    rc = seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    w = len(re.sub(r"\.", "N", pattern))
    pos |= {L - m.start() - w for m in rx.finditer(rc)}  # palindromes give the same site twice; the set removes it
    return sorted(L - p for p in pos)


def enrichment(target_counts, bg_counts, n_perm, rng):
    t, b = np.asarray(target_counts), np.asarray(bg_counts)
    if len(t) == 0:
        return {}
    draws = rng.choice(b, size=(n_perm, len(t)), replace=True).mean(axis=1)
    return {"n_target": int(len(t)), "n_background": int(len(b)),
            "mean_target": round(float(t.mean()), 3), "mean_background": round(float(b.mean()), 3),
            "fold": round(float(t.mean() / b.mean()), 2) if b.mean() > 0 else None,
            "p_perm": float((1 + (draws >= t.mean()).sum()) / (n_perm + 1))}


def main():
    rng = np.random.default_rng(config.SEED)
    rows = []
    for head, seq in read_fasta(DERIVED / "promoters.fna"):
        sp, gene, group, region = head.split("|")
        rec = {"species": sp, "gene": gene, "group": group, "length": len(seq),
               "at_fraction": round(sum(seq.count(x) for x in "AT") / max(len(seq), 1), 3)}
        for label, pat in MOTIFS.items():
            d = hits(seq, pat)
            rec[label] = len(d)
            rec[label + " proximal"] = sum(x <= PROXIMAL for x in d)
            rec[label + " positions"] = d
        rows.append(rec)
    df = pd.DataFrame(rows)
    tor = pd.read_csv(DERIVED / "target_orthologs.csv")
    sym_at = dict(zip(tor.at_gene, tor.symbol))
    sym_focus = tor.dropna(subset=["focus_gene"]).groupby("focus_gene").symbol.agg(lambda s: "/".join(sorted(set(s))))

    out = {"motifs": list(MOTIFS), "proximal_bp": PROXIMAL, "promoter_bp": config.PROMOTER_LEN, "species": {}}
    for sp, sub in df.groupby("species"):
        t, b = sub[sub.group == "target"], sub[sub.group == "background"]
        res = {"at_target": round(float(t.at_fraction.mean()), 3) if len(t) else None,
               "at_background": round(float(b.at_fraction.mean()), 3)}
        for label in MOTIFS:
            full = enrichment(t[label], b[label], config.N_PERMUTATIONS, rng)
            prox = enrichment(t[label + " proximal"], b[label + " proximal"], config.N_PERMUTATIONS, rng)
            has_t, has_b = int((t[label + " proximal"] > 0).sum()), int((b[label + " proximal"] > 0).sum())
            fisher = stats.fisher_exact([[has_t, len(t) - has_t], [has_b, len(b) - has_b]], alternative="greater") if len(t) else None
            # motif density along the promoter, per kb, in 100 bp bins (bin 0 = closest to the gene)
            nb = config.PROMOTER_LEN // BIN

            def density(group):
                counts = np.zeros(nb)
                for d in group[label + " positions"]:
                    for x in d:
                        counts[min(int((x - 1) // BIN), nb - 1)] += 1
                return (counts / max(len(group), 1) * (1000 / BIN)).round(3).tolist()

            res[label] = {"full": full, "proximal": prox,
                          "with_proximal_target": has_t, "with_proximal_background": has_b,
                          "fisher_p": float(fisher.pvalue) if fisher else None,
                          "odds_ratio": float(fisher.statistic) if fisher and np.isfinite(fisher.statistic) else None,
                          "density_target": density(t), "density_background": density(b)}
        sym = sym_at if sp == config.REFERENCE_SPECIES else sym_focus.to_dict()
        res["genes"] = [{"gene": r.gene, "group": r.group, "symbol": sym.get(r.gene, ""),
                         "motif": int(r[list(MOTIFS)[0]]), "proximal": int(r[list(MOTIFS)[0] + " proximal"]),
                         "positions": r[list(MOTIFS)[0] + " positions"]}
                        for _, r in sub[sub.group != "background"].iterrows()]
        out["species"][sp] = res
        m = res[list(MOTIFS)[0]]
        log(f"{sp}: targets {m['full'].get('mean_target')} vs background {m['full'].get('mean_background')} motifs per promoter, p={m['full'].get('p_perm')}")
    df.drop(columns=[c for c in df.columns if c.endswith("positions")]).to_csv(RESULTS / "promoter_motifs.csv", index=False)
    save_json(RESULTS / "promoters.json", out)


if __name__ == "__main__":
    main()
