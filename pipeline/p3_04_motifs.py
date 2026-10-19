"""Phase 3, step 4: do the crop genes that respond to stress carry SOG1's binding motif in their promoters?

Phase 2 found the palindrome CTT(N)7AAG (Ogita et al. 2018) in Arabidopsis SOG1 target promoters. Expression alone
cannot say whether SOG1 itself drives a stress response, because thousands of genes move under stress. This step
adds the regulatory side: take the SOG1 target orthologs that go up in stress contrasts, fetch their promoters from
the reference genome, and compare the motif rate with a background of random genes from the same species.

A promoter is the 1,500 bp upstream of the annotated transcription start, on the gene's own strand (config.PROMOTER_LEN).
"""
import numpy as np
import pandas as pd

import config
from util import RESULTS, ensembl_get, log, save_json

OUT = RESULTS / "phase3"
RNG = np.random.default_rng(config.SEED)
MOTIF = config.SOG1_MOTIF  # CTT.......AAG
BG_PER_SPECIES = 400


def count_motif(seq):
    """CTT(N)7AAG is its own reverse complement, so one pass over the forward strand counts both strands."""
    import re
    return len(re.findall(f"(?=({MOTIF}))", seq.upper()))


def promoter(species, gene):
    g = ensembl_get(f"/lookup/id/{gene}", species=species)
    if not g or "seq_region_name" not in g:
        return None
    strand, start, end = g["strand"], g["start"], g["end"]
    if strand == 1:
        lo, hi = max(1, start - config.PROMOTER_LEN), start - 1
    else:
        lo, hi = end + 1, end + config.PROMOTER_LEN
    if hi <= lo:
        return None
    r = ensembl_get(f"/sequence/region/{species}/{g['seq_region_name']}:{lo}..{hi}:{strand}")
    if not r or "seq" not in r:
        return None
    return r["seq"].upper()


def background_genes(species, n):
    """Random protein-coding genes from the species, as a null set of promoters."""
    ids = ensembl_get(f"/info/assembly/{species}")
    regions = [r["name"] for r in (ids or {}).get("top_level_region", []) if r.get("coord_system") == "chromosome"]
    picked = []
    for region in regions:
        if len(picked) >= n:
            break
        length = next(r["length"] for r in ids["top_level_region"] if r["name"] == region)
        for _ in range(max(1, n // max(1, len(regions)))):
            pos = int(RNG.integers(1, max(2, length - 100000)))
            hits = ensembl_get(f"/overlap/region/{species}/{region}:{pos}..{pos + 100000}",
                               feature="gene", biotype="protein_coding") or []
            for h in hits[:3]:
                picked.append(h["id"])
            if len(picked) >= n:
                break
    return list(dict.fromkeys(picked))[:n]


def main():
    genes = pd.read_csv(OUT / "stress_genes.csv")
    if genes.empty:
        log("no stress gene table; nothing to do")
        return
    up = genes[(genes.padj < 0.05) & (genes.log2fc > 0)]
    sets = (up.groupby(["species", "gene", "symbol"]).size().reset_index(name="contrasts")
            .sort_values("contrasts", ascending=False))
    sets.to_csv(OUT / "stress_up_genes.csv", index=False)
    log(f"{len(sets)} species-gene pairs are up in at least one stress contrast")

    rows, bg_rows = [], []
    for species, grp in sets.groupby("species"):
        log(f"{species}: {len(grp)} responsive target orthologs")
        for _, r in grp.iterrows():
            seq = promoter(species, r.gene)
            if seq:
                rows.append({"species": species, "gene": r.gene, "symbol": r.symbol, "contrasts": r.contrasts,
                             "promoter_len": len(seq), "motifs": count_motif(seq),
                             "gc": round((seq.count("G") + seq.count("C")) / max(1, len(seq)), 3)})
        for gid in background_genes(species, BG_PER_SPECIES):
            seq = promoter(species, gid)
            if seq:
                bg_rows.append({"species": species, "gene": gid, "promoter_len": len(seq),
                                "motifs": count_motif(seq),
                                "gc": round((seq.count("G") + seq.count("C")) / max(1, len(seq)), 3)})

    tar = pd.DataFrame(rows)
    bg = pd.DataFrame(bg_rows)
    tar.to_csv(OUT / "motif_targets.csv", index=False)
    bg.to_csv(OUT / "motif_background.csv", index=False)

    out = []
    for species in sorted(set(tar.species) & set(bg.species)):
        t, b = tar[tar.species == species], bg[bg.species == species]
        t_rate, b_rate = (t.motifs > 0).mean(), (b.motifs > 0).mean()
        # permutation: random background sets of the same size
        k = len(t)
        draws = RNG.integers(0, len(b), size=(10000, k))
        null = (b.motifs.to_numpy()[draws] > 0).mean(axis=1)
        p = float((np.sum(null >= t_rate) + 1) / 10001)
        out.append({"species": species, "n_targets": int(k), "n_background": int(len(b)),
                    "target_rate": float(t_rate), "background_rate": float(b_rate),
                    "ratio": float(t_rate / b_rate) if b_rate else None, "p": p,
                    "target_gc": float(t.gc.mean()), "background_gc": float(b.gc.mean())})
        log(f"  {species}: motif in {t_rate:.0%} of responsive targets vs {b_rate:.0%} of random genes, p = {p:.4f}")
    save_json(OUT / "motif_summary.json", out)


if __name__ == "__main__":
    main()
