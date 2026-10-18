"""Step 5: SQ/TQ phosphorylation motifs, disorder and AlphaFold confidence along each SOG1-clade protein.

ATM and ATR activate SOG1 by phosphorylating serines/threonines followed by
glutamine (SQ/TQ) in its C-terminal tail. Questions:
  1. Does each protein's tail hold more SQ/TQ sites than its own amino-acid
     composition would give by chance? (shuffle the tail 10,000 times)
  2. Does it carry an SQ/TQ cluster domain (>= 3 sites in 100 residues)?
  3. Are the Arabidopsis SOG1 sites at the same alignment positions in banana?
  4. Is the tail disordered (metapredict), and does that agree with AlphaFold pLDDT?
"""
import hashlib
import re

import numpy as np
import pandas as pd
from scipy import stats

from util import DERIVED, RESULTS, config, load_json, log, read_fasta, save_json

SQ = re.compile(r"(?=([ST])Q)")


def sq_positions(seq):
    """0-based positions of the S/T in each S/T-Q pair."""
    return [m.start() for m in SQ.finditer(seq)]


def count_sq_matrix(mat):
    """mat: (n, L) array of single characters; returns the number of [ST]Q pairs in each row."""
    if mat.shape[1] < 2:
        return np.zeros(mat.shape[0], dtype=int)
    st = (mat[:, :-1] == "S") | (mat[:, :-1] == "T")
    return (st & (mat[:, 1:] == "Q")).sum(axis=1)


def shuffle_test(segment, n, rng):
    """One-sided permutation p-value that a segment has as many [ST]Q as observed, given its composition."""
    obs = len(sq_positions(segment))
    if len(segment) < 2:
        return obs, float("nan"), 1.0
    arr = np.array(list(segment))
    null = count_sq_matrix(rng.permuted(np.tile(arr, (n, 1)), axis=1))
    return obs, float(null.mean()), float((1 + (null >= obs).sum()) / (n + 1))


def has_scd(positions, window=config.SCD_WINDOW, k=config.SCD_MIN):
    pos = sorted(positions)
    return any(pos[i + k - 1] - pos[i] < window for i in range(len(pos) - k + 1))


def column_map(aligned):
    """aligned sequence -> list mapping residue index to alignment column."""
    return [c for c, a in enumerate(aligned) if a != "-"]


def disorder_scores(seqs):
    """metapredict disorder (0-1 per residue), cached by sequence so offline reruns work."""
    cache_file = RESULTS / "disorder_cache.json"
    cache = load_json(cache_file) if cache_file.exists() else {}
    try:
        import metapredict as meta
    except Exception:
        meta = None
    out, changed = {}, False
    for name, seq in seqs.items():
        key = hashlib.md5(seq.encode()).hexdigest()
        if key not in cache and meta is not None:
            cache[key] = [round(float(x), 3) for x in meta.predict_disorder(seq)]
            changed = True
        if key in cache:
            out[name] = cache[key]
    if changed:
        save_json(cache_file, cache)
    if meta is None:
        log(f"metapredict not installed; used cached disorder for {len(out)} of {len(seqs)} proteins")
    return out


def main():
    rng = np.random.default_rng(config.SEED)
    nac = pd.read_csv(DERIVED / "nac_genes.csv")
    clade = nac[nac.sog1_clade].copy()
    allseq = dict(read_fasta(DERIVED / "sog1_clade_proteins.faa"))
    seqs = {f"{r.short}|{r.gene}": allseq[f"{r.short}|{r.gene}"] for r in clade.itertuples()}

    rows, positions = [], {}
    for r in clade.itertuples():
        key = f"{r.short}|{r.gene}"
        seq = seqs[key]
        dom_a, dom_b = int(r.env_from) - 1, int(r.env_to)  # 0-based, end exclusive
        tail = seq[dom_b:]
        pos = sq_positions(seq)
        obs, exp, p = shuffle_test(tail, config.N_PERMUTATIONS, rng)
        tail_pos = [x for x in pos if x >= dom_b]
        rows.append({"key": key, "species": r.species, "short": r.short, "gene": r.gene, "label": r.clade_label,
                     "length": len(seq), "domain_start": dom_a + 1, "domain_end": dom_b, "tail_length": len(tail),
                     "sq_total": len(pos), "sq_domain": sum(dom_a <= x < dom_b for x in pos), "sq_tail": obs,
                     "sq_tail_expected": round(exp, 2), "sq_tail_p": p,
                     "sq_per_100_tail": round(100 * obs / max(len(tail), 1), 2),
                     "sq_per_100_domain": round(100 * sum(dom_a <= x < dom_b for x in pos) / max(dom_b - dom_a, 1), 2),
                     "scd": has_scd(tail_pos)})
        positions[key] = [x + 1 for x in pos]
    df = pd.DataFrame(rows)
    # Benjamini-Hochberg across proteins
    order = np.argsort(df.sq_tail_p.values)
    q = df.sq_tail_p.values[order] * len(df) / (np.arange(len(df)) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    df.loc[df.index[order], "sq_tail_q"] = np.minimum(q, 1)

    # tail vs domain density, paired across proteins
    w = stats.wilcoxon(df.sq_per_100_tail, df.sq_per_100_domain, alternative="greater") if len(df) > 5 else None

    # conservation of the Arabidopsis SOG1 sites, using the full alignment from step 4
    cons = []
    tree_sum = RESULTS / "tree" / "summary.json"
    if tree_sum.exists():
        names = load_json(tree_sum)["names"]
        aln = {names[h]: s for h, s in read_fasta(RESULTS / "tree" / "sog1_clade.aln.faa")}
        ref = f"At|{next(g for g, s in config.AT_CLADE.items() if s == 'SOG1')}"
        refmap = column_map(aln[ref])
        colpos = {k: {c: i for i, c in enumerate(column_map(s))} for k, s in aln.items()}
        for site in positions.get(ref, []):
            col = refmap[site - 1]
            entry = {"at_position": site, "motif": aln[ref].replace("-", "")[site - 1:site + 1], "column": col, "present": {}}
            for k, s in aln.items():
                if k not in seqs:
                    continue
                i = colpos[k].get(col)
                plain = seqs[k]
                # a site counts as conserved if S/T sits in the same column and is followed by Q
                entry["present"][k] = bool(i is not None and plain[i] in "ST" and i + 1 < len(plain) and plain[i + 1] == "Q")
            labels = df.set_index("key").label
            sog1_like = [k for k in entry["present"] if labels.get(k) == "SOG1" and k != ref]
            entry["conserved_in_sog1_orthologs"] = round(np.mean([entry["present"][k] for k in sog1_like]), 3) if sog1_like else None
            cons.append(entry)

    # disorder and AlphaFold confidence
    dis = disorder_scores(seqs)
    af = load_json(DERIVED / "alphafold_plddt.json") if (DERIVED / "alphafold_plddt.json").exists() else {}
    agree = []
    idx = df.set_index("key")
    for key in seqs:
        r = idx.loc[key]
        d = np.array(dis[key]) if key in dis else None
        if d is not None:
            df.loc[df.key == key, "disorder_tail"] = round(float((d[int(r.domain_end):] >= 0.5).mean()), 3) if r.tail_length else np.nan
            df.loc[df.key == key, "disorder_domain"] = round(float((d[int(r.domain_start) - 1:int(r.domain_end)] >= 0.5).mean()), 3)
        a = af.get(key)
        if a and a.get("same_sequence") and len(a["plddt"]) == len(seqs[key]):
            pl = np.array(a["plddt"])
            rho = stats.spearmanr(d, pl).statistic if d is not None else np.nan
            agree.append({"key": key, "uniprot": a["uniprot"], "rho": None if np.isnan(rho) else round(float(rho), 3),
                          "tail_low_plddt": round(float((pl[int(r.domain_end):] < 50).mean()), 3),
                          "domain_mean_plddt": round(float(pl[int(r.domain_start) - 1:int(r.domain_end)].mean()), 1)})
    df.to_csv(RESULTS / "sq_motifs.csv", index=False)
    save_json(RESULTS / "proteins.json", {
        "positions": positions, "conservation": cons, "disorder": dis,
        "plddt": {k: v["plddt"] for k, v in af.items() if v.get("same_sequence")},
        "alphafold": agree,
        "tail_vs_domain": {"median_tail": float(df.sq_per_100_tail.median()), "median_domain": float(df.sq_per_100_domain.median()),
                           "wilcoxon_p": float(w.pvalue) if w else None, "n": len(df)}})
    log(f"{int((df.sq_tail_q < 0.05).sum())} of {len(df)} tails enriched for SQ/TQ (q<0.05); SCD in {int(df.scd.sum())}")


if __name__ == "__main__":
    main()
