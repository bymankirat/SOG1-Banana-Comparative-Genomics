"""Phase 3, step 3: is SOG1 (and the genes it switches on) induced by stress, across crops?

For every stress contrast in the Expression Atlas we read the published per-gene log2 fold change and adjusted
p-value, pull out SOG1 and its target orthologs, and ask whether that set is shifted up compared with the rest of
the genome in the same contrast. Two tests, both against the same background:

  * sign test on the median log2 fold change (how far the set sits above the genome-wide median), and
  * a permutation test: 10,000 random gene sets of the same size from the same contrast.

Genotoxic contrasts (bleomycin, gamma rays, UV) act as the positive control: SOG1 targets are known to go up there,
so if the method cannot see it in those, nothing else in this step can be believed.
"""
import numpy as np
import pandas as pd

import config
from util import RAW, RESULTS, log, save_json

OUT = RESULTS / "phase3"
CACHE = RAW / "atlas"
N_PERM = 10000
RNG = np.random.default_rng(config.SEED)
MIN_SET = 3  # a set smaller than this is not worth testing


def analytics_path(acc, array):
    return CACHE / (f"{acc}_{array}-analytics.tsv" if array else f"{acc}-analytics.tsv")


def read_analytics(path):
    """-> {contrast_id: DataFrame(gene, log2fc, padj)} from one Atlas analytics table."""
    df = pd.read_csv(path, sep="\t", low_memory=False)
    gene_col = df.columns[0]
    out = {}
    for col in df.columns:
        if not col.endswith(".log2foldchange"):
            continue
        cid = col[: -len(".log2foldchange")]
        p = f"{cid}.p-value"
        if p not in df.columns:
            continue
        sub = df[[gene_col, col, p]].copy()
        sub.columns = ["gene", "log2fc", "padj"]
        sub["gene"] = sub.gene.astype(str).str.strip().str.upper()
        sub = sub.dropna(subset=["log2fc"])
        if len(sub):
            out[cid] = sub
    return out


def permutation(values, background, n=N_PERM):
    """One-sided p: how often does a random set of the same size have a median at least this high?"""
    obs = float(np.median(values))
    k = len(values)
    idx = RNG.integers(0, len(background), size=(n, k))
    null = np.median(background[idx], axis=1)
    p_up = (np.sum(null >= obs) + 1) / (n + 1)
    p_down = (np.sum(null <= obs) + 1) / (n + 1)
    return obs, float(np.median(null)), p_up, p_down


def main():
    con = pd.read_csv(OUT / "atlas_contrasts.csv")
    orth = pd.read_csv(OUT / "orthologs.csv")
    orth["gene"] = orth.gene.astype(str).str.strip().str.upper()
    by_species = {sp: g for sp, g in orth.groupby(orth.species.str.lower())}

    rows, missing = [], []
    for (acc, array), grp in con.groupby(["accession", con.array.fillna("")]):
        path = analytics_path(acc, array)
        if not path.exists():
            missing.append(str(path.name))
            continue
        try:
            tables = read_analytics(path)
        except Exception as e:  # noqa: BLE001
            missing.append(f"{path.name}: {e}")
            continue
        sp = grp.species.iloc[0].lower().replace(" ", "_")
        genes = by_species.get(sp)
        if genes is None:
            continue
        sog1_ids = set(genes[genes.at_symbol == "SOG1"].gene)
        target_ids = set(genes[genes.at_symbol != "SOG1"].gene)
        for _, c in grp.iterrows():
            tab = tables.get(c.contrast)
            if tab is None or len(tab) < 1000:
                continue
            bg = tab.log2fc.to_numpy(dtype=float)
            row = {"accession": acc, "array": array, "contrast": c.contrast, "name": c["name"],
                   "stress": c.stress, "species": sp, "n_genes": len(tab),
                   "bg_median": float(np.median(bg))}
            s1 = tab[tab.gene.isin(sog1_ids)]
            if len(s1):
                row["sog1_log2fc"] = float(s1.log2fc.max())
                row["sog1_padj"] = float(s1.padj.min())
            tg = tab[tab.gene.isin(target_ids)]
            row["n_targets_found"] = len(tg)
            if len(tg) >= MIN_SET:
                obs, null_med, p_up, p_down = permutation(tg.log2fc.to_numpy(dtype=float), bg)
                row.update({"targets_median_log2fc": obs, "null_median": null_med,
                            "p_up": p_up, "p_down": p_down,
                            "targets_sig_up": int(((tg.padj < 0.05) & (tg.log2fc > 0)).sum()),
                            "targets_sig_down": int(((tg.padj < 0.05) & (tg.log2fc < 0)).sum())})
                for _, t in tg.iterrows():
                    sym = genes[genes.gene == t.gene].at_symbol.iloc[0]
                    rows.append({"kind": "gene", "accession": acc, "contrast": c.contrast, "stress": c.stress,
                                 "species": sp, "symbol": sym, "gene": t.gene,
                                 "log2fc": float(t.log2fc), "padj": float(t.padj)})
            rows.append({"kind": "contrast", **row})

    df = pd.DataFrame(rows)
    contrasts = df[df.kind == "contrast"].drop(columns="kind")
    genes_df = df[df.kind == "gene"].drop(columns="kind")
    contrasts.to_csv(OUT / "stress_contrasts.csv", index=False)
    genes_df.to_csv(OUT / "stress_genes.csv", index=False)

    tested = contrasts.dropna(subset=["p_up"])
    if len(tested):
        # Benjamini-Hochberg across all tested contrasts
        p = tested.p_up.to_numpy()
        order = np.argsort(p)
        q = np.empty_like(p)
        q[order] = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
        tested = tested.assign(q_up=np.clip(q, 0, 1))
        tested.to_csv(OUT / "stress_contrasts.csv", index=False)
        summary = (tested.groupby("stress")
                   .agg(contrasts=("contrast", "size"), species=("species", "nunique"),
                        median_shift=("targets_median_log2fc", "median"),
                        significant=("q_up", lambda s: int((s < 0.05).sum())))
                   .sort_values("contrasts", ascending=False))
        log("\n" + summary.to_string())
        save_json(OUT / "stress_summary.json",
                  {"by_stress": summary.reset_index().to_dict("records"),
                   "n_contrasts": int(len(tested)), "n_experiments": int(tested.accession.nunique()),
                   "n_species": int(tested.species.nunique()), "missing": missing[:50],
                   "n_missing": len(missing)})
    log(f"{len(contrasts)} contrasts analysed, {len(missing)} tables missing")


if __name__ == "__main__":
    main()
