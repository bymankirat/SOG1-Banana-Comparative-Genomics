"""Step 6: how strongly is each part of SOG1 constrained? Pairwise dN/dS (Nei-Gojobori 1986).

Coding sequences are threaded onto the protein alignment from step 4 (codon
alignment), then for every pair of clade genes dN, dS and omega = dN/dS are
computed separately for the NAC domain and the C-terminal tail (columns
defined by Arabidopsis SOG1). Pairs whose dS is below 0.01 or above 2 are
dropped (too similar to estimate, or saturated). Paralog pairs inside the
banana genome also give Ks, which dates their duplication relative to the
known banana whole-genome duplications.
"""
import itertools

import numpy as np
import pandas as pd
from scipy import stats

from util import DERIVED, RESULTS, config, load_json, log, read_fasta, save_json, translate

STOPS = {"TAA", "TAG", "TGA"}


def codon_alignment(prot_aln, cds):
    """Thread a CDS onto an aligned protein. Returns list of codons ('---' for gaps) or None if they disagree."""
    cds = cds.upper()
    plain = prot_aln.replace("-", "")
    if len(cds) >= 3 and cds[-3:] in STOPS:
        cds = cds[:-3]
    if len(cds) != 3 * len(plain) or translate(cds) != plain:
        return None
    out, i = [], 0
    for a in prot_aln:
        if a == "-":
            out.append("---")
        else:
            out.append(cds[3 * i:3 * i + 3])
            i += 1
    return out


def dnds(c1, c2, cols):
    from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds
    a, b = [], []
    for c in cols:
        x, y = c1[c], c2[c]
        if "-" in x or "-" in y or x in STOPS or y in STOPS or set(x + y) - set("ACGT"):
            continue
        a.append(x)
        b.append(y)
    if len(a) < 30:
        return np.nan, np.nan, len(a)
    try:
        dn, ds = cal_dn_ds(CodonSeq("".join(a)), CodonSeq("".join(b)), method="NG86")
    except Exception:
        return np.nan, np.nan, len(a)
    return dn, ds, len(a)


def main():
    out_csv = RESULTS / "selection_pairs.csv"
    try:
        import Bio  # noqa: F401
    except ImportError:
        if out_csv.exists():
            log("Biopython not installed; keeping existing selection results")
            return
        raise SystemExit("Biopython is needed for step 6")

    names = load_json(RESULTS / "tree" / "summary.json")["names"]
    aln = {names[h]: s for h, s in read_fasta(RESULTS / "tree" / "sog1_clade.aln.faa")}
    cds = dict(read_fasta(DERIVED / "sog1_clade_cds.fna"))
    nac = pd.read_csv(DERIVED / "nac_genes.csv")
    clade = nac[nac.sog1_clade].assign(key=lambda d: d.short + "|" + d.gene).set_index("key")

    codons, bad = {}, []
    for k in clade.index:
        c = codon_alignment(aln[k], cds.get(k, "")) if k in aln else None
        if c is None:
            bad.append(k)
        else:
            codons[k] = c
    log(f"{len(codons)} codon alignments, {len(bad)} CDS did not match their protein: {bad[:5]}")

    ref = f"At|{next(g for g, s in config.AT_CLADE.items() if s == 'SOG1')}"
    refcols = [c for c, a in enumerate(aln[ref]) if a != "-"]
    d0, d1 = refcols[int(clade.loc[ref, "env_from"]) - 1], refcols[int(clade.loc[ref, "env_to"]) - 1]
    ncol = len(aln[ref])
    regions = {"full": range(ncol), "domain": range(d0, d1 + 1), "tail": range(d1 + 1, ncol)}

    rows = []
    for a, b in itertools.combinations(sorted(codons), 2):
        row = {"a": a, "b": b, "species_a": clade.loc[a, "species"], "species_b": clade.loc[b, "species"],
               "label_a": clade.loc[a, "clade_label"], "label_b": clade.loc[b, "clade_label"]}
        for reg, cols in regions.items():
            dn, ds, n = dnds(codons[a], codons[b], cols)
            row.update({f"dn_{reg}": dn, f"ds_{reg}": ds, f"codons_{reg}": n,
                        f"omega_{reg}": dn / ds if ds and ds > 0 and np.isfinite(ds) else np.nan})
        rows.append(row)
    df = pd.DataFrame(rows)
    lo, hi = config.DS_RANGE
    ok = lambda reg: df[f"ds_{reg}"].between(lo, hi) & (df[f"dn_{reg}"] >= 0) & df[f"omega_{reg}"].notna()
    df["usable"] = ok("full") & ok("domain") & ok("tail")
    df.round(4).to_csv(out_csv, index=False)

    use = df[df.usable]
    w = stats.wilcoxon(use.omega_tail, use.omega_domain, alternative="greater") if len(use) > 5 else None
    focus = config.FOCUS_SPECIES
    para = df[(df.species_a == focus) & (df.species_b == focus)]
    save_json(RESULTS / "selection.json", {
        "pairs": len(df), "usable_pairs": int(len(use)), "ds_range": [lo, hi],
        "domain_columns": [d0, d1], "codon_alignments": len(codons), "cds_mismatch": bad,
        "median_omega_domain": float(use.omega_domain.median()) if len(use) else None,
        "median_omega_tail": float(use.omega_tail.median()) if len(use) else None,
        "median_omega_full": float(use.omega_full.median()) if len(use) else None,
        "wilcoxon_p": float(w.pvalue) if w else None,
        "tail_higher_fraction": float((use.omega_tail > use.omega_domain).mean()) if len(use) else None,
        "focus_paralogs": [{"a": r.a, "b": r.b, "ks": None if not np.isfinite(r.ds_full) else round(float(r.ds_full), 3),
                            "omega": None if not np.isfinite(r.omega_full) else round(float(r.omega_full), 3),
                            "label_a": r.label_a, "label_b": r.label_b} for r in para.itertuples()],
    })
    if len(use):
        log(f"omega domain {use.omega_domain.median():.3f} vs tail {use.omega_tail.median():.3f} over {len(use)} pairs")


if __name__ == "__main__":
    main()
