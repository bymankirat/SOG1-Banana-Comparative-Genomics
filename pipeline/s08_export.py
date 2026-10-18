"""Step 8: collect every result into docs/results.js for the interactive report."""
import json

import numpy as np
import pandas as pd

from util import DERIVED, RESULTS, ROOT, LABEL, LINEAGE, SHORT, config, load_json, log, read_fasta


def clean(x):
    """Replace NaN/inf anywhere in nested data with None so the output is valid JSON."""
    if isinstance(x, dict):
        return {k: clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, float) and not np.isfinite(x):
        return None
    return x


def records(df):
    return json.loads(df.replace({np.nan: None}).to_json(orient="records"))


def compara_agreement(clade):
    f = DERIVED / "compara_orthologs.csv"
    if not f.exists() or f.stat().st_size < 10:
        return None
    try:
        comp = pd.read_csv(f)
    except pd.errors.EmptyDataError:
        return None
    if comp.empty:
        return None
    rows = []
    for sp in sorted(set(clade.species)):
        if sp == config.REFERENCE_SPECIES:
            continue
        ours = set(clade[clade.species == sp].gene)
        theirs = set(comp[comp.species == sp].gene)
        if not theirs and sp not in set(comp.species):
            continue
        rows.append({"species": sp, "short": SHORT.get(sp, sp), "both": len(ours & theirs),
                     "hmm_only": sorted(ours - theirs), "compara_only": sorted(theirs - ours)})
    return rows


def main():
    prov = load_json(DERIVED / "provenance.json")
    counts = pd.read_csv(RESULTS / "nac_counts.csv")
    nac = pd.read_csv(DERIVED / "nac_genes.csv")
    clade = nac[nac.sog1_clade].copy()
    gs = pd.read_csv(DERIVED / "gene_structure.csv") if (DERIVED / "gene_structure.csv").exists() else pd.DataFrame(columns=["gene"])
    clade = clade.merge(gs[[c for c in gs.columns if c in ("gene", "coding_exons", "introns", "exon_lengths")]], on="gene", how="left")
    sq = pd.read_csv(RESULTS / "sq_motifs.csv")
    clade["key"] = clade.short + "|" + clade.gene
    clade = clade.merge(sq.drop(columns=["species", "short", "gene", "label", "length"]), on="key", how="left")
    clade["lineage"] = clade.species.map(LINEAGE)
    clade["species_label"] = clade.species.map(LABEL)

    tree = {}
    if (RESULTS / "tree" / "sog1_clade.treefile").exists():
        summ = load_json(RESULTS / "tree" / "summary.json")
        tree = {"newick": (RESULTS / "tree" / "sog1_clade.treefile").read_text().strip(), **summ}

    prot = load_json(RESULTS / "proteins.json")
    seqs = dict(read_fasta(DERIVED / "sog1_clade_proteins.faa"))
    sel = load_json(RESULTS / "selection.json") if (RESULTS / "selection.json").exists() else None
    pairs = pd.read_csv(RESULTS / "selection_pairs.csv") if (RESULTS / "selection_pairs.csv").exists() else None
    prom = load_json(RESULTS / "promoters.json")
    tor = pd.read_csv(DERIVED / "target_orthologs.csv")

    focus = [k for k in clade.key if k.startswith(SHORT[config.FOCUS_SPECIES] + "|")]
    out = {
        "meta": {"release": prov.get("ensembl_plants_release"), "downloaded": prov.get("downloaded"),
                 "species": [{"species": s[0], "short": s[1], "label": s[2], "lineage": s[3],
                              "included": s[0] in prov["species"]} for s in config.SPECIES],
                 "focus_species": config.FOCUS_SPECIES, "focus_gene": config.FOCUS_GENE,
                 "at_clade": config.AT_CLADE, "targets": config.SOG1_TARGETS},
        "counts": records(counts),
        "clade": records(clade.drop(columns=["seq_region", "start", "end", "strand"], errors="ignore")),
        "compara": compara_agreement(clade),
        "tree": tree,
        "sq": {"positions": prot["positions"], "conservation": prot["conservation"],
               "tail_vs_domain": prot["tail_vs_domain"], "alphafold": prot["alphafold"]},
        "profiles": {k: {"seq": seqs[k], "disorder": prot["disorder"].get(k), "plddt": prot["plddt"].get(k)}
                     for k in clade.key if k in focus or k.startswith("At|")},
        "selection": sel,
        "pairs": records(pairs[pairs.usable][["a", "b", "species_a", "species_b", "label_a", "label_b",
                                              "ds_full", "omega_full", "omega_domain", "omega_tail"]].round(4)) if pairs is not None else [],
        "targets": records(tor),
        "promoters": prom,
    }
    js = "window.SOG1 = " + json.dumps(clean(out), separators=(",", ":"), allow_nan=False) + ";\n"
    (ROOT / "docs" / "results.js").write_text(js)
    log(f"docs/results.js written ({len(js) / 1024:.0f} kB)")


if __name__ == "__main__":
    main()
