"""Phase 3, step 2: find SOG1 and its target genes in every species that has stress experiments.

SOG1 itself comes from the phase 2 Compara table. The Arabidopsis targets (config.SOG1_TARGETS) are mapped to each
species with Ensembl Compara orthology; one-to-one orthologs are preferred, and many-to-one copies are all kept
because a crop often has several copies of a repair gene.
"""
import pandas as pd

import config
from util import DERIVED, RESULTS, ensembl_get, log

OUT = RESULTS / "phase3"


def species_slug(name):
    return name.strip().lower().replace(" ", "_")


def compara(at_gene, target):
    data = ensembl_get(f"/homology/id/arabidopsis_thaliana/{at_gene}", type="orthologues",
                       target_species=target, format="condensed")
    rows = []
    for h in (data or {}).get("data", []):
        for x in h.get("homologies", []):
            rows.append({"gene": x.get("id"), "type": x.get("type")})
    return rows


def main():
    con = pd.read_csv(OUT / "atlas_contrasts.csv")
    species = sorted({species_slug(s) for s in con.species})
    log("species with stress contrasts:", species)
    sog1 = pd.read_csv(DERIVED / "compara_orthologs.csv")
    sog1 = sog1[sog1.at_symbol == "SOG1"]
    rows = []
    for sp in species:
        if sp == "arabidopsis_thaliana":
            rows.append({"species": sp, "at_gene": "AT1G25580", "at_symbol": "SOG1", "gene": "AT1G25580", "type": "self"})
            for g, sym in config.SOG1_TARGETS.items():
                rows.append({"species": sp, "at_gene": g, "at_symbol": sym, "gene": g, "type": "self"})
            continue
        s1 = sog1[sog1.species == sp]
        if s1.empty:  # not in the phase 2 table: ask Compara directly
            for h in compara("AT1G25580", sp):
                rows.append({"species": sp, "at_gene": "AT1G25580", "at_symbol": "SOG1", **h})
        for _, r in s1.iterrows():
            rows.append({"species": sp, "at_gene": "AT1G25580", "at_symbol": "SOG1", "gene": r.gene, "type": r.type})
        for g, sym in config.SOG1_TARGETS.items():
            hits = compara(g, sp)
            if not hits:
                log(f"  {sp}: no ortholog for {sym}")
            for h in hits:
                rows.append({"species": sp, "at_gene": g, "at_symbol": sym, **h})
    df = pd.DataFrame(rows).drop_duplicates()
    df.to_csv(OUT / "orthologs.csv", index=False)
    log(df.groupby(["species", "at_symbol"]).size().unstack(fill_value=0).T.to_string())


if __name__ == "__main__":
    main()
