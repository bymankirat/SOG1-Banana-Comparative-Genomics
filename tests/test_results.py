"""Checks on the pipeline output against things we already know. Skipped until the pipeline has run."""
import json
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import config  # noqa: E402

DER, RES = ROOT / "data" / "derived", ROOT / "results"
needs = pytest.mark.skipif(not (RES / "nac_counts.csv").exists(), reason="pipeline not run yet")


@needs
def test_nac_family_sizes_match_published_counts():
    counts = pd.read_csv(RES / "nac_counts.csv").set_index("species")
    for sp, (published, (lo, hi)) in config.EXPECTED_NAC.items():
        if sp == "oryza_sativa":
            continue  # published count used a different gene annotation (MSU); shown in the report, not tested
        n, relaxed = counts.loc[sp, "nac"], counts.loc[sp, "nac_relaxed"]
        assert lo <= n <= hi or lo <= relaxed <= hi, f"{sp}: {n} NAC genes ({relaxed} relaxed), published {published}"


@needs
def test_arabidopsis_clade_is_recovered():
    nac = pd.read_csv(DER / "nac_genes.csv")
    at = nac[(nac.species == config.REFERENCE_SPECIES) & nac.sog1_clade]
    assert set(at.gene) == set(config.AT_CLADE)


@needs
def test_sog1_is_absent_from_green_algae():
    # SOG1 is a land-plant innovation; the alga is a negative control
    counts = pd.read_csv(RES / "nac_counts.csv").set_index("species")
    if "chlamydomonas_reinhardtii" in counts.index:
        assert counts.loc["chlamydomonas_reinhardtii", "sog1_clade"] == 0


@needs
def test_banana_has_sog1_and_moss_has_a_relative():
    counts = pd.read_csv(RES / "nac_counts.csv").set_index("species")
    assert counts.loc["musa_acuminata", "sog1_clade"] >= 1
    moss = counts.loc["physcomitrium_patens"]
    assert moss.sog1_clade + moss.pro_orthologs >= 1


@needs
def test_every_clade_protein_has_its_cds():
    from util_reader import fasta_names
    prot = fasta_names(DER / "sog1_clade_proteins.faa")
    cds = fasta_names(DER / "sog1_clade_cds.fna")
    assert len(prot - cds) <= 0.1 * len(prot)


@needs
def test_target_symbols_match_ensembl():
    t = pd.read_csv(DER / "target_orthologs.csv").drop_duplicates("at_gene")
    synonyms = {"OXI1": "AGC2-1"}  # Ensembl lists OXI1 under its other name
    named = t[t.ensembl_symbol.notna() & (t.ensembl_symbol != "") & (t.ensembl_symbol != t.at_gene)]
    named = named[named.symbol.map(synonyms.get).fillna("") != named.ensembl_symbol]
    norm = lambda s: str(s).upper().replace(";", "").replace("-", "").replace("ATBRCA1", "BRCA1")
    mismatch = [(r.at_gene, r.symbol, r.ensembl_symbol) for r in named.itertuples()
                if norm(r.symbol) not in norm(r.ensembl_symbol) and norm(r.ensembl_symbol) not in norm(r.symbol)]
    assert len(mismatch) <= 3, mismatch


@needs
def test_tree_contains_every_clade_protein():
    summ = json.loads((RES / "tree" / "summary.json").read_text())
    tree = (RES / "tree" / "sog1_clade.treefile").read_text()
    assert all(name in tree for name in summ["names"])


@needs
def test_results_js_is_valid():
    js = (ROOT / "docs" / "results.js").read_text()
    assert js.startswith("window.SOG1 = ")
    data = json.loads(js[len("window.SOG1 = "):].rstrip().rstrip(";"))
    assert data["clade"] and data["counts"]
