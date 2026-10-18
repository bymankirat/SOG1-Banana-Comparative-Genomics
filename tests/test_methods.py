"""Unit tests of the building blocks, on small made-up inputs with known answers."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pipeline"))

from util import longest_isoforms, parse_header, revcomp, translate  # noqa: E402
from s05_proteins import count_sq_matrix, has_scd, shuffle_test, sq_positions  # noqa: E402
from s06_selection import codon_alignment  # noqa: E402
from s07_promoters import hits  # noqa: E402
from s04_phylogeny import classify, parse_newick  # noqa: E402

HEAD = ("AT1G25580.1 pep chromosome:TAIR10:1:8979977:8983008:-1 gene:AT1G25580 transcript:AT1G25580.1 "
        "gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:SOG1 description:NAC domain")


def test_parse_header():
    h = parse_header(HEAD)
    assert h["gene"] == "AT1G25580" and h["transcript"] == "AT1G25580.1"
    assert (h["seq_region"], h["start"], h["end"], h["strand"]) == ("1", 8979977, 8983008, -1)
    assert h["gene_symbol"] == "SOG1"


def test_longest_isoform_kept():
    recs = [(HEAD, "M" + "A" * 40), (HEAD.replace("AT1G25580.1", "AT1G25580.2"), "M" + "A" * 60 + "*"),
            (HEAD.replace("protein_coding transcript", "protein_coding transcript").replace("transcript_biotype:protein_coding", "transcript_biotype:nonsense_mediated_decay").replace("AT1G25580.1", "AT1G25580.3"), "M" * 99)]
    best = longest_isoforms(recs)
    assert best["AT1G25580"]["transcript"] == "AT1G25580.2"
    assert len(best["AT1G25580"]["seq"]) == 61  # stop removed, NMD isoform ignored


def test_translate_and_revcomp():
    assert translate("ATGTCTCAATAA") == "MSQ*"
    assert revcomp("AACGTT") == "AACGTT" and revcomp("ATGC") == "GCAT"


def test_sq_positions():
    assert sq_positions("MSQATQSQQ") == [1, 4, 6]
    assert count_sq_matrix(np.array([list("SQTQ"), list("QSQS")])).tolist() == [2, 1]


def test_scd_needs_three_sites_in_window():
    assert has_scd([10, 50, 105])
    assert not has_scd([10, 50, 115])
    assert not has_scd([10, 20])


def test_shuffle_test_detects_clustering():
    rng = np.random.default_rng(1)
    clustered = "SQ" * 6 + "A" * 88          # 6 SQ pairs; shuffling breaks most of them
    obs, exp, p = shuffle_test(clustered, 2000, rng)
    assert obs == 6 and exp < 2 and p < 0.01
    obs, exp, p = shuffle_test("SAQ" * 10, 2000, rng)  # no SQ, plenty of S and Q
    assert obs == 0 and p == 1.0


def test_codon_alignment_threads_gaps_and_rejects_mismatch():
    cds = "ATGTCTCAAGGTTAA"  # M S Q G stop
    assert codon_alignment("MS-QG", cds) == ["ATG", "TCT", "---", "CAA", "GGT"]
    assert codon_alignment("MSQA", cds) is None


def test_motif_hits_palindrome_counted_once_and_positioned():
    L = 100
    site = "CTTACGTACGAAG"  # CTT N7 AAG, a palindrome
    seq = "G" * 40 + site + "G" * (L - 40 - len(site))
    d = hits(seq, "CTT.......AAG")
    assert d == [L - 40]                     # one site, 60 bp upstream of the gene start
    assert hits("G" * 100, "CTT.......AAG") == []
    # non-palindromic control is searched on both strands
    ctrl = "G" * 10 + "CTTAAAAAAACTT" + "G" * 10
    assert len(hits(ctrl, "CTT.......CTT")) == 1
    assert len(hits(revcomp(ctrl), "CTT.......CTT")) == 1


def test_tree_clade_rule():
    # moss gene sits next to the SOG1 group; maize gene sits with a different Arabidopsis NAC
    nwk = ("(At__AT2G46770:0.5,((At__AT1G25580:0.1,Ma__g1:0.1)99:0.1,(At__AT3G01600:0.1,At__AT5G14490:0.1)100:0.1,"
           "Pp__g2:0.4)90:0.2,(At__AT4G28500:0.1,Zm__g3:0.1)100:0.2);")
    names = {k: k.replace("__", "|") for k in ["At__AT2G46770", "At__AT1G25580", "Ma__g1", "At__AT3G01600",
                                              "At__AT5G14490", "Pp__g2", "At__AT4G28500", "Zm__g3"]}
    calls = classify(parse_newick(nwk), names, competitors=["AT4G28500"])
    assert calls["Ma|g1"][:2] == (True, "SOG1") and calls["Ma|g1"][2] == 99
    assert calls["Pp|g2"][0] is True
    assert calls["Zm|g3"][0] is False
    assert calls["At|AT4G28500"][0] is False and calls["At|AT3G01600"] == (True, "NAC044", None)


def test_tree_pro_ortholog():
    # moss gene branches before SOG1 and its sister NAC group split
    nwk = "(At__AT2G46770:0.5,(((At__AT1G25580:0.1,At__AT3G01600:0.1,At__AT5G14490:0.1):0.1,At__AT4G28500:0.2):0.1,Pp__g2:0.4):0.2);"
    names = {k: k.replace("__", "|") for k in ["At__AT2G46770", "At__AT1G25580", "At__AT3G01600", "At__AT5G14490", "At__AT4G28500", "Pp__g2"]}
    assert classify(parse_newick(nwk), names, competitors=["AT4G28500"])["Pp|g2"][0] == "pro-ortholog"
    assert classify(parse_newick(nwk), names, competitors=[])["Pp|g2"][0] is False
