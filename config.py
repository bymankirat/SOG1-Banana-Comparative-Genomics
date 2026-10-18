"""Settings for the SOG1 comparative genomics pipeline.

Everything that decides what the analysis does lives here, so a run can be
repeated or changed without touching the code in pipeline/.
"""

ENSEMBL_REST = "https://rest.ensembl.org"
ENSEMBL_FTP = "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants"
PFAM_NAM = "PF02365"  # NAM / NAC domain
INTERPRO_HMM = "https://www.ebi.ac.uk/interpro/api/entry/pfam/{acc}?annotation=hmm"
UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"

# Ensembl Plants species, from green algae to flowering plants.
# (ensembl name, short code, display name, lineage, required)
# Species that are not in the current Ensembl Plants release are skipped and logged.
SPECIES = [
    ("chlamydomonas_reinhardtii", "Cr", "Chlamydomonas reinhardtii", "Green alga", False),
    ("chara_braunii", "Cb", "Chara braunii", "Charophyte alga", False),
    ("marchantia_polymorpha", "Mp", "Marchantia polymorpha", "Liverwort", False),
    ("physcomitrium_patens", "Pp", "Physcomitrium patens", "Moss", True),
    ("selaginella_moellendorffii", "Sm", "Selaginella moellendorffii", "Lycophyte", False),
    ("amborella_trichopoda", "Atr", "Amborella trichopoda", "Early angiosperm", False),
    ("musa_acuminata", "Ma", "Musa acuminata (banana)", "Monocot", True),
    ("musa_balbisiana", "Mb", "Musa balbisiana (wild banana)", "Monocot", False),
    ("dioscorea_rotundata", "Dr", "Dioscorea rotundata (yam)", "Monocot", False),
    ("ananas_comosus", "Ac", "Ananas comosus (pineapple)", "Monocot", False),
    ("oryza_sativa", "Os", "Oryza sativa (rice)", "Monocot", True),
    ("brachypodium_distachyon", "Bd", "Brachypodium distachyon", "Monocot", False),
    ("sorghum_bicolor", "Sb", "Sorghum bicolor", "Monocot", False),
    ("zea_mays", "Zm", "Zea mays (maize)", "Monocot", False),
    ("vitis_vinifera", "Vv", "Vitis vinifera (grapevine)", "Eudicot", False),
    ("solanum_lycopersicum", "Sl", "Solanum lycopersicum (tomato)", "Eudicot", False),
    ("populus_trichocarpa", "Pt", "Populus trichocarpa (poplar)", "Eudicot", False),
    ("glycine_max", "Gm", "Glycine max (soybean)", "Eudicot", False),
    ("arabidopsis_thaliana", "At", "Arabidopsis thaliana", "Eudicot", True),
]

FOCUS_SPECIES = "musa_acuminata"
REFERENCE_SPECIES = "arabidopsis_thaliana"

# The SOG1 clade in Arabidopsis: SOG1 and its two closest paralogs, which also act in the DNA damage response
AT_CLADE = {"AT1G25580": "SOG1", "AT3G01600": "NAC044", "AT5G14490": "NAC085"}
# Other Arabidopsis NACs used to root the tree
AT_OUTGROUP = {"AT1G52890": "ANAC019", "AT4G27410": "RD26", "AT2G46770": "NST1"}

# Banana gene studied during the 2023 NABI training project (set after the first run; None = not yet known)
FOCUS_GENE = None

# Published NAC family sizes, used as a sanity check of the domain search (not as a hard rule).
# Arabidopsis: Ooka et al. 2003 (105, older annotation; ~115 on TAIR10). Rice: Nuruzzaman et al. 2010 (151).
# Banana: Cenci et al. 2014 (167 on the first genome assembly).
EXPECTED_NAC = {
    "arabidopsis_thaliana": (105, (95, 140)),
    "oryza_sativa": (151, (120, 190)),
    "musa_acuminata": (167, (130, 220)),
}

# Arabidopsis genes reported as SOG1-regulated after DNA damage (most bound directly in ChIP-seq,
# Ogita et al. 2018 Plant J; Bourbousse et al. 2018 PNAS). Symbols are checked against Ensembl at run time.
SOG1_TARGETS = {
    "AT4G21070": "BRCA1", "AT5G20850": "RAD51", "AT5G66130": "RAD17", "AT4G02390": "PARP2",
    "AT1G04020": "BARD1", "AT4G19130": "RPA1E", "AT5G24280": "GMI1", "AT3G27060": "TSO2",
    "AT5G02220": "SMR4", "AT1G07500": "SMR5", "AT3G27630": "SMR7", "AT3G01600": "NAC044",
    "AT5G14490": "NAC085", "AT5G48720": "XRI1", "AT4G37490": "CYCB1;1", "AT1G02970": "WEE1",
    "AT3G25250": "OXI1", "AT1G19250": "FMO1",
}

# SOG1 binds the palindrome CTT(N)7AAG in target promoters (Ogita et al. 2018)
SOG1_MOTIF = "CTT.......AAG"
PROMOTER_LEN = 1500
BACKGROUND_N = 1000
SEED = 20231107

# ATM/ATR phosphorylation sites; an SQ/TQ cluster domain (SCD) = at least 3 sites within 100 residues
SQ_PATTERN = "[ST]Q"
SCD_WINDOW, SCD_MIN = 100, 3
N_PERMUTATIONS = 10000

# Pairs with synonymous divergence outside this range are not used for dN/dS (too close or saturated)
DS_RANGE = (0.01, 2.0)
