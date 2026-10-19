"""Phase 3, step 1: find public plant stress experiments in the EBI Expression Atlas and download their results.

Expression Atlas re-analyses public RNA-seq and microarray studies with one pipeline and publishes, for every
comparison ("contrast"), a log2 fold change and adjusted p-value for every gene. That lets us ask the same question
of many studies without re-processing raw reads: is SOG1 (and the genes it switches on) induced by a stress?
"""
import re
import xml.etree.ElementTree as ET

import pandas as pd

from util import RAW, RESULTS, download, http, log, save_json

ATLAS_LIST = "https://www.ebi.ac.uk/gxa/json/experiments"
ATLAS_FTP = "https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments"
OUT = RESULTS / "phase3"
CACHE = RAW / "atlas"

SPECIES = ["arabidopsis thaliana", "oryza sativa", "zea mays", "solanum lycopersicum", "triticum aestivum",
           "hordeum vulgare", "glycine max", "sorghum bicolor", "musa acuminata", "brachypodium distachyon",
           "physcomitrium patens", "physcomitrella patens", "vitis vinifera", "solanum tuberosum"]

# Stress classes, checked in this order against the contrast name (first match wins)
STRESS = [
    ("DNA damage", r"bleomycin|zeocin|mitomycin|gamma|ionizing|ionising|x-ray|\bx ray|genotoxic|dna damage|"
                   r"methyl methanesulfonate|\bmms\b|cisplatin|hydroxyurea|camptothecin|aphidicolin|uv-b|uvb|\buv\b|ultraviolet"),
    ("Heat", r"\bheat|high temperature|thermo|\b3[7-9] ?°?c|\b4[0-5] ?°?c|warm"),
    ("Cold", r"\bcold|chill|low temperature|freez|\b4 ?°?c\b"),
    ("Drought", r"drought|water deficit|water stress|dehydrat|desiccat|withhold|peg\b|polyethylene glycol|osmotic|mannitol"),
    ("Salt", r"\bsalt|nacl|salin|sodium chloride"),
    ("Oxidative", r"oxidative|hydrogen peroxide|\bh2o2|methyl viologen|paraquat|ozone|high light"),
    ("Flooding / hypoxia", r"flood|submerg|hypoxi|anoxi|waterlog"),
]
KEYWORDS = re.compile("|".join(p for _, p in STRESS) + r"|stress|abiotic", re.I)


def classify(text):
    t = text.lower()
    for name, pat in STRESS:
        if re.search(pat, t):
            return name
    return None


def experiment_list():
    r = http("GET", ATLAS_LIST, headers={"Accept": "application/json"})
    data = r.json()
    exps = data.get("experiments", data if isinstance(data, list) else [])
    log(f"Expression Atlas lists {len(exps)} experiments; keys: {sorted(exps[0].keys()) if exps else '-'}")
    rows = []
    for e in exps:
        species = str(e.get("species", "")).lower()
        etype = str(e.get("experimentType", e.get("rawExperimentType", "")))
        desc = str(e.get("experimentDescription", e.get("description", "")))
        if "DIFFERENTIAL" not in etype.upper():
            continue
        if not any(species.startswith(s) for s in SPECIES):
            continue
        factors = " ".join(map(str, e.get("experimentalFactors", []) or []))
        if not KEYWORDS.search(desc + " " + factors):
            continue
        rows.append({"accession": e.get("experimentAccession"), "species": species, "type": etype,
                     "description": desc, "factors": factors,
                     "arrays": ";".join(e.get("arrayDesigns", []) or []) if isinstance(e.get("arrayDesigns"), list) else ""})
    return pd.DataFrame(rows)


def contrasts(acc):
    """Contrast ids and names from the experiment's configuration file."""
    path = download(f"{ATLAS_FTP}/{acc}/{acc}-configuration.xml", CACHE / f"{acc}-configuration.xml")
    root = ET.parse(path).getroot()
    out = []
    for ad in root.iter("analytics"):
        array = (ad.findtext("array_design") or "").strip()
        for c in ad.iter("contrast"):
            out.append({"accession": acc, "array": array, "contrast": c.get("id"),
                        "name": (c.findtext("name") or "").strip()})
    return out


def analytics(acc, array):
    name = f"{acc}_{array}-analytics.tsv" if array else f"{acc}-analytics.tsv"
    return download(f"{ATLAS_FTP}/{acc}/{name}", CACHE / name)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)
    cat = experiment_list()
    log(f"{len(cat)} differential stress-related experiments in the chosen species")
    log(cat.groupby("species").size().to_string())
    cat.to_csv(OUT / "atlas_experiments.csv", index=False)

    rows, failed = [], []
    for acc in cat.accession:
        try:
            cs = contrasts(acc)
        except Exception as e:  # noqa: BLE001 - one broken experiment must not stop the rest
            failed.append({"accession": acc, "step": "configuration", "error": str(e)[:200]})
            continue
        for c in cs:
            c["stress"] = classify(c["name"])
            rows.append(c)
    con = pd.DataFrame(rows)
    keep = con[con.stress.notna()].copy()
    log(f"{len(con)} contrasts, {len(keep)} classified as a stress: " + keep.stress.value_counts().to_dict().__repr__())

    for (acc, array), _ in keep.groupby(["accession", "array"]):
        try:
            analytics(acc, array)
        except Exception as e:  # noqa: BLE001
            failed.append({"accession": acc, "step": f"analytics {array}", "error": str(e)[:200]})
    con.to_csv(OUT / "atlas_contrasts_all.csv", index=False)
    keep.merge(cat[["accession", "species", "description"]], on="accession").to_csv(OUT / "atlas_contrasts.csv", index=False)
    save_json(OUT / "atlas_failed.json", failed)
    log(f"done; {len(failed)} downloads failed (see atlas_failed.json)")


if __name__ == "__main__":
    main()
