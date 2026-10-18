"""Step 3 (needs internet): gene structure, Ensembl Compara orthologs, AlphaFold models and promoters.

- Exon-intron structure of every SOG1-clade gene (Ensembl REST lookup).
- Ensembl Compara orthologs of SOG1/NAC044/NAC085: an independent check of the clade from step 2.
- AlphaFold DB per-residue confidence (pLDDT) for clade proteins that have a UniProt entry
  with the same sequence; low pLDDT is a good marker of disordered regions.
- 1.5 kb upstream of each target, clade and background gene (Ensembl REST sequence/region).
"""
import difflib
import json

import pandas as pd

from util import DERIVED, config, ensembl_get, ensembl_post, http, log, save_json, write_fasta


def gene_structures(clade):
    out = []
    ids = list(clade.gene)
    for i in range(0, len(ids), 50):
        try:
            res = ensembl_post("/lookup/id", {"ids": ids[i:i + 50], "expand": 1})
        except Exception as e:
            log("lookup failed:", e)
            continue
        for g, obj in (res or {}).items():
            if not obj:
                continue
            tid = clade.set_index("gene").transcript.get(g)
            txs = obj.get("Transcript", [])
            tx = next((t for t in txs if t.get("id") == tid), None) or next((t for t in txs if t.get("is_canonical")), None)
            if not tx:
                continue
            tr = tx.get("Translation") or {}
            exons = tx.get("Exon", [])
            cs, ce = tr.get("start"), tr.get("end")
            coding = [e for e in exons if cs and e["end"] >= cs and e["start"] <= ce]
            out.append({"gene": g, "transcript": tx.get("id"), "strand": obj.get("strand"),
                        "exons": len(exons), "coding_exons": len(coding), "introns": max(len(coding) - 1, 0),
                        "protein_length": tr.get("length"),
                        "exon_lengths": json.dumps([min(e["end"], ce) - max(e["start"], cs) + 1 for e in sorted(coding, key=lambda e: e["start"] * (obj.get("strand") or 1))])})
    return pd.DataFrame(out)


def compara():
    rows = []
    for g, sym in config.AT_CLADE.items():
        data = None
        for path in (f"/homology/id/{config.REFERENCE_SPECIES}/{g}", f"/homology/id/{g}"):
            try:
                data = ensembl_get(path, compara="plants", type="orthologues", format="condensed")
            except Exception as e:
                log(f"  compara {path}: {str(e)[:120]}")
                data = None
            if data and data.get("data"):
                break
        if not data or not data.get("data"):
            log(f"no Compara result for {sym}")
            continue
        for h in data["data"][0].get("homologies", []):
            rows.append({"at_gene": g, "at_symbol": sym, "species": h.get("species"), "gene": h.get("id"),
                         "type": h.get("type")})
    return pd.DataFrame(rows)


def uniprot_match(gene, transcript, seq):
    for q in (gene, transcript, transcript.split(".")[0] if transcript else None):
        if not q:
            continue
        try:
            r = http("GET", config.UNIPROT_SEARCH, params={"query": q, "fields": "accession,sequence", "format": "json", "size": 10})
        except Exception as e:
            log(f"  UniProt {q}: {str(e)[:100]}")
            continue
        best = None
        for hit in r.json().get("results", []):
            useq = hit.get("sequence", {}).get("value", "")
            if not useq or abs(len(useq) - len(seq)) > 0.05 * len(seq):
                continue
            ratio = 1.0 if useq == seq else difflib.SequenceMatcher(None, useq, seq, autojunk=False).ratio()
            if ratio >= 0.95 and (best is None or ratio > best[1]):
                best = (hit["primaryAccession"], ratio, useq)
        if best:
            return best
    return None


def alphafold(acc):
    r = http("GET", config.ALPHAFOLD_API.format(acc=acc), ok404=True)
    if r is None or not r.json():
        return None
    entry = r.json()[0]
    pdb = http("GET", entry["pdbUrl"]).text
    plddt = [float(l[60:66]) for l in pdb.splitlines() if l.startswith("ATOM") and l[12:16].strip() == "CA"]
    return {"model": entry.get("entryId"), "version": entry.get("latestVersion"), "plddt": plddt}


def structures(clade, seqs):
    out = {}
    for r in clade.itertuples():
        key = f"{r.short}|{r.gene}"
        m = uniprot_match(r.gene, r.transcript, seqs[key])
        if not m:
            continue
        acc, ratio, useq = m
        try:
            af = alphafold(acc)
        except Exception as e:
            log(f"  AlphaFold {acc}: {str(e)[:100]}")
            af = None
        if af:
            out[key] = {"uniprot": acc, "identity": round(ratio, 3), "same_sequence": useq == seqs[key], **af}
            log(f"  {key}: {acc} {af['model']} ({len(af['plddt'])} residues)")
    return out


def promoters(loc):
    L = config.PROMOTER_LEN
    recs, missing = [], 0
    for sp, sub in loc.groupby("species"):
        regions, meta = [], {}
        for r in sub.itertuples():
            if r.strand == 1:
                a, b = max(1, r.start - L), r.start - 1
            else:
                a, b = r.end + 1, r.end + L
            if b < a:
                continue
            reg = f"{r.seq_region}:{a}..{b}:{r.strand}"
            regions.append(reg)
            meta[reg] = r
        got = {}
        for i in range(0, len(regions), 50):
            batch = regions[i:i + 50]
            try:
                res = ensembl_post(f"/sequence/region/{sp}", {"regions": batch})
                if len(res) == len(batch):
                    for item, reg in zip(res, batch):
                        got[reg] = item.get("seq", "")
                else:
                    for item in res:
                        got[item.get("query", "")] = item.get("seq", "")
            except Exception as e:
                log(f"  batch failed ({str(e)[:80]}), fetching one by one")
                for reg in batch:
                    try:
                        item = ensembl_get(f"/sequence/region/{sp}/{reg}")
                        if item:
                            got[reg] = item.get("seq", "")
                    except Exception:
                        pass
        for reg in regions:
            seq = got.get(reg, "")
            r = meta[reg]
            if len(seq) < 0.5 * L:
                missing += 1
                continue
            recs.append((f"{sp}|{r.gene}|{r.group}|{reg}", seq.upper()))
        log(f"  {sp}: {sum(1 for h, _ in recs if h.startswith(sp))} promoters")
    write_fasta(DERIVED / "promoters.fna", recs)
    log(f"{len(recs)} promoters written, {missing} skipped (too short or not returned)")
    return len(recs)


def main():
    nac = pd.read_csv(DERIVED / "nac_genes.csv")
    clade = nac[nac.sog1_clade | ((nac.species == config.REFERENCE_SPECIES) & nac.gene.isin(list(config.AT_OUTGROUP)))]
    from util import read_fasta
    seqs = dict(read_fasta(DERIVED / "sog1_clade_proteins.faa"))

    log("gene structures")
    gene_structures(clade).to_csv(DERIVED / "gene_structure.csv", index=False)
    log("Ensembl Compara orthologs")
    compara().to_csv(DERIVED / "compara_orthologs.csv", index=False)
    log("promoters")
    if promoters(pd.read_csv(DERIVED / "promoter_genes.csv")) == 0:
        raise SystemExit("no promoter sequences retrieved")
    log("UniProt + AlphaFold")
    save_json(DERIVED / "alphafold_plddt.json", structures(clade, seqs))


if __name__ == "__main__":
    main()
