"""Step 2: find every NAC gene in every genome, and the SOG1 clade within them.

1. Longest protein per protein-coding gene.
2. hmmsearch with the Pfam NAM domain (PF02365) using Pfam's own gathering
   thresholds, so the cut-off is not tuned by me.
3. Each NAC protein is searched (phmmer) against all Arabidopsis NAC proteins.
   A protein belongs to the SOG1 clade when its best Arabidopsis match is SOG1,
   NAC044 or NAC085. The tree in step 4 checks this assignment independently.
"""
import gzip
import io

import pyhmmer

from util import DERIVED, RAW, RESULTS, LABEL, LINEAGE, SHORT, config, load_json, log, longest_isoforms, read_fasta, write_fasta

ALPHABET = pyhmmer.easel.Alphabet.amino()


def _s(x):
    return x.decode() if isinstance(x, (bytes, bytearray)) else x


def digitize(items):
    """items: list of (name, seq) -> DigitalSequenceBlock"""
    out = []
    for name, seq in items:
        try:
            ts = pyhmmer.easel.TextSequence(name=name, sequence=seq)
        except TypeError:
            ts = pyhmmer.easel.TextSequence(name=name.encode(), sequence=seq)
        out.append(ts.digitize(ALPHABET))
    return pyhmmer.easel.DigitalSequenceBlock(ALPHABET, out)


def load_hmm():
    raw = (RAW / f"{config.PFAM_NAM}.hmm.gz").read_bytes()
    data = gzip.decompress(raw) if raw[:2] == b"\x1f\x8b" else raw
    with pyhmmer.plan7.HMMFile(io.BytesIO(data)) as fh:
        return fh.read()


def proteome(sp, prov):
    return longest_isoforms(read_fasta(RAW / prov["species"][sp]["pep"]["file"]))


def search_domain(hmm, prot, names=None, gathering=True):
    names = sorted(prot) if names is None else sorted(names)
    if not names:
        return {}
    block = digitize([(g, prot[g]["seq"]) for g in names])
    kw = {"bit_cutoffs": "gathering"} if gathering else {"E": 1e-5, "domE": 1e-3}
    found = {}
    for hits in pyhmmer.hmmsearch([hmm], block, **kw):
        for hit in hits:
            if not hit.included:
                continue
            doms = [d for d in hit.domains if d.included] or list(hit.domains)
            d = max(doms, key=lambda d: d.score)
            found[_s(hit.name)] = {"score": round(hit.score, 1), "env_from": d.env_from, "env_to": d.env_to,
                                   "n_domains": len(doms)}
    return found


def best_hits(queries, targets, top=1, E=1e-5):
    """phmmer each query against targets. Returns {query: [(target, bitscore, evalue), ...]} best first."""
    qblock, tblock = digitize(queries), digitize(targets)
    out = {}
    for (qname, _), hits in zip(queries, pyhmmer.phmmer(qblock, tblock, E=E)):
        q = getattr(hits, "query_name", None) or getattr(getattr(hits, "query", None), "name", None)
        if q is not None and _s(q) != qname:
            raise RuntimeError(f"phmmer results out of order: {_s(q)} vs {qname}")
        out[qname] = sorted(((_s(h.name), h.score, h.evalue) for h in hits), key=lambda x: -x[1])[:top]
    return out


def main():
    import pandas as pd
    prov = load_json(DERIVED / "provenance.json")
    hmm = load_hmm()
    species = list(prov["species"])
    at = config.REFERENCE_SPECIES
    prot, nac, counts, extra = {}, {}, [], {}
    for sp in species:
        prot[sp] = proteome(sp, prov)
        nac[sp] = search_domain(hmm, prot[sp])
        relaxed = search_domain(hmm, prot[sp], gathering=False)
        exp = config.EXPECTED_NAC.get(sp)
        counts.append({"species": sp, "short": SHORT[sp], "label": LABEL[sp], "lineage": LINEAGE[sp],
                       "genes": len(prot[sp]), "nac": len(nac[sp]), "nac_relaxed": len(relaxed),
                       "published": exp[0] if exp else None,
                       "in_range": (exp[1][0] <= len(nac[sp]) <= exp[1][1] or exp[1][0] <= len(relaxed) <= exp[1][1]) if exp else None})
        log(f"{SHORT[sp]:>4}: {len(prot[sp]):6d} protein-coding genes, {len(nac[sp]):4d} NAC (Pfam cut-off), {len(relaxed):4d} at E<1e-5")
        extra[sp] = relaxed

    for g in list(config.AT_CLADE) + list(config.AT_OUTGROUP):
        if g not in nac[at]:
            raise SystemExit(f"{g} not detected as NAC in Arabidopsis; check the domain search")

    # Genes the Pfam cut-off misses (e.g. a slightly divergent domain) but that clearly match a clade protein
    clade_q = [(g, prot[at][g]["seq"]) for g in config.AT_CLADE]
    for sp in species:
        hits = best_hits(clade_q, [(g, prot[sp][g]["seq"]) for g in sorted(prot[sp])], top=25, E=1e-30)
        new = {h[0] for q in hits for h in hits[q]} - set(nac[sp])
        for g in new:
            if g in extra[sp]:
                nac[sp][g] = {**extra[sp][g], "pfam_cutoff": False}
                log(f"  {SHORT[sp]} {g}: matches a SOG1-clade protein, NAC domain below Pfam cut-off; added")
    at_nac = [(g, prot[at][g]["seq"]) for g in sorted(nac[at])]

    rows, nac_fa = [], []
    for sp in species:
        queries = [(g, prot[sp][g]["seq"]) for g in sorted(nac[sp])]
        if not queries:
            continue
        hits = best_hits(queries, at_nac, top=8)
        for g, _ in queries:
            p, d = prot[sp][g], nac[sp][g]
            ranked = [h for h in hits.get(g, []) if not (sp == at and h[0] == g)]
            best = ranked[0] if ranked else (None, 0.0, 1.0)
            best_clade = next((h for h in ranked if h[0] in config.AT_CLADE), None)
            best_other = next((h for h in ranked if h[0] not in config.AT_CLADE), None)
            if sp == at:
                candidate = g in config.AT_CLADE
            else:
                # candidate: a clade protein is among the 8 best Arabidopsis matches and scores close to the best;
                # the tree (step 4) makes the final call
                candidate = bool(best_clade) and best_clade[1] >= 0.7 * best[1] and len(p["seq"]) >= config.MIN_LENGTH
            rows.append({"species": sp, "short": SHORT[sp], "gene": g, "transcript": p["transcript"],
                         "symbol": p["symbol"], "length": len(p["seq"]), "pfam_cutoff": d.get("pfam_cutoff", True),
                         "score": d["score"], "env_from": d["env_from"], "env_to": d["env_to"], "n_domains": d["n_domains"],
                         "best_at": best[0], "best_score": round(best[1], 1),
                         "best_clade_at": best_clade[0] if best_clade else None,
                         "clade_score": round(best_clade[1], 1) if best_clade else None,
                         "other_at": best_other[0] if best_other else None,
                         "other_score": round(best_other[1], 1) if best_other else None,
                         "best_hit_rule": bool(g in config.AT_CLADE) if sp == at else bool(best[0] in config.AT_CLADE),
                         "candidate": bool(candidate),
                         "seq_region": p["seq_region"], "start": p["start"], "end": p["end"], "strand": p["strand"]})
            nac_fa.append((f"{SHORT[sp]}|{g}", p["seq"]))
    df = pd.DataFrame(rows)

    # Arabidopsis NACs that compete with the clade for these candidates go into the tree, so the tree can tell them apart
    cand = df[df.candidate & (df.species != at)]
    competitors = [g for g in cand.other_at.value_counts().index if g not in config.AT_CLADE][:config.N_COMPETITORS]
    tree_at = list(config.AT_CLADE) + [g for g in competitors if g not in config.AT_OUTGROUP] + list(config.AT_OUTGROUP)
    df["in_tree"] = df.candidate | ((df.species == at) & df.gene.isin(tree_at))
    df.to_csv(DERIVED / "nac_genes.csv", index=False)
    write_fasta(DERIVED / "nac_proteins.faa", nac_fa)
    seqs = dict(nac_fa)
    write_fasta(DERIVED / "sog1_candidates.faa", [(f"{r.short}|{r.gene}", seqs[f"{r.short}|{r.gene}"]) for r in df[df.in_tree].itertuples()])
    per = df[df.candidate].groupby("species").size()
    for c in counts:
        c["candidates"] = int(per.get(c["species"], 0))
    pd.DataFrame(counts).to_csv(RESULTS / "nac_counts.csv", index=False)
    log("SOG1-clade candidates per species:", {k: int(v) for k, v in per.items()}, "| competitors in tree:", competitors)

    # coding sequences of every protein in the tree (for dN/dS)
    want = {r.transcript: f"{r.short}|{r.gene}" for r in df[df.in_tree].itertuples()}
    cds = []
    for sp in species:
        if "cds" not in prov["species"][sp]:
            continue
        for head, seq in read_fasta(RAW / prov["species"][sp]["cds"]["file"]):
            tid = head.split()[0]
            if tid in want:
                cds.append((want[tid], seq.upper()))
    write_fasta(DERIVED / "sog1_clade_cds.fna", cds)
    log(f"{len(cds)} of {len(want)} coding sequences found")

    # orthologs of the reported SOG1 target genes in the focus species, by reciprocal best hit
    focus = config.FOCUS_SPECIES
    at_all = [(g, prot[at][g]["seq"]) for g in sorted(prot[at])]
    fo_all = [(g, prot[focus][g]["seq"]) for g in sorted(prot[focus])]
    tq = [(g, prot[at][g]["seq"]) for g in config.SOG1_TARGETS if g in prot[at]]
    fwd = best_hits(tq, fo_all, top=6)
    cand = sorted({h[0] for g in fwd for h in fwd[g] if h[2] < 1e-10})
    back = best_hits([(g, prot[focus][g]["seq"]) for g in cand], at_all, top=1)
    trows = []
    for g, sym in config.SOG1_TARGETS.items():
        ens_sym = prot[at].get(g, {}).get("symbol", "")
        for rank, (h, score, ev) in enumerate(fwd.get(g, [])):
            bh = back.get(h, [])
            if ev < 1e-10 and bh and bh[0][0] == g:
                trows.append({"at_gene": g, "symbol": sym, "ensembl_symbol": ens_sym, "focus_gene": h,
                              "score": round(score, 1), "evalue": ev, "rank": rank})
        if not any(r["at_gene"] == g for r in trows):
            trows.append({"at_gene": g, "symbol": sym, "ensembl_symbol": ens_sym, "focus_gene": None,
                          "score": None, "evalue": None, "rank": None})
    pd.DataFrame(trows).to_csv(DERIVED / "target_orthologs.csv", index=False)
    log(f"{sum(r['focus_gene'] is not None for r in trows)} banana orthologs for {len(config.SOG1_TARGETS)} target genes")

    # gene locations for promoter extraction: targets, clade, and a random background
    import random
    rng = random.Random(config.SEED)
    loc = []
    for sp, targets in ((at, set(config.SOG1_TARGETS)), (focus, {r["focus_gene"] for r in trows if r["focus_gene"]})):
        clade = set(df[(df.species == sp) & df.candidate].gene)
        pool = sorted(g for g in prot[sp] if g not in targets and g not in clade and prot[sp][g]["start"])
        bg = rng.sample(pool, min(config.BACKGROUND_N, len(pool)))
        for group, genes in (("target", sorted(targets)), ("sog1_clade", sorted(clade)), ("background", bg)):
            for g in genes:
                p = prot[sp].get(g)
                if p and p["start"]:
                    loc.append({"species": sp, "gene": g, "group": group, "seq_region": p["seq_region"],
                                "start": p["start"], "end": p["end"], "strand": p["strand"]})
    pd.DataFrame(loc).to_csv(DERIVED / "promoter_genes.csv", index=False)


if __name__ == "__main__":
    main()
