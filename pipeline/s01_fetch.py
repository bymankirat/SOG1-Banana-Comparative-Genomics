"""Step 1 (needs internet): download proteomes, coding sequences and the NAC domain model.

Sources: Ensembl Plants FTP (protein and CDS FASTA of every species in config.SPECIES)
and InterPro (Pfam profile HMM PF02365). Every file is fingerprinted (MD5) in
data/derived/provenance.json so a later run can show it used the same data.
"""
import re
import time

from util import RAW, DERIVED, config, download, ensembl_get, http, log, md5, save_json


def ensembl_release():
    try:
        info = ensembl_get("/info/eg_version")
        return str(info["version"]) if info else None
    except Exception as e:  # the REST call is only used to pin the release
        log("could not read Ensembl Genomes release:", e)
        return None


def listing(url):
    r = http("GET", url, ok404=True)
    return None if r is None else re.findall(r'href="([^"?/][^"]*)"', r.text)


def fasta_url(base, species, kind):
    """kind = 'pep' or 'cds'. Returns the URL of <Species>.<assembly>.<kind>.all.fa.gz or None."""
    folder = f"{base}/fasta/{species}/{kind}/"
    files = listing(folder)
    if not files:
        return None
    hits = [f for f in files if f.endswith(f".{kind}.all.fa.gz")]
    return folder + hits[0] if hits else None


def main():
    release = ensembl_release()
    base = f"{config.ENSEMBL_FTP}/release-{release}" if release else f"{config.ENSEMBL_FTP}/current"
    if release and listing(f"{base}/fasta/") is None:
        log(f"release-{release} folder not found, using current/")
        base = f"{config.ENSEMBL_FTP}/current"
    log("Ensembl Plants release:", release, "->", base)

    prov = {"ensembl_plants_release": release, "ftp_base": base, "downloaded": time.strftime("%Y-%m-%d"),
            "species": {}, "skipped": []}
    for sp, short, label, lineage, required in config.SPECIES:
        urls = {k: fasta_url(base, sp, k) for k in ("pep", "cds")}
        if not urls["pep"]:
            msg = f"{sp} not in this Ensembl Plants release"
            if required:
                raise SystemExit(msg)
            log("skip:", msg)
            prov["skipped"].append(sp)
            continue
        entry = {}
        for kind, url in urls.items():
            if not url:
                continue
            dest = RAW / url.rsplit("/", 1)[1]
            log(f"{short}: {dest.name}")
            download(url, dest)
            entry[kind] = {"file": dest.name, "url": url, "md5": md5(dest), "bytes": dest.stat().st_size}
        prov["species"][sp] = entry

    hmm = RAW / f"{config.PFAM_NAM}.hmm.gz"
    download(config.INTERPRO_HMM.format(acc=config.PFAM_NAM), hmm)
    prov["pfam_hmm"] = {"accession": config.PFAM_NAM, "md5": md5(hmm)}
    save_json(DERIVED / "provenance.json", prov)
    log(f"{len(prov['species'])} species downloaded, {len(prov['skipped'])} skipped")


if __name__ == "__main__":
    main()
