"""Small helpers shared by every step: paths, logging, FASTA, HTTP with retries."""
import gzip
import hashlib
import json
import re
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import config  # noqa: E402

RAW = ROOT / "data" / "raw"          # downloads (not committed)
DERIVED = ROOT / "data" / "derived"  # small extracts of public data (committed)
RESULTS = ROOT / "results"
for _d in (RAW, DERIVED, RESULTS):
    _d.mkdir(parents=True, exist_ok=True)

SHORT = {s[0]: s[1] for s in config.SPECIES}
LABEL = {s[0]: s[2] for s in config.SPECIES}
LINEAGE = {s[0]: s[3] for s in config.SPECIES}


def log(*parts):
    print(time.strftime("%H:%M:%S"), *parts, flush=True)


# ---------------------------------------------------------------- FASTA
def read_fasta(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    records, head, buf = [], None, []
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if head is not None:
                    records.append((head, "".join(buf)))
                head, buf = line[1:], []
            else:
                buf.append(line)
    if head is not None:
        records.append((head, "".join(buf)))
    return records


def write_fasta(path, records, width=60):
    with open(path, "w") as fh:
        for head, seq in records:
            fh.write(f">{head}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


def parse_header(head):
    """Ensembl FASTA header -> dict with id, gene, transcript, biotype, location."""
    out = {"id": head.split()[0]}
    for key, val in re.findall(r"(\w+):(\S+)", head):
        out.setdefault(key, val)
    loc = re.search(r"(?:chromosome|scaffold|primary_assembly|supercontig|contig|plasmid|chloroplast|mitochondrion):[^:\s]*:([^:\s]+):(\d+):(\d+):(-?1)", head)
    if loc:
        out["seq_region"], out["start"], out["end"], out["strand"] = loc.group(1), int(loc.group(2)), int(loc.group(3)), int(loc.group(4))
    return out


def longest_isoforms(records):
    """Keep the longest protein of each protein-coding gene. Returns {gene: info dict with seq}."""
    best = {}
    for head, seq in records:
        h = parse_header(head)
        if h.get("transcript_biotype", "protein_coding") != "protein_coding":
            continue
        if h.get("gene_biotype", "protein_coding") != "protein_coding":
            continue
        seq = seq.rstrip("*").upper()
        gene = h.get("gene", h["id"])
        if "*" in seq or len(seq) < 30:
            continue
        cur = best.get(gene)
        if cur is None or len(seq) > len(cur["seq"]) or (len(seq) == len(cur["seq"]) and h["id"] < cur["transcript"]):
            best[gene] = {"gene": gene, "transcript": h.get("transcript", h["id"]), "seq": seq,
                          "symbol": h.get("gene_symbol", ""), "seq_region": h.get("seq_region"),
                          "start": h.get("start"), "end": h.get("end"), "strand": h.get("strand")}
    return best


# ---------------------------------------------------------------- genetic code
BASES = "TCAG"
AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {a + b + c: AMINO[16 * i + 4 * j + k]
               for i, a in enumerate(BASES) for j, b in enumerate(BASES) for k, c in enumerate(BASES)}


def translate(cds):
    cds = cds.upper()
    return "".join(CODON_TABLE.get(cds[i:i + 3], "X") for i in range(0, len(cds) - 2, 3))


def revcomp(seq):
    return seq.upper().translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


# ---------------------------------------------------------------- files
def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def save_json(path, obj):
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=1, sort_keys=True, default=_jsonable)


def load_json(path):
    with open(path) as fh:
        return json.load(fh)


def _jsonable(x):
    try:
        import numpy as np
        if isinstance(x, np.integer):
            return int(x)
        if isinstance(x, np.floating):
            return None if np.isnan(x) else float(x)
        if isinstance(x, np.ndarray):
            return x.tolist()
    except ImportError:
        pass
    raise TypeError(f"not JSON serialisable: {type(x)}")


# ---------------------------------------------------------------- HTTP
_session = None


def session():
    global _session
    if _session is None:
        import requests
        _session = requests.Session()
        _session.headers["User-Agent"] = "SOG1-comparative-genomics (github.com/bymankirat)"
    return _session


def http(method, url, *, ok404=False, tries=7, **kw):
    import requests
    kw.setdefault("timeout", 180)
    for attempt in range(tries):
        try:
            r = session().request(method, url, **kw)
        except requests.RequestException as e:
            log(f"  network error ({e.__class__.__name__}), retrying {url}")
            time.sleep(2 ** attempt)
            continue
        if r.status_code == 404 and ok404:
            return None
        if r.status_code == 429 or r.status_code >= 500:
            wait = float(r.headers.get("Retry-After", 2 ** attempt))
            log(f"  HTTP {r.status_code}, waiting {wait:.0f}s: {url}")
            time.sleep(wait)
            continue
        if r.status_code >= 400:
            raise RuntimeError(f"HTTP {r.status_code} for {url}: {r.text[:300]}")
        return r
    raise RuntimeError(f"gave up after {tries} tries: {url}")


def ensembl_get(path, **params):
    r = http("GET", config.ENSEMBL_REST + path, params=params, headers={"Accept": "application/json"}, ok404=True)
    time.sleep(0.08)  # stay under Ensembl's 15 requests per second
    return None if r is None else r.json()


def ensembl_post(path, body):
    r = http("POST", config.ENSEMBL_REST + path, json=body,
             headers={"Accept": "application/json", "Content-Type": "application/json"})
    time.sleep(0.1)
    return r.json()


def download(url, dest):
    dest = Path(dest)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    tmp = dest.with_suffix(dest.suffix + ".part")
    r = http("GET", url, stream=True)
    with open(tmp, "wb") as fh:
        for chunk in r.iter_content(1 << 20):
            fh.write(chunk)
    tmp.rename(dest)
    return dest
