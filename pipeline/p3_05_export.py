"""Phase 3, step 5: write docs/phase3_results.js so the report page can be opened straight from GitHub Pages."""
import json

import numpy as np
import pandas as pd

from util import ROOT, RESULTS, load_json, log

OUT = RESULTS / "phase3"
DOCS = ROOT / "docs"


def clean(obj):
    """NaN and numpy types are not valid JSON; turn them into null and plain Python."""
    if isinstance(obj, dict):
        return {k: clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [clean(v) for v in obj]
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating, float)):
        return None if (obj is None or np.isnan(obj)) else float(obj)
    if obj is pd.NA or obj is None:
        return None
    return obj


def table(path, **kw):
    if not path.exists():
        return []
    return clean(pd.read_csv(path, **kw).to_dict("records"))


def main():
    DOCS.mkdir(exist_ok=True)
    contrasts = pd.read_csv(OUT / "stress_contrasts.csv") if (OUT / "stress_contrasts.csv").exists() else pd.DataFrame()
    genes = pd.read_csv(OUT / "stress_genes.csv") if (OUT / "stress_genes.csv").exists() else pd.DataFrame()
    data = {
        "generated": pd.Timestamp.utcnow().strftime("%Y-%m-%d %H:%M UTC"),
        "experiments": table(OUT / "atlas_experiments.csv"),
        "contrasts": clean(contrasts.to_dict("records")) if len(contrasts) else [],
        "genes": clean(genes.to_dict("records")) if len(genes) else [],
        "summary": load_json(OUT / "stress_summary.json") if (OUT / "stress_summary.json").exists() else {},
        "motif": load_json(OUT / "motif_summary.json") if (OUT / "motif_summary.json").exists() else [],
        "motif_targets": table(OUT / "motif_targets.csv"),
        "up_genes": table(OUT / "stress_up_genes.csv"),
        "orthologs": table(OUT / "orthologs.csv"),
    }
    js = "window.P3 = " + json.dumps(data, separators=(",", ":")) + ";\n"
    (DOCS / "phase3_results.js").write_text(js)
    log(f"wrote docs/phase3_results.js ({len(js) / 1e6:.2f} MB): "
        f"{len(data['contrasts'])} contrasts, {len(data['genes'])} gene rows")


if __name__ == "__main__":
    main()
