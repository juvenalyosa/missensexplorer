"""Would more information improve the missense model? Gene-split cross-validation on its training set.

    python benchmarks/feature_tests.py

Adds one candidate input at a time to the model's four predictors and measures the change in AUC on
rare training variants (5-fold cross-validation split by gene, three seeds; 95% paired bootstrap
interval of the gain). Candidates:
  * gene constraint from gnomAD v4.1 (missense o/e upper bound, missense z, LOEUF);
  * AlphaFold confidence at the residue (pLDDT) and UniProt functional sites;
  * predictors without clinical training labels from dbNSFP: ESM1b, MutationAssessor, phyloP,
    PROVEAN, SIFT4G, PrimateAI.
A candidate was to be adopted only with a gain of at least 0.003 whose interval excludes zero and a
coefficient of the expected sign. None qualified (see docs/model.md): LOEUF reached +0.003 but with
the wrong sign (it rewards LoF-tolerant genes, i.e. it learns which genes ClinVar's clinical testing
covers, and costs sensitivity in constrained genes); PrimateAI also took the wrong sign.
"""

from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

import httpx
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import clinvar_benchmark as B
import compare_tools as C
import train_model as T

from missensexplorer import PredictOptions, predict
from missensexplorer.cache import default_cache_dir
from missensexplorer.predict import _num

CONSTRAINT = (
    "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv"
)
GENE = {"moeuf": "mis.oe_ci.upper", "mis_z": "mis.z_score", "loeuf": "lof.oe_ci.upper"}
DBNSFP = {
    "esm1b": "ESM1b",
    "mutation_assessor": "MutationAssessor",
    "phylop": "phyloP (100 vertebrates)",
    "provean": "PROVEAN",
    "sift4g": "SIFT4G",
    "primateai": "PrimateAI",
}


def constraint() -> pd.DataFrame:
    f = default_cache_dir() / "validation" / "gnomad.v4.1.constraint_metrics.tsv"
    if not f.exists():
        f.parent.mkdir(parents=True, exist_ok=True)
        print("Downloading gnomAD v4.1 constraint metrics (95 MB) ...", flush=True)
        with httpx.stream("GET", CONSTRAINT, timeout=600) as r, open(f.with_suffix(".part"), "wb") as fh:
            for chunk in r.iter_bytes(1 << 20):
                fh.write(chunk)
        f.with_suffix(".part").rename(f)
    c = pd.read_csv(f, sep="\t", low_memory=False)
    return c[c["mane_select"].astype(str).str.lower() == "true"].drop_duplicates("gene").set_index("gene")


def table() -> pd.DataFrame:
    train = pd.read_csv(T.OUT / "model_training_set.csv.gz", dtype={"clinvar_id": str})
    cands = {c["id"]: c for c in B.load_candidates(B.download())}
    picked = [cands[i] for i in train["clinvar_id"] if i in cands]
    res = predict(
        text="\n".join(p["input"] for p in picked),
        options=PredictOptions(clinvar=False, protein_context=True, sequences=False, workers=4),
        progress=lambda s, k, n: print(f"  {s}: {k}/{n}", end="\r", flush=True),
    )
    print()
    by = defaultdict(list)
    for r in res.table.to_dict("records"):
        by[r["input"]].append(r)
    rows = []
    for p in picked:
        cand = [r for r in by.get(p["input"], []) if "missense_variant" in str(r["consequence"])]
        r = ([r for r in cand if r["gene"] == p["gene"]] or cand or [None])[0]
        if r is not None:
            rows.append({**p, "row": r})
    ext = C.fetch_dbnsfp(rows)
    con = constraint()
    recs = []
    for x in rows:
        r, g = x["row"], (con.loc[x["gene"]] if x["gene"] in con.index else None)
        ks = r.get("key_site")
        recs.append(
            {
                "gene": x["gene"],
                "label": x["label"],
                "grpmax": _num(r.get("gnomad_popmax_af")) or 0.0,
                **{f: _num(r.get(f)) for f in T.FEATURES},
                "plddt": _num(r.get("plddt")),
                "key_site": 1.0 if isinstance(ks, str) and ks else 0.0,
                **{k: ext[x["id"]].get(v) for k, v in DBNSFP.items()},
                **{k: (None if g is None else _num(g[c])) for k, c in GENE.items()},
            }
        )
    df = pd.DataFrame(recs)
    for c in df.columns.drop("gene"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def main() -> int:
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import GroupKFold

    df = table()
    y, rare = df["label"].values, (df["grpmax"] < T.RARE).values
    w = T.class_weights(y)
    print("coverage:", df.drop(columns=["gene", "label", "grpmax"]).notna().mean().round(3).to_dict())

    def oof(feats, seed):
        p = np.zeros(len(df))
        for tr, te in GroupKFold(n_splits=5, shuffle=True, random_state=seed).split(df, y, df["gene"]):
            m = T.pipeline().fit(df[feats].iloc[tr], y[tr], logisticregression__sample_weight=w[tr])
            p[te] = m.predict_proba(df[feats].iloc[te])[:, 1]
        return p

    base = list(T.FEATURES)
    rng = np.random.default_rng(0)
    seeds = (7, 11, 23)
    pb = {s: oof(base, s) for s in seeds}
    print(f"base: rare AUC {np.mean([roc_auc_score(y[rare], pb[s][rare]) for s in seeds]):.4f}")
    for cand in [*list(GENE), "plddt", "key_site", *list(DBNSFP)]:
        feats = [*base, cand]
        pc = {s: oof(feats, s) for s in seeds}
        gains = [roc_auc_score(y[rare], pc[s][rare]) - roc_auc_score(y[rare], pb[s][rare]) for s in seeds]
        a, b, yr = pc[seeds[0]][rare], pb[seeds[0]][rare], y[rare]
        bs = []
        for _ in range(500):
            i = rng.integers(0, len(yr), len(yr))
            bs.append(roc_auc_score(yr[i], a[i]) - roc_auc_score(yr[i], b[i]))
        coef = T.pipeline().fit(df[feats], y, logisticregression__sample_weight=w)[-1].coef_[0][-1]
        print(
            f"+ {cand:18s} gain {np.mean(gains):+.4f} [{np.percentile(bs, 2.5):+.4f}, {np.percentile(bs, 97.5):+.4f}]"
            f"  coefficient {coef:+.3f}",
            flush=True,
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
