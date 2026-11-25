"""Validation on EGFR: predictions against a laboratory transformation assay of rare tumour variants.

    python benchmarks/validate_egfr.py [--bootstrap 2000]

Ikeuchi et al. (Oncogene 2022; doi:10.1038/s41388-022-02263-4; CC BY 4.0) took every EGFR variant
seen in at least two tumours in COSMIC v81 (357 variants, most of them rare) and measured whether
each one transforms cells:

- NIH-3T3 focus formation, scored 1 (no foci) to 4 (severe); the paper calls score 3-4 oncogenic;
- Ba/F3 growth without IL-3 (transformed or not).

Their Supplementary Table 4 gives both results per variant, with OncoKB's annotation. OncoKB's terms
do not allow sharing its content outside publications, so its annotation is used only for one
summary AUC and is not written out per variant. This script downloads the table, runs every variant
through the pipeline as a user would type it ("EGFR G719S", "EGFR p.Glu746_Ala750del") and compares
the calls with the assay.

This is a test outside what the model was built for, on purpose. The model separates disease-causing
germline variants from harmless ones; the EGFR variants that drive tumours are activating (gain of
function), and a variant that breaks the kinase is damaging but does not transform. The question is
how far a "damaging" call carries over to "activating".

Labels for the AUCs: 3T3 score 3-4 = transforming, score 1 = not transforming; score 2 (mild) is
shown but left out. None of the assay's variants is in the model's training set (checked below).
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import httpx
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import compare_tools as C

from missensexplorer import PredictOptions, model, predict
from missensexplorer.aa import to_three_letter
from missensexplorer.cache import default_cache_dir
from missensexplorer.model_params import MODEL
from missensexplorer.predict import _num, prediction_status

TABLE = (
    "https://static-content.springer.com/esm/art%3A10.1038%2Fs41388-022-02263-4/MediaObjects/"
    "41388_2022_2263_MOESM5_ESM.xlsx"
)
PAPER = "Ikeuchi et al., Oncogene 2022, doi:10.1038/s41388-022-02263-4 (Supplementary Table 4)"
OUT = Path(__file__).resolve().parent / "results"
TRAINING = OUT / "model_training_set.csv.gz"
KINASE = (712, 979)  # UniProt P00533 "Protein kinase" domain
SCORES = {  # column -> (name, direction: +1 higher = more damaging)
    "model": ("MissensExplorer v" + MODEL["version"], 1),
    "alphamissense": ("AlphaMissense", 1),
    "revel": ("REVEL", 1),
    "cadd_phred": ("CADD", 1),
    "popeve": ("popEVE", -1),
}
STATUS_ORDER = ["Damaging", "Uncertain", "Benign / silent", "Not assessed", "mixed", "error"]
SHORT = {"Damaging": "damaging", "Uncertain": "uncertain", "Benign / silent": "benign", "Not assessed": "not assessed"}

MISSENSE = re.compile(r"^([A-Z])(\d+)([A-Z])$")
STOP = re.compile(r"^([A-Z])(\d+)_STOP$")
INDEL = re.compile(r"^([A-Z])(\d+)(?:_([A-Z])(\d+))?(del[A-Z]*|ins[A-Z]+|-[A-Z]+)$")


def read_table(path: Path) -> pd.DataFrame:
    """Supplementary Table 4 as one row per variant (the WT control dropped)."""
    raw = pd.read_excel(path, header=2)
    t = raw.iloc[:, :9].copy()
    t.columns = ["variant", "ffa", "type", "baf3", "cosmic", "genie", "oncokb", "clinvar_paper", "uncommon"]
    t = t[t.variant.astype(str).str.strip().ne("WT")].copy()
    t["variant"] = t.variant.astype(str).str.strip()
    t["ffa"] = t.ffa.astype(int)
    t["baf3"] = t.baf3.astype(str).str.contains("○")
    t["uncommon"] = t.uncommon.astype(str).str.contains("○")
    t["oncokb"] = t.oncokb.where(t.oncokb.notna(), "").astype(str).str.strip()
    t["clinvar_paper"] = t.clinvar_paper.where(t.clinvar_paper.notna(), "").astype(str).str.strip()
    return t.reset_index(drop=True)


def three(letters: str) -> str:
    return "".join(to_three_letter(c) or "?" for c in letters)


def to_input(name: str) -> tuple[str, str | None]:
    """(kind, what a user would type) for one of the paper's variant names; None when not one variant.

    The paper writes deletion-insertions as "E746_A750-AP" and single-residue deletions as "E746delE".
    """
    if m := MISSENSE.match(name):
        return "missense", f"EGFR {name}"
    if m := STOP.match(name):
        return "nonsense", f"EGFR {m[1]}{m[2]}*"
    if m := INDEL.match(name):
        a, p, b, q, rest = m.groups()
        span = f"{three(a)}{p}" + (f"_{three(b)}{q}" if b else "")
        if rest.startswith("del"):
            change = "del"
        elif rest.startswith("ins"):
            change = "ins" + three(rest[3:])
        else:
            change = "delins" + three(rest[1:])
        return ("insertion" if change.startswith("ins") else "deletion-insertion" if "ins" in change else "deletion"), (
            f"EGFR p.{span}{change}"
        )
    return "compound", None  # several changes on one allele: not one input


def training_overlap(table: pd.DataFrame, client: httpx.Client) -> list[str]:
    """Assay variants that are also in the model's training set (EGFR rows, matched via ClinVar)."""
    train = pd.read_csv(TRAINING)
    ids = train.loc[train.gene == "EGFR", "clinvar_id"].astype(int).astype(str).tolist()
    if not ids:
        return []
    r = client.get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
        params={"db": "clinvar", "id": ",".join(ids), "retmode": "json"},
    )
    r.raise_for_status()
    res = r.json()["result"]
    names = set()
    for u in res["uids"]:
        if m := re.search(r"\(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)", res[u]["title"]):
            names.add(f"{C_ONE[m[1]]}{m[2]}{C_ONE[m[3]]}")
    return sorted(names & set(table.variant))


C_ONE = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H",
    "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
    "Tyr": "Y", "Val": "V",
}  # fmt: skip


def summarise(item: dict, rows: list[dict], error: str | None) -> dict:
    """One output row per assay variant; several nucleotide changes for one protein change are averaged."""
    out = dict(item)
    if not rows:
        out.update(status="error", error=error or "no row in EGFR", n_dna=0)
        return out
    statuses = {prediction_status(r.get("prediction")) for r in rows}
    out["status"] = statuses.pop() if len(statuses) == 1 else "mixed"
    out["prediction"] = " / ".join(sorted({str(r.get("prediction")) for r in rows}))
    out["n_dna"] = len(rows)
    out["hgvsc"] = ", ".join(str(r.get("hgvsc") or "").split(":")[-1] for r in rows)
    ps = [ms.p for r in rows if (ms := model.score(r)) is not None] if item["kind"] == "missense" else []
    out["model"] = float(np.mean(ps)) if ps else None
    for col in ("alphamissense", "revel", "cadd_phred", "popeve"):
        vals = [v for r in rows if (v := _num(r.get(col))) is not None]
        out[col] = float(np.mean(vals)) if vals else None
    afs = [v for r in rows if (v := _num(r.get("gnomad_af"))) is not None]
    out["gnomad_af"] = max(afs) if afs else None
    out["clinvar"] = "; ".join(sorted({str(r["clinvar"]) for r in rows if r.get("clinvar")}))
    out["notes"] = "; ".join(sorted({str(r["notes"]) for r in rows if r.get("notes")}))
    return out


def pick_rows(item: dict, rows: list[dict]) -> list[dict]:
    """The pipeline's rows for one assay variant: EGFR, and for a substitution the same protein change."""
    rows = [r for r in rows if r.get("gene") == "EGFR"]
    if item["kind"] == "missense":
        same = [r for r in rows if r.get("protein_change") == item["variant"]]
        return [r for r in same if "missense_variant" in str(r.get("consequence"))]
    if item["kind"] == "nonsense":
        return [r for r in rows if "stop_gained" in str(r.get("consequence"))]
    return rows


def auc_table(df: pd.DataFrame, label: pd.Series, rng: np.random.Generator, n_boot: int) -> list[dict]:
    """AUC of every score on the variants where `label` is 0/1 (NaN = left out), with 95% bootstrap intervals."""
    keep = label.notna() & df["model"].notna()
    out = []
    for col, (name, direction) in SCORES.items():
        k = keep & df[col].notna()
        y = label[k].astype(int).to_numpy()
        s = direction * df.loc[k, col].astype(float).to_numpy()
        ref = df.loc[k, "model"].astype(float).to_numpy()
        if len(set(y)) < 2:
            continue
        ci, diff = C.bootstrap(y, s, ref, rng, n_boot)
        out.append(
            {
                "score": name,
                "n": len(y),
                "positives": int(y.sum()),
                "auc": round(C.auc(y, s), 4),
                "ci95": ci,
                "minus_model_ci95": None if col == "model" else diff,
            }
        )
    return out


def calls_table(df: pd.DataFrame) -> dict:
    """Counts of the pipeline's call per 3T3 score."""
    tab: dict = {}
    for ffa in (4, 3, 2, 1):
        c = Counter(df.loc[df.ffa == ffa, "status"])
        tab[str(ffa)] = {s: c.get(s, 0) for s in STATUS_ORDER if c.get(s, 0) or s in SHORT}
    return tab


MOTIFS = {  # catalytic motifs of the EGFR kinase (residues of NP_005219 / P00533)
    "glycine-rich loop GSGAFG (719-724)": (719, 724),
    "catalytic loop HRDLAARN (835-842)": (835, 842),
    "DFG motif (855-857)": (855, 857),
}


def baf3_negative(ms: pd.DataFrame) -> dict:
    """The missense variants that did not transform Ba/F3: their 3T3 scores, calls and kinase motifs."""
    neg, pos = ms[~ms.baf3], ms[ms.baf3]
    in_motif = {name: sorted(neg.loc[neg.position.between(a, b), "variant"]) for name, (a, b) in MOTIFS.items()}
    return {
        "variants": len(neg),
        "3t3_score_1": int((neg.ffa == 1).sum()),
        "damaging": int((neg.status == "Damaging").sum()),
        "kinase_domain": int(neg.position.between(*KINASE).sum()),
        "median_model": round(float(neg.model.median()), 3),
        "median_model_transformed": round(float(pos.model.median()), 3),
        "in_motifs": in_motif,
    }


def rates(df: pd.DataFrame) -> dict:
    """Share called damaging among transforming (3-4) and among non-transforming (1) variants."""
    pos, neg = df[df.ffa >= 3], df[df.ffa == 1]

    def share(d: pd.DataFrame, s: str) -> list:
        return [int((d.status == s).sum()), len(d)]

    return {
        "damaging_if_transforming": share(pos, "Damaging"),
        "benign_if_transforming": share(pos, "Benign / silent"),
        "damaging_if_not_transforming": share(neg, "Damaging"),
        "benign_if_not_transforming": share(neg, "Benign / silent"),
    }


def fmt(v: float | None, nd: int = 2) -> str:
    return "" if v is None or v != v else f"{v:.{nd}f}"


def markdown(rep: dict, df: pd.DataFrame) -> str:
    L = [
        "# EGFR: predictions against a transformation assay",
        "",
        f"Data: {PAPER}, CC BY 4.0. Every EGFR variant seen in at least two tumours in COSMIC v81, tested for",
        "transformation of NIH-3T3 cells (focus formation: 1 = none, 2 = mild, 3 = moderate, 4 = severe) and of",
        "Ba/F3 cells (growth without IL-3). The table lists 356 variants and 39 with 3T3 score 2, where the paper's",
        "text gives 357 and 29; the table is used as it is. Produced by `benchmarks/validate_egfr.py` on "
        f"{str(rep['run'].get('started_utc', ''))[:10]} (Ensembl {rep['run'].get('ensembl_release', '?')}).",
        "",
        "The model scores how likely a missense change is to cause *germline disease*; the EGFR variants that drive",
        "tumours are *activating*. A variant that breaks the kinase is damaging but does not transform, so this",
        "tests how far a damaging call carries over to an activating one. None of these variants is in the model's",
        f"training set (overlap: {rep['training_overlap'] or 'none'}).",
        "",
        "## Missense variants",
        "",
        f"{rep['missense']['variants']} substitutions; {rep['missense']['scored']} received a model score. "
        "Transforming = 3T3 score 3-4; not transforming = score 1; score 2 (mild) is left out of the AUCs.",
        "",
        "| Score | Variants | Transforming | AUC (95% CI) | Difference from the model |",
        "|---|---:|---:|---|---|",
    ]
    for a in rep["missense"]["auc_3t3"]:
        d = a["minus_model_ci95"]
        L.append(
            f"| {a['score']} | {a['n']} | {a['positives']} | {a['auc']:.3f} ({a['ci95'][0]:.3f}-{a['ci95'][1]:.3f}) | "
            + ("" if d is None else f"{d[0]:+.3f} to {d[1]:+.3f}")
            + " |"
        )
    for key, title in (
        ("auc_3t3_kinase", "Inside the kinase domain (residues 712-979) only"),
        ("auc_3t3_outside", "Outside the kinase domain only"),
        ("auc_baf3", "Against the Ba/F3 assay (transformed or not), all scores 1-4"),
        ("auc_oncokb", "Against OncoKB (Oncogenic or Likely Oncogenic vs Likely Neutral)"),
    ):
        L += ["", f"{title}:", "", "| Score | Variants | Positives | AUC (95% CI) |", "|---|---:|---:|---|"]
        for a in rep["missense"][key]:
            ci = f"{a['auc']:.3f} ({a['ci95'][0]:.3f}-{a['ci95'][1]:.3f})"
            L.append(f"| {a['score']} | {a['n']} | {a['positives']} | {ci} |")
    r = rep["missense"]["rates"]

    def pct(x: list) -> str:
        return f"{x[0]} of {x[1]} ({100 * x[0] / x[1]:.0f}%)" if x[1] else "-"

    L += [
        "",
        "Calls of the pipeline per 3T3 score (one row per protein change; *mixed* = its nucleotide changes got",
        "different calls):",
        "",
        "| 3T3 score | " + " | ".join(STATUS_ORDER) + " |",
        "|---|" + "---:|" * len(STATUS_ORDER),
    ]
    for ffa, c in rep["missense"]["calls"].items():
        L.append(f"| {ffa} | " + " | ".join(str(c.get(s, 0)) for s in STATUS_ORDER) + " |")
    L += [
        "",
        f"- Transforming (3-4) called damaging: {pct(r['damaging_if_transforming'])}; called benign: "
        f"{pct(r['benign_if_transforming'])}.",
        f"- Not transforming (1) called damaging: {pct(r['damaging_if_not_transforming'])}; called benign: "
        f"{pct(r['benign_if_not_transforming'])}.",
        "",
    ]
    b = rep["missense"]["baf3_negative"]
    motifs = "; ".join(f"{k}: {', '.join(v)}" for k, v in b["in_motifs"].items() if v)
    L += [
        "## Reading",
        "",
        "- **Activation is not what the scores measure.** Against 3T3 transformation every score is close to chance",
        "  (AUC about 0.5 to 0.55). A damaging call covers most transforming variants, but also most of the",
        "  variants that do not transform.",
        f"- **Loss of function is.** The {b['variants']} substitutions that did not transform Ba/F3 cells all also "
        f"failed in 3T3 ({b['3t3_score_1']} of {b['variants']} with score 1), {b['kinase_domain']} lie in the kinase "
        f"domain, and {b['damaging']} are called damaging (median model score {b['median_model']:.2f}, against "
        f"{b['median_model_transformed']:.2f} for the variants that transformed Ba/F3). Several sit on the kinase's "
        f"catalytic motifs ({motifs}): changes that break the kinase. Hence the Ba/F3 AUCs below 0.5: there, "
        '"does not transform" mostly means "broken", which the predictors detect.',
        "- So, for a somatic EGFR variant, a damaging call says the protein is changed where it matters, not that it",
        "  is activating; whether it drives a tumour needs curated or functional data (see also `--somatic`).",
        f"  A benign call does not rule activation out either: {pct(r['benign_if_transforming'])} of the "
        "transforming variants get one, "
        f"against {pct(r['benign_if_not_transforming'])} of the others.",
        "",
        "## The uncommon mutations of lung cancer",
        "",
        "The variants the paper's table marks as *uncommon EGFR mutations* (the column also marks L858R and T790M).",
        "3T3: focus formation score; Ba/F3: transformed (yes/no).",
        "",
    ]
    L += variant_table(df[df.uncommon])
    L += [
        "",
        "## In-frame deletions and insertions",
        "",
        "Typed in protein notation; the pipeline maps a protein deletion to one DNA change and cannot map a protein",
        "insertion (many DNA sequences give it). In-frame changes get no model score: they are called *likely",
        "damaging* when CADD is at least 20.",
        "",
        f"{rep['indels']['variants']} variants; {rep['indels']['annotated']} annotated, {rep['indels']['errors']} not "
        f"mapped ({rep['indels']['error_kinds']}). Ensembl's CADD scored {rep['indels']['with_cadd']} of the "
        f"{rep['indels']['annotated']} annotated deletions, so they stay *protein-altering in-frame* (uncertain).",
        "",
        "| 3T3 score | " + " | ".join(STATUS_ORDER) + " |",
        "|---|" + "---:|" * len(STATUS_ORDER),
    ]
    for ffa, c in rep["indels"]["calls"].items():
        L.append(f"| {ffa} | " + " | ".join(str(c.get(s, 0)) for s in STATUS_ORDER) + " |")
    L += [
        "",
        f"{rep['compound']} variants with several changes on one allele (such as L858R_T790M) are not single inputs",
        "and are left out.",
        "",
        "## Every missense variant",
        "",
        "<details><summary>Table</summary>",
        "",
    ]
    L += variant_table(df[df.kind.isin(["missense", "nonsense"])])
    L += ["", "</details>", ""]
    return "\n".join(L)


def variant_table(d: pd.DataFrame) -> list[str]:
    L = [
        "| Variant | COSMIC | 3T3 | Ba/F3 | DNA change | Call | Model | AlphaMissense | REVEL | CADD | popEVE |",
        "|---|---:|---:|---|---|---|---:|---:|---:|---:|---:|",
    ]
    for _, x in d.sort_values("position").iterrows():
        call = x.prediction if isinstance(x.get("prediction"), str) else x.status
        if x.status == "error":
            short = "insertion in protein notation" if "insertion written as a protein" in str(x.error) else x.error
            call = f"not annotated ({short})"
        dna = str(x.hgvsc) if isinstance(x.get("hgvsc"), str) else ""
        L.append(
            f"| {x.variant} | {int(x.cosmic) if x.cosmic == x.cosmic else ''} | {x.ffa} | "
            f"{'yes' if x.baf3 else 'no'} | "
            f"{dna} | {call} | {fmt(x.get('model'))} | {fmt(x.get('alphamissense'))} | "
            f"{fmt(x.get('revel'))} | "
            f"{fmt(x.get('cadd_phred'), 1)} | {fmt(x.get('popeve'))} |"
        )
    return L


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bootstrap", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args(argv)

    d = default_cache_dir() / "validation"
    d.mkdir(parents=True, exist_ok=True)
    path = d / "egfr_ikeuchi2022_table4.xlsx"
    with httpx.Client(timeout=60, follow_redirects=True, headers={"User-Agent": "missensexplorer-benchmark"}) as client:
        if not path.exists():
            r = client.get(TABLE)
            r.raise_for_status()
            path.write_bytes(r.content)
        table = read_table(path)
        overlap = training_overlap(table, client)

    items = []
    for x in table.to_dict("records"):
        kind, text = to_input(x["variant"])
        pos = int(re.search(r"\d+", x["variant"])[0])
        items.append({**x, "kind": kind, "input": text, "position": pos})
    inputs = sorted({i["input"] for i in items if i["input"]})
    print(f"{len(table)} assay variants; {len(inputs)} inputs", flush=True)

    res = predict(
        text="\n".join(inputs),
        options=PredictOptions(protein_context=False, sequences=False, gnomad_check=False, workers=6),
        progress=lambda s, k, n: print(f"  {s}: {k}/{n}", end="\r", flush=True),
    )
    print()
    failed = {k: v for k, v in res.sources.items() if str(v).startswith("failed")}
    if failed:
        raise SystemExit(f"annotation incomplete: {failed}. Run again (cached answers are kept).")
    by = defaultdict(list)
    for r in res.table.to_dict("records"):
        by[r["input"]].append(r)
    errors = dict(zip(res.errors["input"], res.errors["error"], strict=True)) if len(res.errors) else {}

    rows = [summarise(i, pick_rows(i, by.get(i["input"], [])), errors.get(i["input"])) for i in items if i["input"]]
    df = pd.DataFrame(rows)
    for i in items:
        if i["kind"] == "compound":
            df = pd.concat([df, pd.DataFrame([{**i, "status": "compound"}])], ignore_index=True)
    df.loc[df.status.isna(), "status"] = "error"

    rng = np.random.default_rng(args.seed)
    ms = df[df.kind == "missense"].reset_index(drop=True)
    y3t3 = ms.ffa.map(lambda f: 1.0 if f >= 3 else 0.0 if f == 1 else np.nan)
    kin = ms.position.between(*KINASE)
    ybaf3 = ms.baf3.astype(float)
    yonc = ms.oncokb.map({"Oncogenic": 1.0, "Likely Oncogenic": 1.0, "Likely Neutral": 0.0})
    indels = df[df.kind.isin(["deletion", "insertion", "deletion-insertion"])]
    err_kinds = Counter(str(e) for e in indels.loc[indels.status == "error", "error"])
    rep = {
        "paper": PAPER,
        "table": TABLE,
        "training_overlap": overlap,
        "missense": {
            "variants": len(ms),
            "scored": int(ms.model.notna().sum()),
            "auc_3t3": auc_table(ms, y3t3, rng, args.bootstrap),
            "auc_3t3_kinase": auc_table(ms, y3t3.where(kin), rng, args.bootstrap),
            "auc_3t3_outside": auc_table(ms, y3t3.where(~kin), rng, args.bootstrap),
            "auc_baf3": auc_table(ms, ybaf3, rng, args.bootstrap),
            "auc_oncokb": auc_table(ms, yonc, rng, args.bootstrap),
            "calls": calls_table(ms),
            "rates": rates(ms),
            "baf3_negative": baf3_negative(ms),
        },
        "indels": {
            "variants": len(indels),
            "annotated": int((indels.status != "error").sum()),
            "errors": int((indels.status == "error").sum()),
            "with_cadd": int(indels.cadd_phred.notna().sum()),
            "error_kinds": "; ".join(f"{n} × {e}" for e, n in err_kinds.most_common()),
            "calls": calls_table(indels),
        },
        "compound": int((df.kind == "compound").sum()),
        "run": res.run,
    }
    OUT.mkdir(exist_ok=True)
    cols = [
        "variant", "kind", "input", "position", "cosmic", "genie", "ffa", "baf3", "clinvar_paper", "uncommon",
        "status", "prediction", "n_dna", "hgvsc", "model", "alphamissense", "revel", "cadd_phred", "popeve",
        "gnomad_af", "clinvar", "notes", "error",
    ]  # fmt: skip
    df.sort_values("position").reindex(columns=cols).to_csv(OUT / "egfr_validation.csv", index=False)
    (OUT / "egfr_validation.json").write_text(json.dumps(rep, indent=1, default=str))
    (OUT / "egfr_validation.md").write_text(markdown(rep, df))
    print((OUT / "egfr_validation.md").read_text().split("## The uncommon")[0])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
