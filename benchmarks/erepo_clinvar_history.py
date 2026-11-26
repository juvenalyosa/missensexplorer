"""PS1/PM5 against the expert panels with the ClinVar that existed when each panel decided.

    python benchmarks/erepo_clinvar_history.py

validate_erepo.py takes PS1/PM5 from today's ClinVar, so a reference variant classified after a panel's
decision can count. Here the tool's rule is applied with the references that existed then:

- ClinVar as it was: for each classification, the last monthly ClinVar release (archived variant_summary
  files) dated before the month of the panel's approval (the lag is reported; about four to nine weeks);
- plus the expert panels' own pathogenic and likely pathogenic classifications (ClinGen Evidence
  Repository) approved on or before that day, which a panel could use before ClinVar showed them (each
  variant's latest classification only, so earlier ones later changed are missing); a sensitivity
  analysis counts only those approved strictly before.

The rule is residue.py's: reference variants are germline Pathogenic, Pathogenic/Likely pathogenic or
Likely pathogenic missense variants with review criteria (one star or more; older releases join several
assertions with ";", read the same way when they are all pathogenic or likely pathogenic), on the same
transcript (the RefSeq accession of the panel's variant name), at the same protein position and reference
amino acid, the variant itself excluded (by ClinVar id, or by c. change in releases without one). A
reference counts at full strength when pathogenic with two or more review stars, one step lower when
pathogenic with one star or likely pathogenic with two or more, and not at all when likely pathogenic with
one star. The same amino-acid change gives PS1 (PS1_Moderate one step lower), another change PM5
(PM5_Supporting); the stronger is kept; none when SpliceAI >= 0.2. The same code on the current release
("now") checks the re-implementation against the tool's own calls ("tool").

Known limits, reported: references are matched on the panel's RefSeq accession, so a reference disappears
from a release in which ClinVar named it on another transcript (flagged when a transcript's references drop
below half of their earlier maximum; results are also given without the flagged classifications); before
2024, ClinVar's ClinicalSignificance could mix germline and somatic assertions.

Each archived file is streamed once and only the pathogenic missense variants of the panels' transcripts
are kept (benchmark/clinvar_history_v2/, a few MB per release); the archives themselves are not stored.
Reads benchmarks/results/erepo_validation_variants.csv.gz; writes benchmarks/results/erepo_clinvar_history.
{json,md} and the per-classification calls (erepo_clinvar_history_variants.csv.gz).
"""

# ruff: noqa: E501  (report lines are kept whole)
from __future__ import annotations

import gzip
import hashlib
import io
import json
import re
import sys
from collections import Counter
from pathlib import Path

import httpx
import numpy as np
import pandas as pd
from scipy.stats import binomtest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import clinvar_benchmark as B  # noqa: E402

from missensexplorer.cache import default_cache_dir  # noqa: E402

OUT = Path(__file__).resolve().parent / "results"
ARCHIVE = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/"
NAME = re.compile(r"^(N[MR]_\d+)(?:\.\d+)?\([^)]*\):(c\.\S+) \(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)$")
LEVELS = {"pathogenic": "P", "likely pathogenic": "LP", "pathogenic/likely pathogenic": "P"}
NEUTRAL = {"not provided"}  # may accompany pathogenic assertions in releases before 2017
STARS = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, single submitter": 1,
}
SPLICE_CAUTION = 0.2
FIRST = (2016, 1)
VERSION = "2"  # of the extract format; part of the cache file names
B_REPS = 10000


def level(sig: str) -> str | None:
    """'P', 'LP' or None for a ClinicalSignificance value; 'Likely pathogenic;Pathogenic' counts as P."""
    parts = {p.strip().lower() for p in re.split(r"[;,]", sig) if p.strip()}
    parts -= NEUTRAL
    if not parts or not parts <= set(LEVELS):
        return None
    return "P" if any(LEVELS[p] == "P" for p in parts) else "LP"


def release_urls(year: int, month: int) -> list[str]:
    name = f"variant_summary_{year}-{month:02d}.txt.gz"
    return [f"{ARCHIVE}{name}", f"{ARCHIVE}{year}/{name}"]  # recent years at the top, older ones by year


def rank(lv: str, stars: int) -> int:
    """residue.rank: 2 full strength, 1 one step lower, 0 not used."""
    if lv == "P" and stars >= 2:
        return 2
    return 1 if lv == "P" or stars >= 2 else 0


def extract(lines, transcripts: set[str]) -> pd.DataFrame:
    """Pathogenic/likely pathogenic reviewed missense variants of `transcripts` from variant_summary lines."""
    header = next(lines).lstrip("#").rstrip("\n").split("\t")
    ix = {k: i for i, k in enumerate(header)}
    vid_col = "VariationID" if "VariationID" in ix else None
    rows, seen = [], set()
    for line in lines:
        if "athogenic" not in line or "(p." not in line:  # fast skip of most lines
            continue
        f = line.rstrip("\n").split("\t")
        lv = level(f[ix["ClinicalSignificance"]])
        stars = STARS.get(f[ix["ReviewStatus"]], 0)
        if not lv or not stars:
            continue
        m = NAME.match(f[ix["Name"]])
        if not m or m.group(1) not in transcripts or m.group(5) in ("Ter", m.group(3)):
            continue
        key = (m.group(1), m.group(2))
        if key in seen:
            continue
        seen.add(key)
        rows.append({"nm": m.group(1), "pos": int(m.group(4)), "ref": m.group(3), "alt": m.group(5),
                     "vid": f[ix[vid_col]] if vid_col else "", "level": lv, "stars": stars, "cdna": m.group(2)})  # fmt: skip
    return pd.DataFrame(rows, columns=["nm", "pos", "ref", "alt", "vid", "level", "stars", "cdna"])


def snapshot(year: int, month: int, transcripts: set[str], http: httpx.Client) -> pd.DataFrame | None:
    """The extract of one archived release (cached as a small CSV), or None when it is not archived."""
    tag = hashlib.sha1((VERSION + ",".join(sorted(transcripts))).encode()).hexdigest()[:10]
    d = default_cache_dir() / "benchmark" / "clinvar_history_v2"
    d.mkdir(parents=True, exist_ok=True)
    f = d / f"{year}-{month:02d}_{tag}.csv.gz"
    if f.exists():
        return pd.read_csv(f, dtype={"vid": str, "cdna": str})
    print(f"  ClinVar {year}-{month:02d} ...", flush=True)
    for url in release_urls(year, month):
        with http.stream("GET", url, timeout=1800) as r:
            if r.status_code == 404:
                continue
            r.raise_for_status()
            raw = io.BytesIO()
            for chunk in r.iter_bytes(1 << 20):
                raw.write(chunk)
        raw.seek(0)
        with gzip.open(raw, "rt", encoding="utf-8", errors="replace") as fh:
            df = extract(iter(fh), transcripts)
        df.to_csv(f, index=False)
        return df
    return None


def current(transcripts: set[str]) -> pd.DataFrame:
    with gzip.open(B.download(), "rt", encoding="utf-8", errors="replace") as fh:
        return extract(iter(fh), transcripts)


Ref = tuple  # (alt, vid, cdna, level, stars)


def index(refs: pd.DataFrame | None) -> dict[tuple[str, int, str], list[Ref]]:
    """Reference variants by (transcript, protein position, reference amino acid), for fast lookups."""
    out: dict[tuple[str, int, str], list[Ref]] = {}
    if refs is None:
        return out
    for r in refs.itertuples(index=False):
        out.setdefault((r.nm, int(r.pos), r.ref), []).append((r.alt, str(r.vid), str(r.cdna), r.level, int(r.stars)))
    return out


def call(v: dict, *indexes: dict[tuple[str, int, str], list[Ref]]) -> str | None:
    """PS1/PM5 code for one eRepo variant from indexed reference variants (residue.py's rule)."""
    if v["spliceai"] is not None and v["spliceai"] >= SPLICE_CAUTION:
        return None
    key = (v["nm"], v["pos"], v["ref"])
    refs = [x for ix in indexes for x in ix.get(key, [])]
    refs = [x for x in refs if x[1] != v["vid"] and x[2] != v["cdna"]]  # never the variant itself
    usable = [(alt, rank(lv, st)) for alt, _, _, lv, st in refs if rank(lv, st) > 0]
    if not usable:
        return None
    same = [r for a, r in usable if a == v["alt"]]
    other = [r for a, r in usable if a != v["alt"]]
    options = []
    if same:
        options.append(("PS1", 4) if max(same) == 2 else ("PS1_Moderate", 2))
    if other:
        options.append(("PM5", 2) if max(other) == 2 else ("PM5_Supporting", 1))
    return max(options, key=lambda o: o[1])[0]  # PS1 first on a tie


def base(c: str | None) -> str | None:
    return None if c is None else c.split("_")[0]


def agreement(d: pd.DataFrame, key: str, code: str) -> dict:
    ev = d[d[f"exp_{code}"].notna()]
    e, g = ev[f"exp_{code}"].astype(bool), ev[key] == code
    tp, fn, fp, tn = (int(x.sum()) for x in (e & g, e & ~g, ~e & g, ~e & ~g))
    pct = lambda a, b: round(a / b, 6) if b else None  # noqa: E731
    return {"tp": tp, "fn": fn, "fp": fp, "tn": tn, "positive_agreement": pct(tp, tp + fn),
            "negative_agreement": pct(tn, tn + fp), "ppv": pct(tp, tp + fp)}  # fmt: skip


def panel_bootstrap(d: pd.DataFrame, stat) -> list[float]:
    """95% percentile interval of stat over panels resampled with replacement (sums of per-panel counts)."""
    groups = [g for _, g in d.groupby("panel")]
    counts = np.array([stat(g) for g in groups], float)  # per panel: (numerators..., denominators...)
    picks = np.random.default_rng(1).integers(0, len(groups), (B_REPS, len(groups)))
    tot = counts[picks].sum(1)
    k = counts.shape[1] // 2
    with np.errstate(invalid="ignore", divide="ignore"):
        v = tot[:, :k] / tot[:, k:]
    v = v[:, 0] - v[:, 1] if k == 2 else v[:, 0]
    v = v[np.isfinite(v)]
    return [round(float(x), 4) for x in np.percentile(v, [2.5, 97.5])]


def main() -> int:
    df = pd.read_csv(OUT / "erepo_validation_variants.csv.gz", dtype={"clinvar_id": str})
    if "spliceai_max" not in df.columns:
        raise SystemExit(
            "this needs the full eRepo table with SpliceAI scores (run validate_erepo.py), not the public one"
        )
    variants = []
    for r in df.to_dict("records"):
        m = NAME.match(str(r["name"]))
        if not m:
            continue
        sp = r.get("spliceai_max")
        tool = str(r["tool_ps1_pm5"]).split(" ")[0]
        variants.append({
            "clinvar_id": r["clinvar_id"], "vid": r["clinvar_id"], "panel": r["panel"], "nm": m.group(1),
            "pos": int(m.group(4)), "ref": m.group(3), "alt": m.group(5), "cdna": m.group(2),
            "approved": r["approved"], "assertion": r["assertion"], "spliceai": None if pd.isna(sp) else float(sp),
            "tool": tool if tool.startswith(("PS1", "PM5")) else None,
            "exp_PS1": None if pd.isna(r["exp_PS1"]) else bool(r["exp_PS1"] > 0),
            "exp_PM5": None if pd.isna(r["exp_PM5"]) else bool(r["exp_PM5"] > 0),
        })  # fmt: skip
    transcripts = {v["nm"] for v in variants}
    print(f"{len(variants)} eRepo variants on {len(transcripts)} transcripts", flush=True)
    # the panels' own pathogenic/likely pathogenic classifications, usable from their approval day
    panel_refs = pd.DataFrame([
        {"nm": v["nm"], "pos": v["pos"], "ref": v["ref"], "alt": v["alt"], "vid": v["vid"], "cdna": v["cdna"],
         "level": "P" if v["assertion"] == "Pathogenic" else "LP", "stars": 3, "approved": v["approved"]}
        for v in variants if v["assertion"] in ("Pathogenic", "Likely Pathogenic")
    ])  # fmt: skip
    panel_index: dict[str, dict[tuple[str, int, str], list[tuple]]] = {}
    for r in panel_refs.itertuples(index=False):
        panel_index.setdefault(r.nm, {}).setdefault((r.nm, int(r.pos), r.ref), []).append(
            (r.alt, str(r.vid), str(r.cdna), r.level, int(r.stars), r.approved)
        )

    now = index(current(transcripts))
    for v in variants:
        v["now"] = call(v, now)
    same = sum(1 for v in variants if (v["now"] or "") == (v["tool"] or ""))
    applied = [v for v in variants if v["now"] or v["tool"]]
    same_applied = sum(1 for v in applied if (v["now"] or "") == (v["tool"] or ""))
    print(f"re-implementation on the current release agrees with the tool's calls for {same} of {len(variants)} "
          f"({same_applied} of {len(applied)} classifications where either applies PS1/PM5)", flush=True)  # fmt: skip

    http = httpx.Client(follow_redirects=True, headers={"User-Agent": "missensexplorer-benchmark"})
    releases: dict[tuple[int, int], dict | None] = {}
    months = [(y, m) for y in range(FIRST[0], 2027) for m in range(1, 13)]
    for v in variants:
        y, mo = int(v["approved"][:4]), int(v["approved"][5:7])
        prior = [x for x in months if x < (y, mo)]
        rel, refs = None, None
        while prior:
            cand = prior.pop()
            if cand < FIRST:
                break
            if cand not in releases:
                snap = snapshot(*cand, transcripts, http)
                releases[cand] = None if snap is None else {"index": index(snap), "counts": Counter(snap["nm"])}
            if releases[cand] is not None:
                rel, refs = cand, releases[cand]
                break
        v["release"] = f"{rel[0]}-{rel[1]:02d}" if rel else None
        clinvar = refs["index"] if refs is not None else {}
        v["then_clinvar"] = call(v, clinvar) if refs is not None else None
        own = {k: [x[:5] for x in xs if x[5] <= v["approved"]] for k, xs in panel_index.get(v["nm"], {}).items()}
        v["then"] = call(v, clinvar, own)
        # sensitivity: only the panels' classifications approved strictly before the decision
        earlier = {k: [x[:5] for x in xs if x[5] < v["approved"]] for k, xs in panel_index.get(v["nm"], {}).items()}
        v["then_strict"] = call(v, clinvar, earlier)
        v["lag_days"] = (pd.Timestamp(v["approved"]) - pd.Timestamp(f"{rel[0]}-{rel[1]:02d}-01")).days if rel else None
    http.close()

    # transcripts whose references drop below half of an earlier maximum: ClinVar named them on another one
    counts = {k: r["counts"] for k, r in releases.items() if r is not None}
    flagged = set()
    for nm in transcripts:
        best = 0
        for k in sorted(counts):
            c = counts[k].get(nm, 0)
            if best >= 10 and c < best / 2:
                flagged.add((nm, f"{k[0]}-{k[1]:02d}"))
            best = max(best, c)
    for v in variants:
        v["transcript_switch"] = (v["nm"], v["release"]) in flagged
        v["now_base"], v["then_base"], v["tool_base"] = base(v["now"]), base(v["then"]), base(v["tool"])
        v["then_clinvar_base"] = base(v["then_clinvar"])
        v["then_strict_base"] = base(v["then_strict"])
    d = pd.DataFrame(variants)

    def pm5_counts(g: pd.DataFrame, a: str, b: str) -> tuple:
        e = g["exp_PM5"] == True  # noqa: E712
        return (int((e & (g[a] == "PM5")).sum()), int((e & (g[b] == "PM5")).sum()), int(e.sum()), int(e.sum()))

    def ppv_counts(g: pd.DataFrame, a: str, b: str) -> tuple:
        ev = g[g["exp_PM5"].notna()]
        e = ev["exp_PM5"] == True  # noqa: E712
        return (int((e & (ev[a] == "PM5")).sum()), int((e & (ev[b] == "PM5")).sum()),
                int((ev[a] == "PM5").sum()), int((ev[b] == "PM5").sum()))  # fmt: skip

    ev = d[d["exp_PM5"].notna()]
    met = ev[ev["exp_PM5"] == True]  # noqa: E712
    lost = int(((met["now_base"] == "PM5") & (met["then_base"] != "PM5")).sum())
    gained = int(((met["now_base"] != "PM5") & (met["then_base"] == "PM5")).sum())
    disc_now = ev[(ev["exp_PM5"] == False) & (ev["now_base"] == "PM5")]  # noqa: E712
    conc_now = met[met["now_base"] == "PM5"]
    disc_removed = int((disc_now["then_base"] != "PM5").sum())
    conc_removed = int((conc_now["then_base"] != "PM5").sum())
    tool_disc = ev[(ev["exp_PM5"] == False) & (ev["tool_base"] == "PM5")]  # noqa: E712
    attribution = {
        "tool_pm5_against_panel": len(tool_disc),
        "gross_no_pm5_reference_then": int((tool_disc["then_base"] != "PM5").sum()),
        "gross_share": round(float((tool_disc["then_base"] != "PM5").mean()), 4),
        "discordant_now": len(disc_now),
        "discordant_removed_then": disc_removed,
        "concordant_now": len(conc_now),
        "concordant_removed_then": conc_removed,
        "excess_share": round((disc_removed - conc_removed / len(conc_now) * len(disc_now)) / len(disc_now), 4),
        "clinvar_only_discordant_removed": int((disc_now["then_clinvar_base"] != "PM5").sum()),
        "clinvar_only_concordant_removed": int((conc_now["then_clinvar_base"] != "PM5").sum()),
    }

    # shares on the re-implementation's basis (now -> then), with intervals over panels and by period
    def shares(g: pd.DataFrame, then: str = "then_base") -> tuple:
        e = g["exp_PM5"]
        disc = g[(e == False) & (g["now_base"] == "PM5")]  # noqa: E712
        conc = g[(e == True) & (g["now_base"] == "PM5")]  # noqa: E712
        return (int((disc[then] != "PM5").sum()), len(disc), int((conc[then] != "PM5").sum()), len(conc))

    def rates(t: tuple) -> dict:
        dr = t[0] / t[1] if t[1] else None
        cr = t[2] / t[3] if t[3] else None
        net = None if dr is None or cr is None else dr - cr
        return {"discordant_lost": [t[0], t[1]], "concordant_lost": [t[2], t[3]],
                "gross": None if dr is None else round(dr, 4), "net": None if net is None else round(net, 4)}  # fmt: skip

    def with_ci(g: pd.DataFrame) -> dict:
        """The shares of `g` with 95% intervals from a bootstrap of its panels."""
        groups = [x for _, x in g.groupby("panel")]
        cells = np.array([shares(x) for x in groups], float)
        picks = np.random.default_rng(1).integers(0, len(groups), (B_REPS, len(groups)))
        tot = cells[picks].sum(1)
        with np.errstate(invalid="ignore", divide="ignore"):
            gross_b, net_b = tot[:, 0] / tot[:, 1], tot[:, 0] / tot[:, 1] - tot[:, 2] / tot[:, 3]
        return {
            **rates(shares(g)),
            "gross_ci": [round(float(x), 4) for x in np.nanpercentile(gross_b, [2.5, 97.5])],
            "net_ci": [round(float(x), 4) for x in np.nanpercentile(net_b, [2.5, 97.5])],
            "panels": len(groups),
        }

    attribution["reimplementation_basis"] = {
        **with_ci(ev),
        "strictly_earlier_panel_references": rates(shares(ev, "then_strict_base")),
    }
    yr = ev["approved"].str[:4].astype(int)
    pku = ev["panel"].str.startswith("Phenylketonuria")
    attribution["by_period"] = {
        name: with_ci(ev[m]) for name, m in
        (("2016-2021", yr <= 2021), ("2022-2023", (yr >= 2022) & (yr <= 2023)), ("2024-2026", yr >= 2024),
         ("2016-2021 without Phenylketonuria", (yr <= 2021) & ~pku))
    }  # fmt: skip
    attribution["phenylketonuria_2016_2021"] = rates(shares(ev[(yr <= 2021) & pku]))
    # sensitivity: the concordant loss rate taken within each panel (panels with >= 5 concordant applications)
    within = [shares(g) for _, g in ev.groupby("panel")]
    kept = [t for t in within if t[3] >= 5 and t[1]]
    attribution["within_panel_net"] = {
        "net": round(sum(t[0] - t[1] * t[2] / t[3] for t in kept) / sum(t[1] for t in kept), 4),
        "discordant_in_panels_used": [sum(t[1] for t in kept), sum(t[1] for t in within)],
    }
    lags = pd.Series([v["lag_days"] for v in variants if v["lag_days"] is not None])
    attribution["release_lag_days"] = {"median": float(lags.median()), "min": int(lags.min()), "max": int(lags.max())}
    attribution["clinvar_only_excess_share"] = round(
        (attribution["clinvar_only_discordant_removed"]
         - attribution["clinvar_only_concordant_removed"] / len(conc_now) * len(disc_now)) / len(disc_now), 4)  # fmt: skip
    rep = {
        "variants": len(variants),
        "transcripts": len(transcripts),
        "releases_archived": sum(r is not None for r in releases.values()),
        "releases_used": dict(sorted(Counter(v["release"] for v in variants if v["release"]).items())),
        "reimplementation_agrees_with_tool": same,
        "reimplementation_agrees_with_tool_where_applied": [same_applied, len(applied)],
        "tool_calls": dict(Counter(v["tool"] for v in variants if v["tool"])),
        "calls_now": dict(Counter(v["now"] for v in variants if v["now"])),
        "calls_at_approval": dict(Counter(v["then"] for v in variants if v["then"])),
        "tool_call_without_reference_at_approval": sum(1 for v in variants if v["tool_base"] and not v["then_base"]),
        "tool_calls_total": sum(1 for v in variants if v["tool_base"]),
        "agreement": {
            k: {c: agreement(d, f"{k}_base", c) for c in ("PS1", "PM5")}
            for k in ("tool", "now", "then", "then_strict", "then_clinvar")
        },
        "pm5_paired": {
            "panel_pm5_met": len(met),
            "lost_then": lost,
            "gained_then": gained,
            "mcnemar_p_unclustered": float(binomtest(lost, lost + gained).pvalue) if lost + gained else None,
            "positive_agreement_change_ci": panel_bootstrap(d, lambda g: pm5_counts(g, "then_base", "now_base")),
            "ppv_change_ci": panel_bootstrap(d, lambda g: ppv_counts(g, "then_base", "now_base")),
        },
        "pm5_attribution": attribution,
        "without_transcript_switches": {
            "flagged_classifications": int(d["transcript_switch"].sum()),
            "flagged_without_reference_then": int(
                (d["transcript_switch"] & d["tool_base"].notna() & d["then_base"].isna()).sum()
            ),
            "then": {c: agreement(d[~d["transcript_switch"]], "then_base", c) for c in ("PS1", "PM5")},
        },
    }
    OUT.mkdir(exist_ok=True)
    keep = ["clinvar_id", "panel", "approved", "release", "lag_days", "tool", "now", "then", "then_strict",
            "then_clinvar", "exp_PS1", "exp_PM5", "transcript_switch"]  # fmt: skip
    d[keep].to_csv(OUT / "erepo_clinvar_history_variants.csv.gz", index=False, compression="gzip")
    (OUT / "erepo_clinvar_history.json").write_text(json.dumps(rep, indent=1))
    md = ["# PS1/PM5 against the expert panels, with the references of each panel's decision", "",
          "Produced by `benchmarks/erepo_clinvar_history.py`. `tool`: the tool's calls (current ClinVar, MANE "
          "transcript); `now`: this re-implementation on the 22 Nov 2026 release (panel's transcript); `then`: the "
          "same rule with the last monthly release before the panel's approval month plus the panels' own "
          "classifications approved by that day; `then_clinvar`: the archived release alone.", "",
          "```", json.dumps(rep, indent=1), "```", ""]  # fmt: skip
    (OUT / "erepo_clinvar_history.md").write_text("\n".join(md))
    print(json.dumps({k: rep[k] for k in ("reimplementation_agrees_with_tool_where_applied", "tool_call_without_reference_at_approval", "pm5_paired", "pm5_attribution")}, indent=1))  # fmt: skip
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
