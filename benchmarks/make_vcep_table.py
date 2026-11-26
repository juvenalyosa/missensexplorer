"""Build src/missensexplorer/vcep_params.py from ClinGen's Criteria Specification Registry (CSpec).

    python benchmarks/make_vcep_table.py [--texts DIR | --from DIR]

ClinGen Variant Curation Expert Panels (VCEPs) specify the ACMG/AMP criteria for their genes: which
strengths of each criterion they allow, and their own thresholds. CSpec (cspec.genome.network; ClinGen
curated content, CC0) records, for every criterion and strength, whether it is applicable, plus free
text with the thresholds ("Applicable", "Applicable as originally described" and "Applicable with VCEP
specification" all count as applicable). For each released specification and each of its genes, this
script keeps:

- the strengths it allows for the criteria the tool computes (PP3, BP4, PS1, PM5, PM2, BA1, BS1);
- the criteria it does not use for missense variants although CSpec marks them applicable (the text
  says so: "Do not apply to missense variants", "Missense: Do not use", ...), with the quote;
- its conditions on PS1/PM5/PP3/BP4/PM2 that the tool cannot check (the reference variant must be
  classified under the panel's specification, combinations that are capped, ...), quoted;
- its BA1 threshold, read from the text (allele frequency, operator, and whether the panel applies it to
  the filtering allele frequency), with the corrections in BA1_FIX for texts the reader gets wrong;
- the diseases of the gene's rule set (MONDO labels, from EBI's Ontology Lookup Service).

The app uses it to name the specification and to note where its evidence differs from what the panel
allows. Thresholds other than BA1 are free text and are not read. --texts saves the full JSON of each rule
set; --from builds the table from such a folder instead of the registry (no network except MONDO).
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import pprint
import re
import time
from pathlib import Path
from typing import Any

import httpx

API = "https://cspec.genome.network/cspec/api"
OLS = "https://www.ebi.ac.uk/ols4/api/ontologies/mondo/terms"
ROOT = Path(__file__).resolve().parents[1]
PARAMS = ROOT / "src" / "missensexplorer" / "vcep_params.py"
CRITERIA = ("PP3", "BP4", "PS1", "PM5", "PM2", "BA1", "BS1")
ORDER = ("Supporting", "Moderate", "Strong", "Very Strong", "Stand Alone")

# PP3/BP4 marked applicable but not used for missense variants (checked against every rule set's text)
NOT_MISSENSE = re.compile(
    r"[^.\n]*(do not (apply|use)[^.\n]{0,40}\bmissense|\bmissense[^.\n]{0,60}(do not use|not applicable|should not "
    r"be used|not (be )?appl(ied|y))|not appl(ied|y)[^.\n]{0,20}(for|to) missense|missense predictors should not be "
    r"used|do not use computational prediction)[^.\n]*",
    re.I,
)
# PP3/BP4 that a specification restricts to splicing or to non-coding variants in words the regex above does
# not catch (checked by hand against the texts, November 2026)
NOT_MISSENSE_FIX: dict[str, dict[str, str]] = {
    "GN007": {  # CDH1 (Gastric Cancer VCEP): both are splicing predictions only (quoted verbatim)
        "PP3": "At least three in silico splicing predictors in agreement (SpliceAI, MaxEntScan, SSF, GeneSplicer, "
        "HSF, TraP, varSEAK).",
        "BP4": "Splicing predictions only. At least three in silico splicing predictors in agreement (SpliceAI, "
        "MaxEntScan, SSF, GeneSplicer, HSF, TraP, varSEAK).",
    },
    "GN018": {  # AKT3, MTOR, PIK3CA, PIK3R2 (Brain Malformations VCEP)
        "BP4": "Not applicable for any variant type except for synonymous, intronic positions (except canonical "
        "splice sites) and non-coding variants in the UTRs",
    },
}
# conditions the tool cannot check, quoted (from the matching clause) when the tool gives the criterion
CONDITION = re.compile(
    r"cannot be (used|combined|applied)|may not be combined|must (also )?(be|reach|have been)[^.;]{0,80}classif|"
    r"must reach a classification|classifi\w*[^.;]{0,40}(by|using|under|according to|with)[^.;]{0,60}"
    r"(specification|guideline|criteria|rules|VCEP)|likely pathogenic[^.;]{0,40}(do|does)? ?not apply|"
    r"do not apply (if|for|when)|without PM5|Grantham|splicing impact for either variant|"
    r"not (be )?(use|appl)\w* (if|when|for)|meet criteria for pathogenicity|"
    r"according to[^.;]{0,60}(specification|VCEP|guideline|criteria)|"
    r"established (by this VCEP|using|by interpretation)|per criteria established|"
    r"based on[^.;]{0,60}\b(specifications?|VCEP)\b|\bby (the )?[^.;]{0,60}\b(VCEP|Expert Panel)\b|"
    r"reach (a|an) (likely )?pathogenic classification",
    re.I,
)
NUM = r"(\d+(?:\.\d+)?)"
# BA1 texts the reader gets wrong (checked by hand against each text, November 2026): value, operator, FAF
BA1_FIX: dict[str, tuple[float | None, str, bool]] = {
    "GN005": (None, "", False),  # 0.5% (autosomal recessive genes) or 0.1% (autosomal dominant): by gene
    "GN020": (0.005, ">", True),  # "Grpmax Filtering AF >.5%"
    "GN021": (0.007, ">=", False),  # "MAF >=0.007 (0.7%)" (the 1:30,000 is the prevalence)
    "GN032": (0.000083, ">=", False),  # "present at >=0.000083 (0.0083%) in any sub-population"
    "GN034": (0.000083, ">=", False),
    "GN067": (0.0002, ">", False),  # "above 0.02%" (the >= is about the number of alleles)
    "GN068": (0.0001, ">", False),
    "GN069": (0.0001, ">", False),
    "GN070": (0.0001, ">", False),
    "GN071": (0.000333, ">=", False),  # "greater or equal to 0.0333% (or 0.000333)"
    "GN076": (0.003, ">", False),
    "GN077": (0.001, ">", True),  # "Grpmax Filtering AF >.1%"
    "GN088": (0.004, ">", True),  # "popmax filtering allele frequency >0.00400" (1:N is the prevalence)
    "GN112": (0.004, ">=", False),  # gnomAD: "greater than or equal to 0.004 (0.4%)" (5% is for ESP/ExAC)
    # MAF thresholds written as fractions, without a % sign: VWF "Popmax MAF of >0.1" (10%), ABCA4 "greater than
    # 0.163" (16.3%); read as fractions, as the other MAF thresholds of these texts are
    "GN081": (0.1, ">", False),
    "GN090": (0.1, ">", False),
    "GN164": (0.163, ">", False),
    "GN156": (0.01, ">", True),  # "above 0.010000 (1.0%), Grpmax Filtering Allele Frequency"
}


def get(client: httpx.Client, url: str) -> dict:
    for attempt in range(4):
        try:
            r = client.get(url, timeout=120)
            if r.status_code == 200:
                return r.json()
        except httpx.HTTPError:
            pass
        time.sleep(2 * (attempt + 1))
    raise SystemExit(f"service did not answer: {url}")


def _clean(text: str) -> str:
    """Plain text: no citation links, HTML tags, markdown emphasis or repeated spaces (2x10<sup>-4</sup> is 2x10^-4)."""
    t = re.sub(r"\[<sup>[^<]*</sup>\]\([^)]*\)", "", text or "")  # citation links such as [<sup>16</sup>](#pmid…)
    t = re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", t)  # other markdown links: their text
    t = re.sub(r"<sup>([^<]*)</sup>", r"^\1", t)
    t = re.sub(r"<[^>]+>", "", t)
    return re.sub(r"\s+", " ", re.sub(r"[*_\\]", "", t).replace("\xa0", " ")).strip()


def _applicable(e: dict[str, Any]) -> bool:
    return str(e.get("applicability", "")).lower().startswith("applicable")


# BA1 texts with conditions beyond the threshold: allele counts, sex, other datasets, exclusions, founders
BA1_CONDITIONS = re.compile(
    r"alleles|allele count|\bmales?\b|hemizyg|homozyg|Exome Sequencing Project|\bESP\b|ExAC|1000 Genomes|"
    r"exclu|founder|except|at least|observation|non-cancer|quality control|NOT MET|fail|caution|allele balance",
    re.I,
)
HGVS = re.compile(r"c\.[0-9_+\-*]+[A-Za-z0-9>_+\-*]*|p\.\(?[A-Z][a-z]{2}\d+(?:[A-Z][a-z]{2}|Ter|fs|\*)[A-Za-z0-9*]*\)?")


def ba1_conditions(raw: str) -> tuple[list[str], list[str]]:
    """The clauses of a BA1 text with further conditions, and the variants it lists as BA1 exclusions."""
    excluded: list[str] = []
    if m := re.search(r"exclusion variants?(.*)$", raw or "", re.I | re.S):  # HGVS keeps its underscores here
        excluded = list(dict.fromkeys(x.strip("()") for x in HGVS.findall(m.group(1).replace("\\", ""))))
    out = []
    for c in clauses(raw):
        if re.search(r"exclusion variants", c, re.I) and excluded:
            continue  # replaced by the list itself, below
        if HGVS.fullmatch(c.split(" ")[-1].strip("()")) and excluded:
            continue  # a line of the exclusion list
        if BA1_CONDITIONS.search(c):
            out.append(c if len(c) <= 220 else c[:219].rsplit(" ", 1)[0] + "…")
    if excluded:  # 'c.14524G>A (p.Val4842Met)' when the protein change follows the DNA change
        shown: list[str] = []
        for x in excluded:
            if x.startswith("p.") and shown and shown[-1].startswith("c.") and "(" not in shown[-1]:
                shown[-1] += f" ({x})"
            else:
                shown.append(x)
        out.append("BA1 exclusion variants: " + ", ".join(shown))
    return out[:4], excluded


def ba1_threshold(sid: str, text: str, raw: str = "") -> dict[str, Any] | None:
    """The BA1 allele-frequency threshold in a specification's text, or None when there is no text.

    Args:
        sid: The specification's id (for the hand-checked corrections in BA1_FIX).
        text: The applicable BA1 texts, cleaned and joined.
        raw: The same texts as written (lines kept), for the conditions and the exclusion list.
    """
    t = _clean(text)
    if not t:
        return None
    if sid in BA1_FIX:
        value, op, faf = BA1_FIX[sid]
    else:
        value = None
        if m := re.search(r"1\s*:\s*([\d,]+)", t):
            value = 1 / float(m.group(1).replace(",", ""))
        elif m := re.search(NUM + r"\s*[x×]\s*10\s*(?:<sup>|\^)?\s*-\s*(\d+)", t):
            value = float(m.group(1)) * 10 ** -int(m.group(2))
        else:
            for m in re.finditer(
                r"(?:≥|>=|>|above|greater than(?: or equal to)?|cutoff of|cut-?off|MAF of|frequency of|is|at)\s*"
                r"(?:than\s*)?" + NUM + r"\s*(%)?",
                t,
                re.I,
            ):
                v = float(m.group(1))
                value = v / 100 if m.group(2) or v >= 1 else v
                break
        op = ">=" if re.search(r"≥|>=|greater than or equal|or greater", t, re.I) else ">"
        faf = bool(re.search(r"filter|FAF|lower bound|95%|99\.99% CI", t, re.I))
    conditions, excluded = ba1_conditions(raw or text)
    return {"value": value, "op": op, "faf": faf, "quote": t, "conditions": conditions, "excluded": excluded}


def clauses(raw: str) -> list[str]:
    """The clauses of a specification text: its lines, bullets and sentences, cleaned."""
    parts = re.split(r"\n+|(?:^|\s)\*\s+|(?<=[.;])\s+(?=[A-Z*])", raw or "")
    return [c for c in (_clean(x) for x in parts) if c]


def condition_quotes(raw_texts: list[str]) -> list[str]:
    """The clauses that state a condition the tool cannot check, each cut to 200 characters around the condition."""
    out = []
    for raw in raw_texts:
        for c in clauses(raw):
            if re.search(r"splic\w* variant|nonsense|frameshift|synonymous|intronic", c, re.I):
                continue  # conditions for other variant types
            if m := CONDITION.search(c):
                start = max(0, m.start() - 100)  # enough to keep what the reference must be ("classified as ...")
                start = c.rfind(" ", 0, start) + 1 if start else 0
                q = ("…" if start else "") + c[start:]
                out.append(q if len(q) <= 200 else q[:199].rsplit(" ", 1)[0] + "…")
    return list(dict.fromkeys(out))


def rule_set_entry(svi: dict[str, Any], rs: dict[str, Any], label: str) -> dict[str, Any]:
    """The table entry of one rule set (shared by its genes; diseases are added per gene)."""
    sid = svi["@id"].rsplit("/", 1)[-1]
    allowed: dict[str, list[str]] = {}
    not_missense: dict[str, str] = {}
    conditions: dict[str, list[str]] = {}
    ba1 = None
    for c in rs.get("criteriaCodes", []):
        crit = c.get("label")
        if crit not in CRITERIA:
            continue
        strengths = c.get("evidenceStrengths", [])
        ok = [e["label"] for e in strengths if _applicable(e) and e.get("label") in ORDER]
        allowed[crit] = sorted(ok, key=ORDER.index)
        texts = [_clean(e.get("description") or "") for e in strengths if _applicable(e)]
        if crit in ("PP3", "BP4"):
            for t in texts:
                if m := NOT_MISSENSE.search(t):
                    not_missense[crit] = m.group(0).strip()[:200]
                    break
            if crit in NOT_MISSENSE_FIX.get(sid, {}):
                not_missense[crit] = NOT_MISSENSE_FIX[sid][crit]
        if crit in ("PP3", "BP4", "PS1", "PM5", "PM2"):
            raw = [e.get("description") or "" for e in strengths if _applicable(e)]
            if found := condition_quotes(raw):
                conditions[crit] = found
        if crit == "BA1":
            raw_ba1 = "\n".join(e.get("description") or "" for e in strengths if _applicable(e))
            ba1 = ba1_threshold(sid, " / ".join(texts), raw_ba1)
    return {
        "panel": label,
        "id": sid,
        "version": svi.get("version"),
        "url": svi.get("url"),
        "allowed": allowed,
        "not_missense": not_missense,
        "conditions": conditions,
        "ba1": ba1,
        "inheritance": structured_inheritance(rs) or inheritance(ba1["quote"] if ba1 else ""),
    }


def structured_inheritance(rs: dict[str, Any]) -> str | None:
    """The mode of inheritance of the rule set's diseases (CSpec's modeOfInheritance), when there is one."""
    found = {
        re.sub(r"\s*inheritance$", "", str(m.get("@label") or "").strip(), flags=re.I).lower().replace("x-", "X-")
        for g in rs.get("genes") or []
        for d in g.get("diseases") or []
        for m in d.get("modeOfInheritance") or []
        if m.get("@label") and not str(m.get("@label")).lower().startswith(("other", "undetermined"))
    }
    found.discard("")
    return " or ".join(sorted(found)) or None


def _all_text(rs: dict[str, Any]) -> str:
    return " ".join(_clean(e.get("description") or "") for c in rs.get("criteriaCodes", []) for e in
                    c.get("evidenceStrengths", []))  # fmt: skip


def inheritance(text: str) -> str | None:
    """'autosomal dominant' or 'autosomal recessive' when the text names only one of them."""
    found = {
        "autosomal dominant" if m.group(0).lower() in ("ad variants", "autosomal dominant") else "autosomal recessive"
        for m in re.finditer(r"\b(?:AD|AR) variants\b|autosomal (?:dominant|recessive)", text, re.I)
    }
    return found.pop() if len(found) == 1 else None


def mondo_labels(client: httpx.Client, ids: set[str]) -> dict[str, str]:
    out = {}
    for k, mid in enumerate(sorted(ids)):
        d = get(client, f"{OLS}?obo_id={mid}")
        terms = (d.get("_embedded") or {}).get("terms") or []
        if terms:
            out[mid] = terms[0]["label"]
        time.sleep(0.3)
        if (k + 1) % 25 == 0:
            print(f"  MONDO {k + 1}/{len(ids)}", flush=True)
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--texts", type=Path, help="also save each rule set's full JSON here")
    g.add_argument("--from", dest="src", type=Path, help="build from rule sets saved with --texts")
    ap.add_argument("--labels", type=Path, help="JSON cache of MONDO labels (read, and updated with new ones)")
    args = ap.parse_args(argv)
    client = httpx.Client(headers={"User-Agent": "missensexplorer-benchmark"}, follow_redirects=True)
    pairs: list[tuple[dict, dict]] = []
    if args.src:
        for f in sorted(args.src.glob("*.json")):
            d = json.loads(f.read_text())
            pairs.append((d["svi"], d["ruleset"]))
        accessed = dt.date.fromtimestamp(max(f.stat().st_mtime for f in args.src.glob("*.json"))).isoformat()
    else:
        svis = [s for s in get(client, f"{API}/svis")["data"] if s.get("status") == "Released"]
        for k, s in enumerate(svis):
            for rs in s["ruleSets"]:
                rid = rs["@id"].rsplit("/", 1)[-1]
                d = get(client, f"{API}/RuleSet/id/{rid}?detail=high")
                d = d.get("data", d)
                if args.texts:
                    args.texts.mkdir(parents=True, exist_ok=True)
                    (args.texts / f"{s['@id'].rsplit('/', 1)[-1]}_{rid}.json").write_text(
                        json.dumps({"svi": s, "ruleset": d})
                    )
                pairs.append((s, d))
            print(f"  {k + 1}/{len(svis)} {s['affiliation']['label']}", flush=True)
        accessed = dt.date.today().isoformat()
    mondo = {x["label"] for _, rs in pairs for g in rs.get("genes") or [] for x in g.get("diseases") or []}
    known = json.loads(args.labels.read_text()) if args.labels and args.labels.exists() else {}
    labels = {**known, **mondo_labels(client, mondo - set(known))}
    if args.labels:
        args.labels.write_text(json.dumps(labels, indent=1, sort_keys=True))
    genes: dict[str, list[dict]] = {}
    for svi, rs in pairs:
        base = rule_set_entry(svi, rs, svi["affiliation"]["label"])
        for gene in rs.get("genes") or []:
            dis = [re.sub(r"^obsolete ", "", labels.get(x["label"], x["label"])) for x in gene.get("diseases") or []]
            genes.setdefault(gene["label"], []).append({**base, "diseases": dis})
    ids = {e["id"] for v in genes.values() for e in v}
    data = {
        "source": "ClinGen Criteria Specification Registry (CSpec), https://cspec.genome.network; ClinGen curated "
        "content, CC0 1.0; disease names from MONDO (CC BY 4.0) via EBI OLS",
        "accessed": accessed,
        "specifications": len({s["@id"] for s, _ in pairs}),
        "specifications_with_genes": len(ids),
        "genes": dict(sorted(genes.items())),
    }
    PARAMS.write_text(
        '"""ClinGen VCEP specifications by gene: allowed strengths, exclusions, conditions and BA1.\n\n'
        "Generated by benchmarks/make_vcep_table.py from ClinGen's CSpec registry (CC0); do not edit by hand.\n"
        '"""\n\n'
        f"VCEP = {pprint.pformat(data, width=110, sort_dicts=False)}\n"
    )
    print(f"wrote {PARAMS.relative_to(ROOT)}: {len(genes)} genes, {data['specifications']} released specifications "
          f"({len(ids)} with genes)")  # fmt: skip
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
