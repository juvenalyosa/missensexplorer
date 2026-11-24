"""Exports (CSV/TSV/Excel/JSON/FASTA) and terminal rendering."""

from __future__ import annotations

import io
import json
import textwrap
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

import pandas as pd

from .predict import PredictionResult, is_damaging

if TYPE_CHECKING:
    from openpyxl.cell.cell import Cell
    from rich.console import Console

# Compact view used by the terminal and the web table.
SUMMARY_COLUMNS = {
    "input": "Input",
    "sample": "Sample",
    "zygosity": "Zygosity",
    "gene": "Gene",
    "protein_change": "Protein change",
    "hgvsp": "HGVS p.",
    "consequence": "Consequence",
    "prediction": "Prediction",
    "effect_score": "Effect score",
    "score_with_frequency": "Score with frequency",
    "alphamissense": "AlphaMissense",
    "revel": "REVEL",
    "cadd_phred": "CADD",
    "eve": "EVE",
    "grantham": "Grantham",
    "plddt": "pLDDT",
    "key_site": "Functional site",
    "protein_features": "Protein region",
    "clinvar": "ClinVar",
    "clinvar_review": "ClinVar review",
    "gnomad_af": "gnomAD AF",
    "rsid": "dbSNP",
}

# Columns of the terminal table: (column, header).
TERMINAL_COLUMNS = (
    ("input", "Input"),
    ("gene", "Gene"),
    ("protein_change", "Change"),
    ("consequence", "Consequence"),
    ("prediction", "Prediction"),
    ("alphamissense", "AM"),
    ("revel", "REVEL"),
    ("cadd_phred", "CADD"),
    ("plddt", "pLDDT"),
    ("clinvar", "ClinVar"),
    ("gnomad_af", "gnomAD"),
)

# Lines of the terminal detail card: (label, column).
DETAIL_FIELDS = (
    ("Prediction", "prediction_basis"),
    ("Score with frequency", "score_with_frequency"),
    ("ACMG criteria", "acmg_codes"),
    ("ACMG PP3/BP4", "acmg_pp3_bp4"),
    ("ACMG PS1/PM5", "acmg_ps1_pm5"),
    ("ACMG PM2", "acmg_pm2"),
    ("Population frequency", "population_frequency"),
    ("Protein", "protein_name"),
    ("Sequence", "protein_window"),
    ("Mutant protein", "mutant_protein"),
    ("Physico-chemistry", "physchem_notes"),
    ("Grantham", "grantham_class"),
    ("Region", "protein_features"),
    ("Functional site", "key_site"),
    ("AlphaFold pLDDT", "plddt_class"),
    ("UniProt variants", "uniprot_variants_here"),
    ("Mutagenesis", "uniprot_mutagenesis_here"),
    ("ClinVar", "clinvar"),
    ("ClinVar review", "clinvar_review"),
    ("Conditions", "clinvar_conditions"),
    ("Oncogenicity", "clinvar_oncogenicity"),
    ("Transcript", "transcript"),
    ("HGVS", "hgvsc"),
    ("NMD", "nmd"),
    ("Notes", "notes"),
)

TERMINAL_ROWS = 200  # rows shown in the terminal table; the output file has them all


def _clean(df: pd.DataFrame) -> pd.DataFrame:
    """Drop the columns that are empty in every row (a table without rows is returned as is)."""
    return df.dropna(axis=1, how="all") if len(df) else df


def to_bytes(res: PredictionResult, fmt: str) -> bytes:
    """Serialise a result in one of the export formats.

    Args:
        res: The prediction result.
        fmt: ``"csv"``, ``"tsv"``, ``"json"`` (the table with the run information, warnings, statistics
            and the inputs not analysed), ``"xlsx"`` (sheets Summary, All columns, Not analysed and Run
            info) or ``"fasta"``/``"fa"`` (protein sequences, see :func:`to_fasta`).

    Returns:
        The file contents (UTF-8 for the text formats).

    Raises:
        ValueError: If ``fmt`` is not one of these.
    """
    df = res.table
    if fmt == "csv":
        return df.to_csv(index=False).encode()
    if fmt == "tsv":
        return df.to_csv(index=False, sep="\t").encode()
    if fmt == "json":
        recs = json.loads(df.to_json(orient="records"))
        return json.dumps(
            {
                "run": res.run,
                "assembly": res.assembly,
                "warnings": res.warnings,
                "stats": res.stats,
                "variants": recs,
                "errors": res.errors.to_dict("records"),
            },
            indent=2,
            default=str,
        ).encode()
    if fmt == "xlsx":
        return _xlsx_bytes(res)
    if fmt in ("fasta", "fa"):
        return to_fasta(res).encode()
    raise ValueError(f"unknown format {fmt}")


def _xlsx_bytes(res: PredictionResult) -> bytes:
    """Return the Excel workbook: summary, all columns, inputs not analysed and run information."""
    df = res.table
    buf = io.BytesIO()
    with pd.ExcelWriter(buf, engine="openpyxl") as xw:
        view = df[[c for c in SUMMARY_COLUMNS if c in df.columns]].rename(columns=SUMMARY_COLUMNS)
        _clean(view).to_excel(xw, sheet_name="Summary", index=False)
        df.to_excel(xw, sheet_name="All columns", index=False)
        if len(res.errors):
            res.errors.to_excel(xw, sheet_name="Not analysed", index=False)
        run_info(res).to_excel(xw, sheet_name="Run info", index=False)
        for ws in xw.book.worksheets:
            ws.freeze_panes = "A2"
            for col in ws.columns:
                width = max((len(str(c.value)) for c in col[:200] if c.value is not None), default=8)
                first = cast("Cell", col[0])  # these sheets have no merged cells
                ws.column_dimensions[first.column_letter].width = min(60, max(8, width + 2))
    return buf.getvalue()


def run_info(res: PredictionResult) -> pd.DataFrame:
    """Provenance as two columns (for the Excel 'Run info' sheet).

    Returns:
        Columns ``item`` and ``value`` (text): tool, versions, options, input files, sources, the
        missense model and the warnings.
    """
    run = res.run or {}
    rows = [(k, run.get(k)) for k in ("tool", "started_utc", "seconds", "assembly", "annotation", "python")]
    for k, v in (run.get("options") or {}).items():
        rows.append((f"option: {k}", v))
    for f in (run.get("inputs") or {}).get("files", []):
        rows.append((f"input file: {f['name']}", f"{f['bytes']} bytes, sha256 {f['sha256']}"))
    for k, v in (run.get("sources") or {}).items():
        rows.append((f"source: {k}", v))
    for k, v in (run.get("missense_model") or {}).items():
        rows.append((f"missense model: {k}", v))
    for w in res.warnings:
        rows.append(("warning", w))
    return pd.DataFrame([(k, "" if v is None else str(v)) for k, v in rows], columns=["item", "value"])


def to_fasta(res: PredictionResult, include_wildtype: bool = True) -> str:
    """Return the mutant protein sequences as FASTA, 60 residues per line.

    Args:
        res: The prediction result.
        include_wildtype: Also write each protein's wild-type sequence (once, before its first mutant).

    Returns:
        Records ``>WT|gene|protein`` and ``>MUT|header|note``; an empty string when there are none.
    """
    lines: list[str] = []
    seen_wt: set[Any] = set()
    for i, s in res.sequences.items():
        row = res.table.iloc[i]
        pid = row.get("protein_id") or "protein"
        if include_wildtype and pid not in seen_wt:
            seen_wt.add(pid)
            lines.append(f">WT|{row.get('gene')}|{pid}")
            lines += textwrap.wrap(s["wt"], 60)
        lines.append(f">MUT|{s['header']}|{s['note']}")
        lines += textwrap.wrap(s["mut"], 60) or [""]
    return "\n".join(lines) + ("\n" if lines else "")


def write(res: PredictionResult, path: str | Path) -> Path:
    """Write a result to a file, in the format given by its extension.

    Args:
        res: The prediction result.
        path: Output file: .csv, .tsv (or .txt), .json, .xlsx (or .xls) or .fasta (or .fa).

    Returns:
        The path written.

    Raises:
        ValueError: If the extension is not one of these.
    """
    path = Path(path)
    fmt = path.suffix.lower().lstrip(".")
    fmt = {"txt": "tsv", "fa": "fasta", "xls": "xlsx"}.get(fmt, fmt)
    path.write_bytes(to_bytes(res, fmt))
    return path


# --------------------------------------------------------------------------- terminal


def _fmt(v: Any, nd: int = 2) -> str:
    """Format a table value: '' if missing, ``nd`` decimals for floats (scientific below 0.001)."""
    if v is None or (isinstance(v, float) and v != v):
        return ""
    if isinstance(v, float):
        if 0 < abs(v) < 1e-3:
            return f"{v:.1e}"
        return f"{v:.{nd}f}"
    return str(v)


def print_result(res: PredictionResult, console: Console | None = None, details: bool | None = None) -> None:
    """Print a result to the terminal: summary line, warnings, table, detail cards and errors.

    Args:
        res: The prediction result.
        console: Rich console to print to; None for a new one (standard output).
        details: Print a detail card for every row; None to print them only for up to 4 rows.
    """
    from rich.console import Console

    con = console or Console()
    df = res.table
    _print_header(con, res)

    if len(df):
        _print_table(con, df)

    if details is None:
        details = len(df) <= 4
    if details:
        _print_details(con, df)

    if len(res.errors):
        con.print(f"[red]{len(res.errors)} input(s) could not be analysed:[/]")
        for _, e in res.errors.head(50).iterrows():
            con.print(f"  [red]-[/] {e['input']} ({e['source']}): {e['error']}")


def _print_header(con: Console, res: PredictionResult) -> None:
    """Print the summary line, the warnings and the count of variants without a protein effect."""
    st = res.stats
    con.print(
        f"[bold]MissensExplorer[/] · {res.assembly} · {st.get('inputs', 0)} input(s) · "
        f"{st.get('rows', 0)} protein-level result(s) · {st.get('likely_damaging', 0)} likely damaging · "
        f"{st.get('seconds', 0)} s"
    )
    for w in res.warnings:
        con.print(f"[yellow]! {w}[/]")
    if st.get("no_protein_effect"):
        con.print(
            f"[dim]{st['no_protein_effect']} variant(s) have no coding/splice consequence (use --all to list them).[/]"
        )
    skipped = [
        f"{st[key]} {what}"
        for key, what in (
            ("file_skipped_not_carried", "not carried by the sample(s)"),
            ("file_skipped_outside_genes", "outside the selected genes"),
            ("file_skipped_symbolic", "symbolic/structural alleles"),
            ("file_skipped_filter", "not passing FILTER"),
        )
        if st.get(key)
    ]
    if skipped:
        con.print(f"[dim]Skipped file records: {' · '.join(skipped)}.[/]")


def _print_table(con: Console, df: pd.DataFrame) -> None:
    """Print the results table (the first TERMINAL_ROWS rows; fewer columns on a narrow terminal)."""
    from rich.table import Table

    t = Table(show_lines=False, header_style="bold", expand=False)
    cols = list(TERMINAL_COLUMNS)
    if con.width < 150:  # narrow terminal: keep the essentials
        cols = [c for c in cols if c[0] not in ("consequence", "plddt", "gnomad_af")]
    if df["zygosity"].notna().any() and con.width >= 120:
        cols.insert(1, ("zygosity", "Zyg."))
    for c, label in cols:
        t.add_column(
            label,
            overflow="fold",
            max_width=34,
            no_wrap=c in ("gene", "protein_change", "alphamissense", "revel", "cadd_phred", "plddt"),
        )
    for _, r in df.head(TERMINAL_ROWS).iterrows():
        pred = r.get("prediction") or ""
        style = (
            "red"
            if is_damaging(pred)
            else "green"
            if "benign" in pred.lower() or pred.startswith("Silent")
            else "yellow"
        )
        vals = []
        for c, _ in cols:
            v = _table_value(r, c)
            vals.append(f"[{style}]{pred}[/]" if c == "prediction" else _fmt(v, 1 if c == "cadd_phred" else 2))
        t.add_row(*vals)
    con.print(t)
    if len(df) > TERMINAL_ROWS:
        con.print(f"[dim]... {len(df) - TERMINAL_ROWS} more rows in the output file[/]")


def _table_value(r: pd.Series, c: str) -> Any:
    """Return column ``c`` of row ``r`` shortened for the terminal table (consequence, ClinVar, zygosity)."""
    v = r.get(c)
    if c == "consequence" and isinstance(v, str):
        v = v.replace("_variant", "").replace("_", " ")
    if c == "clinvar" and isinstance(v, str):
        v = v.split(",")[0].replace("_", " ").replace(" (submissions)", "")
        if r.get("clinvar_stars") is not None and r.get("clinvar_stars") == r.get("clinvar_stars"):
            v += " " + "★" * int(r["clinvar_stars"])
    if c == "zygosity" and isinstance(v, str):
        v = {"heterozygous": "het", "homozygous": "hom", "hemizygous": "hemi"}.get(v, v)
    return v


def _print_details(con: Console, df: pd.DataFrame) -> None:
    """Print a detail card (DETAIL_FIELDS) for every row."""
    from rich.panel import Panel

    for _, r in df.iterrows():
        body = []
        for label, key in DETAIL_FIELDS:
            v = r.get(key)
            if v is None or (isinstance(v, float) and v != v) or v == "":
                continue
            if key == "plddt_class":
                v = f"{r.get('plddt')} ({v})"
            if key == "score_with_frequency":
                v = f"{v:.2f} (model score and gnomAD frequency together; for ranking, not ACMG evidence)"
            if key == "grantham_class":
                v = f"{int(r['grantham'])} ({v}); BLOSUM62 {int(r['blosum62'])}"
            s = str(v)
            body.append(f"[bold]{label}:[/] {s if len(s) < 400 else s[:400] + '…'}")
        title = f"{r.get('gene')} {r.get('hgvsp') or r.get('protein_change') or ''} · {r.get('input')}"
        con.print(Panel("\n".join(body), title=title, expand=False))
