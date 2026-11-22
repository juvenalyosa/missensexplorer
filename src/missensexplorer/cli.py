"""Command line interface.

    missensexplorer "BRAF V600E" "TP53: R175H, R248Q" rs28934578
    missensexplorer -f sample.vcf.gz --genes TP53,BRCA1 -o results.xlsx
    missensexplorer --gene TP53 R175H R248Q R273H --fasta mutants.fa
    missensexplorer web
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from . import __version__


def _predict_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="missensexplorer",
        description="Predict how DNA variants change the encoded human protein "
                    "(Ensembl VEP + AlphaMissense, REVEL, CADD, EVE, SpliceAI, UniProt, AlphaFold).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("variants", nargs="*", help="variants as text (use '-' to read from stdin)")
    p.add_argument("-f", "--file", action="append", default=[], metavar="FILE",
                   help="VCF/VCF.gz/CSV/TSV/XLSX/TXT input (repeatable)")
    p.add_argument("-g", "--gene", help="default gene for bare changes such as 'V600E' or 'c.1799T>A'")
    p.add_argument("--genes", help="comma-separated genes: only report variants in these genes")
    p.add_argument("--sample", help="VCF sample to report (default: all carriers)")
    p.add_argument("--pass-only", action="store_true", help="VCF: keep only FILTER=PASS records")
    p.add_argument("-a", "--assembly", default="auto", choices=["auto", "GRCh38", "GRCh37"],
                   help="genome build for coordinates (default: auto-detect)")
    p.add_argument("--all", action="store_true", help="also report non-coding consequences")
    p.add_argument("--no-predictors", action="store_true", help="skip AlphaMissense/REVEL/CADD/EVE/SpliceAI (faster)")
    p.add_argument("--no-protein-context", action="store_true", help="skip UniProt/AlphaFold annotation")
    p.add_argument("-o", "--output", action="append", default=[], metavar="FILE",
                   help="write results (.csv, .tsv, .xlsx, .json); repeatable")
    p.add_argument("--fasta", metavar="FILE", help="write wild-type and mutant protein sequences")
    p.add_argument("--details", action="store_true", help="print a detail card for every variant")
    p.add_argument("-w", "--workers", type=int, default=4, help="parallel API requests (default 4)")
    p.add_argument("--no-cache", action="store_true", help="do not use the local response cache")
    p.add_argument("-q", "--quiet", action="store_true", help="no terminal table")
    p.add_argument("-v", "--verbose", action="store_true")
    p.add_argument("--version", action="version", version=f"missensexplorer {__version__}")
    return p


def run_predict(argv: list[str]) -> int:
    args = _predict_parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO if args.verbose else logging.WARNING, format="%(levelname)s %(message)s")
    from rich.console import Console
    from rich.progress import BarColumn, Progress, TextColumn, TimeElapsedColumn

    from .predict import PredictOptions, predict
    from .report import print_result, to_fasta, write

    err = Console(stderr=True)
    texts = []
    for v in args.variants:
        texts.append(sys.stdin.read() if v == "-" else v)
    for f in args.file:
        if not Path(f).exists():
            err.print(f"[red]File not found: {f}[/]")
            return 2
    if not texts and not args.file:
        _predict_parser().print_help()
        return 1

    opts = PredictOptions(
        assembly=args.assembly, default_gene=args.gene,
        genes=[g for g in (args.genes or "").replace(";", ",").split(",") if g.strip()] or None,
        sample=args.sample, pass_only=args.pass_only, include_noncoding=args.all,
        predictors=not args.no_predictors, protein_context=not args.no_protein_context,
        workers=args.workers, cache=not args.no_cache,
    )
    with Progress(TextColumn("{task.description}"), BarColumn(), TextColumn("{task.completed}/{task.total}"),
                  TimeElapsedColumn(), console=err, transient=True, disable=args.quiet) as prog:
        tasks = {}

        def cb(stage, done, total):
            if stage not in tasks:
                tasks[stage] = prog.add_task(stage, total=max(total, 1))
            prog.update(tasks[stage], completed=done, total=max(total, 1))

        try:
            res = predict(text="\n".join(texts), files=args.file, options=opts, progress=cb)
        except ValueError as e:
            err.print(f"[red]{e}[/]")
            return 2

    if not args.quiet:
        print_result(res, Console(), details=True if args.details else None)
    outputs = list(args.output)
    if not outputs and args.file:
        outputs = [f"{Path(args.file[0]).name.split('.')[0]}.missensexplorer.csv"]
    for out in outputs:
        p = write(res, out)
        err.print(f"[green]Saved[/] {p} ({len(res.table)} rows)")
    if args.fasta:
        Path(args.fasta).write_text(to_fasta(res))
        err.print(f"[green]Saved[/] {args.fasta} ({len(res.sequences)} mutant sequences)")
    return 0


def run_web(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(prog="missensexplorer web", description="Start the web interface")
    ap.add_argument("--port", type=int, default=8501)
    ap.add_argument("--host", default="localhost")
    args = ap.parse_args(argv)
    try:
        from streamlit.web import cli as stcli
    except ImportError:
        print("The web interface needs Streamlit:  pip install 'missensexplorer[web]'", file=sys.stderr)
        return 1
    app = Path(__file__).with_name("webapp.py")
    sys.argv = ["streamlit", "run", str(app), "--server.port", str(args.port), "--server.address", args.host,
                "--browser.gatherUsageStats", "false"]
    return stcli.main()


def run_cache(argv: list[str]) -> int:
    from .cache import Cache

    c = Cache()
    if argv[:1] == ["clear"]:
        c.clear()
        print(f"Cache cleared ({c.path})")
    else:
        print(c.path)
    return 0


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] == ["web"]:
        return run_web(argv[1:])
    if argv[:1] == ["cache"]:
        return run_cache(argv[1:])
    if argv[:1] == ["predict"]:
        argv = argv[1:]
    return run_predict(argv)


if __name__ == "__main__":
    sys.exit(main())
