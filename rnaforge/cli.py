"""Komut satırı girişi."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from rnaforge import __version__
from rnaforge.config import ConfigError, load_config
from rnaforge.gates import GateFailure
from rnaforge.metadata import MetadataError
from rnaforge.modules.m01_validate import run_validation
from rnaforge.modules.m02_qc import run_qc
from rnaforge.modules.m03_trim import run_trim
from rnaforge.modules.m04_quant import run_quant
from rnaforge.modules.m05_counts import run_counts
from rnaforge.modules.m06_de import run_de
from rnaforge.platform import UnsupportedPlatformError
from rnaforge.quality import load_profile
from rnaforge.report.confidence import write_confidence_card
from rnaforge.state import resolve_run_dir


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="rnaforge", description="Bulk RNA-seq pipeline")
    parser.add_argument("--version", action="version", version=__version__)
    sub = parser.add_subparsers(dest="command")

    validate = sub.add_parser("validate", help="validate config, metadata and FASTQ inputs")
    validate.add_argument("--config", required=True, type=Path)
    validate.add_argument("--metadata", required=True, type=Path)
    validate.add_argument("--runs-dir", type=Path, default=Path("runs"))
    validate.add_argument("--run-id", default="run")
    validate.add_argument(
        "--force", action="store_true",
        help="re-run even if this module already completed in this run directory",
    )

    qc = sub.add_parser("qc", help="run FastQC on raw reads (m02)")
    qc.add_argument("--config", required=True, type=Path)
    qc.add_argument("--metadata", required=True, type=Path)
    qc.add_argument("--runs-dir", type=Path, default=Path("runs"))
    qc.add_argument("--run-id", default="run")
    qc.add_argument(
        "--force", action="store_true",
        help="re-run even if m02 already completed in this run directory",
    )

    trim = sub.add_parser("trim", help="gently trim reads with fastp (m03)")
    trim.add_argument("--config", required=True, type=Path)
    trim.add_argument("--metadata", required=True, type=Path)
    trim.add_argument("--runs-dir", type=Path, default=Path("runs"))
    trim.add_argument("--run-id", default="run")
    trim.add_argument(
        "--force", action="store_true",
        help="re-run even if m03 already completed in this run directory",
    )

    quant = sub.add_parser("quant", help="align reads to reference (m04)")
    quant.add_argument("--config", required=True, type=Path)
    quant.add_argument("--metadata", required=True, type=Path)
    quant.add_argument("--runs-dir", type=Path, default=Path("runs"))
    quant.add_argument("--run-id", default="run")
    quant.add_argument(
        "--force", action="store_true",
        help="re-run even if m04 already completed in this run directory",
    )

    counts = sub.add_parser("counts", help="build gene x sample count matrix (m05)")
    counts.add_argument("--config", required=True, type=Path)
    counts.add_argument("--metadata", required=True, type=Path)
    counts.add_argument("--runs-dir", type=Path, default=Path("runs"))
    counts.add_argument("--run-id", default="run")
    counts.add_argument(
        "--force", action="store_true",
        help="re-run even if m05 already completed in this run directory",
    )

    de = sub.add_parser("de", help="differential expression with DESeq2 (m06)")
    de.add_argument("--config", required=True, type=Path)
    de.add_argument("--metadata", required=True, type=Path)
    de.add_argument("--runs-dir", type=Path, default=Path("runs"))
    de.add_argument("--run-id", default="run")
    de.add_argument(
        "--force", action="store_true",
        help="re-run even if m06 already completed in this run directory",
    )

    figures = sub.add_parser("figures", help="render figures from DE results (m07)")
    figures.add_argument("--config", required=True, type=Path)
    figures.add_argument("--metadata", required=True, type=Path)
    figures.add_argument("--runs-dir", type=Path, default=Path("runs"))
    figures.add_argument("--run-id", default="run")
    figures.add_argument(
        "--force", action="store_true",
        help="re-run even if m07 already completed in this run directory",
    )

    enrich = sub.add_parser("enrich", help="GO over-representation on DEGs (m09)")
    enrich.add_argument("--config", required=True, type=Path)
    enrich.add_argument("--metadata", required=True, type=Path)
    enrich.add_argument("--runs-dir", type=Path, default=Path("runs"))
    enrich.add_argument("--run-id", default="run")
    enrich.add_argument(
        "--force", action="store_true",
        help="re-run even if m09 already completed in this run directory",
    )

    kegg = sub.add_parser("kegg", help="KEGG pathway over-representation on DEGs (m10)")
    kegg.add_argument("--config", required=True, type=Path)
    kegg.add_argument("--metadata", required=True, type=Path)
    kegg.add_argument("--runs-dir", type=Path, default=Path("runs"))
    kegg.add_argument("--run-id", default="run")
    kegg.add_argument(
        "--force", action="store_true",
        help="re-run even if m10 already completed in this run directory",
    )

    gsea = sub.add_parser("gsea", help="GSEA on ranked genes with fgsea (m11)")
    gsea.add_argument("--config", required=True, type=Path)
    gsea.add_argument("--metadata", required=True, type=Path)
    gsea.add_argument("--runs-dir", type=Path, default=Path("runs"))
    gsea.add_argument("--run-id", default="run")
    gsea.add_argument(
        "--force", action="store_true",
        help="re-run even if m11 already completed in this run directory",
    )

    semantic = sub.add_parser("semantic", help="REVIGO-like semantic reduction of GO terms (m12)")
    semantic.add_argument("--config", required=True, type=Path)
    semantic.add_argument("--metadata", required=True, type=Path)
    semantic.add_argument("--runs-dir", type=Path, default=Path("runs"))
    semantic.add_argument("--run-id", default="run")
    semantic.add_argument(
        "--force", action="store_true",
        help="re-run even if m12 already completed in this run directory",
    )

    amr = sub.add_parser("amr", help="AMR + virulence gene overlay with abricate (m13)")
    amr.add_argument("--config", required=True, type=Path)
    amr.add_argument("--metadata", required=True, type=Path)
    amr.add_argument("--runs-dir", type=Path, default=Path("runs"))
    amr.add_argument("--run-id", default="run")
    amr.add_argument(
        "--force", action="store_true",
        help="re-run even if m13 already completed in this run directory",
    )

    report = sub.add_parser("report", help="assemble self-contained HTML report (m08)")
    report.add_argument("--config", required=True, type=Path)
    report.add_argument("--metadata", required=True, type=Path)
    report.add_argument("--runs-dir", type=Path, default=Path("runs"))
    report.add_argument("--run-id", default="run")
    report.add_argument(
        "--force", action="store_true",
        help="re-run even if m08 already completed in this run directory",
    )
    return parser


def _cmd_validate(args) -> int:
    config = load_config(args.config)
    run_dir = resolve_run_dir(args.runs_dir, args.run_id)
    profile = load_profile(config.organism_type, config.quality)
    try:
        summary = run_validation(config, args.metadata, run_dir, force=args.force)
    except GateFailure:
        # Kart FAIL yolunda da yazilmali: teshis raporu kosu basarisiz oldugunda
        # tam da bu veriye ihtiyac duyar (m01_validate kapilari enforce etmeden
        # ONCE gates.json'a yazar, o yuzden burada card icin veri hazirdir).
        write_confidence_card(run_dir, profile)
        raise
    if summary.get("resumed"):
        # Atlanan iş görünür olmalı: kullanıcı "koştu" sanıp eski sonuca bakmasın.
        print("m01_validate already completed in this run directory — reusing its result "
              "(use --force to re-run).")
    print(
        f"validation OK: {summary['n_samples']} sample(s), "
        f"platform={summary['platform']}, "
        f"organism_type={summary['organism_type']}, "
        f"conditions={summary['conditions']}"
    )
    print(f"run directory: {run_dir}")

    card_path = write_confidence_card(run_dir, profile)
    card = json.loads(card_path.read_text())
    print(f"quality verdict: {card['verdict']} "
          f"(PASS={card['counts']['PASS']} WARN={card['counts']['WARN']} "
          f"FAIL={card['counts']['FAIL']}, profile={profile.name})")
    if profile.permissive:
        print("note: this profile is permissive — thresholds are deliberately loose "
              "and results should be read with that in mind.")
    for gate, value in profile.overrides().items():
        print(f"note: threshold {gate} was overridden to {value} in the config.")
    return 0


def _cmd_qc(args) -> int:
    config = load_config(args.config)
    run_dir = resolve_run_dir(args.runs_dir, args.run_id)
    profile = load_profile(config.organism_type, config.quality)
    summary = run_qc(config, args.metadata, run_dir, force=args.force)
    if summary.get("resumed"):
        print("m02_qc already completed in this run directory — reusing its result "
              "(use --force to re-run).")
    print(f"FastQC OK: {summary['n_samples']} sample(s)")
    print(f"run directory: {run_dir}")
    card_path = write_confidence_card(run_dir, profile)
    card = json.loads(card_path.read_text())
    print(f"quality verdict: {card['verdict']} "
          f"(PASS={card['counts']['PASS']} WARN={card['counts']['WARN']} "
          f"FAIL={card['counts']['FAIL']}, profile={profile.name})")
    return 0


def _cmd_trim(args) -> int:
    config = load_config(args.config)
    run_dir = resolve_run_dir(args.runs_dir, args.run_id)
    profile = load_profile(config.organism_type, config.quality)
    try:
        summary = run_trim(config, args.metadata, run_dir, force=args.force)
    except GateFailure:
        # FAIL'de de güvence kartı yaz (gates.json diskte; verdict INVALID olur).
        write_confidence_card(run_dir, profile)
        raise
    if summary.get("resumed"):
        print("m03_trim already completed in this run directory — reusing its result "
              "(use --force to re-run).")
    print(f"trimming OK: {summary['n_samples']} sample(s)")
    print(f"run directory: {run_dir}")
    card_path = write_confidence_card(run_dir, profile)
    card = json.loads(card_path.read_text())
    print(f"quality verdict: {card['verdict']} "
          f"(PASS={card['counts']['PASS']} WARN={card['counts']['WARN']} "
          f"FAIL={card['counts']['FAIL']}, profile={profile.name})")
    return 0


def _cmd_quant(args) -> int:
    config = load_config(args.config)
    run_dir = resolve_run_dir(args.runs_dir, args.run_id)
    profile = load_profile(config.organism_type, config.quality)
    try:
        summary = run_quant(config, args.metadata, run_dir, force=args.force)
    except GateFailure:
        write_confidence_card(run_dir, profile)
        raise
    if summary.get("resumed"):
        print("m04_quant already completed in this run directory — reusing its result "
              "(use --force to re-run).")
    print(f"alignment OK: {summary['n_samples']} sample(s)")
    print(f"run directory: {run_dir}")
    card_path = write_confidence_card(run_dir, profile)
    card = json.loads(card_path.read_text())
    print(f"quality verdict: {card['verdict']} "
          f"(PASS={card['counts']['PASS']} WARN={card['counts']['WARN']} "
          f"FAIL={card['counts']['FAIL']}, profile={profile.name})")
    return 0


def _cmd_counts(args) -> int:
    config = load_config(args.config)
    run_dir = resolve_run_dir(args.runs_dir, args.run_id)
    profile = load_profile(config.organism_type, config.quality)
    try:
        summary = run_counts(config, args.metadata, run_dir, force=args.force)
    except GateFailure:
        write_confidence_card(run_dir, profile)
        raise
    if summary.get("resumed"):
        print("m05_counts already completed in this run directory — reusing its result "
              "(use --force to re-run).")
    print(f"count matrix OK: {summary['n_genes']} genes x {summary['n_samples']} sample(s)")
    print(f"run directory: {run_dir}")
    card_path = write_confidence_card(run_dir, profile)
    card = json.loads(card_path.read_text())
    print(f"quality verdict: {card['verdict']} "
          f"(PASS={card['counts']['PASS']} WARN={card['counts']['WARN']} "
          f"FAIL={card['counts']['FAIL']}, profile={profile.name})")
    return 0


def _cmd_de(args) -> int:
    config = load_config(args.config)
    run_dir = resolve_run_dir(args.runs_dir, args.run_id)
    profile = load_profile(config.organism_type, config.quality)
    summary = run_de(config, args.metadata, run_dir, force=args.force)
    if summary.get("resumed"):
        print("m06_de already completed in this run directory — reusing its result "
              "(use --force to re-run).")
    print(f"DE OK: {summary['n_significant']} significant / {summary['n_genes']} genes "
          f"({summary['contrast']})")
    print(f"run directory: {run_dir}")
    card_path = write_confidence_card(run_dir, profile)
    card = json.loads(card_path.read_text())
    print(f"quality verdict: {card['verdict']} "
          f"(PASS={card['counts']['PASS']} WARN={card['counts']['WARN']} "
          f"FAIL={card['counts']['FAIL']}, profile={profile.name})")
    return 0


def _cmd_figures(args) -> int:
    from rnaforge.modules.m07_figures import run_figures
    config = load_config(args.config)
    run_dir = resolve_run_dir(args.runs_dir, args.run_id)
    profile = load_profile(config.organism_type, config.quality)
    summary = run_figures(config, args.metadata, run_dir, force=args.force)
    if summary.get("resumed"):
        print("m07_figures already completed in this run directory — reusing its result "
              "(use --force to re-run).")
    print(f"figures OK: {summary['n_figures']} figure(s)")
    print(f"run directory: {run_dir}")
    card_path = write_confidence_card(run_dir, profile)
    card = json.loads(card_path.read_text())
    print(f"quality verdict: {card['verdict']} "
          f"(PASS={card['counts']['PASS']} WARN={card['counts']['WARN']} "
          f"FAIL={card['counts']['FAIL']}, profile={profile.name})")
    return 0


def _cmd_enrich(args) -> int:
    from rnaforge.modules.m09_enrichment import run_enrichment
    config = load_config(args.config)
    run_dir = resolve_run_dir(args.runs_dir, args.run_id)
    profile = load_profile(config.organism_type, config.quality)
    summary = run_enrichment(config, args.metadata, run_dir, force=args.force)
    if summary.get("resumed"):
        print("m09_enrichment already completed in this run directory — reusing its result "
              "(use --force to re-run).")
    print(f"enrichment OK: {summary['n_sig_up']} up + {summary['n_sig_down']} down "
          f"significant GO term(s); background={summary['background_size']} "
          f"(annotated={summary['n_annotated']})")
    print(f"run directory: {run_dir}")
    # GATE YOK: güvence kartı yalnız m06/m07 verdict'ini taşır (m09 değiştirmez).
    card_path = write_confidence_card(run_dir, profile)
    card = json.loads(card_path.read_text())
    print(f"quality verdict: {card['verdict']} "
          f"(PASS={card['counts']['PASS']} WARN={card['counts']['WARN']} "
          f"FAIL={card['counts']['FAIL']}, profile={profile.name})")
    return 0


def _cmd_kegg(args) -> int:
    from rnaforge.modules.m10_kegg import run_kegg
    config = load_config(args.config)
    run_dir = resolve_run_dir(args.runs_dir, args.run_id)
    profile = load_profile(config.organism_type, config.quality)
    summary = run_kegg(config, args.metadata, run_dir, force=args.force)
    if summary.get("resumed"):
        print("m10_kegg already completed in this run directory — reusing its result "
              "(use --force to re-run).")
    print(f"KEGG OK: {summary['n_sig_up']} up + {summary['n_sig_down']} down "
          f"significant pathway(s); background={summary['background_size']} "
          f"(annotated={summary['n_annotated']}, org={summary['organism']})")
    print(f"run directory: {run_dir}")
    # GATE YOK: güvence kartı yalnız m06/m07 verdict'ini taşır.
    card_path = write_confidence_card(run_dir, profile)
    card = json.loads(card_path.read_text())
    print(f"quality verdict: {card['verdict']} "
          f"(PASS={card['counts']['PASS']} WARN={card['counts']['WARN']} "
          f"FAIL={card['counts']['FAIL']}, profile={profile.name})")
    return 0


def _cmd_gsea(args) -> int:
    from rnaforge.modules.m11_gsea import run_gsea
    config = load_config(args.config)
    run_dir = resolve_run_dir(args.runs_dir, args.run_id)
    profile = load_profile(config.organism_type, config.quality)
    summary = run_gsea(config, args.metadata, run_dir, force=args.force)
    if summary.get("resumed"):
        print("m11_gsea already completed in this run directory — reusing its result "
              "(use --force to re-run).")
    colls = ", ".join(f"{k}: +{v['n_sig_pos']}/-{v['n_sig_neg']}"
                      for k, v in summary["collections"].items())
    print(f"GSEA OK: {summary['n_ranked']} ranked genes; significant NES ({colls})")
    print(f"run directory: {run_dir}")
    # GATE YOK: güvence kartı yalnız m06/m07 verdict'ini taşır.
    card_path = write_confidence_card(run_dir, profile)
    card = json.loads(card_path.read_text())
    print(f"quality verdict: {card['verdict']} "
          f"(PASS={card['counts']['PASS']} WARN={card['counts']['WARN']} "
          f"FAIL={card['counts']['FAIL']}, profile={profile.name})")
    return 0


def _cmd_semantic(args) -> int:
    from rnaforge.modules.m12_semantic import run_semantic
    config = load_config(args.config)
    run_dir = resolve_run_dir(args.runs_dir, args.run_id)
    profile = load_profile(config.organism_type, config.quality)
    summary = run_semantic(config, args.metadata, run_dir, force=args.force)
    if summary.get("resumed"):
        print("m12_semantic already completed in this run directory — reusing its result "
              "(use --force to re-run).")
    colls = ", ".join(f"{k}: {v['n_terms']}→{v['n_representatives']}"
                      for k, v in summary["collections"].items())
    print(f"semantic OK: GO terms reduced ({colls})")
    print(f"run directory: {run_dir}")
    # GATE YOK: güvence kartı yalnız m06/m07 verdict'ini taşır.
    card_path = write_confidence_card(run_dir, profile)
    card = json.loads(card_path.read_text())
    print(f"quality verdict: {card['verdict']} "
          f"(PASS={card['counts']['PASS']} WARN={card['counts']['WARN']} "
          f"FAIL={card['counts']['FAIL']}, profile={profile.name})")
    return 0


def _cmd_amr(args) -> int:
    from rnaforge.modules.m13_amr import run_amr
    config = load_config(args.config)
    run_dir = resolve_run_dir(args.runs_dir, args.run_id)
    profile = load_profile(config.organism_type, config.quality)
    summary = run_amr(config, args.metadata, run_dir, force=args.force)
    if summary.get("resumed"):
        print("m13_amr already completed in this run directory — reusing its result "
              "(use --force to re-run).")
    print(f"AMR OK: {summary['n_amr_genes']} AMR gene(s) ({summary['n_amr_de']} DE, "
          f"{summary['amr_db']}); {summary['n_vir_genes']} virulence gene(s) "
          f"({summary['n_vir_de']} DE, {summary['virulence_db']})")
    print(f"run directory: {run_dir}")
    # GATE YOK: güvence kartı yalnız m06/m07 verdict'ini taşır.
    card_path = write_confidence_card(run_dir, profile)
    card = json.loads(card_path.read_text())
    print(f"quality verdict: {card['verdict']} "
          f"(PASS={card['counts']['PASS']} WARN={card['counts']['WARN']} "
          f"FAIL={card['counts']['FAIL']}, profile={profile.name})")
    return 0


def _cmd_report(args) -> int:
    from rnaforge.modules.m08_report import run_report
    config = load_config(args.config)
    run_dir = resolve_run_dir(args.runs_dir, args.run_id)
    profile = load_profile(config.organism_type, config.quality)
    summary = run_report(config, args.metadata, run_dir, force=args.force)
    if summary.get("resumed"):
        print("m08_report already completed in this run directory — reusing its result "
              "(use --force to re-run).")
    print(f"report OK: {run_dir / summary['report']}")
    print(f"run directory: {run_dir}")
    card_path = write_confidence_card(run_dir, profile)
    card = json.loads(card_path.read_text())
    print(f"quality verdict: {card['verdict']} "
          f"(PASS={card['counts']['PASS']} WARN={card['counts']['WARN']} "
          f"FAIL={card['counts']['FAIL']}, profile={profile.name})")
    return 0


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command is None:
        print("error: no command given (try: rnaforge validate --help)", file=sys.stderr)
        return 2
    try:
        if args.command == "qc":
            return _cmd_qc(args)
        if args.command == "trim":
            return _cmd_trim(args)
        if args.command == "quant":
            return _cmd_quant(args)
        if args.command == "counts":
            return _cmd_counts(args)
        if args.command == "de":
            return _cmd_de(args)
        if args.command == "figures":
            return _cmd_figures(args)
        if args.command == "enrich":
            return _cmd_enrich(args)
        if args.command == "kegg":
            return _cmd_kegg(args)
        if args.command == "gsea":
            return _cmd_gsea(args)
        if args.command == "semantic":
            return _cmd_semantic(args)
        if args.command == "amr":
            return _cmd_amr(args)
        if args.command == "report":
            return _cmd_report(args)
        return _cmd_validate(args)
    except GateFailure as exc:
        print(f"error: {exc}", file=sys.stderr)
        print("no results were produced: the data did not pass the quality gates.",
              file=sys.stderr)
        return 1
    except (ConfigError, MetadataError, UnsupportedPlatformError,
            FileNotFoundError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
