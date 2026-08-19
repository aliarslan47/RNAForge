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
from rnaforge.modules.m00_basecall import run_basecall
from rnaforge.modules.m01_validate import run_validation
from rnaforge.modules.m02_qc import run_qc
from rnaforge.modules.m03_trim import run_trim
from rnaforge.modules.m04_quant import run_quant
from rnaforge.modules.m05_counts import run_counts
from rnaforge.modules.m06_de import run_de
from rnaforge.platform import UnsupportedPlatformError
from rnaforge.basecall import basecalled_metadata_path
from rnaforge.quality import load_profile, profile_name_for
from rnaforge.report.confidence import write_confidence_card
from rnaforge.routing import resolve_read_type
from rnaforge.state import resolve_run_dir


def _effective_metadata(metadata, run_dir):
    """m00 basecall ham sinyali FASTQ'ya çevirdiyse çözülmüş metadata'yı kullan —
    TÜM downstream aşamalar (m01–m18) için, yalnız m01 değil. basecall'ın KENDİSİ
    orijinali kullanır (çözülmüş metadata'yı O üretir)."""
    resolved = basecalled_metadata_path(run_dir)
    return resolved if resolved.exists() else metadata


def _load_run_profile(config, run_dir):
    """Koşunun read_type'ına göre kalite profilini çöz (long → `<organism>_long`,
    permissive/damgalı). raw_statistics yoksa (m01 öncesi) kısa profile düşer — güven
    kartı m01'den sonra zaten yeniden yazılır. Kartın doğru profili (ör. prokaryote_long)
    damgalaması için TEK yer."""
    try:
        read_type = resolve_read_type(run_dir)
    except ValueError:
        read_type = "short"
    return load_profile(profile_name_for(config.organism_type, read_type), config.quality)


# Orkestratör aşama sırası. Çekirdek zincir m01→m08 (validate→report); opsiyoneller
# (m09-m18) rapordan ÖNCE eklenir ki rapor onların çıktısını içerebilsin.
_CORE_ORDER = ["validate", "qc", "trim", "quant", "counts", "de", "figures", "report"]
_OPTIONAL_STAGES = ["seqqc", "alignqc", "enrich", "kegg", "gsea", "semantic",
                    "amr", "operon", "ppi", "multiqc"]


def build_run_sequence(start: str | None = None, end: str | None = None,
                       include=None) -> list[str]:
    """`rnaforge run`'ın çalıştıracağı sıralı aşama listesini üretir (saf; I/O yok).

    start/end çekirdek zinciri dilimler; include opsiyonel aşamaları (kanonik sırada)
    rapordan önce (rapor dilimde değilse en sona) yerleştirir."""
    include = list(include or [])
    unknown = [s for s in include if s not in _OPTIONAL_STAGES]
    if unknown:
        raise ValueError(
            f"unknown --include stage(s): {', '.join(unknown)}; "
            f"available: {', '.join(_OPTIONAL_STAGES)}"
        )
    for label, stage in (("--from", start), ("--to", end)):
        if stage is not None and stage not in _CORE_ORDER:
            raise ValueError(
                f"{label} {stage!r} is not a core stage; core: {', '.join(_CORE_ORDER)}"
            )
    i0 = _CORE_ORDER.index(start) if start else 0
    i1 = _CORE_ORDER.index(end) if end else len(_CORE_ORDER) - 1
    if i0 > i1:
        raise ValueError(f"--from {start!r} comes after --to {end!r}")
    sliced = _CORE_ORDER[i0:i1 + 1]
    ordered_incl = [s for s in _OPTIONAL_STAGES if s in include]
    if "report" in sliced:
        idx = sliced.index("report")
        return sliced[:idx] + ordered_incl + sliced[idx:]
    return sliced + ordered_incl


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="rnaforge", description="Bulk RNA-seq pipeline")
    parser.add_argument("--version", action="version", version=__version__)
    sub = parser.add_subparsers(dest="command")

    run = sub.add_parser(
        "run", help="run the whole pipeline in order, stop-on-FAIL, resumable (m01→m08)")
    run.add_argument("--config", required=True, type=Path)
    run.add_argument("--metadata", required=True, type=Path)
    run.add_argument("--runs-dir", type=Path, default=Path("runs"))
    run.add_argument("--run-id", default="run")
    run.add_argument("--force", action="store_true",
                     help="re-run every stage even if already completed")
    run.add_argument("--from", dest="from_stage", default=None,
                     help=f"start at this core stage (one of: {', '.join(_CORE_ORDER)})")
    run.add_argument("--to", dest="to_stage", default=None,
                     help="stop after this core stage")
    run.add_argument("--include", default=None,
                     help="comma-separated optional stages to run before report "
                          f"(any of: {', '.join(_OPTIONAL_STAGES)})")

    basecall = sub.add_parser(
        "basecall", help="basecall raw signal (FAST5/POD5) to FASTQ with dorado GPU (m00)")
    basecall.add_argument("--config", required=True, type=Path)
    basecall.add_argument("--metadata", required=True, type=Path)
    basecall.add_argument("--runs-dir", type=Path, default=Path("runs"))
    basecall.add_argument("--run-id", default="run")
    basecall.add_argument(
        "--force", action="store_true",
        help="re-run even if this module already completed in this run directory",
    )

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

    operon = sub.add_parser("operon", help="operon prediction + DE coordination (m14)")
    operon.add_argument("--config", required=True, type=Path)
    operon.add_argument("--metadata", required=True, type=Path)
    operon.add_argument("--runs-dir", type=Path, default=Path("runs"))
    operon.add_argument("--run-id", default="run")
    operon.add_argument(
        "--force", action="store_true",
        help="re-run even if m14 already completed in this run directory",
    )

    ppi = sub.add_parser("ppi", help="STRING PPI network + community detection (m15)")
    ppi.add_argument("--config", required=True, type=Path)
    ppi.add_argument("--metadata", required=True, type=Path)
    ppi.add_argument("--runs-dir", type=Path, default=Path("runs"))
    ppi.add_argument("--run-id", default="run")
    ppi.add_argument(
        "--force", action="store_true",
        help="re-run even if m15 already completed in this run directory",
    )

    seqqc = sub.add_parser("seqqc", help="rRNA%% (SortMeRNA) + strandedness (RSeQC) QC gates (m16)")
    seqqc.add_argument("--config", required=True, type=Path)
    seqqc.add_argument("--metadata", required=True, type=Path)
    seqqc.add_argument("--runs-dir", type=Path, default=Path("runs"))
    seqqc.add_argument("--run-id", default="run")
    seqqc.add_argument(
        "--force", action="store_true",
        help="re-run even if m16 already completed in this run directory",
    )

    alignqc = sub.add_parser("alignqc",
                             help="insert-size + coverage + read-distribution QC (m17, diagnostik)")
    alignqc.add_argument("--config", required=True, type=Path)
    alignqc.add_argument("--metadata", required=True, type=Path)
    alignqc.add_argument("--runs-dir", type=Path, default=Path("runs"))
    alignqc.add_argument("--run-id", default="run")
    alignqc.add_argument(
        "--force", action="store_true",
        help="re-run even if m17 already completed in this run directory",
    )

    multiqc = sub.add_parser("multiqc",
                             help="MultiQC toplu görünüm — tüm araç çıktılarını birleştirir (m18)")
    multiqc.add_argument("--config", required=True, type=Path)
    multiqc.add_argument("--metadata", required=True, type=Path)
    multiqc.add_argument("--runs-dir", type=Path, default=Path("runs"))
    multiqc.add_argument("--run-id", default="run")
    multiqc.add_argument(
        "--force", action="store_true",
        help="re-run even if m18 already completed in this run directory",
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


def _cmd_basecall(args) -> int:
    config = load_config(args.config)
    run_dir = resolve_run_dir(args.runs_dir, args.run_id)
    summary = run_basecall(config, args.metadata, run_dir, force=args.force)
    if summary.get("resumed"):
        print("m00_basecall already completed in this run directory — reusing its result "
              "(use --force to re-run).")
    n_signal = sum(1 for v in summary["samples"].values() if v["input_kind"] != "fastq")
    print(f"basecall OK: {summary['n_samples']} sample(s), {n_signal} basecalled from raw signal "
          f"(model={summary['model']}, device={summary['device']})")
    print(f"run directory: {run_dir}")
    return 0


def _cmd_validate(args) -> int:
    config = load_config(args.config)
    run_dir = resolve_run_dir(args.runs_dir, args.run_id)
    profile = _load_run_profile(config, run_dir)
    try:
        summary = run_validation(config, _effective_metadata(args.metadata, run_dir), run_dir, force=args.force)
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
    profile = _load_run_profile(config, run_dir)
    summary = run_qc(config, _effective_metadata(args.metadata, run_dir), run_dir, force=args.force)
    if summary.get("resumed"):
        print("m02_qc already completed in this run directory — reusing its result "
              "(use --force to re-run).")
    _qc_tool = "NanoPlot" if summary.get("read_type") == "long" else "FastQC"
    print(f"{_qc_tool} OK: {summary['n_samples']} sample(s)")
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
    profile = _load_run_profile(config, run_dir)
    try:
        summary = run_trim(config, _effective_metadata(args.metadata, run_dir), run_dir, force=args.force)
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
    profile = _load_run_profile(config, run_dir)
    try:
        summary = run_quant(config, _effective_metadata(args.metadata, run_dir), run_dir, force=args.force)
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
    profile = _load_run_profile(config, run_dir)
    try:
        summary = run_counts(config, _effective_metadata(args.metadata, run_dir), run_dir, force=args.force)
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
    profile = _load_run_profile(config, run_dir)
    summary = run_de(config, _effective_metadata(args.metadata, run_dir), run_dir, force=args.force)
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
    profile = _load_run_profile(config, run_dir)
    summary = run_figures(config, _effective_metadata(args.metadata, run_dir), run_dir, force=args.force)
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
    profile = _load_run_profile(config, run_dir)
    summary = run_enrichment(config, _effective_metadata(args.metadata, run_dir), run_dir, force=args.force)
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
    profile = _load_run_profile(config, run_dir)
    summary = run_kegg(config, _effective_metadata(args.metadata, run_dir), run_dir, force=args.force)
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
    profile = _load_run_profile(config, run_dir)
    summary = run_gsea(config, _effective_metadata(args.metadata, run_dir), run_dir, force=args.force)
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
    profile = _load_run_profile(config, run_dir)
    summary = run_semantic(config, _effective_metadata(args.metadata, run_dir), run_dir, force=args.force)
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
    profile = _load_run_profile(config, run_dir)
    summary = run_amr(config, _effective_metadata(args.metadata, run_dir), run_dir, force=args.force)
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


def _cmd_operon(args) -> int:
    from rnaforge.modules.m14_operon import run_operon
    config = load_config(args.config)
    run_dir = resolve_run_dir(args.runs_dir, args.run_id)
    profile = _load_run_profile(config, run_dir)
    summary = run_operon(config, _effective_metadata(args.metadata, run_dir), run_dir, force=args.force)
    if summary.get("resumed"):
        print("m14_operon already completed in this run directory — reusing its result "
              "(use --force to re-run).")
    print(f"operon OK: {summary['n_operons']} operon(s) predicted "
          f"({summary['n_multi_gene']} multi-gene, {summary['n_coordinated']} coordinated DE)")
    print(f"run directory: {run_dir}")
    # GATE YOK: güvence kartı yalnız m06/m07 verdict'ini taşır.
    card_path = write_confidence_card(run_dir, profile)
    card = json.loads(card_path.read_text())
    print(f"quality verdict: {card['verdict']} "
          f"(PASS={card['counts']['PASS']} WARN={card['counts']['WARN']} "
          f"FAIL={card['counts']['FAIL']}, profile={profile.name})")
    return 0


def _cmd_ppi(args) -> int:
    from rnaforge.modules.m15_ppi import run_ppi
    config = load_config(args.config)
    run_dir = resolve_run_dir(args.runs_dir, args.run_id)
    profile = _load_run_profile(config, run_dir)
    summary = run_ppi(config, _effective_metadata(args.metadata, run_dir), run_dir, force=args.force)
    if summary.get("resumed"):
        print("m15_ppi already completed in this run directory — reusing its result "
              "(use --force to re-run).")
    print(f"PPI OK: {summary['n_deg_in_network']}/{summary['n_deg']} DEGs in network, "
          f"{summary['n_edges']} edges, {summary['n_communities']} community/-ies "
          f"(STRING {summary['taxid']}, score>={summary['min_score']})")
    print(f"run directory: {run_dir}")
    # GATE YOK: güvence kartı yalnız m06/m07 verdict'ini taşır.
    card_path = write_confidence_card(run_dir, profile)
    card = json.loads(card_path.read_text())
    print(f"quality verdict: {card['verdict']} "
          f"(PASS={card['counts']['PASS']} WARN={card['counts']['WARN']} "
          f"FAIL={card['counts']['FAIL']}, profile={profile.name})")
    return 0


def _cmd_alignqc(args) -> int:
    from rnaforge.modules.m17_alignqc import run_alignqc
    config = load_config(args.config)
    run_dir = resolve_run_dir(args.runs_dir, args.run_id)
    summary = run_alignqc(config, _effective_metadata(args.metadata, run_dir), run_dir, force=args.force)
    if summary.get("resumed"):
        print("m17_alignqc already completed in this run directory — reusing its result "
              "(use --force to re-run).")
    ins = (f"{summary['insert_size_mean']:.0f} bç" if summary.get("paired") else "yok (single-end)")
    print(f"alignqc OK: insert-size {ins}, genome mean depth "
          f"{summary.get('genome_mean_depth')}×, read-dist {summary.get('read_distribution')}")
    print(f"run directory: {run_dir}")
    return 0


def _cmd_multiqc(args) -> int:
    from rnaforge.modules.m18_multiqc import run_multiqc
    config = load_config(args.config)
    run_dir = resolve_run_dir(args.runs_dir, args.run_id)
    summary = run_multiqc(config, _effective_metadata(args.metadata, run_dir), run_dir, force=args.force)
    if summary.get("resumed"):
        print("m18_multiqc already completed in this run directory — reusing its result "
              "(use --force to re-run).")
    print(f"multiqc OK: {summary.get('n_modules')} modül, rapor {summary.get('report_relpath')}")
    print(f"run directory: {run_dir}")
    return 0


def _cmd_seqqc(args) -> int:
    from rnaforge.modules.m16_seqqc import run_seqqc
    config = load_config(args.config)
    run_dir = resolve_run_dir(args.runs_dir, args.run_id)
    profile = _load_run_profile(config, run_dir)
    summary = run_seqqc(config, _effective_metadata(args.metadata, run_dir), run_dir, force=args.force)
    if summary.get("resumed"):
        print("m16_seqqc already completed in this run directory — reusing its result "
              "(use --force to re-run).")
    match = "match" if summary["strandedness_match"] else "MISMATCH"
    print(f"seqqc OK: mean rRNA {summary['mean_rrna_fraction']:.1%}, "
          f"strandedness inferred={summary['inferred_strandedness']} vs "
          f"declared={summary['declared_strandedness']} ({match})")
    print(f"run directory: {run_dir}")
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
    profile = _load_run_profile(config, run_dir)
    summary = run_report(config, _effective_metadata(args.metadata, run_dir), run_dir, force=args.force)
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


# Aşama adı → komut fonksiyonu. Orkestratör buradan çağırır; testler monkeypatch'ler.
_STAGE_DISPATCH = {
    "validate": _cmd_validate, "qc": _cmd_qc, "trim": _cmd_trim, "quant": _cmd_quant,
    "counts": _cmd_counts, "de": _cmd_de, "figures": _cmd_figures, "report": _cmd_report,
    "seqqc": _cmd_seqqc, "alignqc": _cmd_alignqc, "enrich": _cmd_enrich, "kegg": _cmd_kegg,
    "gsea": _cmd_gsea, "semantic": _cmd_semantic, "amr": _cmd_amr, "operon": _cmd_operon,
    "ppi": _cmd_ppi, "multiqc": _cmd_multiqc,
}


def _cmd_run(args) -> int:
    include = [s.strip() for s in (args.include or "").split(",") if s.strip()]
    sequence = build_run_sequence(args.from_stage, args.to_stage, include)
    print(f"pipeline: {' → '.join(sequence)}")
    for i, name in enumerate(sequence, start=1):
        print(f"=== [{i}/{len(sequence)}] {name} ===")
        # GateFailure BURADA yakalanmaz: main'in handler'ına yükselir → pipeline durur,
        # sonraki aşamalar KOŞMAZ, exit 1. Bu stop-on-FAIL'in ta kendisi.
        rc = _STAGE_DISPATCH[name](args)
        if rc != 0:
            print(f"stage {name} returned non-zero ({rc}); stopping pipeline", file=sys.stderr)
            return rc
    print(f"pipeline complete: {len(sequence)} stage(s) — {sequence[-1]} was last")
    return 0


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command is None:
        print("error: no command given (try: rnaforge validate --help)", file=sys.stderr)
        return 2
    try:
        if args.command == "run":
            return _cmd_run(args)
        if args.command == "basecall":
            return _cmd_basecall(args)
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
        if args.command == "operon":
            return _cmd_operon(args)
        if args.command == "ppi":
            return _cmd_ppi(args)
        if args.command == "seqqc":
            return _cmd_seqqc(args)
        if args.command == "alignqc":
            return _cmd_alignqc(args)
        if args.command == "multiqc":
            return _cmd_multiqc(args)
        if args.command == "report":
            return _cmd_report(args)
        return _cmd_validate(args)
    except GateFailure as exc:
        print(f"error: {exc}", file=sys.stderr)
        print("no results were produced: the data did not pass the quality gates.",
              file=sys.stderr)
        return 1
    except (ConfigError, MetadataError, UnsupportedPlatformError,
            OSError, ValueError) as exc:
        # OSError (FileNotFoundError, gzip.BadGzipFile, izin hataları…) traceback yerine
        # temiz kullanıcı-hatası olarak yüzeye çıkmalı — Kural 7 (sessiz/çirkin hata yok).
        print(f"error: {exc}", file=sys.stderr)
        return 1
