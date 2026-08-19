"""m05 — Count Matrix (prokaryot: featureCounts).

m04 BAM'lerini anotasyona göre sayıp gen×örnek count matrisi (ortak sözleşme,
PLAN §5) üretir. Veri kapısı `assignment_rate`: featureCounts'un gene atadığı
okuma oranı profil eşiğinin altındaysa FAIL — çok düşük atama yanlış anotasyon/
tür demektir, sayımlar güvenilmez."""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

from rnaforge.config import Config
from rnaforge.featurecounts import run_featurecounts, tpm_fpkm
from rnaforge.gates import FAIL, PASS, WARN, GateResult, raise_if_failed, write_gate_results
from rnaforge.metadata import load_metadata
from rnaforge.quality import Profile, load_profile, profile_name_for
from rnaforge.minimap2 import count_primary_alignments
from rnaforge.nanocount import run_nanocount
from rnaforge.routing import resolve_platform, resolve_read_type
from rnaforge.state import RunState
from rnaforge.tximport import parse_tx2gene, run_tximport

MODULE_NAME = "m05_counts"
_GATE = "assignment_rate"


def build_count_gates(assignment_rates: dict[str, float],
                      profile: Profile, warn_only: bool = False) -> list[GateResult]:
    """assignment_rate kapısı. warn_only=True → eşiğin altı WARN (uzun-okuma: ONT
    CDS-only sayımda düşük atama şüpheli ama geçersiz değil)."""
    threshold = profile.threshold(_GATE)
    offenders = sorted(sid for sid, r in assignment_rates.items() if r < threshold)
    lowest = min(assignment_rates.values(), default=1.0)
    overridden = _GATE in profile.overrides()
    if offenders:
        status = WARN if warn_only else FAIL
        message = (
            f"gene atama oranı eşiğin altında ({len(offenders)} örnek: "
            f"{', '.join(offenders)}); en düşük {lowest:.2f} < {threshold:.2f}. "
            "Düşük atama yanlış anotasyon/tür → güvenilmez sayımlar."
        )
    else:
        status = PASS
        message = f"tüm örnekler assignment ≥ {threshold:.2f} (en düşük {lowest:.2f})."
    return [GateResult(
        name=_GATE, module=MODULE_NAME, status=status, message=message,
        remedy=("Anotasyon (GFF/GTF) ile referans genomun eşleştiğini ve feature_type/"
                "attribute config'inin anotasyon formatına uyduğunu doğrulayın."),
        measured=lowest, threshold=threshold, overridden=overridden,
        samples=tuple(offenders),
    )]


def run_counts(config: Config, metadata_path: Path, run_dir: Path,
               force: bool = False) -> dict:
    run_dir = Path(run_dir)
    quant_dir = run_dir / "quantification"
    stats_dir = run_dir / "statistics"
    logs_dir = run_dir / "logs"
    for d in (quant_dir, stats_dir, logs_dir):
        d.mkdir(parents=True, exist_ok=True)
    state = RunState(run_dir)
    stats_path = stats_dir / "count_statistics.json"

    if not force and state.is_done(MODULE_NAME) and stats_path.exists():
        summary = json.loads(stats_path.read_text())
        summary["resumed"] = True
        return summary

    if not state.is_done("m04_quant"):
        raise ValueError(
            "m05 (counts) requires m04 (quant) to have completed in this run directory "
            f"first: {run_dir}. Run `rnaforge quant` with the same --run-id, then re-run counts."
        )
    # ROUTER: önce organism_type (ökaryot → tximport), sonra read_type (prokaryot: kısa
    # featureCounts kapılı / uzun featureCounts -L diagnostik). Step-1'in `require_short_read`
    # muhafızı read_type dispatch'le değişti — uzun-okuma kolunun SON Step-1 muhafızıydı.
    if config.organism_type == "eukaryote":
        read_type = resolve_read_type(run_dir)
        if read_type == "long":
            summary = _counts_euk_long(config, metadata_path, run_dir,
                                       quant_dir, stats_dir, logs_dir, state)
        else:
            summary = _counts_euk(config, metadata_path, run_dir,
                                  quant_dir, stats_dir, logs_dir, state)
    else:
        read_type = resolve_read_type(run_dir)
        if read_type == "long":
            summary = _counts_long(config, metadata_path, run_dir,
                                  quant_dir, stats_dir, logs_dir, state)
        else:
            summary = _counts_short(config, metadata_path, run_dir,
                                  quant_dir, stats_dir, logs_dir, state)

    state.mark_done(MODULE_NAME, [str(stats_path), str(logs_dir / "counts.log")])
    return summary


def _write_count_outputs(result, sample_ids: list[str], quant_dir: Path,
                         log, feature_type: str, attribute: str) -> dict[str, float]:
    """counts.tsv (gene\\t<sample_id...>) + TPM/FPKM yaz; sütun→sample_id KONUMLA.
    Boş matris → yüksek sesle hata. Ortak (short + long). assignment_by_sample döner."""
    if not result.gene_ids:
        raise ValueError(
            "featureCounts assigned reads to no genes (empty matrix). Likely the "
            f"feature_type ({feature_type!r}) or attribute ({attribute!r}) does not "
            "match the annotation."
        )
    columns = list(result.counts.keys())
    assignment_by_sample = {
        sid: result.assignment_rates[col] for sid, col in zip(sample_ids, columns)
    }
    matrix_path = quant_dir / "counts.tsv"
    with matrix_path.open("w") as fh:
        fh.write("gene\t" + "\t".join(sample_ids) + "\n")
        for i, gene in enumerate(result.gene_ids):
            row = [str(result.counts[col][i]) for col in columns]
            fh.write(gene + "\t" + "\t".join(row) + "\n")
    log(f"count matrix written: {matrix_path} ({len(result.gene_ids)} genes)")

    # TPM / FPKM (gen uzunluğuyla normalize). Uzunluk yoksa (ör. mock) atla.
    if result.lengths:
        _cols, tpm, fpkm = tpm_fpkm(result.gene_ids, result.counts, result.lengths)
        for name, mat in (("tpm.tsv", tpm), ("fpkm.tsv", fpkm)):
            with (quant_dir / name).open("w") as fh:
                fh.write("gene\t" + "\t".join(sample_ids) + "\n")
                for i, gene in enumerate(result.gene_ids):
                    fh.write(gene + "\t" + "\t".join(f'{mat[c][i]:g}' for c in _cols) + "\n")
        log("expression matrices written: tpm.tsv, fpkm.tsv")
    return assignment_by_sample


def _counts_euk(config: Config, metadata_path: Path, run_dir: Path,
                quant_dir: Path, stats_dir: Path, logs_dir: Path,
                state: RunState) -> dict:
    """Ökaryot sayım (tximport, lengthScaledTPM). counts.tsv sözleşmesi.
    Salmon zaten hizalamada atadı → assignment FAIL kapısı yok (diagnostik)."""
    stats_path = stats_dir / "count_statistics.json"
    log_path = logs_dir / "counts.log"
    with log_path.open("w") as log_file:
        def log(msg: str) -> None:
            log_file.write(msg + "\n")
            log_file.flush()

        samples = load_metadata(metadata_path)
        quant_sfs = {s.sample_id: quant_dir / s.sample_id / "quant.sf" for s in samples}
        missing = [sid for sid, p in quant_sfs.items() if not p.exists()]
        if missing:
            raise ValueError(
                f"m05 eukaryote: quant.sf eksik örnek(ler): {missing} "
                "(m04 salmon koştu mu?)")
        state.heartbeat()
        res = run_tximport(quant_sfs, config.reference.tx2gene, quant_dir)
        if not res.gene_ids:
            raise ValueError("tximport 0 gen döndürdü (tx2gene eşleşmedi mi?)")
        sample_ids = [s.sample_id for s in samples]
        matrix_path = quant_dir / "counts.tsv"
        with matrix_path.open("w") as fh:
            fh.write("gene\t" + "\t".join(sample_ids) + "\n")
            for i, gene in enumerate(res.gene_ids):
                row = [f"{int(round(res.counts[sid][i]))}" for sid in sample_ids]
                fh.write(gene + "\t" + "\t".join(row) + "\n")
        log(f"count matrix written: {matrix_path} ({len(res.gene_ids)} genes)")
        summary = {
            "read_type": "short", "organism_type": "eukaryote",
            "n_samples": len(samples), "n_genes": len(res.gene_ids),
            "gate_counts": {},
        }
        stats_path.write_text(json.dumps(summary, indent=2))
    return summary


def _counts_euk_long(config: Config, metadata_path: Path, run_dir: Path,
                     quant_dir: Path, stats_dir: Path, logs_dir: Path,
                     state: RunState) -> dict:
    """Ökaryot uzun-okuma sayımı: transkriptom-hizalı BAM'den primer-hizalama sayımı →
    tx2gene ile gen'e topla (gen-içi izoform çoklu-eşleşmesi aynı gene toplanır).
    counts.tsv ortak sözleşme. Diagnostik (kapı yok — hizalama zaten eledi)."""
    stats_path = stats_dir / "count_statistics.json"
    log_path = logs_dir / "counts.log"
    with log_path.open("w") as log_file:
        def log(msg: str) -> None:
            log_file.write(msg + "\n")
            log_file.flush()

        samples = load_metadata(metadata_path)
        tx2gene = parse_tx2gene(config.reference.tx2gene)
        per_sample_gene: dict[str, dict[str, int]] = {}
        genes_seen: set[str] = set()
        for sample in samples:
            state.heartbeat()
            bam = quant_dir / sample.sample_id / "aligned.sorted.bam"
            if not bam.exists():
                raise ValueError(
                    f"m05 eukaryote-long: BAM eksik: {bam} (m04 salmon değil minimap2 koştu mu?)")
            txc = count_primary_alignments(bam)
            gc: dict[str, int] = {}
            for tx, c in txc.items():
                g = tx2gene.get(tx) or tx2gene.get(tx.split(".")[0])   # versiyonsuz yedek
                if g is None:
                    continue
                gc[g] = gc.get(g, 0) + c
            per_sample_gene[sample.sample_id] = gc
            genes_seen.update(gc)
            log(f"{sample.sample_id}: {sum(txc.values())} primer okuma → {len(gc)} gen")
        if not genes_seen:
            raise ValueError(
                "m05 eukaryote-long: 0 gen (tx2gene eşleşmedi / hizalama boş?). "
                "tx2gene transkript ID'leri BAM referans adlarıyla (transkriptom FASTA) eşleşmeli.")
        sample_ids = [s.sample_id for s in samples]
        genes = sorted(genes_seen)
        matrix_path = quant_dir / "counts.tsv"
        with matrix_path.open("w") as fh:
            fh.write("gene\t" + "\t".join(sample_ids) + "\n")
            for g in genes:
                fh.write(g + "\t" + "\t".join(
                    str(per_sample_gene[s].get(g, 0)) for s in sample_ids) + "\n")
        log(f"count matrix written: {matrix_path} ({len(genes)} genes)")
        summary = {
            "read_type": "long", "organism_type": "eukaryote",
            "n_samples": len(samples), "n_genes": len(genes), "gate_counts": {},
        }
        # İzoform-düzeyi (NanoCount EM) matrisi — ADDITIVE, best-effort. Gen matrisi (üstte)
        # değişmez; NanoCount yoksa/başarısızsa gen-düzeyi korunur (yüksek sesle log).
        n_tx = _write_transcript_matrix(samples, quant_dir, sample_ids, log)
        if n_tx is not None:
            summary["n_transcripts"] = n_tx
        stats_path.write_text(json.dumps(summary, indent=2))
    return summary


def _write_transcript_matrix(samples, quant_dir: Path, sample_ids: list[str],
                             log) -> int | None:
    """NanoCount ile izoform-düzeyi sayım → counts_transcript.tsv (transkript × örnek).
    est_count (kesirli EM tahmini) DESeq2 için yuvarlanır. Best-effort: NanoCount yoksa/
    herhangi bir örnek başarısızsa None döner (gen-düzeyi yol bozulmaz)."""
    try:
        per_sample_tx: dict[str, dict[str, float]] = {}
        tx_seen: set[str] = set()
        for sample in samples:
            bam = quant_dir / sample.sample_id / "aligned.sorted.bam"
            nc_tsv = quant_dir / sample.sample_id / "nanocount.tsv"
            est = run_nanocount(bam, nc_tsv)
            per_sample_tx[sample.sample_id] = est
            tx_seen.update(est)
            log(f"{sample.sample_id}: NanoCount {len(est)} transkript (izoform EM)")
    except Exception as exc:  # NanoCount yok/başarısız → izoform atla, gen-düzeyi korunur
        log(f"izoform niceleme ATLANDI (NanoCount yok/başarısız; gen-düzeyi korunur): {exc}")
        return None
    if not tx_seen:
        log("izoform niceleme: 0 transkript (NanoCount boş) — izoform matrisi yazılmadı")
        return None
    tx_ids = sorted(tx_seen)
    tx_matrix = quant_dir / "counts_transcript.tsv"
    with tx_matrix.open("w") as fh:
        fh.write("transcript\t" + "\t".join(sample_ids) + "\n")
        for t in tx_ids:
            fh.write(t + "\t" + "\t".join(
                str(round(per_sample_tx[s].get(t, 0.0))) for s in sample_ids) + "\n")
    log(f"transcript (isoform) matrix written: {tx_matrix} ({len(tx_ids)} transcripts)")
    return len(tx_ids)


def _counts_short(config: Config, metadata_path: Path, run_dir: Path,
                  quant_dir: Path, stats_dir: Path, logs_dir: Path,
                  state: RunState) -> dict:
    """Kısa-okuma sayım (featureCounts). assignment_rate FAIL kapısı korunur."""
    stats_path = stats_dir / "count_statistics.json"
    profile = load_profile(config.organism_type, config.quality)
    log_path = logs_dir / "counts.log"
    with log_path.open("w") as log_file:
        def log(msg: str) -> None:
            log_file.write(msg + "\n")
            log_file.flush()

        samples = load_metadata(metadata_path)
        bams = [quant_dir / s.sample_id / "aligned.sorted.bam" for s in samples]
        paired = any(s.fastq_2 is not None for s in samples)
        log(f"m05 featureCounts: {len(samples)} sample(s), "
            f"feature_type={config.quantification.feature_type}, "
            f"attribute={config.quantification.attribute}, paired={paired}")
        result = run_featurecounts(
            bams, config.reference.annotation_gff, quant_dir / "_featurecounts",
            feature_type=config.quantification.feature_type,
            attribute=config.quantification.attribute,
            paired=paired, threads=config.resources.threads,
        )
        state.heartbeat()
        sample_ids = [s.sample_id for s in samples]
        assignment_by_sample = _write_count_outputs(
            result, sample_ids, quant_dir, log,
            config.quantification.feature_type, config.quantification.attribute)

        gates = build_count_gates(assignment_by_sample, profile)
        summary = {
            "read_type": "short",
            "n_samples": len(samples), "n_genes": len(result.gene_ids),
            "samples": {sid: {"assignment_rate": assignment_by_sample[sid]} for sid in sample_ids},
            "gate_counts": dict(Counter(g.status for g in gates)),
            "expression_values": ["tpm.tsv", "fpkm.tsv"],
        }
        stats_path.write_text(json.dumps(summary, indent=2))
        write_gate_results(run_dir, gates)
        for g in gates:
            log(f"gate {g.name}: {g.status} — {g.message}")
        raise_if_failed(gates)
    return summary


def _counts_long(config: Config, metadata_path: Path, run_dir: Path,
                 quant_dir: Path, stats_dir: Path, logs_dir: Path,
                 state: RunState) -> dict:
    """Uzun-okuma sayım (featureCounts -L). Diagnostik — FAIL kapısı yok (long profil
    Step 6); assignment_rate yalnız istatistik. Aynı counts.tsv sözleşmesi → m06+ aynen."""
    platform = resolve_platform(run_dir)
    stats_path = stats_dir / "count_statistics.json"
    log_path = logs_dir / "counts.log"
    with log_path.open("w") as log_file:
        def log(msg: str) -> None:
            log_file.write(msg + "\n")
            log_file.flush()

        samples = load_metadata(metadata_path)
        bams = [quant_dir / s.sample_id / "aligned.sorted.bam" for s in samples]
        log(f"m05 featureCounts -L ({platform}): {len(samples)} sample(s), "
            f"feature_type={config.quantification.feature_type}, "
            f"attribute={config.quantification.attribute}")
        result = run_featurecounts(
            bams, config.reference.annotation_gff, quant_dir / "_featurecounts",
            feature_type=config.quantification.feature_type,
            attribute=config.quantification.attribute,
            paired=False, threads=config.resources.threads, long_read=True,
        )
        state.heartbeat()
        sample_ids = [s.sample_id for s in samples]
        assignment_by_sample = _write_count_outputs(
            result, sample_ids, quant_dir, log,
            config.quantification.feature_type, config.quantification.attribute)

        # Step 6: uzun-okuma assignment WARN kapısı (prokaryote_long; asla FAIL — ONT
        # CDS-only sayımda düşük atama şüpheli ama geçersiz değil).
        profile = load_profile(profile_name_for(config.organism_type, "long"),
                               config.quality)
        gates = build_count_gates(assignment_by_sample, profile, warn_only=True)
        summary = {
            "read_type": "long",
            "platform": platform,
            "n_samples": len(samples), "n_genes": len(result.gene_ids),
            "samples": {sid: {"assignment_rate": assignment_by_sample[sid]} for sid in sample_ids},
            "gate_counts": dict(Counter(g.status for g in gates)),
            "expression_values": ["tpm.tsv", "fpkm.tsv"],
        }
        stats_path.write_text(json.dumps(summary, indent=2))
        write_gate_results(run_dir, gates)
        for g in gates:
            log(f"gate {g.name}: {g.status} — {g.message}")
        log(f"count statistics written: {stats_path}")
    return summary
