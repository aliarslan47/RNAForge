"""m05 metatranscriptome dalı testi: katalog anotasyonuna (config.reference.
catalog_annotation) featureCounts sayımı, m04'ün ürettiği BAM'ler üzerinde
(quantification/<sid>/aligned.sorted.bam). run_featurecounts monkeypatch'lenir
(gerçek subread gerektirmez). Odak: assignment_rate DİAGNOSTİK — permissive
metatranscriptome profilinde düşük/sıfır atama GateFailure fırlatmaz; counts.tsv
+ tpm/fpkm ortak sözleşme yoluna (m06 girdisi) yazılır."""
from __future__ import annotations

import json
import textwrap
from pathlib import Path

import pytest

from rnaforge.config import load_config
from rnaforge.featurecounts import FeatureCountsResult
from rnaforge.gates import GateFailure
from rnaforge.metadata import load_metadata
from rnaforge.modules import m05_counts
from rnaforge.modules.m05_counts import run_counts
from rnaforge.state import RunState
from tests.conftest import write_fastq


def _setup(tmp_path):
    cat = tmp_path / "catalog.fa"
    cat.write_text(">g1\nACGT\n")
    ann = tmp_path / "catalog.gff"
    ann.write_text("g1\t.\tgene\t1\t4\t.\t+\t.\tID=g1\n")
    for n in ("c1.fastq", "c2.fastq", "t1.fastq", "t2.fastq"):
        write_fastq(tmp_path / n, 50, 100, "I")
    config_path = tmp_path / "config.yaml"
    config_path.write_text(textwrap.dedent(f"""
        organism: "gut community"
        organism_type: "metatranscriptome"
        platform: "illumina"
        reference:
          gene_catalog_fasta: "{cat}"
          catalog_annotation: "{ann}"
        rrna:
          db_fasta: "{cat}"
          env: "rnaforge-seqqc"
    """))
    metadata_path = tmp_path / "samples.tsv"
    metadata_path.write_text(
        "sample_id\tcondition\tfastq_1\n"
        "s1\tcontrol\tc1.fastq\n" "s2\tcontrol\tc2.fastq\n"
        "s3\ttreated\tt1.fastq\n" "s4\ttreated\tt2.fastq\n"
    )
    return config_path, metadata_path


def _seed_m04_done(tmp_path, run_dir):
    """m01+m03+m_rrna_deplete+m04_quant done state + sözleşme yolunda fake BAM'ler
    (Task 7'nin _quant_meta'sının ürettiği path: quantification/<sid>/aligned.sorted.bam).
    m04'ün kendi mantığını gerçek çalıştırmak bu testin odağı değil (o test_m04_meta.py'de);
    burada yalnız sözleşme çıktısı gerekiyor."""
    _, metadata_path = _setup(tmp_path)
    samples = load_metadata(metadata_path)
    state = RunState(run_dir)
    for m in ("m01_validate", "m03_trim", "m_rrna_deplete", "m04_quant"):
        state.mark_done(m, [])
    for sample in samples:
        d = run_dir / "quantification" / sample.sample_id
        d.mkdir(parents=True, exist_ok=True)
        (d / "aligned.sorted.bam").write_bytes(b"BAM")
    return samples


def _fake_featurecounts(monkeypatch, rate=0.02, n_genes=3):
    """rate defaults LOW (0.02) — permissive metatranscriptome eşiğine (0.02) eşit/
    tipik gerçek katalog atama oranını taklit eder; asla FAIL üretmemeli."""
    def fake_run(bams, gff, out_dir, feature_type, attribute, paired=False, threads=4,
                 env="rnaforge-quant-prok", long_read=False):
        Path(out_dir).mkdir(parents=True, exist_ok=True)
        genes = [f"g{i}" for i in range(n_genes)]
        counts = {str(b): [10 + i for i in range(n_genes)] for b in bams}
        rates = {str(b): rate for b in bams}
        return FeatureCountsResult(gene_ids=genes, counts=counts, assignment_rates=rates)
    monkeypatch.setattr(m05_counts, "run_featurecounts", fake_run)


def test_counts_meta_requires_m04_done(tmp_path, monkeypatch):
    _fake_featurecounts(monkeypatch)
    config_path, metadata_path = _setup(tmp_path)
    run_dir = tmp_path / "run"
    with pytest.raises(ValueError, match="m04"):
        run_counts(load_config(config_path), metadata_path, run_dir)


def test_counts_meta_writes_matrix_and_never_raises_on_low_assignment(tmp_path, monkeypatch):
    """KRİTİK: assignment_rate çok düşük (0.001 < 0.02 eşik) olsa bile GateFailure
    YÜKSELMEMELİ — metatranscriptome profili permissive/diagnostik."""
    config_path, metadata_path = _setup(tmp_path)
    run_dir = tmp_path / "run"
    _seed_m04_done(tmp_path, run_dir)
    _fake_featurecounts(monkeypatch, rate=0.001, n_genes=3)

    summary = run_counts(load_config(config_path), metadata_path, run_dir)   # must NOT raise

    assert summary["read_type"] == "short"
    assert summary["organism_type"] == "metatranscriptome"
    assert summary["n_samples"] == 4
    assert summary["n_genes"] == 3

    matrix = (run_dir / "quantification" / "counts.tsv").read_text().splitlines()
    assert matrix[0] == "gene\ts1\ts2\ts3\ts4"
    assert matrix[1].split("\t")[0] == "g0"

    stats = json.loads((run_dir / "statistics" / "count_statistics.json").read_text())
    assert stats["organism_type"] == "metatranscriptome"
    assert stats["samples"]["s1"]["assignment_rate"] == 0.001

    gates = json.loads((run_dir / "quality" / "gates.json").read_text())["gates"]
    assign_gates = [g for g in gates
                    if g["module"] == "m05_counts" and g["name"] == "assignment_rate"]
    assert assign_gates
    assert all(g["status"] != "FAIL" for g in assign_gates)   # diagnostik: asla FAIL


def test_counts_meta_zero_assignment_never_raises(tmp_path, monkeypatch):
    config_path, metadata_path = _setup(tmp_path)
    run_dir = tmp_path / "run"
    _seed_m04_done(tmp_path, run_dir)
    _fake_featurecounts(monkeypatch, rate=0.0, n_genes=3)
    try:
        run_counts(load_config(config_path), metadata_path, run_dir)
    except GateFailure:
        pytest.fail("metatranscriptome assignment_rate gate must never FAIL (diagnostic only)")


def test_counts_meta_uses_catalog_annotation_not_gff(tmp_path, monkeypatch):
    config_path, metadata_path = _setup(tmp_path)
    run_dir = tmp_path / "run"
    _seed_m04_done(tmp_path, run_dir)
    seen_ann = []

    def fake_run(bams, gff, out_dir, feature_type, attribute, paired=False, threads=4,
                 env="rnaforge-quant-prok", long_read=False):
        seen_ann.append(Path(gff))
        Path(out_dir).mkdir(parents=True, exist_ok=True)
        return FeatureCountsResult(
            gene_ids=["g0"], counts={str(b): [5] for b in bams},
            assignment_rates={str(b): 0.02 for b in bams},
        )
    monkeypatch.setattr(m05_counts, "run_featurecounts", fake_run)

    config = load_config(config_path)
    run_counts(config, metadata_path, run_dir)
    assert seen_ann == [config.reference.catalog_annotation]


def test_counts_meta_empty_matrix_raises(tmp_path, monkeypatch):
    config_path, metadata_path = _setup(tmp_path)
    run_dir = tmp_path / "run"
    _seed_m04_done(tmp_path, run_dir)
    _fake_featurecounts(monkeypatch, n_genes=0)
    with pytest.raises(ValueError, match="no genes"):
        run_counts(load_config(config_path), metadata_path, run_dir)


def test_counts_meta_resumes(tmp_path, monkeypatch):
    config_path, metadata_path = _setup(tmp_path)
    run_dir = tmp_path / "run"
    _seed_m04_done(tmp_path, run_dir)
    _fake_featurecounts(monkeypatch, rate=0.02)
    run_counts(load_config(config_path), metadata_path, run_dir)

    calls = []
    monkeypatch.setattr(m05_counts, "run_featurecounts", lambda *a, **k: calls.append(1))
    summary = run_counts(load_config(config_path), metadata_path, run_dir)
    assert summary.get("resumed") is True
    assert calls == []
