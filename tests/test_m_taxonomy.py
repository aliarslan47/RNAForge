"""m_taxonomy testi: run_kraken2/run_bracken/parse_bracken monkeypatch'lenir
(Kraken2/Bracken gerektirmez). Odak: abundance_matrix birleştirme doğruluğu
(taxa union, eksik=0.0), per-sample .bracken yazımı, state/resume, m_rrna_deplete
bağımlılık guard'ı."""
from __future__ import annotations

import json
import textwrap
from pathlib import Path

import pytest

from rnaforge.config import load_config
from rnaforge.metadata import load_metadata
from rnaforge.modules import m_taxonomy
from rnaforge.modules.m01_validate import run_validation
from rnaforge.modules.m_rrna_deplete import rrna_depleted_reads
from rnaforge.modules.m_taxonomy import (
    bracken_path, build_abundance_matrix, run_taxonomy, write_abundance_matrix,
)
from rnaforge.state import RunState
from tests.conftest import write_fastq


def test_build_abundance_matrix_union_missing_is_zero():
    per_sample = {
        "s1": {"E.coli": 0.6, "B.fragilis": 0.4},
        "s2": {"E.coli": 0.9},
    }
    taxa, matrix = build_abundance_matrix(per_sample)
    assert taxa == ["B.fragilis", "E.coli"]
    assert matrix["E.coli"] == {"s1": 0.6, "s2": 0.9}
    assert matrix["B.fragilis"] == {"s1": 0.4, "s2": 0.0}   # s2'de görülmedi -> 0.0


def test_write_abundance_matrix_tsv_shape(tmp_path):
    taxa, matrix = build_abundance_matrix({
        "s1": {"E.coli": 0.6}, "s2": {"E.coli": 0.9, "B.fragilis": 0.1},
    })
    out = write_abundance_matrix(tmp_path / "abundance_matrix.tsv", ["s1", "s2"], taxa, matrix)
    lines = out.read_text().strip().split("\n")
    assert lines[0] == "taxon\ts1\ts2"
    rows = {line.split("\t")[0]: line.split("\t")[1:] for line in lines[1:]}
    assert rows["E.coli"] == ["0.6", "0.9"]
    assert rows["B.fragilis"] == ["0.0", "0.1"]


def _setup(tmp_path):
    cat = tmp_path / "catalog.fa"
    cat.write_text(">g1\nACGT\n")
    ann = tmp_path / "catalog.gff"
    ann.write_text("g1\t.\tgene\t1\t4\t.\t+\t.\tID=g1\n")
    for n in ("c1.fastq", "c2.fastq", "t1.fastq", "t2.fastq"):
        write_fastq(tmp_path / n, 50, 100, "I")
    kraken_db = tmp_path / "krakendb"
    kraken_db.mkdir()
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
        taxonomy:
          kraken2_db: "{kraken_db}"
          bracken_read_len: 100
          bracken_level: "S"
          env: "rnaforge-meta"
    """))
    metadata_path = tmp_path / "samples.tsv"
    metadata_path.write_text(
        "sample_id\tcondition\tfastq_1\n"
        "s1\tcontrol\tc1.fastq\n" "s2\tcontrol\tc2.fastq\n"
        "s3\ttreated\tt1.fastq\n" "s4\ttreated\tt2.fastq\n"
    )
    return config_path, metadata_path


def _mark_m_rrna_deplete_done(run_dir, metadata_path):
    """m_taxonomy m_rrna_deplete tamamlanmışlığını VE gerçek rrna_depleted_reads()
    dosyalarını ister; gerçek SortMeRNA/fastp zincirini çalıştırmak yerine sözleşme
    yoluna doğrudan sahte rRNA'sız FASTQ yazıp state'i işaretliyoruz — m_rrna_deplete'nin
    kendi test dosyasındaki desenle aynı."""
    samples = load_metadata(metadata_path)
    state = RunState(run_dir)
    for sample in samples:
        out_dir = run_dir / "rrna_depleted" / sample.sample_id
        out_dir.mkdir(parents=True, exist_ok=True)
        out1 = out_dir / "other_1.fastq.gz"
        out1.write_bytes(b"")
        state.mark_item_done("m_rrna_deplete", sample.sample_id, {"depletion_rate": 0.5})
    state.mark_done("m_rrna_deplete", [])


def _run_prereqs(config_path, metadata_path, run_dir):
    run_validation(load_config(config_path), metadata_path, run_dir)
    _mark_m_rrna_deplete_done(run_dir, metadata_path)


def _fake_kraken_bracken(monkeypatch, bracken_data: dict[str, dict[str, float]]):
    """run_kraken2 sahte bir report dosyası yazar; run_bracken sahte bir bracken dosyası
    yazar (gerçek biçimde: header + taxon satırları, fraction_total_reads son sütun);
    parse_bracken zaten gerçek (rnaforge.kraken2) — sahte dosyayı gerçekten parse eder."""
    def fake_kraken2(reads, db, out_prefix, paired=False, threads=4, env="rnaforge-meta"):
        out_prefix = Path(out_prefix)
        out_prefix.parent.mkdir(parents=True, exist_ok=True)
        report = Path(str(out_prefix) + ".report")
        report.write_text("100.00\t1\t1\tS\t1\troot\n")
        return report

    def fake_bracken(kraken_report, db, out_path, read_len=100, level="S", env="rnaforge-meta"):
        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        sid = out_path.stem
        fractions = bracken_data[sid]
        lines = ["name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\t"
                 "new_est_reads\tfraction_total_reads"]
        for i, (taxon, frac) in enumerate(fractions.items(), start=1):
            lines.append(f"{taxon}\t{i}\tS\t100\t0\t100\t{frac}")
        out_path.write_text("\n".join(lines) + "\n")
        return out_path

    monkeypatch.setattr(m_taxonomy, "run_kraken2", fake_kraken2)
    monkeypatch.setattr(m_taxonomy, "run_bracken", fake_bracken)


def test_run_taxonomy_requires_m_rrna_deplete(tmp_path, monkeypatch):
    _fake_kraken_bracken(monkeypatch, {"s1": {"E.coli": 1.0}, "s2": {"E.coli": 1.0}})
    config_path, metadata_path = _setup(tmp_path)
    with pytest.raises(ValueError, match="m_rrna_deplete"):
        run_taxonomy(load_config(config_path), metadata_path, tmp_path / "run")


def test_run_taxonomy_requires_m_rrna_deplete_even_after_m01(tmp_path, monkeypatch):
    _fake_kraken_bracken(monkeypatch, {"s1": {"E.coli": 1.0}, "s2": {"E.coli": 1.0}})
    config_path, metadata_path = _setup(tmp_path)
    run_dir = tmp_path / "run"
    run_validation(load_config(config_path), metadata_path, run_dir)   # yalnız m01
    with pytest.raises(ValueError, match="m_rrna_deplete"):
        run_taxonomy(load_config(config_path), metadata_path, run_dir)


def test_run_taxonomy_writes_bracken_and_abundance_matrix(tmp_path, monkeypatch):
    _fake_kraken_bracken(monkeypatch, {
        "s1": {"E.coli": 0.6, "B.fragilis": 0.4},
        "s2": {"E.coli": 0.9},
        "s3": {"E.coli": 0.5, "B.fragilis": 0.5},
        "s4": {"B.fragilis": 1.0},
    })
    config_path, metadata_path = _setup(tmp_path)
    run_dir = tmp_path / "run"
    _run_prereqs(config_path, metadata_path, run_dir)
    config = load_config(config_path)
    summary = run_taxonomy(config, metadata_path, run_dir)

    assert summary["n_samples"] == 4
    assert summary["n_taxa"] == 2

    samples = load_metadata(metadata_path)
    for sample in samples:
        p = bracken_path(run_dir, sample)
        assert p.exists()
        assert p.name == f"{sample.sample_id}.bracken"

    matrix_path = run_dir / "taxonomy" / "abundance_matrix.tsv"
    assert matrix_path.exists()
    lines = matrix_path.read_text().strip().split("\n")
    header = lines[0].split("\t")
    assert header == ["taxon", "s1", "s2", "s3", "s4"]
    rows = {line.split("\t")[0]: line.split("\t")[1:] for line in lines[1:]}
    assert set(rows) == {"E.coli", "B.fragilis"}
    assert rows["E.coli"] == ["0.6", "0.9", "0.5", "0.0"]     # s4'te yok -> 0.0
    assert rows["B.fragilis"] == ["0.4", "0.0", "0.5", "1.0"]  # s2'de yok -> 0.0

    state = RunState(run_dir)
    assert state.is_done("m_taxonomy")
    for sid in ("s1", "s2", "s3", "s4"):
        assert state.is_item_done("m_taxonomy", sid)


def test_run_taxonomy_resumes_without_rerunning(tmp_path, monkeypatch):
    _fake_kraken_bracken(monkeypatch, {
        "s1": {"E.coli": 0.6, "B.fragilis": 0.4},
        "s2": {"E.coli": 0.9},
        "s3": {"E.coli": 0.5, "B.fragilis": 0.5},
        "s4": {"B.fragilis": 1.0},
    })
    config_path, metadata_path = _setup(tmp_path)
    run_dir = tmp_path / "run"
    _run_prereqs(config_path, metadata_path, run_dir)
    config = load_config(config_path)
    run_taxonomy(config, metadata_path, run_dir)

    calls = []
    def boom_kraken(*a, **k):
        calls.append("kraken2")
        raise AssertionError("run_kraken2 should not be called again on resume")
    def boom_bracken(*a, **k):
        calls.append("bracken")
        raise AssertionError("run_bracken should not be called again on resume")
    monkeypatch.setattr(m_taxonomy, "run_kraken2", boom_kraken)
    monkeypatch.setattr(m_taxonomy, "run_bracken", boom_bracken)

    summary = run_taxonomy(config, metadata_path, run_dir)
    assert summary["resumed"] is True
    assert summary["n_samples"] == 4
    assert not calls


def test_run_taxonomy_per_item_resume_skips_completed_sample(tmp_path, monkeypatch):
    """Kısmi durumda (bir örnek zaten bitmiş, run_dir'e state ile işaretli) yalnız
    kalan örnekler işlenmeli — örnek-başı resume (force olmadan tam-modül resume yolu
    devreye girmesin diye modül state'i hiç `mark_done` ile tamamlanmamış olmalı)."""
    bracken_data = {
        "s1": {"E.coli": 0.6}, "s2": {"E.coli": 0.9},
        "s3": {"E.coli": 0.5}, "s4": {"E.coli": 0.3},
    }
    _fake_kraken_bracken(monkeypatch, bracken_data)
    config_path, metadata_path = _setup(tmp_path)
    run_dir = tmp_path / "run"
    _run_prereqs(config_path, metadata_path, run_dir)
    config = load_config(config_path)

    samples = load_metadata(metadata_path)
    s1 = next(s for s in samples if s.sample_id == "s1")
    out1 = bracken_path(run_dir, s1)
    out1.parent.mkdir(parents=True, exist_ok=True)
    out1.write_text(
        "name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\t"
        "new_est_reads\tfraction_total_reads\nE.coli\t1\tS\t100\t0\t100\t0.6\n"
    )
    state = RunState(run_dir)
    state.mark_item_done("m_taxonomy", "s1", {"n_taxa": 1})

    called = []
    def tracking_kraken2(reads, db, out_prefix, paired=False, threads=4, env="rnaforge-meta"):
        called.append(Path(out_prefix).parent.name)
        out_prefix = Path(out_prefix)
        out_prefix.parent.mkdir(parents=True, exist_ok=True)
        report = Path(str(out_prefix) + ".report")
        report.write_text("100.00\t1\t1\tS\t1\troot\n")
        return report
    monkeypatch.setattr(m_taxonomy, "run_kraken2", tracking_kraken2)

    summary = run_taxonomy(config, metadata_path, run_dir)
    assert called == ["s2", "s3", "s4"]   # s1 atlandı (zaten bitmiş), gerisi işlendi
    assert summary["n_samples"] == 4


def test_run_taxonomy_raises_loudly_when_no_depleted_reads(tmp_path, monkeypatch):
    """rRNA'sız okuma bulunamazsa (m_rrna_deplete boş üretmiş) sessizce atlamak yerine
    yüksek sesle raise etmeli (Kural 7: gürültülü hata)."""
    _fake_kraken_bracken(monkeypatch, {"s1": {"E.coli": 1.0}, "s2": {"E.coli": 1.0}})
    config_path, metadata_path = _setup(tmp_path)
    run_dir = tmp_path / "run"
    run_validation(load_config(config_path), metadata_path, run_dir)
    state = RunState(run_dir)
    state.mark_done("m_rrna_deplete", [])   # işaretli ama HİÇ rrna_depleted/ çıktısı yok
    config = load_config(config_path)
    with pytest.raises(ValueError, match="no rRNA-depleted reads"):
        run_taxonomy(config, metadata_path, run_dir)
