from __future__ import annotations

import json
import textwrap
from pathlib import Path

import pytest

from rnaforge.config import load_config
from rnaforge.deseq2 import DeseqResult
from rnaforge.gates import FAIL, PASS, WARN
from rnaforge.modules import m06_de
from rnaforge.modules.m06_de import build_de_gates, count_up_down, run_de
from rnaforge.quality import load_profile
from rnaforge.state import RunState


def test_write_coldata_includes_subject_when_present(tmp_path):
    """~subject + condition tasarımı için coldata subject sütununu içermeli;
    aksi halde DESeq2 'subject not found' ile derin çöker (aracın kendi önerisi)."""
    from rnaforge.metadata import Sample
    samples = [
        Sample("s1", "control", tmp_path / "a.fq", None, batch="b1", subject="p1"),
        Sample("s2", "treated", tmp_path / "b.fq", None, batch="b1", subject="p1"),
        Sample("s3", "control", tmp_path / "c.fq", None, batch="b2", subject="p2"),
        Sample("s4", "treated", tmp_path / "d.fq", None, batch="b2", subject="p2"),
    ]
    out = tmp_path / "coldata.tsv"
    m06_de._write_coldata(samples, out)
    lines = out.read_text().splitlines()
    assert lines[0] == "sample\tcondition\tbatch\tsubject"
    assert lines[1] == "s1\tcontrol\tb1\tp1"


def test_write_coldata_omits_subject_when_absent(tmp_path):
    """subject yoksa sütun eklenmemeli (geriye uyumlu — doğrulanmış koşular değişmez)."""
    from rnaforge.metadata import Sample
    samples = [
        Sample("s1", "control", tmp_path / "a.fq", None),
        Sample("s2", "treated", tmp_path / "b.fq", None),
    ]
    out = tmp_path / "coldata.tsv"
    m06_de._write_coldata(samples, out)
    assert out.read_text().splitlines()[0] == "sample\tcondition"


def test_count_up_down():
    res = [
        {"padj": 1e-8, "log2FoldChange": 3.0},   # up
        {"padj": 1e-6, "log2FoldChange": -2.0},  # down
        {"padj": 1e-9, "log2FoldChange": 0.2},   # |lfc|<1 -> neither
        {"padj": 0.9,  "log2FoldChange": 5.0},   # padj high -> neither
        {"padj": None, "log2FoldChange": 4.0},   # NA -> neither
    ]
    up, down = count_up_down(res, fdr=0.05, lfc=1.0)
    assert (up, down) == (1, 1)


def test_high_correlation_passes():
    profile = load_profile("prokaryote")  # replicate_correlation = 0.85
    gates = build_de_gates(0.95, profile)
    assert len(gates) == 1
    assert gates[0].name == "replicate_correlation"
    assert gates[0].module == "m06_de"
    assert gates[0].status == PASS


def test_low_correlation_warns_never_fails():
    profile = load_profile("prokaryote")
    gates = build_de_gates(0.40, profile)
    assert gates[0].status == WARN          # FAIL degil!
    assert gates[0].measured == 0.40
    assert gates[0].threshold == 0.85
    assert all(g.status != FAIL for g in gates)


def test_override_marks_overridden():
    profile = load_profile("prokaryote", {"replicate_correlation": 0.20})
    gates = build_de_gates(0.40, profile)
    assert gates[0].status == PASS
    assert gates[0].overridden is True


def _setup(tmp_path):
    (tmp_path / "ref").mkdir()
    (tmp_path / "ref" / "genome.fa").write_text(">c1\nACGT\n")
    (tmp_path / "ref" / "genes.gtf").write_text('c1\ts\texon\t1\t80\t.\t+\t.\tgene_id "g1";\n')
    # m06 FASTQ kullanmaz ama load_metadata dosya varlığı ister → dummy oluştur.
    for n in ("c1.fastq", "c2.fastq", "t1.fastq", "t2.fastq"):
        (tmp_path / n).write_text("@r\nACGT\n+\nIIII\n")
    config_path = tmp_path / "config.yaml"
    config_path.write_text(textwrap.dedent(f"""
        organism: "E. coli"
        organism_type: "prokaryote"
        reference:
          genome_fasta: "{tmp_path / 'ref' / 'genome.fa'}"
          annotation_gff: "{tmp_path / 'ref' / 'genes.gtf'}"
        de:
          design: "~condition"
          reference: control
    """))
    metadata_path = tmp_path / "samples.tsv"
    metadata_path.write_text(
        "sample_id\tcondition\tfastq_1\n"
        "s1\tcontrol\tc1.fastq\n" "s2\tcontrol\tc2.fastq\n"
        "s3\ttreated\tt1.fastq\n" "s4\ttreated\tt2.fastq\n"
    )
    return config_path, metadata_path


def _mark_m05_done(run_dir):
    """m06 yalnız m05'in done+counts.tsv olmasını ister; zinciri gerçekten koşmadan
    minimal ortam kurulur."""
    run_dir = Path(run_dir)
    (run_dir / "quantification").mkdir(parents=True, exist_ok=True)
    counts = run_dir / "quantification" / "counts.tsv"
    counts.write_text("gene\ts1\ts2\ts3\ts4\ng0\t100\t110\t400\t420\ng1\t50\t55\t52\t53\n")
    (run_dir / "statistics").mkdir(exist_ok=True)
    (run_dir / "statistics" / "count_statistics.json").write_text("{}")
    RunState(run_dir).mark_done("m05_counts", [str(counts)])


def _fake_deseq2(monkeypatch, min_corr=0.95, sig_gene=True):
    def fake_run(counts_tsv, coldata_tsv, design, out_dir, reference=None,
                 contrasts=None, env="rnaforge-de"):
        out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
        rp = out_dir / "deseq2_results.tsv"
        padj = "1e-6" if sig_gene else "0.9"
        rp.write_text("gene\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\n"
                      f"g0\t250\t2.5\t0.3\t8\t1e-8\t{padj}\n"
                      "g1\t52\t0.05\t0.3\t0.1\t0.9\tNA\n")
        np = out_dir / "normalized_counts.tsv"; np.write_text("gene\ts1\ng0\t100\n")
        mp = out_dir / "de_metrics.tsv"
        mp.write_text(f"min_replicate_correlation\t{min_corr}\ncontrast\ttreated vs control\nn_genes\t2\n")
        return DeseqResult(
            results=[{"gene": "g0", "log2FoldChange": 2.5, "padj": float(padj)},
                     {"gene": "g1", "log2FoldChange": 0.05, "padj": None}],
            metrics={"min_replicate_correlation": min_corr, "contrast": "treated vs control", "n_genes": 2.0},
            results_path=rp, normalized_path=np)
    monkeypatch.setattr(m06_de, "run_deseq2", fake_run)


def test_run_de_requires_m05_done(tmp_path, monkeypatch):
    _fake_deseq2(monkeypatch)
    config_path, metadata_path = _setup(tmp_path)
    with pytest.raises(ValueError, match="m05"):
        run_de(load_config(config_path), metadata_path, tmp_path / "run")


def test_run_de_writes_results_and_coldata(tmp_path, monkeypatch):
    config_path, metadata_path = _setup(tmp_path)
    run_dir = tmp_path / "run"; _mark_m05_done(run_dir)
    _fake_deseq2(monkeypatch, min_corr=0.95, sig_gene=True)
    summary = run_de(load_config(config_path), metadata_path, run_dir)

    assert summary["n_significant"] == 1               # g0 padj<0.05 & |lfc|>=1
    assert summary["contrast"] == "treated vs control"
    coldata = (run_dir / "differential_expression" / "coldata.tsv").read_text()
    assert "s1\tcontrol" in coldata and "s3\ttreated" in coldata
    assert (run_dir / "differential_expression" / "deseq2_results.tsv").exists()
    gates = json.loads((run_dir / "quality" / "gates.json").read_text())["gates"]
    assert any(g["module"] == "m06_de" and g["status"] == "PASS" for g in gates)


def test_run_de_low_correlation_warns_but_completes(tmp_path, monkeypatch):
    config_path, metadata_path = _setup(tmp_path)
    run_dir = tmp_path / "run"; _mark_m05_done(run_dir)
    _fake_deseq2(monkeypatch, min_corr=0.40)
    summary = run_de(load_config(config_path), metadata_path, run_dir)  # DURMAZ
    assert summary["gate_counts"].get("WARN") == 1
    assert RunState(run_dir).is_done("m06_de") is True
    gates = json.loads((run_dir / "quality" / "gates.json").read_text())["gates"]
    assert any(g["module"] == "m06_de" and g["status"] == "WARN" for g in gates)


def test_run_de_resumes(tmp_path, monkeypatch):
    config_path, metadata_path = _setup(tmp_path)
    run_dir = tmp_path / "run"; _mark_m05_done(run_dir)
    _fake_deseq2(monkeypatch)
    run_de(load_config(config_path), metadata_path, run_dir)
    calls = []
    monkeypatch.setattr(m06_de, "run_deseq2", lambda *a, **k: calls.append(1))
    summary = run_de(load_config(config_path), metadata_path, run_dir)
    assert summary.get("resumed") is True
    assert calls == []


def test_cli_de_returns_zero_and_prints_verdict(tmp_path, monkeypatch, capsys):
    from rnaforge.cli import main
    from rnaforge.state import resolve_run_dir
    config_path, metadata_path = _setup(tmp_path)
    run_dir_base = tmp_path / "runs"
    _fake_deseq2(monkeypatch, min_corr=0.95)
    common = ["--config", str(config_path), "--metadata", str(metadata_path),
              "--runs-dir", str(run_dir_base), "--run-id", "demo"]
    assert main(["validate", *common]) == 0
    run_dir = resolve_run_dir(run_dir_base, "demo")
    _mark_m05_done(run_dir)
    capsys.readouterr()
    assert main(["de", *common]) == 0
    assert "quality verdict" in capsys.readouterr().out
