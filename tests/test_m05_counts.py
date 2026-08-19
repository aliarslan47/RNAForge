from __future__ import annotations

import json
import textwrap
from pathlib import Path

import pytest

from rnaforge.config import load_config
from rnaforge.featurecounts import FeatureCountsResult
from rnaforge.gates import FAIL, PASS, WARN, GateFailure
from rnaforge.modules import m05_counts
from rnaforge.modules.m05_counts import build_count_gates, run_counts
from rnaforge.quality import load_profile
from tests.conftest import write_fastq


def test_all_above_threshold_passes():
    profile = load_profile("prokaryote")  # assignment_rate = 0.50
    gates = build_count_gates({"s1": 0.95, "s2": 0.80}, profile)
    assert len(gates) == 1
    assert gates[0].name == "assignment_rate"
    assert gates[0].module == "m05_counts"
    assert gates[0].status == PASS


def test_below_threshold_fails():
    profile = load_profile("prokaryote")
    gates = build_count_gates({"s1": 0.95, "s2": 0.10}, profile)
    g = gates[0]
    assert g.status == FAIL
    assert g.samples == ("s2",)
    assert g.measured == 0.10
    assert g.threshold == 0.50


def test_build_count_gates_warn_only_is_warn_not_fail():
    """Uzun-okuma: düşük atama ŞÜPHELİ (WARN), geçersiz DEĞİL (ONT CDS-only doğal düşük)."""
    profile = load_profile("prokaryote_long")   # assignment_rate = 0.05
    gates = build_count_gates({"s1": 0.02}, profile, warn_only=True)
    assert gates[0].status == WARN
    assert gates[0].samples == ("s1",)


def test_override_marks_overridden():
    profile = load_profile("prokaryote", {"assignment_rate": 0.05})
    gates = build_count_gates({"s1": 0.10}, profile)
    assert gates[0].status == PASS
    assert gates[0].overridden is True


def _setup(tmp_path):
    (tmp_path / "ref").mkdir()
    (tmp_path / "ref" / "genome.fa").write_text(">c1\n" + "ACGT" * 25 + "\n")
    (tmp_path / "ref" / "genes.gtf").write_text('c1\ts\texon\t1\t80\t.\t+\t.\tgene_id "g1";\n')
    for n in ("c1.fastq", "c2.fastq", "t1.fastq", "t2.fastq"):
        write_fastq(tmp_path / n, 200, 150, "I")
    config_path = tmp_path / "config.yaml"
    config_path.write_text(textwrap.dedent(f"""
        organism: "E. coli"
        organism_type: "prokaryote"
        reference:
          genome_fasta: "{tmp_path / 'ref' / 'genome.fa'}"
          annotation_gff: "{tmp_path / 'ref' / 'genes.gtf'}"
    """))
    metadata_path = tmp_path / "samples.tsv"
    metadata_path.write_text(
        "sample_id\tcondition\tfastq_1\n"
        "s1\tcontrol\tc1.fastq\n" "s2\tcontrol\tc2.fastq\n"
        "s3\ttreated\tt1.fastq\n" "s4\ttreated\tt2.fastq\n"
    )
    return config_path, metadata_path


def _prep_through_m04(config_path, metadata_path, run_dir, monkeypatch):
    """m01(gerçek)+m03(fake fastp)+m04(fake bowtie2) done state + BAM üretir."""
    from rnaforge.modules.m01_validate import run_validation
    from rnaforge.modules import m03_trim, m04_quant
    from rnaforge.modules.m03_trim import run_trim, trimmed_name
    from rnaforge.modules.m04_quant import run_quant
    from rnaforge.fastp import FastpResult
    from rnaforge.bowtie2 import AlignmentResult
    run_validation(load_config(config_path), metadata_path, run_dir)

    def fake_fastp(fastq_1, out_dir, min_length, fastq_2=None, aggressive_quality=False, env="rnaforge-qc"):
        out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
        out1 = out_dir / trimmed_name(Path(fastq_1)); out1.write_text("@r\nACGT\n+\nIIII\n")
        (out_dir / "fastp.json").write_text("{}")
        return FastpResult(200, 196, 0.98, out1=out1)
    monkeypatch.setattr(m03_trim, "run_fastp", fake_fastp)
    run_trim(load_config(config_path), metadata_path, run_dir)

    monkeypatch.setattr(m04_quant, "build_index", lambda g, i, env="rnaforge-quant-prok": Path(i) / "genome")
    def fake_align(index_prefix, out_dir, fastq_1, fastq_2=None, threads=4, env="rnaforge-quant-prok"):
        out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
        bam = out_dir / "aligned.sorted.bam"; bam.write_bytes(b"BAM")
        return AlignmentResult(bam=bam, alignment_rate=0.95)
    monkeypatch.setattr(m04_quant, "run_bowtie2", fake_align)
    run_quant(load_config(config_path), metadata_path, run_dir)


def _fake_featurecounts(monkeypatch, rate=0.9, n_genes=3):
    def fake_run(bams, gff, out_dir, feature_type, attribute, paired=False, threads=4, env="rnaforge-quant-prok"):
        Path(out_dir).mkdir(parents=True, exist_ok=True)
        genes = [f"g{i}" for i in range(n_genes)]
        counts = {str(b): [10 + i for i in range(n_genes)] for b in bams}
        rates = {str(b): rate for b in bams}
        return FeatureCountsResult(gene_ids=genes, counts=counts, assignment_rates=rates)
    monkeypatch.setattr(m05_counts, "run_featurecounts", fake_run)


def test_run_counts_requires_m04_done(tmp_path, monkeypatch):
    _fake_featurecounts(monkeypatch)
    config_path, metadata_path = _setup(tmp_path)
    with pytest.raises(ValueError, match="m04"):
        run_counts(load_config(config_path), metadata_path, tmp_path / "run")


def test_run_counts_writes_matrix_and_passes(tmp_path, monkeypatch):
    config_path, metadata_path = _setup(tmp_path)
    run_dir = tmp_path / "run"
    _prep_through_m04(config_path, metadata_path, run_dir, monkeypatch)
    _fake_featurecounts(monkeypatch, rate=0.9, n_genes=3)
    summary = run_counts(load_config(config_path), metadata_path, run_dir)

    assert summary["n_samples"] == 4
    assert summary["read_type"] == "short"
    assert summary["n_genes"] == 3
    matrix = (run_dir / "quantification" / "counts.tsv").read_text().splitlines()
    assert matrix[0] == "gene\ts1\ts2\ts3\ts4"        # sample_id basliklari (BAM yollari degil)
    assert matrix[1].split("\t")[0] == "g0"
    gates = json.loads((run_dir / "quality" / "gates.json").read_text())["gates"]
    assert any(g["module"] == "m05_counts" and g["status"] == "PASS" for g in gates)
    assert any(g["module"] == "m04_quant" for g in gates)   # onceki kapilar korundu


def test_run_counts_empty_matrix_raises(tmp_path, monkeypatch):
    config_path, metadata_path = _setup(tmp_path)
    run_dir = tmp_path / "run"
    _prep_through_m04(config_path, metadata_path, run_dir, monkeypatch)
    _fake_featurecounts(monkeypatch, n_genes=0)   # yanlis feature_type senaryosu
    with pytest.raises(ValueError, match="no genes"):
        run_counts(load_config(config_path), metadata_path, run_dir)


def test_run_counts_low_assignment_fails(tmp_path, monkeypatch):
    config_path, metadata_path = _setup(tmp_path)
    run_dir = tmp_path / "run"
    _prep_through_m04(config_path, metadata_path, run_dir, monkeypatch)
    _fake_featurecounts(monkeypatch, rate=0.10)
    with pytest.raises(GateFailure):
        run_counts(load_config(config_path), metadata_path, run_dir)
    gates = json.loads((run_dir / "quality" / "gates.json").read_text())["gates"]
    assert any(g["module"] == "m05_counts" and g["status"] == "FAIL" for g in gates)


def test_run_counts_resumes(tmp_path, monkeypatch):
    config_path, metadata_path = _setup(tmp_path)
    run_dir = tmp_path / "run"
    _prep_through_m04(config_path, metadata_path, run_dir, monkeypatch)
    _fake_featurecounts(monkeypatch, rate=0.9)
    run_counts(load_config(config_path), metadata_path, run_dir)
    calls = []
    monkeypatch.setattr(m05_counts, "run_featurecounts", lambda *a, **k: calls.append(1))
    summary = run_counts(load_config(config_path), metadata_path, run_dir)
    assert summary.get("resumed") is True
    assert calls == []


def _seed_long_m04(tmp_path, platform="ont"):
    """Uzun-okuma run dizini: m01..m04 done + fake BAM'ler (m05 sözleşme yolu)."""
    from rnaforge.state import RunState
    run_dir = tmp_path / "run"
    (run_dir / "statistics").mkdir(parents=True)
    (run_dir / "statistics" / "raw_statistics.json").write_text(
        json.dumps({"platform": platform, "read_type": "long", "chemistry": "cdna"})
    )
    st = RunState(run_dir)
    for m in ("m01_validate", "m03_trim", "m04_quant"):
        st.mark_done(m, [])
    for sid in ("ctrl1", "ctrl2", "trt1", "trt2"):
        d = run_dir / "quantification" / sid; d.mkdir(parents=True)
        (d / "aligned.sorted.bam").write_bytes(b"BAM")
    return run_dir


def _long_config_meta(tmp_path):
    (tmp_path / "ref").mkdir(exist_ok=True)
    (tmp_path / "ref" / "genome.fa").write_text(">c1\n" + "ACGT" * 25 + "\n")
    (tmp_path / "ref" / "genes.gtf").write_text(
        'c1\ts\texon\t1\t80\t.\t+\t.\tgene_id "g1";\n')
    for n in ("c1.fastq", "c2.fastq", "t1.fastq", "t2.fastq"):
        write_fastq(tmp_path / n, 20, 600, "I")
    cfg = tmp_path / "config.yaml"
    cfg.write_text(textwrap.dedent(f"""
        organism: "E. coli"
        organism_type: "prokaryote"
        library:
          chemistry: "cdna"
        reference:
          genome_fasta: "{tmp_path / 'ref' / 'genome.fa'}"
          annotation_gff: "{tmp_path / 'ref' / 'genes.gtf'}"
    """))
    meta = tmp_path / "samples.tsv"
    meta.write_text(
        "sample_id\tcondition\tfastq_1\n"
        "ctrl1\tcontrol\tc1.fastq\n" "ctrl2\tcontrol\tc2.fastq\n"
        "trt1\ttreated\tt1.fastq\n" "trt2\ttreated\tt2.fastq\n"
    )
    return cfg, meta


def _fake_fc_capture(monkeypatch, rate=0.85, n_genes=2):
    seen: dict = {}

    def fake_run(bams, gff, out_dir, feature_type, attribute, paired=False,
                 threads=4, env="rnaforge-quant-prok", long_read=False):
        seen["long_read"] = long_read
        seen["paired"] = paired
        Path(out_dir).mkdir(parents=True, exist_ok=True)
        genes = [f"g{i}" for i in range(n_genes)]
        counts = {str(b): [10 + i for i in range(n_genes)] for b in bams}
        rates = {str(b): rate for b in bams}
        return FeatureCountsResult(gene_ids=genes, counts=counts, assignment_rates=rates)

    monkeypatch.setattr(m05_counts, "run_featurecounts", fake_run)
    return seen


def test_run_counts_long_dispatches_L_diagnostic(tmp_path, monkeypatch):
    cfg, meta = _long_config_meta(tmp_path)
    run_dir = _seed_long_m04(tmp_path, "ont")
    seen = _fake_fc_capture(monkeypatch, rate=0.85, n_genes=2)
    summary = run_counts(load_config(cfg), meta, run_dir)
    assert seen["long_read"] is True          # -L long-read mode
    assert seen["paired"] is False            # long reads single-molecule
    assert summary["read_type"] == "long"
    assert summary["n_genes"] == 2
    matrix = (run_dir / "quantification" / "counts.tsv").read_text().splitlines()
    assert matrix[0] == "gene\tctrl1\tctrl2\ttrt1\ttrt2"   # common m06 contract
    # Step 6: uzun-okuma assignment WARN kapısı (prokaryote_long, eşik 0.05); asla FAIL.
    # 0.85 > 0.05 → PASS.
    gates = json.loads((run_dir / "quality" / "gates.json").read_text())["gates"]
    m05_gates = [g for g in gates if g["module"] == "m05_counts"]
    assert any(g["name"] == "assignment_rate" for g in m05_gates)
    assert all(g["status"] != "FAIL" for g in m05_gates)   # long assignment asla FAIL değil
    stats = json.loads((run_dir / "statistics" / "count_statistics.json").read_text())
    assert stats["read_type"] == "long"


def test_run_counts_long_low_assignment_warns_not_fails(tmp_path, monkeypatch):
    """ONT'de çok düşük atama ŞÜPHELİ (WARN) — koşuyu geçersiz kılmaz (FAIL değil)."""
    cfg, meta = _long_config_meta(tmp_path)
    run_dir = _seed_long_m04(tmp_path, "ont")
    _fake_fc_capture(monkeypatch, rate=0.02, n_genes=2)   # 0.02 < 0.05 long floor
    summary = run_counts(load_config(cfg), meta, run_dir)   # must NOT raise
    assert summary["read_type"] == "long"
    gates = json.loads((run_dir / "quality" / "gates.json").read_text())["gates"]
    assert any(g["module"] == "m05_counts" and g["status"] == "WARN" for g in gates)
    assert not any(g["module"] == "m05_counts" and g["status"] == "FAIL" for g in gates)


def test_run_counts_long_empty_matrix_raises(tmp_path, monkeypatch):
    cfg, meta = _long_config_meta(tmp_path)
    run_dir = _seed_long_m04(tmp_path, "ont")
    _fake_fc_capture(monkeypatch, n_genes=0)
    with pytest.raises(ValueError, match="no genes"):
        run_counts(load_config(cfg), meta, run_dir)


def test_cli_counts_returns_zero_and_prints_verdict(tmp_path, monkeypatch, capsys):
    from rnaforge.cli import main
    from rnaforge.modules import m03_trim, m04_quant
    from rnaforge.modules.m03_trim import trimmed_name
    from rnaforge.fastp import FastpResult
    from rnaforge.bowtie2 import AlignmentResult
    config_path, metadata_path = _setup(tmp_path)
    common = ["--config", str(config_path), "--metadata", str(metadata_path),
              "--runs-dir", str(tmp_path / "runs"), "--run-id", "demo"]

    def fake_fastp(fastq_1, out_dir, min_length, fastq_2=None, aggressive_quality=False, env="rnaforge-qc"):
        out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
        out1 = out_dir / trimmed_name(Path(fastq_1)); out1.write_text("@r\nACGT\n+\nIIII\n")
        (out_dir / "fastp.json").write_text("{}")
        return FastpResult(200, 196, 0.98, out1=out1)
    monkeypatch.setattr(m03_trim, "run_fastp", fake_fastp)
    monkeypatch.setattr(m04_quant, "build_index", lambda g, i, env="rnaforge-quant-prok": Path(i) / "genome")
    def fake_align(index_prefix, out_dir, fastq_1, fastq_2=None, threads=4, env="rnaforge-quant-prok"):
        out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
        bam = out_dir / "aligned.sorted.bam"; bam.write_bytes(b"BAM")
        return AlignmentResult(bam=bam, alignment_rate=0.95)
    monkeypatch.setattr(m04_quant, "run_bowtie2", fake_align)
    _fake_featurecounts(monkeypatch, rate=0.9, n_genes=3)

    assert main(["validate", *common]) == 0
    assert main(["trim", *common]) == 0
    assert main(["quant", *common]) == 0
    capsys.readouterr()
    assert main(["counts", *common]) == 0
    assert "quality verdict" in capsys.readouterr().out


def test_run_counts_eukaryote_writes_counts_tsv(tmp_path, monkeypatch):
    import rnaforge.modules.m05_counts as m05
    from rnaforge.tximport import TximportResult
    from rnaforge.state import RunState
    run_dir = tmp_path / "run"; (run_dir / "statistics").mkdir(parents=True)
    (run_dir / "statistics" / "raw_statistics.json").write_text(
        '{"read_type":"short","platform":"illumina"}')
    st = RunState(run_dir); st.mark_done("m04_quant", [])
    qd = run_dir / "quantification"
    for sid in ("s1", "s2"):
        (qd / sid).mkdir(parents=True); (qd / sid / "quant.sf").write_text("Name\n")
    a = tmp_path / "a.fq"; a.write_text("@r\nACGT\n+\nIIII\n")
    b = tmp_path / "b.fq"; b.write_text("@r\nACGT\n+\nIIII\n")
    meta = tmp_path / "m.tsv"
    meta.write_text(f"sample_id\tcondition\tfastq_1\ns1\tctrl\t{a}\ns2\ttrt\t{b}\n")
    monkeypatch.setattr(m05, "run_tximport",
        lambda *a, **k: TximportResult(gene_ids=["g1", "g2"],
                                       counts={"s1": [10.0, 0.0], "s2": [20.0, 5.0]}))
    from rnaforge.config import (Config, Reference, Library, Trimming, DE, Report, Resources)
    cfg = Config(organism="human", organism_type="eukaryote", platform="illumina",
        reference=Reference(transcriptome_fasta=tmp_path / "tx.fa", tx2gene=tmp_path / "t2g.tsv"),
        library=Library(), trimming=Trimming(), de=DE(), report=Report(), resources=Resources())
    summary = m05.run_counts(cfg, meta, run_dir)
    matrix = (qd / "counts.tsv").read_text().splitlines()
    assert matrix[0] == "gene\ts1\ts2"
    assert matrix[1] == "g1\t10\t20"
    assert summary["organism_type"] == "eukaryote"
    assert summary["n_genes"] == 2


def test_run_counts_eukaryote_long_aggregates_tx_to_gene(tmp_path, monkeypatch):
    import rnaforge.modules.m05_counts as m05
    from rnaforge.state import RunState
    run_dir = tmp_path / "run"; (run_dir / "statistics").mkdir(parents=True)
    (run_dir / "statistics" / "raw_statistics.json").write_text('{"read_type":"long","platform":"ont"}')
    st = RunState(run_dir); st.mark_done("m04_quant", [])
    qd = run_dir / "quantification"
    for sid in ("s1", "s2"):
        (qd / sid).mkdir(parents=True); (qd / sid / "aligned.sorted.bam").write_text("")
    t2g = tmp_path / "t2g.tsv"; t2g.write_text("ENST1.1\tENSG1\nENST2.2\tENSG1\nENST3.1\tENSG2\n")
    a = tmp_path / "a.fq"; a.write_text("@r\nA\n+\nI\n"); b = tmp_path / "b.fq"; b.write_text("@r\nA\n+\nI\n")
    meta = tmp_path / "m.tsv"; meta.write_text(f"sample_id\tcondition\tfastq_1\ns1\tctrl\t{a}\ns2\ttrt\t{b}\n")
    cnts = {"s1": {"ENST1.1": 3, "ENST3.1": 1}, "s2": {"ENST2.2": 5}}
    monkeypatch.setattr(m05, "count_primary_alignments",
                        lambda bam, **k: cnts["s1" if "s1" in str(bam) else "s2"])
    # İzoform (NanoCount) best-effort; bu test gen-yoluna odaklı → NanoCount'u atlat.
    monkeypatch.setattr(m05, "run_nanocount",
                        lambda bam, out, **k: (_ for _ in ()).throw(RuntimeError("no NanoCount")))
    from rnaforge.config import (Config, Reference, Library, Trimming, DE, Report, Resources)
    cfg = Config(organism="human", organism_type="eukaryote", platform="ont",
        reference=Reference(transcriptome_fasta=tmp_path / "tx.fa", tx2gene=t2g),
        library=Library(), trimming=Trimming(), de=DE(), report=Report(), resources=Resources())
    summary = m05.run_counts(cfg, meta, run_dir)
    lines = (qd / "counts.tsv").read_text().splitlines()
    assert lines[0] == "gene\ts1\ts2"
    d = {l.split("\t")[0]: l.split("\t")[1:] for l in lines[1:]}
    assert d["ENSG1"] == ["3", "5"] and d["ENSG2"] == ["1", "0"]
    assert summary["read_type"] == "long" and summary["organism_type"] == "eukaryote"
    assert summary["n_genes"] == 2
    # best-effort: NanoCount yoksa izoform matrisi yok, gen-düzeyi bozulmaz.
    assert not (qd / "counts_transcript.tsv").exists()
    assert "n_transcripts" not in summary


def _setup_euk_long(tmp_path, monkeypatch):
    """Ökaryot uzun-okuma euk-long ortak kurulum (gen sayımı mock'lu)."""
    import rnaforge.modules.m05_counts as m05
    from rnaforge.state import RunState
    run_dir = tmp_path / "run"; (run_dir / "statistics").mkdir(parents=True)
    (run_dir / "statistics" / "raw_statistics.json").write_text('{"read_type":"long","platform":"ont"}')
    RunState(run_dir).mark_done("m04_quant", [])
    qd = run_dir / "quantification"
    for sid in ("s1", "s2"):
        (qd / sid).mkdir(parents=True); (qd / sid / "aligned.sorted.bam").write_text("")
    t2g = tmp_path / "t2g.tsv"; t2g.write_text("ENST1.1\tENSG1\nENST2.2\tENSG1\nENST3.1\tENSG2\n")
    a = tmp_path / "a.fq"; a.write_text("@r\nA\n+\nI\n"); b = tmp_path / "b.fq"; b.write_text("@r\nA\n+\nI\n")
    meta = tmp_path / "m.tsv"; meta.write_text(f"sample_id\tcondition\tfastq_1\ns1\tctrl\t{a}\ns2\ttrt\t{b}\n")
    monkeypatch.setattr(m05, "count_primary_alignments",
                        lambda bam, **k: ({"ENST1.1": 3, "ENST3.1": 1} if "s1" in str(bam) else {"ENST2.2": 5}))
    from rnaforge.config import (Config, Reference, Library, Trimming, DE, Report, Resources)
    cfg = Config(organism="human", organism_type="eukaryote", platform="ont",
        reference=Reference(transcriptome_fasta=tmp_path / "tx.fa", tx2gene=t2g),
        library=Library(), trimming=Trimming(), de=DE(), report=Report(), resources=Resources())
    return m05, cfg, meta, run_dir, qd


def test_run_counts_eukaryote_long_writes_transcript_matrix(tmp_path, monkeypatch):
    """İzoform: NanoCount est_count → counts_transcript.tsv (yuvarlı); gen counts.tsv DEĞİŞMEZ."""
    m05, cfg, meta, run_dir, qd = _setup_euk_long(tmp_path, monkeypatch)
    nc = {"s1": {"ENST1.1": 3.4, "ENST3.1": 1.0}, "s2": {"ENST2.2": 5.6}}
    monkeypatch.setattr(m05, "run_nanocount",
                        lambda bam, out, **k: nc["s1" if "s1" in str(bam) else "s2"])
    summary = m05.run_counts(cfg, meta, run_dir)
    # gen matrisi (primer-sayım) değişmedi
    assert (qd / "counts.tsv").read_text().splitlines()[0] == "gene\ts1\ts2"
    assert summary["n_genes"] == 2
    # YENİ izoform matrisi
    tlines = (qd / "counts_transcript.tsv").read_text().splitlines()
    assert tlines[0] == "transcript\ts1\ts2"
    td = {l.split("\t")[0]: l.split("\t")[1:] for l in tlines[1:]}
    assert td["ENST1.1"] == ["3", "0"]   # round(3.4)=3, s2'de yok=0
    assert td["ENST2.2"] == ["0", "6"]   # round(5.6)=6
    assert td["ENST3.1"] == ["1", "0"]
    assert summary["n_transcripts"] == 3
