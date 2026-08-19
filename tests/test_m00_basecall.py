from __future__ import annotations

import json
import textwrap
from pathlib import Path

from rnaforge.config import load_config
from rnaforge.basecall import basecalled_metadata_path
from rnaforge.metadata import load_metadata
from rnaforge.modules import m00_basecall
from rnaforge.modules.m00_basecall import run_basecall


def _config(tmp_path):
    (tmp_path / "ref").mkdir(exist_ok=True)
    (tmp_path / "ref" / "genome.fa").write_text(">c1\nACGT\n")
    (tmp_path / "ref" / "genes.gff").write_text("##gff-version 3\n")
    cfg = tmp_path / "config.yaml"
    cfg.write_text(textwrap.dedent(f"""
        organism: "E. coli"
        organism_type: "prokaryote"
        library:
          chemistry: "cdna"
        reference:
          genome_fasta: "{tmp_path / 'ref' / 'genome.fa'}"
          annotation_gff: "{tmp_path / 'ref' / 'genes.gff'}"
    """))
    return load_config(cfg)


def _fake_dorado(monkeypatch, reads=42):
    calls = []

    def fake_run_dorado(pod5, out_fastq, dorado_bin="dorado", model="hac",
                        device="cuda:all", models_dir=None):
        calls.append((Path(pod5).name, model, device))
        Path(out_fastq).parent.mkdir(parents=True, exist_ok=True)
        Path(out_fastq).write_text("@r\nACGT\n+\nIIII\n")
        return reads

    monkeypatch.setattr(m00_basecall, "run_dorado", fake_run_dorado)
    return calls


def test_basecall_pod5_sample_produces_fastq_and_resolved_metadata(tmp_path, monkeypatch):
    (tmp_path / "sig").mkdir()
    (tmp_path / "sig" / "reads.pod5").write_bytes(b"\x00")   # signal input
    meta = tmp_path / "samples.tsv"
    meta.write_text("sample_id\tcondition\tfastq_1\n"
                    f"s1\tctrl\t{tmp_path / 'sig'}\n")
    run_dir = tmp_path / "run"
    calls = _fake_dorado(monkeypatch, reads=42)

    summary = run_basecall(_config(tmp_path), meta, run_dir)
    assert len(calls) == 1                              # dorado ran once (pod5)
    assert summary["samples"]["s1"]["reads"] == 42
    assert summary["samples"]["s1"]["input_kind"] == "pod5"
    # resolved metadata repoints fastq_1 to the basecalled FASTQ
    rm = basecalled_metadata_path(run_dir)
    assert rm.exists()
    s = load_metadata(rm)[0]
    assert s.fastq_1.suffix == ".fastq" and s.fastq_1.exists()
    assert s.fastq_1.is_absolute()   # göreli yol run_dir göreliyken ikilenirdi
    # diagnostic: no gate written
    assert not (run_dir / "quality" / "gates.json").exists()


def test_basecall_fastq_sample_is_passthrough(tmp_path, monkeypatch):
    fq = tmp_path / "s1.fastq"; fq.write_text("@r\nACGT\n+\nIIII\n")
    meta = tmp_path / "samples.tsv"
    meta.write_text("sample_id\tcondition\tfastq_1\n" f"s1\tctrl\t{fq}\n")
    run_dir = tmp_path / "run"
    calls = _fake_dorado(monkeypatch)

    summary = run_basecall(_config(tmp_path), meta, run_dir)
    assert calls == []                                 # no basecalling for FASTQ
    assert summary["samples"]["s1"]["input_kind"] == "fastq"
    rm = basecalled_metadata_path(run_dir)
    assert load_metadata(rm)[0].fastq_1 == fq          # unchanged


def test_basecall_fast5_converts_then_basecalls(tmp_path, monkeypatch):
    (tmp_path / "sig").mkdir()
    (tmp_path / "sig" / "reads.fast5").write_bytes(b"\x00")
    meta = tmp_path / "samples.tsv"
    meta.write_text("sample_id\tcondition\tfastq_1\n" f"s1\tctrl\t{tmp_path / 'sig'}\n")
    run_dir = tmp_path / "run"
    calls = _fake_dorado(monkeypatch)
    conv = []

    def fake_convert(fast5, out_pod5, env="rnaforge-basecall"):
        conv.append(Path(fast5).name)
        Path(out_pod5).parent.mkdir(parents=True, exist_ok=True)
        Path(out_pod5).write_bytes(b"\x00")
        return Path(out_pod5)

    monkeypatch.setattr(m00_basecall, "convert_fast5_to_pod5", fake_convert)
    summary = run_basecall(_config(tmp_path), meta, run_dir)
    assert conv and len(calls) == 1                     # converted then basecalled
    assert summary["samples"]["s1"]["input_kind"] == "fast5"


def test_basecall_resumes(tmp_path, monkeypatch):
    (tmp_path / "sig").mkdir()
    (tmp_path / "sig" / "reads.pod5").write_bytes(b"\x00")
    meta = tmp_path / "samples.tsv"
    meta.write_text("sample_id\tcondition\tfastq_1\n" f"s1\tctrl\t{tmp_path / 'sig'}\n")
    run_dir = tmp_path / "run"
    _fake_dorado(monkeypatch)
    run_basecall(_config(tmp_path), meta, run_dir)
    calls2 = _fake_dorado(monkeypatch)
    summary = run_basecall(_config(tmp_path), meta, run_dir)
    assert summary.get("resumed") is True
    assert calls2 == []                                # did not re-basecall


def test_config_basecall_defaults(tmp_path):
    cfg = _config(tmp_path)
    assert cfg.basecall.model == "hac"
    assert cfg.basecall.device == "cuda:all"
    assert cfg.basecall.env == "rnaforge-basecall"


def test_resolved_metadata_preserves_covariates(tmp_path):
    """Faz 3: basecall çözülmüş metadata'yı yeniden yazarken keyfi kovaryatları (sex,
    genotype...) düşürmemeli — yoksa ONT yolunda kovaryat design'ları sessizce bozulur."""
    from rnaforge.metadata import Sample
    (tmp_path / "a.fq").write_text("@r\nACGT\n+\nIIII\n")
    (tmp_path / "b.fq").write_text("@r\nACGT\n+\nIIII\n")
    rows = [
        Sample("s1", "control", tmp_path / "a.fq", None,
               batch="b1", covariates={"sex": "M", "lane": "L1"}),
        Sample("s2", "treated", tmp_path / "b.fq", None,
               batch="b1", covariates={"sex": "F", "lane": "L2"}),
    ]
    out = tmp_path / "resolved.tsv"
    m00_basecall._write_resolved_metadata(out, rows)
    reloaded = load_metadata(out)
    assert [s.covariates.get("sex") for s in reloaded] == ["M", "F"]
    assert [s.covariates.get("lane") for s in reloaded] == ["L1", "L2"]
