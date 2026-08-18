import json
import subprocess
from pathlib import Path

import pytest

from rnaforge import salmon as salmon_mod
from rnaforge.salmon import parse_salmon_meta


def test_parse_salmon_meta_percent_mapped(tmp_path):
    mi = tmp_path / "meta_info.json"
    mi.write_text(json.dumps({"num_processed": 1000, "num_mapped": 853,
                              "percent_mapped": 85.3}))
    assert parse_salmon_meta(mi) == 0.853


def test_build_index_decoy_builds_gentrome_and_decoys(tmp_path, monkeypatch):
    tx = tmp_path / "tx.fa"; tx.write_text(">t1\nACGT\n")
    genome = tmp_path / "genome.fa"; genome.write_text(">chr1\nAAAA\n>chr2\nTTTT\n")
    idx = tmp_path / "idx"
    calls = {}

    def fake_run(cmd, **k):
        calls["cmd"] = cmd
        Path(idx).mkdir(parents=True, exist_ok=True)
        class R: returncode = 0; stdout = ""; stderr = ""
        return R()

    monkeypatch.setattr(salmon_mod.subprocess, "run", fake_run)
    out = salmon_mod.build_salmon_index(tx, idx, genome_fasta=genome)
    assert out == idx
    decoys = (idx.parent / "decoys.txt").read_text().split()
    assert decoys == ["chr1", "chr2"]
    assert "-d" in calls["cmd"]


def test_build_index_no_decoy_transcriptome_only(tmp_path, monkeypatch):
    tx = tmp_path / "tx.fa"; tx.write_text(">t1\nACGT\n")
    idx = tmp_path / "idx"
    calls = {}

    def fake_run(cmd, **k):
        calls["cmd"] = cmd
        Path(idx).mkdir(parents=True, exist_ok=True)
        class R: returncode = 0; stdout = ""; stderr = ""
        return R()

    monkeypatch.setattr(salmon_mod.subprocess, "run", fake_run)
    salmon_mod.build_salmon_index(tx, idx, genome_fasta=None)
    assert "-d" not in calls["cmd"]


def test_run_salmon_quant_parses_mapping_rate(tmp_path, monkeypatch):
    idx = tmp_path / "idx"; idx.mkdir()
    r1 = tmp_path / "s_R1.fastq"; r1.write_text("@r\nACGT\n+\nIIII\n")
    out = tmp_path / "s1"

    def fake_run(cmd, **k):
        assert "-1" not in cmd and "-r" in cmd
        aux = out / "aux_info"; aux.mkdir(parents=True, exist_ok=True)
        (aux / "meta_info.json").write_text('{"percent_mapped": 77.0}')
        (out / "quant.sf").write_text("Name\tLength\tTPM\tNumReads\n")
        class R: returncode = 0; stdout = ""; stderr = ""
        return R()

    monkeypatch.setattr(salmon_mod.subprocess, "run", fake_run)
    q = salmon_mod.run_salmon_quant(idx, out, r1)
    assert q.mapping_rate == 0.77
    assert q.quant_sf == out / "quant.sf"


def test_run_salmon_quant_paired_uses_1_2(tmp_path, monkeypatch):
    idx = tmp_path / "idx"; idx.mkdir()
    r1 = tmp_path / "s_R1.fastq"; r1.write_text("@r\nACGT\n+\nIIII\n")
    r2 = tmp_path / "s_R2.fastq"; r2.write_text("@r\nACGT\n+\nIIII\n")
    out = tmp_path / "s1"

    def fake_run(cmd, **k):
        assert "-1" in cmd and "-2" in cmd and "-r" not in cmd
        aux = out / "aux_info"; aux.mkdir(parents=True, exist_ok=True)
        (aux / "meta_info.json").write_text('{"percent_mapped": 80.0}')
        (out / "quant.sf").write_text("Name\tLength\tTPM\tNumReads\n")
        class R: returncode = 0; stdout = ""; stderr = ""
        return R()

    monkeypatch.setattr(salmon_mod.subprocess, "run", fake_run)
    q = salmon_mod.run_salmon_quant(idx, out, r1, fastq_2=r2)
    assert q.mapping_rate == 0.80


def _has_env(name):
    r = subprocess.run(["conda", "run", "-n", name, "salmon", "--version"],
                       capture_output=True, text=True)
    return r.returncode == 0


@pytest.mark.skipif(not _has_env("rnaforge-quant-euk"), reason="salmon env yok")
def test_salmon_index_and_quant_real(tmp_path):
    tx = tmp_path / "tx.fa"
    tx.write_text(">t1\n" + "ACGT" * 60 + "\n>t2\n" + "GGCC" * 60 + "\n")
    idx = tmp_path / "idx"
    from rnaforge.salmon import build_salmon_index, run_salmon_quant
    build_salmon_index(tx, idx)
    r1 = tmp_path / "s.fastq"
    r1.write_text("".join(f"@r{i}\n{'ACGT'*60}\n+\n{'I'*240}\n" for i in range(50)))
    q = run_salmon_quant(idx, tmp_path / "out", r1)
    assert q.quant_sf.exists() and 0.0 <= q.mapping_rate <= 1.0
