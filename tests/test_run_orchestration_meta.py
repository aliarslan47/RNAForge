"""`rnaforge run` orkestrasyonu — metatranscriptome dalı: rrna-deplete + taxonomy
trim'den SONRA, quant'tan ÖNCE koşmalı (m04 _quant_meta'nın m_rrna_deplete state
guard'ını karşılamak için). prokaryot/ökaryot sırası DEĞİŞMEMELİ (Task 9)."""
from __future__ import annotations

import textwrap

import pytest

from rnaforge import cli
from rnaforge.cli import build_run_sequence, main
from rnaforge.gates import FAIL, GateFailure, GateResult


# ---------------------------------------------------------------------------
# build_run_sequence: saf fonksiyon, organism_type parametresi
# ---------------------------------------------------------------------------

def test_metatranscriptome_sequence_inserts_rrna_deplete_and_taxonomy_between_trim_and_quant():
    seq = build_run_sequence(organism_type="metatranscriptome")
    assert seq == [
        "validate", "qc", "trim", "rrna-deplete", "taxonomy",
        "quant", "counts", "de", "figures", "report",
    ]


def test_metatranscriptome_rrna_deplete_taxonomy_between_trim_and_quant_indices():
    seq = build_run_sequence(organism_type="metatranscriptome")
    assert seq.index("trim") < seq.index("rrna-deplete") < seq.index("taxonomy") < seq.index("quant")


@pytest.mark.parametrize("organism_type", [None, "prokaryote", "eukaryote"])
def test_prokaryote_eukaryote_sequence_unchanged(organism_type):
    seq = build_run_sequence(organism_type=organism_type)
    assert seq == ["validate", "qc", "trim", "quant", "counts", "de", "figures", "report"]
    assert "rrna-deplete" not in seq
    assert "taxonomy" not in seq


def test_metatranscriptome_from_to_can_slice_across_new_stages():
    seq = build_run_sequence(start="rrna-deplete", end="quant", organism_type="metatranscriptome")
    assert seq == ["rrna-deplete", "taxonomy", "quant"]


def test_metatranscriptome_from_trim_to_taxonomy():
    seq = build_run_sequence(start="trim", end="taxonomy", organism_type="metatranscriptome")
    assert seq == ["trim", "rrna-deplete", "taxonomy"]


def test_unknown_from_stage_still_rejected_for_metatranscriptome():
    with pytest.raises(ValueError, match="not a core stage"):
        build_run_sequence(start="bogus", organism_type="metatranscriptome")


def test_rrna_deplete_not_a_valid_from_stage_for_prokaryote():
    with pytest.raises(ValueError, match="not a core stage"):
        build_run_sequence(start="rrna-deplete", organism_type="prokaryote")


# ---------------------------------------------------------------------------
# _STAGE_DISPATCH: rrna-deplete/taxonomy artık kayıtlı
# ---------------------------------------------------------------------------

def test_stage_dispatch_has_rrna_deplete_and_taxonomy():
    assert "rrna-deplete" in cli._STAGE_DISPATCH
    assert "taxonomy" in cli._STAGE_DISPATCH
    assert cli._STAGE_DISPATCH["rrna-deplete"] is cli._cmd_rrna_deplete
    assert cli._STAGE_DISPATCH["taxonomy"] is cli._cmd_taxonomy


# ---------------------------------------------------------------------------
# `rnaforge run` uçtan uca: gerçek metatranscriptome config, sahte dispatch
# ---------------------------------------------------------------------------

def _fake_dispatch(monkeypatch, fail_at=None):
    calls = []

    def make(name):
        def fn(args):
            calls.append(name)
            if name == fail_at:
                raise GateFailure([GateResult(
                    name="demo_gate", module=name, status=FAIL,
                    message=f"{name} failed", remedy="fix it")])
            return 0
        return fn

    monkeypatch.setattr(cli, "_STAGE_DISPATCH",
                        {name: make(name) for name in cli._STAGE_DISPATCH})
    return calls


def _write_meta_config(tmp_path):
    cat = tmp_path / "catalog.fa"
    cat.write_text(">g1\nACGT\n")
    ann = tmp_path / "catalog.gff"
    ann.write_text("g1\t.\tgene\t1\t4\t.\t+\t.\tID=g1\n")
    config_path = tmp_path / "c.yaml"
    config_path.write_text(textwrap.dedent(f"""\
        organism: gut community
        organism_type: metatranscriptome
        platform: illumina
        reference:
          gene_catalog_fasta: {cat}
          catalog_annotation: {ann}
        taxonomy:
          kraken2_db: {tmp_path}
          bracken_read_len: 150
        rrna:
          db_fasta: {cat}
        """))
    return config_path


def _run_args(config_path, tmp_path, *extra):
    return ["run", "--config", str(config_path),
            "--metadata", str(tmp_path / "m.tsv"),
            "--runs-dir", str(tmp_path / "runs"), "--run-id", "demo", *extra]


def test_run_executes_metatranscriptome_sequence_in_order(tmp_path, monkeypatch):
    config_path = _write_meta_config(tmp_path)
    calls = _fake_dispatch(monkeypatch)
    assert main(_run_args(config_path, tmp_path)) == 0
    assert calls == [
        "validate", "qc", "trim", "rrna-deplete", "taxonomy",
        "quant", "counts", "de", "figures", "report",
    ]


def test_run_prokaryote_config_sequence_still_unchanged(tmp_path, monkeypatch):
    """Config dosyası hiç yokken (mevcut davranış) sıra prokaryot varsayılanına düşer —
    gerçek hata yine de ilk aşamanın (validate) kendi config yüklemesinde fırlar."""
    calls = _fake_dispatch(monkeypatch)
    args = ["run", "--config", str(tmp_path / "missing.yaml"),
            "--metadata", str(tmp_path / "m.tsv"),
            "--runs-dir", str(tmp_path / "runs"), "--run-id", "demo"]
    assert main(args) == 0
    assert calls == ["validate", "qc", "trim", "quant", "counts", "de", "figures", "report"]


def test_run_stops_on_first_gate_failure_metatranscriptome(tmp_path, monkeypatch):
    config_path = _write_meta_config(tmp_path)
    calls = _fake_dispatch(monkeypatch, fail_at="taxonomy")
    assert main(_run_args(config_path, tmp_path)) == 1
    assert calls == ["validate", "qc", "trim", "rrna-deplete", "taxonomy"]
