from __future__ import annotations

import json

import pytest

from rnaforge.gates import (
    FAIL,
    PASS,
    WARN,
    GateFailure,
    GateResult,
    raise_if_failed,
    write_gate_results,
)


def _result(name="alignment_rate", status=PASS, **kw):
    defaults = dict(
        module="m04", status=status, measured=0.9, threshold=0.7,
        message="alignment rate is the share of reads mapped to the reference",
        remedy="check that the reference genome matches the organism",
        overridden=False,
    )
    defaults.update(kw)
    return GateResult(name=name, **defaults)


def test_remedy_must_not_be_empty():
    """Kapı mesajı NE YAPILACAGINI söylemezse müşteriye faydası yok."""
    with pytest.raises(ValueError, match="remedy"):
        _result(remedy="")


def test_raise_if_failed_is_quiet_when_all_pass():
    raise_if_failed([_result(status=PASS), _result(status=WARN)])


def test_raise_if_failed_raises_on_fail():
    failing = _result(status=FAIL, measured=0.42)
    with pytest.raises(GateFailure) as exc:
        raise_if_failed([_result(status=PASS), failing])
    assert exc.value.failures == [failing]
    assert "alignment_rate" in str(exc.value)


def test_write_gate_results_appends_across_modules(tmp_path):
    """Resume ile uyumlu olmali: m02'nin sonuclari m01'inkileri EZMEMELI."""
    write_gate_results(tmp_path, [_result(name="design_rank", module="m01")])
    write_gate_results(tmp_path, [_result(name="alignment_rate", module="m04")])
    data = json.loads((tmp_path / "quality" / "gates.json").read_text())
    assert [g["name"] for g in data["gates"]] == ["design_rank", "alignment_rate"]


def test_rerunning_a_module_replaces_its_own_results(tmp_path):
    """--force ile yeniden kosulan modul kendi eski sonucunu birakmamali."""
    write_gate_results(tmp_path, [_result(name="design_rank", module="m01", status=FAIL)])
    write_gate_results(tmp_path, [_result(name="design_rank", module="m01", status=PASS)])
    data = json.loads((tmp_path / "quality" / "gates.json").read_text())
    assert len(data["gates"]) == 1
    assert data["gates"][0]["status"] == PASS


def test_corrupted_gates_file_is_preserved_not_silently_discarded(tmp_path):
    """Bozuk gates.json sessizce yok sayilmamali: kenara alinmali, yuksek sesle
    bildirilmeli ve yeni modulun yazimi yine de basarili olmali (crash-survival)."""
    write_gate_results(tmp_path, [_result(name="design_rank", module="m01")])
    gates_path = tmp_path / "quality" / "gates.json"
    # Yariya kesilmis (truncated) yazimi simule et: cokme sirasinda olabilecek durum.
    gates_path.write_text('{"gates": [{"name": "design_rank", "module": "m01"')

    with pytest.warns(UserWarning, match="corrupt"):
        result_path = write_gate_results(tmp_path, [_result(name="alignment_rate", module="m04")])

    assert result_path == gates_path
    # Bozuk dosya silinmemis, kenara alinmis olmali (adli/forensic iz).
    corrupt_candidates = list((tmp_path / "quality").glob("gates.json.corrupt*"))
    assert len(corrupt_candidates) == 1
    assert "design_rank" in corrupt_candidates[0].read_text()

    # Yeni yazim basarili: m04 sonucu yerinde, dosya yine gecerli JSON.
    data = json.loads(gates_path.read_text())
    assert [g["name"] for g in data["gates"]] == ["alignment_rate"]

    # Bozuk dosyadan sonra artik yariya kesik gecici (.tmp) dosya kalmamali.
    leftovers = list((tmp_path / "quality").glob("gates.json.tmp*"))
    assert leftovers == []


def test_samples_field_cannot_be_mutated_after_construction():
    """samples frozen dataclass uzerinde list olursa yerinde mutasyona acik kalir;
    tuple olmali ki .samples uzerinde append/remove calismasin."""
    result = _result(samples=["S1", "S2"])
    assert result.samples == ("S1", "S2")
    with pytest.raises(AttributeError):
        result.samples.append("S3")
