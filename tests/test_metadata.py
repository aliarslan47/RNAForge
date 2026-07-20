from __future__ import annotations

import pytest

from rnaforge.gates import FAIL, PASS
from rnaforge.metadata import (
    MetadataError,
    design_variables,
    load_metadata,
    looks_paired,
    validate_design,
)


def _make_fastqs(tmp_path, *names):
    for n in names:
        (tmp_path / n).write_text("@r\nACGT\n+\nIIII\n")


def _write_meta(tmp_path, body: str):
    path = tmp_path / "samples.tsv"
    path.write_text(body)
    return path


def _gate(gates, name):
    matching = [g for g in gates if g.name == name]
    assert matching, f"gate {name} not reported; got {[g.name for g in gates]}"
    return matching[0]


def test_loads_paired_end_samples(tmp_path):
    _make_fastqs(tmp_path, "a_R1.fastq", "a_R2.fastq", "b_R1.fastq", "b_R2.fastq")
    path = _write_meta(tmp_path, (
        "sample_id\tcondition\tfastq_1\tfastq_2\n"
        "s1\tcontrol\ta_R1.fastq\ta_R2.fastq\n"
        "s2\ttreated\tb_R1.fastq\tb_R2.fastq\n"
    ))
    samples = load_metadata(path)
    assert [s.sample_id for s in samples] == ["s1", "s2"]
    assert samples[0].fastq_2.name == "a_R2.fastq"
    assert samples[0].batch is None


def test_loads_single_end_samples(tmp_path):
    _make_fastqs(tmp_path, "a.fastq", "b.fastq")
    path = _write_meta(tmp_path, (
        "sample_id\tcondition\tfastq_1\n"
        "s1\tcontrol\ta.fastq\n"
        "s2\ttreated\tb.fastq\n"
    ))
    assert load_metadata(path)[0].fastq_2 is None


def test_missing_fastq_file_raises(tmp_path):
    path = _write_meta(tmp_path, (
        "sample_id\tcondition\tfastq_1\n"
        "s1\tcontrol\tyok.fastq\n"
    ))
    with pytest.raises(MetadataError, match="yok.fastq"):
        load_metadata(path)


def test_duplicate_sample_id_raises(tmp_path):
    _make_fastqs(tmp_path, "a.fastq", "b.fastq")
    path = _write_meta(tmp_path, (
        "sample_id\tcondition\tfastq_1\n"
        "s1\tcontrol\ta.fastq\n"
        "s1\ttreated\tb.fastq\n"
    ))
    with pytest.raises(MetadataError, match="s1"):
        load_metadata(path)


def test_missing_required_column_raises(tmp_path):
    path = _write_meta(tmp_path, "sample_id\tfastq_1\ns1\ta.fastq\n")
    with pytest.raises(MetadataError, match="condition"):
        load_metadata(path)


def test_design_variables_parses_formula():
    assert design_variables("~condition") == ["condition"]
    assert design_variables("~batch + condition") == ["batch", "condition"]


def test_design_variable_missing_from_metadata_raises(tmp_path):
    _make_fastqs(tmp_path, "a.fastq", "b.fastq")
    path = _write_meta(tmp_path, (
        "sample_id\tcondition\tfastq_1\n"
        "s1\tcontrol\ta.fastq\n"
        "s2\ttreated\tb.fastq\n"
    ))
    samples = load_metadata(path)
    with pytest.raises(MetadataError, match="batch"):
        validate_design(samples, "~batch + condition")


def test_single_condition_level_fails_the_replication_gate(tmp_path):
    _make_fastqs(tmp_path, "a.fastq", "b.fastq")
    path = _write_meta(tmp_path, (
        "sample_id\tcondition\tfastq_1\n"
        "s1\tcontrol\ta.fastq\n"
        "s2\tcontrol\tb.fastq\n"
    ))
    gates = validate_design(load_metadata(path), "~condition")
    assert _gate(gates, "replication").status == FAIL


def test_condition_without_replicate_fails_the_replication_gate(tmp_path):
    _make_fastqs(tmp_path, "a.fastq", "b.fastq", "c.fastq")
    path = _write_meta(tmp_path, (
        "sample_id\tcondition\tfastq_1\n"
        "s1\tcontrol\ta.fastq\n"
        "s2\tcontrol\tb.fastq\n"
        "s3\ttreated\tc.fastq\n"
    ))
    gate = _gate(validate_design(load_metadata(path), "~condition"), "replication")
    assert gate.status == FAIL
    assert "treated" in gate.message


def test_malformed_formula_still_raises(tmp_path):
    """Bozuk FORMUL kapi degil, gecersiz girdidir — MetadataError kalir."""
    _make_fastqs(tmp_path, "a.fastq", "b.fastq")
    path = _write_meta(tmp_path, (
        "sample_id\tcondition\tfastq_1\n"
        "s1\tcontrol\ta.fastq\n"
        "s2\ttreated\tb.fastq\n"
    ))
    samples = load_metadata(path)
    with pytest.raises(MetadataError, match="no variables"):
        validate_design(samples, "~")
    with pytest.raises(MetadataError, match="unknown variable"):
        validate_design(samples, "~temperature")


def test_valid_design_passes(tmp_path):
    _make_fastqs(tmp_path, "a.fastq", "b.fastq", "c.fastq", "d.fastq")
    path = _write_meta(tmp_path, (
        "sample_id\tcondition\tfastq_1\n"
        "s1\tcontrol\ta.fastq\n"
        "s2\tcontrol\tb.fastq\n"
        "s3\ttreated\tc.fastq\n"
        "s4\ttreated\td.fastq\n"
    ))
    validate_design(load_metadata(path), "~condition")  # raise etmemeli


def test_empty_fastq_1_raises(tmp_path):
    path = _write_meta(tmp_path, (
        "sample_id\tcondition\tfastq_1\n"
        "s1\tcontrol\t\n"
    ))
    with pytest.raises(MetadataError, match="fastq_1"):
        load_metadata(path)


def test_batch_confounded_with_condition_fails_the_rank_gate(tmp_path):
    """Batch condition'la tam confounded ise etkiler ayristirilamaz. DESeq2 bunu
    kriptik 'not full rank' hatasiyla soyler; biz kapiyi dusurup ne yapilacagini soyleriz."""
    _make_fastqs(tmp_path, "a.fastq", "b.fastq", "c.fastq", "d.fastq")
    path = _write_meta(tmp_path, (
        "sample_id\tcondition\tbatch\tfastq_1\n"
        "s1\tcontrol\tb1\ta.fastq\n"
        "s2\tcontrol\tb1\tb.fastq\n"
        "s3\ttreated\tb2\tc.fastq\n"
        "s4\ttreated\tb2\td.fastq\n"
    ))
    gates = validate_design(load_metadata(path), "~batch + condition")
    gate = _gate(gates, "design_rank")
    assert gate.status == FAIL
    assert "confounded" in gate.message
    assert gate.remedy


def test_single_level_batch_fails_the_rank_gate(tmp_path):
    _make_fastqs(tmp_path, "a.fastq", "b.fastq", "c.fastq", "d.fastq")
    path = _write_meta(tmp_path, (
        "sample_id\tcondition\tbatch\tfastq_1\n"
        "s1\tcontrol\tb1\ta.fastq\n"
        "s2\tcontrol\tb1\tb.fastq\n"
        "s3\ttreated\tb1\tc.fastq\n"
        "s4\ttreated\tb1\td.fastq\n"
    ))
    gates = validate_design(load_metadata(path), "~batch + condition")
    assert _gate(gates, "design_rank").status == FAIL


def test_balanced_batch_design_passes_every_gate(tmp_path):
    """Dengeli tasarim GECMELI — kapi sistemi yanlis pozitif uretirse musteri guvenmez."""
    _make_fastqs(tmp_path, "a.fastq", "b.fastq", "c.fastq", "d.fastq")
    path = _write_meta(tmp_path, (
        "sample_id\tcondition\tbatch\tfastq_1\n"
        "s1\tcontrol\tb1\ta.fastq\n"
        "s2\ttreated\tb1\tb.fastq\n"
        "s3\tcontrol\tb2\tc.fastq\n"
        "s4\ttreated\tb2\td.fastq\n"
    ))
    gates = validate_design(load_metadata(path), "~batch + condition")
    assert all(g.status == PASS for g in gates)


def test_subject_column_is_loaded(tmp_path):
    _make_fastqs(tmp_path, "a.fastq", "b.fastq")
    path = _write_meta(tmp_path, (
        "sample_id\tcondition\tsubject\tfastq_1\n"
        "s1\tbefore\tp1\ta.fastq\n"
        "s2\tafter\tp1\tb.fastq\n"
    ))
    samples = load_metadata(path)
    assert [s.subject for s in samples] == ["p1", "p1"]


def test_subject_is_none_when_column_absent(tmp_path):
    _make_fastqs(tmp_path, "a.fastq", "b.fastq")
    path = _write_meta(tmp_path, (
        "sample_id\tcondition\tfastq_1\n"
        "s1\tcontrol\ta.fastq\n"
        "s2\ttreated\tb.fastq\n"
    ))
    assert all(s.subject is None for s in load_metadata(path))


def test_looks_paired_detects_repeated_subject_across_conditions(tmp_path):
    _make_fastqs(tmp_path, "a.fastq", "b.fastq", "c.fastq", "d.fastq")
    path = _write_meta(tmp_path, (
        "sample_id\tcondition\tsubject\tfastq_1\n"
        "s1\tbefore\tp1\ta.fastq\n"
        "s2\tafter\tp1\tb.fastq\n"
        "s3\tbefore\tp2\tc.fastq\n"
        "s4\tafter\tp2\td.fastq\n"
    ))
    assert looks_paired(load_metadata(path)) is True


def test_looks_paired_false_when_each_subject_has_one_condition(tmp_path):
    """Her subject tek condition'da -> eslesmis degil, sadece etiketlenmis."""
    _make_fastqs(tmp_path, "a.fastq", "b.fastq", "c.fastq", "d.fastq")
    path = _write_meta(tmp_path, (
        "sample_id\tcondition\tsubject\tfastq_1\n"
        "s1\tcontrol\tp1\ta.fastq\n"
        "s2\tcontrol\tp2\tb.fastq\n"
        "s3\ttreated\tp3\tc.fastq\n"
        "s4\ttreated\tp4\td.fastq\n"
    ))
    assert looks_paired(load_metadata(path)) is False


def test_looks_paired_ignores_empty_subject_values(tmp_path):
    """Boş veya whitespace-only subject değerleri None'a dönüşür ve yoksayılır.
    Detector sadece paired measurements görmez ve False döner."""
    _make_fastqs(tmp_path, "a.fastq", "b.fastq", "c.fastq")
    path = _write_meta(tmp_path, (
        "sample_id\tcondition\tsubject\tfastq_1\n"
        "s1\tcontrol\t\ta.fastq\n"
        "s2\ttreated\t   \tb.fastq\n"
        "s3\tbefore\tp1\tc.fastq\n"
    ))
    assert looks_paired(load_metadata(path)) is False


def test_looks_paired_same_subject_same_condition_is_not_pairing(tmp_path):
    """Aynı subject aynı condition'da iki kez görülüyorsa bu tekrarlı ölçümdür,
    pairing değil. Gelecekteki optimizasyon (set -> list/counter) bunu kırabileceğinden,
    bu davranışı açık ve kasten test ediyoruz."""
    _make_fastqs(tmp_path, "a.fastq", "b.fastq", "c.fastq")
    path = _write_meta(tmp_path, (
        "sample_id\tcondition\tsubject\tfastq_1\n"
        "s1\tcontrol\tp1\ta.fastq\n"
        "s2\tcontrol\tp1\tb.fastq\n"
        "s3\ttreated\tp2\tc.fastq\n"
    ))
    assert looks_paired(load_metadata(path)) is False


def test_looks_paired_mixed_with_partial_subject_coverage(tmp_path):
    """Bazı satırlar subject'e sahip değildir, bazısı vardır. Yoksayılan satırlar
    (subject=None) işleme alınmaz; geri kalan satırlar kontrol edilir.
    Bu durumda p1 hem before hem de after'da görünüyor -> True."""
    _make_fastqs(tmp_path, "a.fastq", "b.fastq", "c.fastq", "d.fastq")
    path = _write_meta(tmp_path, (
        "sample_id\tcondition\tsubject\tfastq_1\n"
        "s1\tbefore\tp1\ta.fastq\n"
        "s2\tbefore\t\tb.fastq\n"
        "s3\tafter\tp1\tc.fastq\n"
        "s4\tafter\t\td.fastq\n"
    ))
    assert looks_paired(load_metadata(path)) is True
