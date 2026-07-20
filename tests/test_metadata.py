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
    gates = validate_design(load_metadata(path), "~condition")
    # Sadece raise etmemesi yetmiyor: yeni sözleşmede all-FAIL de raise etmez,
    # bu yüzden kapıların gerçekten PASS olduğunu doğrulamalıyız.
    assert all(g.status == PASS for g in gates)


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


def test_single_level_batch_rank_fail_names_the_offending_samples(tmp_path):
    """Finding 4: tek-seviye batch FAIL'i mesajda batch'i adlandırıyor ama
    samples alanını boş bırakıyordu; teşhis raporu bu listeyi render eder."""
    _make_fastqs(tmp_path, "a.fastq", "b.fastq", "c.fastq", "d.fastq")
    path = _write_meta(tmp_path, (
        "sample_id\tcondition\tbatch\tfastq_1\n"
        "s1\tcontrol\tb1\ta.fastq\n"
        "s2\tcontrol\tb1\tb.fastq\n"
        "s3\ttreated\tb1\tc.fastq\n"
        "s4\ttreated\tb1\td.fastq\n"
    ))
    gate = _gate(validate_design(load_metadata(path), "~batch + condition"), "design_rank")
    assert gate.status == FAIL
    assert set(gate.samples) == {"s1", "s2", "s3", "s4"}


def test_saturated_unique_subject_fails_the_rank_gate(tmp_path):
    """CRITICAL (Finding 1): 4 örnek, 4 BENZERSİZ subject, 2 condition seviyesi.
    subject örnek sayısı kadar seviyeye sahip -> doygun (saturated) tasarım,
    residual serbestlik derecesi kalmıyor. DESeq2 bunu kriptik 'model matrix is
    not full rank' hatasıyla söyler; eskiden bu senaryo tüm kapılardan PASS
    alıyordu çünkü _rank_gate yalnızca 'batch'e bakıyordu."""
    _make_fastqs(tmp_path, "a.fastq", "b.fastq", "c.fastq", "d.fastq")
    path = _write_meta(tmp_path, (
        "sample_id\tcondition\tsubject\tfastq_1\n"
        "s1\tcontrol\tp1\ta.fastq\n"
        "s2\tcontrol\tp2\tb.fastq\n"
        "s3\ttreated\tp3\tc.fastq\n"
        "s4\ttreated\tp4\td.fastq\n"
    ))
    gate = _gate(validate_design(load_metadata(path), "~subject + condition"), "design_rank")
    assert gate.status == FAIL
    assert "subject" in gate.message
    assert gate.remedy
    assert set(gate.samples) == {"s1", "s2", "s3", "s4"}


def test_genuinely_paired_subject_design_passes_the_rank_gate(tmp_path):
    """Doğru bir eşleşmiş (paired) tasarımda her subject BİRDEN FAZLA condition'da
    görünür, bu yüzden confounded/saturated DEĞİLDİR ve rank gate PASS vermelidir.
    Burada yanlış-pozitif projenin ana kullanım senaryosunu (paired design) kırar."""
    _make_fastqs(tmp_path, "a.fastq", "b.fastq", "c.fastq", "d.fastq")
    path = _write_meta(tmp_path, (
        "sample_id\tcondition\tsubject\tfastq_1\n"
        "s1\tbefore\tp1\ta.fastq\n"
        "s2\tafter\tp1\tb.fastq\n"
        "s3\tbefore\tp2\tc.fastq\n"
        "s4\tafter\tp2\td.fastq\n"
    ))
    gates = validate_design(load_metadata(path), "~subject + condition")
    assert all(g.status == PASS for g in gates)


def test_single_level_subject_fails_the_rank_gate(tmp_path):
    _make_fastqs(tmp_path, "a.fastq", "b.fastq", "c.fastq", "d.fastq")
    path = _write_meta(tmp_path, (
        "sample_id\tcondition\tsubject\tfastq_1\n"
        "s1\tcontrol\tp1\ta.fastq\n"
        "s2\tcontrol\tp1\tb.fastq\n"
        "s3\ttreated\tp1\tc.fastq\n"
        "s4\ttreated\tp1\td.fastq\n"
    ))
    gate = _gate(validate_design(load_metadata(path), "~subject + condition"), "design_rank")
    assert gate.status == FAIL
    assert "subject" in gate.message
    assert set(gate.samples) == {"s1", "s2", "s3", "s4"}


def test_subject_confounded_but_not_saturated_fails_the_rank_gate(tmp_path):
    """6 örnek, 2 subject seviyesi (her biri 3x tekrar), her subject tek bir
    condition'da -> confounded (doygun değil, çünkü 2 seviye != 6 örnek)."""
    _make_fastqs(tmp_path, "a.fastq", "b.fastq", "c.fastq", "d.fastq", "e.fastq", "f.fastq")
    path = _write_meta(tmp_path, (
        "sample_id\tcondition\tsubject\tfastq_1\n"
        "s1\tcontrol\tp1\ta.fastq\n"
        "s2\tcontrol\tp1\tb.fastq\n"
        "s3\tcontrol\tp1\tc.fastq\n"
        "s4\ttreated\tp2\td.fastq\n"
        "s5\ttreated\tp2\te.fastq\n"
        "s6\ttreated\tp2\tf.fastq\n"
    ))
    gate = _gate(validate_design(load_metadata(path), "~subject + condition"), "design_rank")
    assert gate.status == FAIL
    assert "confounded" in gate.message
    assert "subject" in gate.message


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


def _paired_looking_metadata(tmp_path):
    _make_fastqs(tmp_path, "a.fastq", "b.fastq", "c.fastq", "d.fastq")
    return _write_meta(tmp_path, (
        "sample_id\tcondition\tsubject\tfastq_1\n"
        "s1\tbefore\tp1\ta.fastq\n"
        "s2\tafter\tp1\tb.fastq\n"
        "s3\tbefore\tp2\tc.fastq\n"
        "s4\tafter\tp2\td.fastq\n"
    ))


def test_paired_declared_gate_fails_when_pairing_undeclared(tmp_path):
    """Finding 2 + 3: veri PAIRED görünüyor (p1 ve p2 iki condition'da), design
    'subject' kullanmıyor, paired= de belirtilmemiş -> FAIL. Mesaj hangi
    subject(ler)in eşleşmiş göründüğünü ADLANDIRMALI ve samples alanı bu
    subject'lere ait örnek id'leriyle DOLMALI (teşhis raporu bunu render eder)."""
    path = _paired_looking_metadata(tmp_path)
    gate = _gate(validate_design(load_metadata(path), "~condition"), "paired_declared")
    assert gate.status == FAIL
    assert "p1" in gate.message
    assert "p2" in gate.message
    assert gate.remedy
    assert set(gate.samples) == {"s1", "s2", "s3", "s4"}


def test_paired_declared_gate_passes_when_subject_in_design(tmp_path):
    path = _paired_looking_metadata(tmp_path)
    gate = _gate(validate_design(load_metadata(path), "~subject + condition"), "paired_declared")
    assert gate.status == PASS


def test_paired_declared_gate_passes_with_explicit_paired_false(tmp_path):
    """Kullanıcı bilerek unpaired koşmak istiyorsa paired=False deklarasyonu
    kapının FAIL vermesini engellemeli."""
    path = _paired_looking_metadata(tmp_path)
    gate = _gate(
        validate_design(load_metadata(path), "~condition", paired=False), "paired_declared"
    )
    assert gate.status == PASS


def test_paired_declared_gate_passes_when_data_is_not_paired(tmp_path):
    _make_fastqs(tmp_path, "a.fastq", "b.fastq", "c.fastq", "d.fastq")
    path = _write_meta(tmp_path, (
        "sample_id\tcondition\tsubject\tfastq_1\n"
        "s1\tcontrol\tp1\ta.fastq\n"
        "s2\tcontrol\tp2\tb.fastq\n"
        "s3\ttreated\tp3\tc.fastq\n"
        "s4\ttreated\tp4\td.fastq\n"
    ))
    gate = _gate(validate_design(load_metadata(path), "~condition"), "paired_declared")
    assert gate.status == PASS


def test_subject_in_design_missing_for_some_sample_raises(tmp_path):
    """Finding 3: design 'subject' kullanıyor ama bazı örneklerin subject
    değeri yok -> MetadataError (formül hâlâ bozuk girdi, kapı değil)."""
    _make_fastqs(tmp_path, "a.fastq", "b.fastq", "c.fastq")
    path = _write_meta(tmp_path, (
        "sample_id\tcondition\tsubject\tfastq_1\n"
        "s1\tbefore\tp1\ta.fastq\n"
        "s2\tafter\tp1\tb.fastq\n"
        "s3\tbefore\t\tc.fastq\n"
    ))
    samples = load_metadata(path)
    with pytest.raises(MetadataError, match="subject"):
        validate_design(samples, "~subject + condition")
