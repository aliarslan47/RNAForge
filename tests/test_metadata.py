from __future__ import annotations

import pytest

from rnaforge.metadata import (
    MetadataError,
    design_variables,
    load_metadata,
    validate_design,
)


def _make_fastqs(tmp_path, *names):
    for n in names:
        (tmp_path / n).write_text("@r\nACGT\n+\nIIII\n")


def _write_meta(tmp_path, body: str):
    path = tmp_path / "samples.tsv"
    path.write_text(body)
    return path


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


def test_design_requires_two_condition_levels(tmp_path):
    _make_fastqs(tmp_path, "a.fastq", "b.fastq")
    path = _write_meta(tmp_path, (
        "sample_id\tcondition\tfastq_1\n"
        "s1\tcontrol\ta.fastq\n"
        "s2\tcontrol\tb.fastq\n"
    ))
    with pytest.raises(MetadataError, match="at least 2 levels"):
        validate_design(load_metadata(path), "~condition")


def test_design_requires_replicates(tmp_path):
    """DESeq2 replika olmadan dispersiyon tahmin edemez — erken ve net uyar."""
    _make_fastqs(tmp_path, "a.fastq", "b.fastq")
    path = _write_meta(tmp_path, (
        "sample_id\tcondition\tfastq_1\n"
        "s1\tcontrol\ta.fastq\n"
        "s2\ttreated\tb.fastq\n"
    ))
    with pytest.raises(MetadataError, match="replicate"):
        validate_design(load_metadata(path), "~condition")


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
