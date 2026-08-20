from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from rnaforge.kraken2 import (
    Kraken2ParseError,
    Kraken2RunError,
    parse_bracken,
    parse_kraken2_report,
    run_bracken,
    run_kraken2,
)


def test_parse_kraken2_report(tmp_path):
    """Kraken2 report parse: fraction, clade_reads, taxon_reads, rank, taxid, name."""
    r = tmp_path / "k.report"
    r.write_text(" 50.00\t500\t500\tS\t562\tEscherichia coli\n"
                 " 30.00\t300\t300\tS\t1280\tStaphylococcus aureus\n")
    rows = parse_kraken2_report(r)
    assert rows[0]["name"] == "Escherichia coli"
    assert abs(rows[0]["fraction"] - 0.50) < 1e-6
    assert rows[0]["rank"] == "S"
    assert rows[0]["taxid"] == "562"
    assert rows[0]["reads"] == 500
    assert len(rows) == 2
    assert rows[1]["name"] == "Staphylococcus aureus"


def test_parse_bracken(tmp_path):
    """Bracken output parse: header + name, taxonomy_id, taxonomy_lvl,
    kraken_assigned_reads, added_reads, new_est_reads, fraction_total_reads."""
    b = tmp_path / "b.bracken"
    b.write_text("name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\t"
                 "added_reads\tnew_est_reads\tfraction_total_reads\n"
                 "Escherichia coli\t562\tS\t400\t100\t500\t0.55\n"
                 "Staphylococcus aureus\t1280\tS\t300\t50\t350\t0.45\n")
    d = parse_bracken(b)
    assert abs(d["Escherichia coli"] - 0.55) < 1e-6
    assert abs(d["Staphylococcus aureus"] - 0.45) < 1e-6
    assert len(d) == 2


def test_parse_kraken2_report_rejects_malformed():
    """parse_kraken2_report raises error on missing fields."""
    with pytest.raises(Kraken2ParseError):
        parse_kraken2_report(Path("/nonexistent/file.report"))


def test_parse_bracken_rejects_malformed():
    """parse_bracken raises error on missing file."""
    with pytest.raises(Kraken2ParseError):
        parse_bracken(Path("/nonexistent/file.bracken"))


@pytest.mark.skipif(shutil.which("conda") is None, reason="conda yok")
def test_run_kraken2_requires_db(tmp_path):
    """run_kraken2 raises RunError if db does not exist."""
    reads = tmp_path / "reads.fastq"
    reads.write_text("@r1\nACGT\n+\nIIII\n")
    with pytest.raises(Kraken2RunError, match="database"):
        run_kraken2(reads, tmp_path / "fake_db", tmp_path / "out", paired=False, threads=1)


@pytest.mark.skipif(shutil.which("conda") is None, reason="conda yok")
def test_run_bracken_requires_files(tmp_path):
    """run_bracken raises RunError if input files do not exist."""
    with pytest.raises(Kraken2RunError):
        run_bracken(tmp_path / "fake.report", tmp_path / "fake_db", tmp_path / "out.bracken",
                    read_len=100, level="S")
