from __future__ import annotations

import textwrap

import pytest

from rnaforge.config import ConfigError, load_config


def _write(tmp_path, body: str):
    path = tmp_path / "config.yaml"
    path.write_text(textwrap.dedent(body))
    return path


# DİKKAT: girintisiz tanımlı — testler buna satır EKLİYOR (PROK_BODY + '...').
# Girintili olsaydı eklenen girintisiz satır textwrap.dedent'i bozardı ve
# YAML parse hatası alırdık, beklediğimiz ConfigError'ı değil.
PROK_BODY = """organism: "Escherichia coli"
organism_type: "prokaryote"
reference:
  genome_fasta: "ref/genome.fa"
  annotation_gff: "ref/genes.gff"
de:
  design: "~condition"
"""


def test_valid_prokaryote_config_loads(tmp_path):
    cfg = load_config(_write(tmp_path, PROK_BODY))
    assert cfg.organism_type == "prokaryote"
    assert cfg.reference.genome_fasta.name == "genome.fa"


def test_missing_organism_type_raises(tmp_path):
    path = _write(tmp_path, """
        organism: "Escherichia coli"
        reference:
          genome_fasta: "ref/genome.fa"
          annotation_gff: "ref/genes.gff"
    """)
    with pytest.raises(ConfigError, match="organism_type"):
        load_config(path)


def test_invalid_organism_type_raises(tmp_path):
    path = _write(tmp_path, """
        organism: "X"
        organism_type: "virus"
        reference:
          genome_fasta: "a"
          annotation_gff: "b"
    """)
    with pytest.raises(ConfigError, match="prokaryote"):
        load_config(path)


def test_prokaryote_requires_genome_and_annotation(tmp_path):
    path = _write(tmp_path, """
        organism: "Escherichia coli"
        organism_type: "prokaryote"
        reference:
          transcriptome_fasta: "ref/tx.fa"
    """)
    with pytest.raises(ConfigError, match="genome_fasta"):
        load_config(path)


def test_eukaryote_requires_transcriptome_and_tx2gene(tmp_path):
    path = _write(tmp_path, """
        organism: "Homo sapiens"
        organism_type: "eukaryote"
        reference:
          genome_fasta: "ref/genome.fa"
    """)
    with pytest.raises(ConfigError, match="transcriptome_fasta"):
        load_config(path)


def test_trimming_defaults_are_gentle(tmp_path):
    """PLAN §4.2: agresif trimming ekspresyon tahminlerini bozar (Williams 2016).
    Varsayılan NAZİK olmalı; bu test o kararı sabitler."""
    cfg = load_config(_write(tmp_path, PROK_BODY))
    assert cfg.trimming.aggressive_quality is False
    assert cfg.trimming.min_length >= 1


def test_platform_defaults_to_auto(tmp_path):
    assert load_config(_write(tmp_path, PROK_BODY)).platform == "auto"


def test_invalid_platform_raises(tmp_path):
    path = _write(tmp_path, PROK_BODY + '\nplatform: "ont"\n')
    with pytest.raises(ConfigError, match="platform"):
        load_config(path)


def test_invalid_strandedness_raises(tmp_path):
    path = _write(tmp_path, PROK_BODY + '\nlibrary:\n  strandedness: "sideways"\n')
    with pytest.raises(ConfigError, match="strandedness"):
        load_config(path)
