from pathlib import Path
import pytest
from rnaforge.config import parse_config, ConfigError, ORGANISM_TYPES


def _base(tmp_path):
    cat = tmp_path / "catalog.fa"; cat.write_text(">g1\nACGT\n")
    ann = tmp_path / "catalog.gff"; ann.write_text("g1\t.\tgene\t1\t4\t.\t+\t.\tID=g1\n")
    return cat, ann


def test_metatranscriptome_is_allowed(tmp_path):
    cat, ann = _base(tmp_path)
    cfg = parse_config({
        "organism": "gut community", "organism_type": "metatranscriptome",
        "platform": "illumina",
        "reference": {"gene_catalog_fasta": str(cat), "catalog_annotation": str(ann)},
        "taxonomy": {"kraken2_db": str(tmp_path), "bracken_read_len": 150},
        "rrna": {"db_fasta": str(cat)},
    })
    assert cfg.organism_type == "metatranscriptome"
    assert cfg.reference.gene_catalog_fasta == cat
    assert cfg.taxonomy.bracken_read_len == 150
    assert cfg.taxonomy.bracken_level == "S"
    assert cfg.rrna.db_fasta == cat


def test_metatranscriptome_requires_catalog(tmp_path):
    with pytest.raises(ConfigError, match="gene_catalog_fasta"):
        parse_config({"organism": "x", "organism_type": "metatranscriptome",
                      "platform": "illumina", "reference": {}})


def test_metatranscriptome_in_organism_types():
    assert "metatranscriptome" in ORGANISM_TYPES
