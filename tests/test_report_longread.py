from __future__ import annotations

import textwrap

from rnaforge.config import load_config
from rnaforge.report_html import (
    LABELS,
    section_dataset,
    section_methods,
    section_references,
    section_software,
)


def _config(tmp_path):
    (tmp_path / "ref").mkdir()
    (tmp_path / "ref" / "genome.fa").write_text(">c1\nACGT\n")
    (tmp_path / "ref" / "genes.gff").write_text("##gff-version 3\n")
    cfg = tmp_path / "config.yaml"
    cfg.write_text(textwrap.dedent(f"""
        organism: "E. coli"
        organism_type: "prokaryote"
        reference:
          genome_fasta: "{tmp_path / 'ref' / 'genome.fa'}"
          annotation_gff: "{tmp_path / 'ref' / 'genes.gff'}"
    """))
    return load_config(cfg)


def _raw(read_type):
    return {"organism": "E. coli", "platform": "ont", "read_type": read_type,
            "design": "~condition", "conditions": {"control": 2, "treated": 2},
            "samples": []}


def test_dataset_shows_long_read_type_badge():
    html = section_dataset(_raw("long"), LABELS["tr"])
    assert "Okuma tipi" in html            # read_type rozeti etiketi
    assert "uzun" in html                  # TR gloss: long → uzun
    html_en = section_dataset(_raw("long"), LABELS["en"])
    assert "Read type" in html_en and "long" in html_en


def test_software_long_lists_minimap2_nanoplot_not_bowtie(tmp_path):
    cfg = _config(tmp_path)
    flags = {"short": False, "long": True}
    html = section_software(cfg, LABELS["en"], {}, flags)
    assert "minimap2" in html
    assert "NanoPlot" in html
    assert "Pychopper" in html
    assert "Bowtie2" not in html           # short-only tool hidden on long runs
    assert "FastQC" not in html


def test_software_short_lists_bowtie2_not_minimap2(tmp_path):
    cfg = _config(tmp_path)
    flags = {"short": True, "long": False}
    html = section_software(cfg, LABELS["en"], {}, flags)
    assert "Bowtie2" in html
    assert "FastQC" in html
    assert "minimap2" not in html
    assert "NanoPlot" not in html


def test_methods_long_describes_minimap2(tmp_path):
    cfg = _config(tmp_path)
    html = section_methods(cfg, LABELS["en"], read_type="long")
    assert "minimap2" in html
    assert "featureCounts" in html
    assert "Bowtie2" not in html


def test_methods_short_describes_bowtie2(tmp_path):
    cfg = _config(tmp_path)
    html = section_methods(cfg, LABELS["en"], read_type="short")
    assert "Bowtie2" in html


def test_references_long_cites_long_read_tools_not_short(tmp_path):
    html = section_references(LABELS["en"], read_type="long")
    assert "minimap2" in html.lower()
    assert "NanoPack" in html                 # NanoPlot/chopper atfı
    assert "Bowtie 2" not in html
    assert "FastQC" not in html


def test_references_short_cites_short_read_tools(tmp_path):
    html = section_references(LABELS["en"], read_type="short")
    assert "FastQC" in html and "Bowtie 2" in html
    assert "minimap2" not in html.lower()
