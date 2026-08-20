"""Task 10 — rapor: metatranskriptom taksonomi/topluluk-kompozisyonu bölümü + rozet +
yöntem/atıf. section_taxonomy artık taxonomy/abundance_matrix.tsv'den top-N takson
tablosu üretir; render_report yalnız organism_type=metatranscriptome'da bu bölümü
ekler (prokaryote raporunda YOK — koşullu doğruluk). Kraken2/Bracken/SortMeRNA
atıfları yalnız bu kolda basılır (kullanılmayan aracı atıflamama dürüstlüğü)."""
from __future__ import annotations

from rnaforge.config import load_config
from rnaforge.report_html import (
    LABELS,
    parse_abundance_matrix,
    render_report,
    section_dataset,
    section_methods,
    section_references,
    section_software,
    section_taxonomy,
)


def _cfg(tmp_path, organism_type="metatranscriptome", lang="tr"):
    lines = [
        'organism: "gut community"',
        f'organism_type: "{organism_type}"',
        'platform: "illumina"',
    ]
    if organism_type == "metatranscriptome":
        cat = tmp_path / "catalog.fa"; cat.write_text(">g1\nACGT\n")
        ann = tmp_path / "catalog.gff"; ann.write_text("g1\t.\tgene\t1\t4\t.\t+\t.\tID=g1\n")
        lines += [
            "reference:", f'  gene_catalog_fasta: "{cat}"', f'  catalog_annotation: "{ann}"',
            "taxonomy:", '  kraken2_db: "/db/kraken2"', "  bracken_read_len: 100",
            '  bracken_level: "S"', '  env: "rnaforge-meta"',
            "rrna:", '  db_fasta: "/db/rrna.fa"', '  env: "rnaforge-seqqc"',
        ]
    else:
        genome = tmp_path / "genome.fa"; genome.write_text(">c1\nACGT\n")
        gff = tmp_path / "genes.gff"; gff.write_text("##gff-version 3\n")
        lines += ["reference:", f'  genome_fasta: "{genome}"', f'  annotation_gff: "{gff}"']
    lines += [
        "de:", '  design: "~condition"', "  fdr_threshold: 0.05", "  log2fc_threshold: 1.0",
        "report:", f"  language: {lang}",
    ]
    p = tmp_path / f"config_{organism_type}_{lang}.yaml"
    p.write_text("\n".join(lines) + "\n")
    return load_config(p)


def _abundance_matrix(tmp_path):
    p = tmp_path / "abundance_matrix.tsv"
    p.write_text(
        "taxon\tc1\tt1\n"
        "Escherichia coli\t0.62\t0.55\n"
        "Bacteroides fragilis\t0.30\t0.10\n"
        "Faecalibacterium prausnitzii\t0.05\t0.30\n"
    )
    return p


def _raw(organism_type="metatranscriptome"):
    return {"organism": "gut community", "platform": "illumina", "read_type": "short",
            "organism_type": organism_type, "design": "~condition",
            "conditions": {"control": 1, "treated": 1}, "samples": []}


# ---------------------------------------------------------------------------
# section_dataset — organism_type/community rozeti
# ---------------------------------------------------------------------------

def test_dataset_shows_metatranscriptome_badge():
    html = section_dataset(_raw("metatranscriptome"), LABELS["tr"])
    assert "Organizma tipi" in html and "Topluluk" in html
    html_en = section_dataset(_raw("metatranscriptome"), LABELS["en"])
    assert "Organism type" in html_en and "Community" in html_en
    assert "metatranscriptome" in html_en.lower() or "Community" in html_en


def test_dataset_no_community_badge_for_prokaryote():
    html = section_dataset(_raw("prokaryote"), LABELS["en"])
    assert "Community" not in html


# ---------------------------------------------------------------------------
# parse_abundance_matrix
# ---------------------------------------------------------------------------

def test_parse_abundance_matrix_shape(tmp_path):
    p = _abundance_matrix(tmp_path)
    sample_ids, rows = parse_abundance_matrix(p)
    assert sample_ids == ["c1", "t1"]
    by_taxon = {r["taxon"]: r for r in rows}
    assert by_taxon["Escherichia coli"]["c1"] == 0.62
    assert by_taxon["Bacteroides fragilis"]["t1"] == 0.10


def test_parse_abundance_matrix_missing_is_empty(tmp_path):
    assert parse_abundance_matrix(tmp_path / "nope.tsv") == ([], [])


# ---------------------------------------------------------------------------
# section_taxonomy
# ---------------------------------------------------------------------------

def test_section_taxonomy_present_with_top_taxa(tmp_path):
    _, rows = parse_abundance_matrix(_abundance_matrix(tmp_path))
    inputs = {"taxonomy_samples": ["c1", "t1"], "taxonomy_rows": rows,
              "rrna_depletion": {"c1": {"depletion_rate": 0.92}, "t1": {"depletion_rate": 0.88}},
              "figures_dir": None}
    html = section_taxonomy(inputs, LABELS["tr"], "tr", top_n=15)
    assert 'id="taxonomy"' in html
    assert "Escherichia coli" in html and "Bacteroides fragilis" in html
    assert "90.0%" in html            # ortalama rRNA depletion (0.92+0.88)/2 = 0.90


def test_section_taxonomy_not_run_note():
    html = section_taxonomy({"taxonomy_rows": None}, LABELS["en"], "en")
    assert "was not run" in html
    assert 'id="taxonomy"' in html


def test_section_taxonomy_figure_failure_never_crashes(tmp_path):
    _, rows = parse_abundance_matrix(_abundance_matrix(tmp_path))
    # figures_dir bilinçli olarak var-olmayan bir dosyaya işaret ediyor gibi davranmaz;
    # gerçek best-effort davranışı, conda/env yoksa da section_taxonomy'nin
    # asla exception fırlatmamasıdır (figür üretimi try/except içinde).
    inputs = {"taxonomy_samples": ["c1", "t1"], "taxonomy_rows": rows,
              "rrna_depletion": None, "figures_dir": tmp_path / "figures"}
    html = section_taxonomy(inputs, LABELS["tr"], "tr")   # raise etmemeli
    assert 'id="taxonomy"' in html


def test_section_taxonomy_rows_sorted_by_mean_fraction(tmp_path):
    # Ortalama fraksiyona göre AZALAN sıra: E.coli(0.585) > Bacteroides(0.20) > Faecali(0.175).
    _, rows = parse_abundance_matrix(_abundance_matrix(tmp_path))
    inputs = {"taxonomy_samples": ["c1", "t1"], "taxonomy_rows": rows, "figures_dir": None}
    html = section_taxonomy(inputs, LABELS["tr"], "tr")
    i_ec = html.index("Escherichia coli")
    i_bf = html.index("Bacteroides fragilis")
    i_fp = html.index("Faecalibacterium prausnitzii")
    assert i_ec < i_bf < i_fp


def test_section_taxonomy_truncates_to_top_n(tmp_path):
    # 20 takson, azalan bolluk (Taxon00 en yüksek … Taxon19 en düşük); top_n=15 → yalnız 00..14.
    p = tmp_path / "wide_matrix.tsv"
    lines = ["taxon\tc1\tt1"]
    for i in range(20):
        v = (20 - i) / 100.0            # Taxon00=0.20 … Taxon19=0.01 (hepsi farklı)
        lines.append(f"Taxon{i:02d}\t{v}\t{v}")
    p.write_text("\n".join(lines) + "\n")
    sample_ids, rows = parse_abundance_matrix(p)
    inputs = {"taxonomy_samples": sample_ids, "taxonomy_rows": rows, "figures_dir": None}
    html = section_taxonomy(inputs, LABELS["en"], "en", top_n=15)
    assert "Taxon00" in html and "Taxon14" in html      # top-15 içinde
    assert "Taxon15" not in html and "Taxon19" not in html  # kesildi
    assert html.count("<tr") == 16                      # 1 başlık + 15 veri satırı


def test_section_taxonomy_ran_but_empty_matrix_no_crash():
    # rows=[] ("çalıştı ama boş"): None'dan (çalışmadı) farklı; çökmeden bölüm+not render eder.
    html = section_taxonomy({"taxonomy_rows": [], "taxonomy_samples": ["c1"]}, LABELS["en"], "en")
    assert 'id="taxonomy"' in html
    assert "was not run" not in html                    # "çalışmadı" DEĞİL — çalıştı, veri boş


# ---------------------------------------------------------------------------
# software / methods / references — kraken2/bracken/sortmerna yalnız meta kolunda
# ---------------------------------------------------------------------------

def test_software_meta_lists_kraken2_bracken_sortmerna(tmp_path):
    cfg = _cfg(tmp_path, "metatranscriptome")
    flags = {"short": True, "long": False, "meta": True}
    html = section_software(cfg, LABELS["en"], {}, flags)
    assert "Kraken2" in html and "Bracken" in html and "SortMeRNA" in html


def test_software_prokaryote_omits_kraken2_bracken(tmp_path):
    cfg = _cfg(tmp_path, "prokaryote")
    flags = {"short": True, "long": False, "meta": False}
    html = section_software(cfg, LABELS["en"], {}, flags)
    assert "Kraken2" not in html and "Bracken" not in html


def test_methods_meta_describes_full_flow(tmp_path):
    cfg = _cfg(tmp_path, "metatranscriptome")
    html = section_methods(cfg, LABELS["en"], read_type="short", organism_type="metatranscriptome")
    assert "SortMeRNA" in html
    assert "Kraken2" in html and "Bracken" in html
    assert "Bowtie2" in html and "featureCounts" in html and "DESeq2" in html


def test_references_meta_cites_kraken2_bracken_sortmerna_not_others(tmp_path):
    refs = section_references(LABELS["en"], read_type="short", meta_ran=True)
    assert "10.1093/bioinformatics/bts611" in refs         # SortMeRNA (Kopylova 2012)
    assert "10.1186/s13059-019-1891-0" in refs             # Kraken2 (Wood 2019, Genome Biology)
    assert "10.7717/peerj-cs.104" in refs                  # Bracken (Lu 2017, PeerJ CS)
    refs_no = section_references(LABELS["en"], read_type="short", meta_ran=False)
    assert "10.1093/bioinformatics/bts611" not in refs_no
    assert "10.1186/s13059-019-1891-0" not in refs_no
    assert "10.7717/peerj-cs.104" not in refs_no


# ---------------------------------------------------------------------------
# render_report end-to-end — metatranscriptome vs prokaryote
# ---------------------------------------------------------------------------

def _meta_inputs(tmp_path):
    fig = tmp_path / "figures"; fig.mkdir()
    matrix_path = _abundance_matrix(tmp_path)
    sample_ids, rows = parse_abundance_matrix(matrix_path)
    return {
        "raw": _raw("metatranscriptome"),
        "qc": {"samples": {}}, "trimming": {"samples": {}},
        "alignment": {"samples": {"c1": {"alignment_rate": 0.20}}, "read_type": "short"},
        "count": {"samples": {"c1": {"assignment_rate": 0.10}}, "n_genes": 500},
        "de": {"contrast": "treated vs control", "n_genes": 500, "n_significant": 1,
               "fdr_threshold": 0.05, "log2fc_threshold": 1.0, "min_replicate_correlation": 0.90},
        "confidence": {"verdict": "TRUSTWORTHY", "counts": {"PASS": 5, "WARN": 1, "FAIL": 0},
                       "profile": {"name": "metatranscriptome", "overrides": {},
                                   "permissive": True}, "gates": []},
        "figures": {"figures": []},
        "de_results": [{"gene": "g1", "baseMean": 100.0, "log2FoldChange": 2.5, "padj": 1e-4}],
        "gene_map": {}, "figures_dir": fig,
        "norm_counts": {"g1": {"c1": 10.0, "t1": 40.0}},
        "coldata": [("c1", "control"), ("t1", "treated")],
        "taxonomy_samples": sample_ids, "taxonomy_rows": rows,
        "rrna_depletion": {"c1": {"depletion_rate": 0.91}, "t1": {"depletion_rate": 0.87}},
    }


def _prok_inputs(tmp_path):
    fig = tmp_path / "figures"; fig.mkdir()
    return {
        "raw": _raw("prokaryote"),
        "qc": {"samples": {}}, "trimming": {"samples": {}},
        "alignment": {"samples": {"c1": {"alignment_rate": 0.99}}, "read_type": "short"},
        "count": {"samples": {"c1": {"assignment_rate": 0.85}}, "n_genes": 4000},
        "de": {"contrast": "treated vs control", "n_genes": 4000, "n_significant": 1,
               "fdr_threshold": 0.05, "log2fc_threshold": 1.0, "min_replicate_correlation": 0.98},
        "confidence": {"verdict": "TRUSTWORTHY", "counts": {"PASS": 11, "WARN": 0, "FAIL": 0},
                       "profile": {"name": "prokaryote", "overrides": {}}, "gates": []},
        "figures": {"figures": []},
        "de_results": [{"gene": "g1", "baseMean": 100.0, "log2FoldChange": 2.5, "padj": 1e-4}],
        "gene_map": {}, "figures_dir": fig,
        "norm_counts": {"g1": {"c1": 10.0, "t1": 40.0}},
        "coldata": [("c1", "control"), ("t1", "treated")],
    }


def test_render_report_meta_has_taxonomy_and_citations_and_badge(tmp_path):
    cfg = _cfg(tmp_path, "metatranscriptome", "tr")
    doc = render_report(_meta_inputs(tmp_path), cfg, version="0.1.0")
    assert 'id="taxonomy"' in doc
    assert "Escherichia coli" in doc
    assert "10.1093/bioinformatics/bts611" in doc      # SortMeRNA
    assert "10.1186/s13059-019-1891-0" in doc          # Kraken2
    assert "10.7717/peerj-cs.104" in doc               # Bracken
    assert "Topluluk" in doc                            # organism_type rozeti
    assert "metatranscriptome" in doc                   # permissive-profil damgası (confidence)


def test_render_report_prokaryote_has_no_taxonomy_section(tmp_path):
    cfg = _cfg(tmp_path, "prokaryote", "tr")
    doc = render_report(_prok_inputs(tmp_path), cfg, version="0.1.0")
    assert 'id="taxonomy"' not in doc
    assert "10.1186/s13059-019-1891-0" not in doc       # Kraken2 atfı yok
    assert "10.7717/peerj-cs.104" not in doc            # Bracken atfı yok


def test_render_report_meta_has_exactly_one_more_section_than_prokaryote(tmp_path):
    # n_sections doğruluğu (Task 10 review): taksonomi metatranskriptomda TAM +1 bölümdür;
    # diğer tüm opsiyonel bölümler her iki kolda da <section> olarak render edilir.
    md = tmp_path / "m"; md.mkdir()
    pd = tmp_path / "p"; pd.mkdir()
    meta_doc = render_report(_meta_inputs(md), _cfg(tmp_path, "metatranscriptome", "tr"),
                             version="0.1.0")
    prok_doc = render_report(_prok_inputs(pd), _cfg(tmp_path, "prokaryote", "tr"),
                             version="0.1.0")
    assert meta_doc.count("<section") == prok_doc.count("<section") + 1
