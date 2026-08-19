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
    # `>= 1` vacuous olurdu: config.py zaten min_length < 1'de ConfigError atıyor,
    # yani assertion asla düşemezdi. Varsayılanı sabitlemek literatür gerekçesini
    # (Williams 2016) gerçekten korur — biri 36'yı 2 yaparsa test yakalar.
    assert cfg.trimming.min_length == 36


def test_platform_defaults_to_auto(tmp_path):
    assert load_config(_write(tmp_path, PROK_BODY)).platform == "auto"


def test_invalid_platform_raises(tmp_path):
    # ont/pacbio_hifi artık geçerli (uzun-okuma kolu); gerçekten bilinmeyen bir değer:
    path = _write(tmp_path, PROK_BODY + '\nplatform: "solid"\n')
    with pytest.raises(ConfigError, match="platform"):
        load_config(path)


def test_ont_platform_is_valid(tmp_path):
    # Kullanıcı ONT platformunu açıkça bildirebilir (kısa cDNA yanlış-tespitini ezmek için).
    assert load_config(_write(tmp_path, PROK_BODY + '\nplatform: "ont"\n')).platform == "ont"


def test_invalid_strandedness_raises(tmp_path):
    path = _write(tmp_path, PROK_BODY + '\nlibrary:\n  strandedness: "sideways"\n')
    with pytest.raises(ConfigError, match="strandedness"):
        load_config(path)


def test_non_mapping_section_raises_config_error(tmp_path):
    """`library: "foo"` ham AttributeError sızdırmamalı — ConfigError sözleşmesi
    kapıda tutar; ticari üründe traceback kullanıcıya ne yapacağını söylemez."""
    path = _write(tmp_path, PROK_BODY + '\nlibrary: "foo"\n')
    with pytest.raises(ConfigError, match="library must be a mapping"):
        load_config(path)


def test_non_numeric_threads_raises_config_error(tmp_path):
    path = _write(tmp_path, PROK_BODY + '\nresources:\n  threads: "sekiz"\n')
    with pytest.raises(ConfigError, match="resources.threads"):
        load_config(path)


def test_non_numeric_min_length_raises_config_error(tmp_path):
    path = _write(tmp_path, PROK_BODY + '\ntrimming:\n  min_length: "otuzalti"\n')
    with pytest.raises(ConfigError, match="trimming.min_length"):
        load_config(path)


def test_paired_defaults_to_undeclared(tmp_path):
    assert load_config(_write(tmp_path, PROK_BODY)).paired is None


def test_paired_can_be_declared_false(tmp_path):
    cfg = load_config(_write(tmp_path, PROK_BODY + "\npaired: false\n"))
    assert cfg.paired is False


def test_quality_overrides_are_loaded(tmp_path):
    cfg = load_config(_write(tmp_path, PROK_BODY + "\nquality:\n  alignment_rate: 0.4\n"))
    assert cfg.quality == {"alignment_rate": 0.4}


def test_de_reference_defaults_none(tmp_path):
    assert load_config(_write(tmp_path, PROK_BODY)).de.reference is None


def test_de_reference_loaded(tmp_path):
    body = PROK_BODY.replace('de:\n  design: "~condition"\n',
                             'de:\n  design: "~condition"\n  reference: control\n')
    cfg = load_config(_write(tmp_path, body))
    assert cfg.de.reference == "control"
    assert cfg.de.design == "~condition"


def test_de_contrasts_default_empty(tmp_path):
    assert load_config(_write(tmp_path, PROK_BODY)).de.contrasts == ()


def test_de_contrasts_loaded_as_pairs(tmp_path):
    body = PROK_BODY + (
        "  contrasts:\n"
        "    - [high, control]\n"
        "    - [low, control]\n"
    )
    cfg = load_config(_write(tmp_path, body))
    assert cfg.de.contrasts == (("high", "control"), ("low", "control"))


def test_de_contrasts_rejects_non_pair(tmp_path):
    body = PROK_BODY + "  contrasts:\n    - [only_one]\n"
    with pytest.raises(ConfigError, match="contrasts"):
        load_config(_write(tmp_path, body))


def test_de_contrasts_rejects_identical_levels(tmp_path):
    body = PROK_BODY + "  contrasts:\n    - [same, same]\n"
    with pytest.raises(ConfigError, match="identical"):
        load_config(_write(tmp_path, body))


def test_de_contrasts_rejects_delimiter_in_level(tmp_path):
    body = PROK_BODY + "  contrasts:\n    - [\"a:b\", control]\n"
    with pytest.raises(ConfigError, match="';'|':'"):
        load_config(_write(tmp_path, body))


def test_quantification_defaults(tmp_path):
    cfg = load_config(_write(tmp_path, PROK_BODY))
    assert cfg.quantification.feature_type == "exon"
    assert cfg.quantification.attribute == "gene_id"


def test_quantification_overrides(tmp_path):
    cfg = load_config(_write(tmp_path, PROK_BODY +
                             "\nquantification:\n  feature_type: CDS\n  attribute: locus_tag\n"))
    assert cfg.quantification.feature_type == "CDS"
    assert cfg.quantification.attribute == "locus_tag"


def test_quantification_is_known_top_level_key(tmp_path):
    # 'quantification' reddedilmemeli (KNOWN_TOP_LEVEL_KEYS'te)
    load_config(_write(tmp_path, PROK_BODY + "\nquantification:\n  feature_type: gene\n"))


def test_unknown_top_level_key_raises(tmp_path):
    """`design:` üst seviyede SESSİZCE yok sayılıyordu (doğru şema `de.design`).
    Kullanıcı tasarımını değiştirdiğini sanıp eski `~condition` ile koşardı —
    makul görünen SAHTE sonuç. Bilinmeyen üst anahtar reddedilmeli; mesaj
    doğru yeri (`de.design`) göstermeli değil ama en azından anahtarı adlandırmalı."""
    path = _write(tmp_path, PROK_BODY + '\ndesign: "~subject + condition"\n')
    with pytest.raises(ConfigError, match="design"):
        load_config(path)


def test_typo_in_top_level_key_raises(tmp_path):
    """`organismtype` gibi yazım hatası sessizce yutulup `organism_type` zorunlu
    hatası vermek yerine — burada organism_type zaten var, yani yanlış anahtar
    tamamen görünmez olurdu. Reddedilsin ki kullanıcı hatasını görsün."""
    path = _write(tmp_path, PROK_BODY + '\nrefernce:\n  genome_fasta: "x"\n')
    with pytest.raises(ConfigError, match="refernce"):
        load_config(path)


def test_enrichment_defaults(tmp_path):
    cfg = load_config(_write(tmp_path, PROK_BODY))
    assert cfg.enrichment.min_term_size == 3
    assert cfg.enrichment.top_n == 15
    assert cfg.enrichment.obo is None
    assert cfg.enrichment.gaf is None


def test_enrichment_parsed(tmp_path):
    body = PROK_BODY + (
        "\nenrichment:\n  min_term_size: 5\n  top_n: 10\n"
        "  obo: \"references/go/go-basic.obo\"\n  gaf: \"references/ecoli/ecoli.gaf\"\n"
    )
    cfg = load_config(_write(tmp_path, body))
    assert cfg.enrichment.min_term_size == 5
    assert cfg.enrichment.top_n == 10
    assert cfg.enrichment.obo.name == "go-basic.obo"
    assert cfg.enrichment.gaf.name == "ecoli.gaf"


def test_enrichment_bad_int_raises(tmp_path):
    body = PROK_BODY + "\nenrichment:\n  min_term_size: \"three\"\n"
    with pytest.raises(ConfigError, match="enrichment.min_term_size"):
        load_config(_write(tmp_path, body))


def test_enrichment_kegg_fields(tmp_path):
    body = PROK_BODY + "\nenrichment:\n  kegg_organism: eco\n  kegg_dir: references/kegg/eco\n"
    cfg = load_config(_write(tmp_path, body))
    assert cfg.enrichment.kegg_organism == "eco"
    assert cfg.enrichment.kegg_dir.name == "eco"


def test_enrichment_kegg_defaults(tmp_path):
    cfg = load_config(_write(tmp_path, PROK_BODY))
    assert cfg.enrichment.kegg_organism is None and cfg.enrichment.kegg_dir is None


def test_enrichment_gsea_sizes(tmp_path):
    body = PROK_BODY + "\nenrichment:\n  gsea_min_size: 10\n  gsea_max_size: 300\n"
    cfg = load_config(_write(tmp_path, body))
    assert cfg.enrichment.gsea_min_size == 10 and cfg.enrichment.gsea_max_size == 300


def test_enrichment_gsea_defaults(tmp_path):
    cfg = load_config(_write(tmp_path, PROK_BODY))
    assert cfg.enrichment.gsea_min_size == 15 and cfg.enrichment.gsea_max_size == 500


def test_enrichment_revigo_similarity(tmp_path):
    cfg = load_config(_write(tmp_path, PROK_BODY + "\nenrichment:\n  revigo_similarity: 0.5\n"))
    assert cfg.enrichment.revigo_similarity == 0.5


def test_enrichment_revigo_default(tmp_path):
    assert load_config(_write(tmp_path, PROK_BODY)).enrichment.revigo_similarity == 0.7


def test_amr_config_defaults(tmp_path):
    cfg = load_config(_write(tmp_path, PROK_BODY))
    assert cfg.amr.amr_db == "card" and cfg.amr.virulence_db == "vfdb"
    assert cfg.amr.min_identity == 80.0 and cfg.amr.env == "rnaforge-amr"


def test_amr_config_parsed(tmp_path):
    cfg = load_config(_write(tmp_path, PROK_BODY + "\namr:\n  amr_db: ncbi\n  min_identity: 90\n"))
    assert cfg.amr.amr_db == "ncbi" and cfg.amr.min_identity == 90.0


def test_operon_config_default(tmp_path):
    assert load_config(_write(tmp_path, PROK_BODY)).operon.max_gap == 50


def test_operon_config_parsed(tmp_path):
    assert load_config(_write(tmp_path, PROK_BODY + "\noperon:\n  max_gap: 100\n")).operon.max_gap == 100


def test_ppi_config_defaults(tmp_path):
    cfg = load_config(_write(tmp_path, PROK_BODY))
    assert cfg.ppi.taxid is None and cfg.ppi.min_score == 700 and cfg.ppi.min_community_size == 3


def test_ppi_config_parsed(tmp_path):
    cfg = load_config(_write(tmp_path, PROK_BODY + "\nppi:\n  taxid: '511145'\n  min_score: 400\n"))
    assert cfg.ppi.taxid == "511145" and cfg.ppi.min_score == 400


def test_amr_amrfinder_fields(tmp_path):
    cfg = load_config(_write(tmp_path, PROK_BODY + "\namr:\n  amrfinder_organism: Escherichia\n"))
    assert cfg.amr.amrfinder_organism == "Escherichia" and cfg.amr.amrfinder_env == "ali-amrfinder"


def test_amr_amrfinder_default_none(tmp_path):
    assert load_config(_write(tmp_path, PROK_BODY)).amr.amrfinder_organism is None


def test_library_chemistry_parsed(tmp_path):
    cfg = load_config(_write(tmp_path, PROK_BODY + '\nlibrary:\n  chemistry: "direct_rna"\n'))
    assert cfg.library.chemistry == "direct_rna"


def test_library_chemistry_defaults_to_none(tmp_path):
    assert load_config(_write(tmp_path, PROK_BODY)).library.chemistry is None


def test_library_chemistry_invalid_rejected(tmp_path):
    with pytest.raises(ConfigError):
        load_config(_write(tmp_path, PROK_BODY + '\nlibrary:\n  chemistry: "nanopore"\n'))


def test_library_full_length_cdna_defaults_true(tmp_path):
    assert load_config(_write(tmp_path, PROK_BODY)).library.full_length_cdna is True


def test_library_full_length_cdna_parsed_false(tmp_path):
    cfg = load_config(_write(
        tmp_path, PROK_BODY + '\nlibrary:\n  chemistry: "cdna"\n  full_length_cdna: false\n'))
    assert cfg.library.full_length_cdna is False


def test_library_full_length_cdna_invalid_rejected(tmp_path):
    with pytest.raises(ConfigError, match="full_length_cdna"):
        load_config(_write(
            tmp_path, PROK_BODY + '\nlibrary:\n  full_length_cdna: "maybe"\n'))


EUK_BODY = """organism: "human"
organism_type: "eukaryote"
reference:
  transcriptome_fasta: "ref/tx.fa"
  tx2gene: "ref/t2g.tsv"
de:
  design: "~condition"
"""


def test_eukaryote_genome_fasta_optional(tmp_path):
    cfg = load_config(_write(tmp_path, EUK_BODY))
    assert cfg.reference.genome_fasta is None
    assert cfg.reference.transcriptome_fasta is not None


def test_eukaryote_genome_fasta_parsed_when_present(tmp_path):
    body = EUK_BODY + '  genome_fasta: "ref/genome.fa"\n'
    # not: son iki satır reference bloğunun altına gelmeli — reference alanları bitişik
    body = ('organism: "human"\norganism_type: "eukaryote"\n'
            'reference:\n  transcriptome_fasta: "ref/tx.fa"\n  tx2gene: "ref/t2g.tsv"\n'
            '  genome_fasta: "ref/genome.fa"\nde:\n  design: "~condition"\n')
    cfg = load_config(_write(tmp_path, body))
    assert str(cfg.reference.genome_fasta) == "ref/genome.fa"
