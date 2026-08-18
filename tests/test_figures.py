import json
from pathlib import Path
from rnaforge.figures import gene_name_map
from rnaforge.figures import FIGURE_SPECS, write_gene_map, build_manifest, write_manifest


def test_figure_specs_has_eight_in_order():
    from rnaforge.figures import FIGURE_SPECS
    ids = [s[0] for s in FIGURE_SPECS]
    assert ids == ["pca", "sample_correlation", "expression_dist", "dispersion",
                   "pval_histogram", "volcano", "ma", "heatmap"]


def test_gene_name_map_extracts_named_cds(tmp_path):
    gff = tmp_path / "g.gff"
    gff.write_text(
        "##gff-version 3\n"
        "chr\tx\tCDS\t1\t9\t.\t+\t0\tID=cds1;locus_tag=LT_1;gene=pspA;product=p\n"
        "chr\tx\tCDS\t20\t29\t.\t-\t0\tID=cds2;locus_tag=LT_2;product=hypothetical\n"
        "chr\tx\tgene\t1\t9\t.\t+\t.\tID=gene1;locus_tag=LT_1;gene=pspA\n"
    )
    m = gene_name_map(gff)
    assert m == {"LT_1": "pspA"}   # LT_2 has no gene= -> absent; gene-feature row ignored


def test_write_gene_map_tsv(tmp_path):
    gff = tmp_path / "g.gff"
    gff.write_text("chr\tx\tCDS\t1\t9\t.\t+\t0\tlocus_tag=LT_1;gene=pspA\n")
    out = tmp_path / "gene_map.tsv"
    write_gene_map(gff, out)
    lines = out.read_text().splitlines()
    assert lines[0] == "locus_tag\tgene"
    assert "LT_1\tpspA" in lines


def test_build_manifest_ok_and_missing(tmp_path):
    fig = tmp_path / "figures"; fig.mkdir()
    for _id, base, _title in FIGURE_SPECS:
        (fig / f"{base}.png").write_bytes(b"x")
        (fig / f"{base}.svg").write_text("<svg/>")
    man = build_manifest(fig)
    assert [f["id"] for f in man["figures"]] == [s[0] for s in FIGURE_SPECS]
    assert man["figures"][0]["png"] == "01_pca.png"
    p = write_manifest(fig)
    assert json.loads(p.read_text())["figures"][1]["id"] == "sample_correlation"
    # eksik PNG -> yuksek sesle
    (fig / "01_pca.png").unlink()
    import pytest
    with pytest.raises(FileNotFoundError):
        build_manifest(fig)


def test_gene_name_map_none_gff_returns_empty(tmp_path):
    # Ökaryot: annotation_gff yok → boş map, çökme yok (figürler ID'yle etiketler)
    assert gene_name_map(None) == {}


def test_write_gene_map_none_writes_header_only(tmp_path):
    out = tmp_path / "gene_map.tsv"
    write_gene_map(None, out)
    assert out.read_text() == "locus_tag\tgene\n"
