from pathlib import Path
from rnaforge.figures import gene_name_map


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
