import shutil, subprocess
from pathlib import Path
import pytest
from rnaforge.figures import run_figures_r, build_manifest, write_gene_map


def _has_de_env():
    return subprocess.run(["conda","run","-n","rnaforge-de","Rscript","-e",
        'cat(requireNamespace("ggplot2",quietly=TRUE))'], capture_output=True, text=True
        ).stdout.strip().endswith("TRUE")


@pytest.mark.skipif(not _has_de_env(), reason="rnaforge-de env/ggplot2 yok")
def test_figures_r_renders_all(tmp_path):
    de = tmp_path / "de"; de.mkdir()
    samples = ["c1","c2","t1","t2"]
    with (de/"normalized_counts.tsv").open("w") as f:
        f.write("gene\t"+"\t".join(samples)+"\n")
        for i in range(1,61):
            base = [100,110,90,105] if i>15 else [20,22,400,420]
            f.write(f"LT_{i}\t"+"\t".join(str(b) for b in base)+"\n")
    with (de/"deseq2_results.tsv").open("w") as f:
        f.write("gene\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\n")
        for i in range(1,61):
            lfc = 4.0 if i<=15 else 0.0; padj = 1e-20 if i<=15 else 0.9
            f.write(f"LT_{i}\t200\t{lfc}\t0.2\t5\t1e-22\t{padj}\n")
    (de/"coldata.tsv").write_text("sample\tcondition\nc1\tcontrol\nc2\tcontrol\nt1\ttreated\nt2\ttreated\n")
    gm = tmp_path/"gm.tsv"; (tmp_path/"g.gff").write_text("chr\tx\tCDS\t1\t9\t.\t+\t0\tlocus_tag=LT_1;gene=pspA\n")
    write_gene_map(tmp_path/"g.gff", gm)
    out = tmp_path/"figures"
    run_figures_r(de, gm, 0.05, 1.0, out)
    man = build_manifest(out)
    assert len(man["figures"]) == 4
    for fig in man["figures"]:
        assert (out/fig["png"]).stat().st_size > 1000
