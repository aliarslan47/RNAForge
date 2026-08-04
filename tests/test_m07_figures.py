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


import json as _json
from rnaforge.modules import m07_figures
from rnaforge.state import RunState


def _min_config(tmp_path):
    from rnaforge.config import load_config
    (tmp_path/"c.yaml").write_text(
        "organism: E\norganism_type: prokaryote\nplatform: auto\n"
        "reference:\n  genome_fasta: g.fa\n  annotation_gff: "+str(tmp_path/'g.gff')+"\n"
        "de:\n  design: '~condition'\n  fdr_threshold: 0.05\n  log2fc_threshold: 1.0\n")
    (tmp_path/"g.gff").write_text("chr\tx\tCDS\t1\t9\t.\t+\t0\tlocus_tag=LT_1;gene=pspA\n")
    return load_config(tmp_path/"c.yaml")


def test_run_figures_requires_m06(tmp_path, monkeypatch):
    rd = tmp_path/"run"; (rd/"differential_expression").mkdir(parents=True)
    cfg = _min_config(tmp_path); md = tmp_path/"m.tsv"
    md.write_text("sample_id\tcondition\tfastq_1\nc1\tcontrol\t/x/a.fq\n")
    with __import__("pytest").raises(ValueError):
        m07_figures.run_figures(cfg, md, rd)


def test_run_figures_orchestrates(tmp_path, monkeypatch):
    rd = tmp_path/"run"; de=(rd/"differential_expression"); de.mkdir(parents=True)
    (rd/"statistics").mkdir(); (rd/"logs").mkdir()
    RunState(rd).mark_done("m06_de", [])
    cfg = _min_config(tmp_path); md = tmp_path/"m.tsv"
    md.write_text("sample_id\tcondition\tfastq_1\nc1\tcontrol\t/x/a.fq\n")
    def fake_r(de_dir, gene_map, fdr, lfc, out_dir, env="rnaforge-de"):
        from rnaforge.figures import FIGURE_SPECS
        Path(out_dir).mkdir(parents=True, exist_ok=True)
        for _id,base,_t in FIGURE_SPECS:
            (Path(out_dir)/f"{base}.png").write_bytes(b"x"*2000)
            (Path(out_dir)/f"{base}.svg").write_text("<svg/>")
    monkeypatch.setattr(m07_figures, "run_figures_r", fake_r)
    s = m07_figures.run_figures(cfg, md, rd)
    assert s["n_figures"] == 4
    assert (rd/"figures"/"manifest.json").exists()
    assert RunState(rd).is_done("m07_figures")
    # resume
    s2 = m07_figures.run_figures(cfg, md, rd)
    assert s2.get("resumed") is True
