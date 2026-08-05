"""m11 GSEA girdi hazırlığı testleri: ranked list, GMT inverter, (env-gated) fgsea runner."""
from __future__ import annotations

import subprocess

import pytest

from rnaforge.gsea import invert_to_gmt, run_gsea_r, write_rnk

DESEQ = (
    "gene\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\n"
    "g1\t100\t2.0\t0.3\t6.5\t1e-9\t1e-8\n"
    "g2\t100\t-1.0\t0.2\t-5.0\t1e-6\t1e-5\n"
    "g3\t0\tNA\tNA\tNA\tNA\tNA\n"          # NA stat -> atılır
)


def test_write_rnk_drops_na_and_keeps_stat(tmp_path):
    (tmp_path / "de.tsv").write_text(DESEQ)
    n = write_rnk(tmp_path / "de.tsv", tmp_path / "r.rnk")
    assert n == 2
    lines = (tmp_path / "r.rnk").read_text().splitlines()
    assert lines[0] == "g1\t6.5"
    assert lines[1] == "g2\t-5.0"


def test_write_rnk_missing_stat_column_raises(tmp_path):
    (tmp_path / "de.tsv").write_text("gene\tbaseMean\tlog2FoldChange\tpadj\ng1\t1\t2\t0.01\n")
    with pytest.raises(ValueError, match="stat"):
        write_rnk(tmp_path / "de.tsv", tmp_path / "r.rnk")


def test_invert_to_gmt(tmp_path):
    g2s = {"g1": {"GO:1"}, "g2": {"GO:1", "GO:2"}}
    meta = {"GO:1": ("BP", "alpha"), "GO:2": ("BP", "beta")}
    n = invert_to_gmt(g2s, meta, tmp_path / "x.gmt")
    assert n == 2
    rows = dict(l.split("\t", 1) for l in (tmp_path / "x.gmt").read_text().splitlines())
    assert rows["GO:1"].split("\t")[0] == "alpha"                 # açıklama = set adı
    assert set(rows["GO:1"].split("\t")[1:]) == {"g1", "g2"}      # üyeler
    assert rows["GO:2"].split("\t")[1:] == ["g2"]


def test_invert_to_gmt_name_fallback_to_id(tmp_path):
    n = invert_to_gmt({"g1": {"P9"}}, {}, tmp_path / "x.gmt")     # meta yok -> id ad olur
    assert (tmp_path / "x.gmt").read_text().split("\t")[1] == "P9"


# --- env-gated fgsea entegrasyonu ---
def _has_fgsea():
    return subprocess.run(["conda", "run", "-n", "rnaforge-de", "Rscript", "-e",
        'cat(requireNamespace("fgsea",quietly=TRUE))'], capture_output=True, text=True
        ).stdout.strip().endswith("TRUE")


@pytest.mark.skipif(not _has_fgsea(), reason="rnaforge-de/fgsea yok")
def test_gsea_r_produces_nes(tmp_path):
    # 40 gen: ilk 10'u yüksek pozitif stat, bir set bu 10'u içerir -> pozitif NES beklenir.
    rnk = tmp_path / "r.rnk"
    with rnk.open("w") as f:
        for i in range(1, 41):
            stat = 5.0 - (i * 0.2)          # azalan; ilk genler en yüksek
            f.write(f"g{i}\t{stat}\n")
    gmt = tmp_path / "s.gmt"
    gmt.write_text("SET1\tup set\t" + "\t".join(f"g{i}" for i in range(1, 11)) + "\n")
    gene_map = tmp_path / "gm.tsv"
    gene_map.write_text("locus_tag\tgene\n" + "".join(f"g{i}\tsym{i}\n" for i in range(1, 11)))
    run_gsea_r(rnk, gmt, gene_map, tmp_path, "test", 5, 500, "GSEA test")
    out = (tmp_path / "gsea_test.tsv").read_text().splitlines()
    assert out[0].split("\t")[:2] == ["pathway_id", "name"]
    assert "NES" in out[0]
    row = out[1].split("\t")
    nes = float(row[out[0].split("\t").index("NES")])
    assert nes > 0                          # üst-sıra seti -> pozitif NES
    assert (tmp_path / "gsea_test.png").exists()
