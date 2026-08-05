"""m09 ORA motoru testleri: hipergeometrik, BH, deg_sets, run_ora, TSV."""
from __future__ import annotations

from math import comb

import pytest

from rnaforge.enrichment import (
    bh_fdr, deg_sets, hypergeometric_pvalue, run_ora, write_ora_tsv,
)


def test_hypergeometric_known_value():
    # N=10,K=5,n=4,k=4 -> yalnız i=4: C(5,4)C(5,0)/C(10,4) = 5/210
    assert hypergeometric_pvalue(4, 4, 5, 10) == pytest.approx(5 / 210)


def test_hypergeometric_upper_tail_sums():
    # k=3, n=4, K=5, N=10 -> i=3,4
    expected = (comb(5, 3) * comb(5, 1) + comb(5, 4) * comb(5, 0)) / comb(10, 4)
    assert hypergeometric_pvalue(3, 4, 5, 10) == pytest.approx(expected)


def test_hypergeometric_degenerate_returns_one():
    assert hypergeometric_pvalue(0, 0, 5, 10) == 1.0
    assert hypergeometric_pvalue(1, 4, 0, 10) == 1.0


def test_bh_monotone_and_bounded():
    adj = bh_fdr([0.01, 0.02, 0.03, 0.04])
    assert all(0 <= a <= 1 for a in adj)
    # BH: p*m/rank -> 0.04,0.04,0.04,0.04 (monoton, hepsi 0.04)
    assert adj == pytest.approx([0.04, 0.04, 0.04, 0.04])


def test_bh_preserves_order():
    adj = bh_fdr([0.5, 0.001, 0.3])
    assert adj[1] < adj[0] and adj[1] < adj[2]   # en küçük p en küçük padj


def test_bh_empty():
    assert bh_fdr([]) == []


DESEQ = (
    "gene\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\n"
    "gUp\t100\t2.5\t0.3\t8\t1e-9\t1e-8\n"
    "gDown\t100\t-3.0\t0.3\t-9\t1e-10\t1e-9\n"
    "gNS\t100\t0.1\t0.3\t0.3\t0.7\t0.8\n"
    "gLowFC\t100\t0.5\t0.3\t1\t0.01\t0.02\n"       # anlamlı ama |lfc|<1
    "gNA\t0\tNA\tNA\tNA\tNA\tNA\n"
)


def test_deg_sets_split_by_direction(tmp_path):
    tsv = tmp_path / "de.tsv"
    tsv.write_text(DESEQ)
    up, down = deg_sets(tsv, fdr=0.05, lfc=1.0)
    assert up == ["gUp"]
    assert down == ["gDown"]           # gLowFC lfc eşiğini geçmez, gNA atlanır


def _annotation():
    # 8 genli evren; termT 4 gene, termOther 3 gene.
    gene2go = {
        "g1": {"GO:T"}, "g2": {"GO:T"}, "g3": {"GO:T"}, "g4": {"GO:T"},
        "g5": {"GO:O"}, "g6": {"GO:O"}, "g7": {"GO:O"}, "g8": {"GO:T"},
    }
    go_meta = {"GO:T": ("BP", "target process"), "GO:O": ("MF", "other")}
    gene_symbol = {f"g{i}": f"sym{i}" for i in range(1, 9)}
    return gene2go, go_meta, gene_symbol


def test_run_ora_enriched_term_significant():
    gene2go, go_meta, sym = _annotation()
    background = list(gene2go)
    gene_set = ["g1", "g2", "g3"]          # hepsi GO:T -> aşırı temsil
    rows = run_ora(gene_set, background, gene2go, go_meta, sym, min_term_size=3)
    top = rows[0]
    assert top["go_id"] == "GO:T"
    assert top["study_count"] == 3
    assert top["fold_enrichment"] > 1
    assert "sym1" in top["genes"]
    assert 0 <= top["p_adj"] <= 1


def test_run_ora_respects_min_term_size():
    gene2go, go_meta, sym = _annotation()
    background = list(gene2go)
    # GO:O arka planda 3 genli; min_term_size=4 -> elenmeli
    rows = run_ora(["g5", "g6"], background, gene2go, go_meta, sym, min_term_size=4)
    assert all(r["go_id"] != "GO:O" for r in rows)


def test_run_ora_empty_study_returns_empty():
    gene2go, go_meta, sym = _annotation()
    assert run_ora([], list(gene2go), gene2go, go_meta, sym) == []


def test_write_ora_tsv_empty_header_only(tmp_path):
    out = tmp_path / "e.tsv"
    write_ora_tsv([], out)
    lines = out.read_text().splitlines()
    assert len(lines) == 1
    assert lines[0].split("\t")[0] == "go_id"


def test_write_ora_tsv_roundtrip(tmp_path):
    gene2go, go_meta, sym = _annotation()
    rows = run_ora(["g1", "g2", "g3"], list(gene2go), gene2go, go_meta, sym, min_term_size=3)
    out = tmp_path / "e.tsv"
    write_ora_tsv(rows, out)
    body = out.read_text().splitlines()
    assert body[1].startswith("GO:T\tBP\ttarget process")
