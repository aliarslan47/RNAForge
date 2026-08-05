"""m12 semantic testleri: IC, Lin benzerliği, REVIGO-benzeri indirgeme."""
from __future__ import annotations

from rnaforge.semantic import compute_ic, lin_similarity, reduce_terms

# Küçük DAG: GO:1 kök; GO:2, GO:3 kardeş (GO:1 çocuğu); GO:9 farklı namespace (CC).
OBO = {
    "GO:1": {"parents": set(), "namespace": "BP", "obsolete": False},
    "GO:2": {"parents": {"GO:1"}, "namespace": "BP", "obsolete": False},
    "GO:3": {"parents": {"GO:1"}, "namespace": "BP", "obsolete": False},
    "GO:9": {"parents": set(), "namespace": "CC", "obsolete": False},
}


def test_compute_ic_rarer_term_higher():
    g2go = {"g1": {"GO:1", "GO:2"}, "g2": {"GO:1", "GO:2"}, "g3": {"GO:1"}}
    ic = compute_ic(g2go)
    assert ic["GO:1"] < ic["GO:2"]      # GO:1 tüm 3 gende (yaygın), GO:2 2 gende
    assert ic["GO:1"] == 0.0            # kök: tüm genlerde -> -log(1)=0


def test_lin_identity():
    assert lin_similarity("GO:2", "GO:2", OBO, {"GO:2": 2.0}, {}) == 1.0


def test_lin_siblings_share_parent():
    ic = {"GO:1": 0.5, "GO:2": 2.0, "GO:3": 2.0}
    s = lin_similarity("GO:2", "GO:3", OBO, ic, {})
    # MICA = GO:1 (tek ortak ata); Lin = 2*0.5/(2+2) = 0.25
    assert abs(s - 0.25) < 1e-9


def test_lin_zero_when_no_common_ancestor():
    ic = {"GO:2": 2.0, "GO:9": 2.0}
    assert lin_similarity("GO:2", "GO:9", OBO, ic, {}) == 0.0   # BP vs CC, ortak ata yok


def _terms():
    return [
        {"go_id": "GO:2", "namespace": "BP", "term": "alpha", "padj": 1e-8},
        {"go_id": "GO:3", "namespace": "BP", "term": "beta", "padj": 1e-4},
        {"go_id": "GO:9", "namespace": "CC", "term": "gamma", "padj": 1e-6},
    ]


def test_reduce_collapses_similar_within_namespace():
    # Yüksek IC + düşük eşik -> GO:2, GO:3 tek temsilcide (BP); GO:9 ayrı (CC)
    ic = {"GO:1": 0.5, "GO:2": 2.0, "GO:3": 2.0, "GO:9": 2.0}
    out = reduce_terms(_terms(), OBO, ic, threshold=0.2)   # Lin(GO2,GO3)=0.25 >= 0.2
    bp = [r for r in out if r["namespace"] == "BP"]
    cc = [r for r in out if r["namespace"] == "CC"]
    assert len(bp) == 1 and bp[0]["go_id"] == "GO:2"       # en iyi padj temsilci
    assert bp[0]["n_collapsed"] == 2 and set(bp[0]["members"]) == {"GO:2", "GO:3"}
    assert len(cc) == 1                                     # farklı namespace ayrı


def test_reduce_keeps_dissimilar_separate():
    ic = {"GO:1": 0.5, "GO:2": 2.0, "GO:3": 2.0, "GO:9": 2.0}
    out = reduce_terms(_terms(), OBO, ic, threshold=0.9)   # 0.25 < 0.9 -> ayrı kalır
    assert len([r for r in out if r["namespace"] == "BP"]) == 2


def test_reduce_empty():
    assert reduce_terms([], OBO, {}, 0.7) == []
