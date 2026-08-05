"""m12 — GO terim semantik benzerliği + REVIGO-benzeri indirgeme. Saf Python stdlib.

IC arka plan anotasyonundan (m09 build_gene2go), benzerlik Lin (obo ataları), indirgeme greedy.
Fazlalık GO terimlerini (parent/child, benzer süreç) temsilcilere indirir. numpy gerekmez.
"""
from __future__ import annotations

from math import log

from rnaforge.go_annotation import _ancestors


def compute_ic(gene2go: dict[str, set[str]]) -> dict[str, float]:
    """Information content: IC(t) = -log(count(t)/N). gene2go propagate edilmiş olmalı.
    count(t) = t'ye anotlı gen sayısı, N = toplam anotlı gen. Kök terim IC≈0, nadir terim yüksek."""
    n = len(gene2go)
    if n == 0:
        return {}
    count: dict[str, int] = {}
    for gos in gene2go.values():
        for t in gos:
            count[t] = count.get(t, 0) + 1
    return {t: -log(c / n) for t, c in count.items() if c > 0}


def lin_similarity(a: str, b: str, obo: dict, ic: dict[str, float],
                   cache: dict[str, set[str]]) -> float:
    """Lin benzerliği: 2·IC(MICA)/(IC(a)+IC(b)). [0,1]. a==b→1; ortak ata yok/payda 0→0."""
    if a == b:
        return 1.0
    ica, icb = ic.get(a, 0.0), ic.get(b, 0.0)
    denom = ica + icb
    if denom <= 0:
        return 0.0
    common = ({a} | _ancestors(a, obo, cache)) & ({b} | _ancestors(b, obo, cache))
    if not common:
        return 0.0
    mica_ic = max((ic.get(t, 0.0) for t in common), default=0.0)
    return (2.0 * mica_ic) / denom


def reduce_terms(terms: list[dict], obo: dict, ic: dict[str, float],
                 threshold: float = 0.7) -> list[dict]:
    """REVIGO-benzeri greedy indirgeme. terms: [{go_id, namespace, term, padj}].
    Namespace başına: padj artan sırala; her terim temsilcilere max Lin≥threshold ise o kümeye,
    değilse yeni temsilci. Dönüş: temsilci + n_collapsed + members (go_id listesi)."""
    cache: dict[str, set[str]] = {}
    by_ns: dict[str, list[dict]] = {}
    for t in terms:
        by_ns.setdefault(t.get("namespace", "?"), []).append(t)

    reps: list[dict] = []
    for ns, group in by_ns.items():
        group = sorted(group, key=lambda r: (r.get("padj") if r.get("padj") is not None else 1.0))
        ns_reps: list[dict] = []
        for t in group:
            best_sim, best_rep = 0.0, None
            for rep in ns_reps:
                s = lin_similarity(t["go_id"], rep["go_id"], obo, ic, cache)
                if s > best_sim:
                    best_sim, best_rep = s, rep
            if best_rep is not None and best_sim >= threshold:
                best_rep["_members"].append(t["go_id"])
            else:
                t = dict(t); t["_members"] = [t["go_id"]]
                ns_reps.append(t)
        reps.extend(ns_reps)

    out = []
    for r in reps:
        members = r.pop("_members")
        out.append({
            "go_id": r["go_id"], "namespace": r.get("namespace", "?"),
            "term": r.get("term", ""), "padj": r.get("padj"),
            "n_collapsed": len(members), "members": members,
        })
    return out
