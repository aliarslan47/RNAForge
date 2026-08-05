"""m12 — GO terim semantik benzerliği + REVIGO-benzeri indirgeme. Saf Python stdlib.

IC arka plan anotasyonundan (m09 build_gene2go), benzerlik Lin (obo ataları), indirgeme greedy.
Fazlalık GO terimlerini (parent/child, benzer süreç) temsilcilere indirir. numpy gerekmez.
"""
from __future__ import annotations

import subprocess
from math import log
from pathlib import Path

from rnaforge.go_annotation import _ancestors

_SCRIPT = Path(__file__).parent / "scripts" / "semantic.R"


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


def lin_distance_matrix(go_ids: list[str], obo: dict, ic: dict[str, float]):
    """Temsilci terimler için simetrik Lin uzaklık matrisi (1 − Lin). MDS (cmdscale) girdisi."""
    cache: dict[str, set[str]] = {}
    n = len(go_ids)
    mat = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            d = 1.0 - lin_similarity(go_ids[i], go_ids[j], obo, ic, cache)
            mat[i][j] = mat[j][i] = d
    return mat


def write_distance_matrix(go_ids: list[str], mat: list[list[float]], out_path: Path) -> None:
    """Kare uzaklık matrisi TSV: ilk sütun go_id, sonra mesafeler (semantic.R cmdscale okur)."""
    with Path(out_path).open("w") as f:
        f.write("go_id\t" + "\t".join(go_ids) + "\n")
        for gid, row in zip(go_ids, mat):
            f.write(gid + "\t" + "\t".join(f"{d:.6f}" for d in row) + "\n")


def run_semantic_r(dist_tsv: Path, reduced_tsv: Path, out_dir: Path, basename: str,
                   title: str, env: str = "rnaforge-de") -> str:
    """semantic.R (MDS scatter, cmdscale + ggplot). stdout/stderr döndür, hatada gürültülü yüksel."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = ["conda", "run", "-n", env, "Rscript", str(_SCRIPT),
           str(dist_tsv), str(reduced_tsv), str(out_dir), basename, title]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"semantic.R failed (exit {r.returncode}):\n{r.stderr}")
    return (r.stdout or "") + (r.stderr or "")
