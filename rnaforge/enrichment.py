"""m09 — GO over-representation analysis (ORA). Saf Python + stdlib (math.comb).

Yön ayrı setler (UP/DOWN), arka plan = anotasyonlu test edilen genler, tek-yönlü
hipergeometrik test, BH FDR her namespace içinde. YENİ veri-kapısı YOK (yorumlayıcı).
"""
from __future__ import annotations

import json
import subprocess
from math import comb
from pathlib import Path

_SCRIPT = Path(__file__).parent / "scripts" / "enrichment.R"

_TSV_HEADER = [
    "go_id", "namespace", "term", "study_count", "study_n", "bg_count", "bg_n",
    "expected", "fold_enrichment", "p_value", "p_adj", "genes",
]


def hypergeometric_pvalue(k: int, n: int, K: int, N: int) -> float:
    """Over-representation üst-kuyruk p: P(X >= k), X ~ Hypergeom(N, K, n).

    k=terimdeki set geni, n=set büyüklüğü, K=terimdeki arka plan geni, N=arka plan."""
    if n == 0 or K == 0:
        return 1.0
    top = min(K, n)
    total = comb(N, n)
    if total == 0:
        return 1.0
    p = sum(comb(K, i) * comb(N - K, n - i) for i in range(k, top + 1)) / total
    return min(1.0, p)


def bh_fdr(pvalues: list[float]) -> list[float]:
    """Benjamini–Hochberg düzeltilmiş p (orijinal sırada). Monoton + [0,1]."""
    m = len(pvalues)
    if m == 0:
        return []
    order = sorted(range(m), key=lambda i: pvalues[i])
    adj = [0.0] * m
    prev = 1.0
    for rank in range(m, 0, -1):     # büyükten küçüğe: kümülatif min
        idx = order[rank - 1]
        val = min(prev, pvalues[idx] * m / rank)
        adj[idx] = val
        prev = val
    return [min(1.0, a) for a in adj]


def deg_sets(deseq_tsv: Path, fdr: float, lfc: float) -> tuple[list[str], list[str]]:
    """deseq2_results.tsv -> (up, down) locus_tag listeleri. NA padj atlanır."""
    up: list[str] = []
    down: list[str] = []
    lines = Path(deseq_tsv).read_text().splitlines()
    if not lines:
        return up, down
    header = lines[0].split("\t")
    gi, li, pi = header.index("gene"), header.index("log2FoldChange"), header.index("padj")
    for line in lines[1:]:
        cols = line.split("\t")
        if len(cols) <= max(gi, li, pi):
            continue
        try:
            padj = float(cols[pi])
            l2fc = float(cols[li])
        except ValueError:            # NA padj/log2FC
            continue
        if padj >= fdr:
            continue
        if l2fc >= lfc:
            up.append(cols[gi])
        elif l2fc <= -lfc:
            down.append(cols[gi])
    return up, down


def all_tested_genes(deseq_tsv: Path) -> list[str]:
    """deseq2_results.tsv'deki tüm gen id'leri (ORA arka plan evreni)."""
    lines = Path(deseq_tsv).read_text().splitlines()
    if not lines:
        return []
    gi = lines[0].split("\t").index("gene")
    return [c.split("\t")[gi] for c in lines[1:] if c.split("\t")[gi:gi + 1]]


def run_ora(gene_set: list[str], background: list[str], gene2go: dict[str, set[str]],
            go_meta: dict[str, tuple[str, str]], gene_symbol: dict[str, str],
            min_term_size: int = 3) -> list[dict]:
    """Bir gen seti için ORA. Namespace başına BH. padj artan sıralı satırlar döndürür."""
    annotated_bg = [g for g in background if gene2go.get(g)]
    bg_set = set(annotated_bg)
    N = len(annotated_bg)
    study = [g for g in gene_set if g in bg_set]
    n = len(study)
    if N == 0 or n == 0:
        return []

    # term -> arka plan / set genleri
    term_bg: dict[str, set[str]] = {}
    for g in annotated_bg:
        for go_id in gene2go[g]:
            term_bg.setdefault(go_id, set()).add(g)
    term_study: dict[str, set[str]] = {}
    for g in study:
        for go_id in gene2go[g]:
            term_study.setdefault(go_id, set()).add(g)

    rows: list[dict] = []
    for go_id, bg_genes in term_bg.items():
        K = len(bg_genes)
        if K < min_term_size:
            continue
        study_genes = term_study.get(go_id, set())
        k = len(study_genes)
        if k == 0:
            continue
        ns, name = go_meta.get(go_id, ("?", ""))
        expected = n * K / N
        fold = (k / n) / (K / N)
        p = hypergeometric_pvalue(k, n, K, N)
        genes = sorted(gene_symbol.get(g, g) for g in study_genes)
        rows.append({
            "go_id": go_id, "namespace": ns, "term": name,
            "study_count": k, "study_n": n, "bg_count": K, "bg_n": N,
            "expected": expected, "fold_enrichment": fold,
            "p_value": p, "genes": genes,
        })

    # BH her namespace grubu içinde
    for ns in {r["namespace"] for r in rows}:
        grp = [r for r in rows if r["namespace"] == ns]
        adj = bh_fdr([r["p_value"] for r in grp])
        for r, a in zip(grp, adj):
            r["p_adj"] = a

    rows.sort(key=lambda r: (r["p_adj"], r["p_value"]))
    return rows


def write_ora_tsv(rows: list[dict], path: Path) -> None:
    """ORA satırlarını TSV'ye yaz. Boş set -> yalnız başlık (çökme yok)."""
    with Path(path).open("w") as f:
        f.write("\t".join(_TSV_HEADER) + "\n")
        for r in rows:
            f.write("\t".join([
                r["go_id"], r["namespace"], r["term"],
                str(r["study_count"]), str(r["study_n"]),
                str(r["bg_count"]), str(r["bg_n"]),
                f'{r["expected"]:.4f}', f'{r["fold_enrichment"]:.4f}',
                f'{r["p_value"]:.6e}', f'{r["p_adj"]:.6e}',
                ";".join(r["genes"]),
            ]) + "\n")


def run_enrichment_r(up_tsv: Path, down_tsv: Path, out_dir: Path, top_n: int,
                     title_prefix: str = "GO zenginleştirme",
                     basename_prefix: str = "enrichment", env: str = "rnaforge-de") -> str:
    """enrichment.R'i çalıştır (dot-plot PNG+SVG). GO (varsayılan) ve KEGG (prefix'li) ortak kullanır.
    stdout/stderr döndür, hatada gürültülü yüksel."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = ["conda", "run", "-n", env, "Rscript", str(_SCRIPT),
           str(up_tsv), str(down_tsv), str(out_dir), str(top_n),
           title_prefix, basename_prefix]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"enrichment.R failed (exit {r.returncode}):\n{r.stderr}")
    return (r.stdout or "") + (r.stderr or "")


def build_enrichment_manifest(fig_dir: Path, basename_prefix: str = "enrichment") -> dict:
    """Var olan enrichment figürlerini manifest'e topla. Boş-durumda eksik PNG olmaz
    (R her zaman panel üretir), ama m07'nin aksine eksikte yükselmez — dürüst boş liste.
    basename_prefix ile GO (enrichment) ve KEGG (kegg) için yeniden kullanılır."""
    fig_dir = Path(fig_dir)
    titles = {"up": "Artan", "down": "Azalan"}
    figures = []
    for direction, title in titles.items():
        base = f"{basename_prefix}_{direction}"
        png, svg = fig_dir / f"{base}.png", fig_dir / f"{base}.svg"
        if not png.exists():
            continue
        figures.append({"id": base, "title": title, "png": png.name,
                        "svg": svg.name if svg.exists() else None})
    return {"figures": figures}


def write_enrichment_manifest(fig_dir: Path, basename_prefix: str = "enrichment") -> Path:
    p = Path(fig_dir) / "manifest.json"
    p.write_text(json.dumps(build_enrichment_manifest(fig_dir, basename_prefix), indent=2))
    return p
