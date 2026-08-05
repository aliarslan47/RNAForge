"""m11 — GSEA girdi hazırlığı + fgsea runner. Ranked liste (DESeq2 stat) + GMT gen setleri.

Gen-seti kurucuları m09/m10'dan yeniden kullanılır (build_gene2go / build_gene2pathway); burada
yalnız ters çevirme (gen→set → set→genler) + fgsea çağrısı. Motor fgsea (rnaforge-de, altın standart).
"""
from __future__ import annotations

import subprocess
from pathlib import Path

_SCRIPT = Path(__file__).parent / "scripts" / "gsea.R"


def write_rnk(deseq_tsv: Path, out_path: Path) -> int:
    """deseq2_results.tsv -> ranked .rnk (gene\\tstat). NA/boş stat atılır. Yazılan gen sayısını döndürür."""
    lines = Path(deseq_tsv).read_text().splitlines()
    if not lines:
        Path(out_path).write_text("")
        return 0
    header = lines[0].split("\t")
    if "stat" not in header:
        raise ValueError(
            f"m11 (gsea): deseq2_results.tsv has no 'stat' column (found: {header}). "
            "GSEA needs the Wald statistic as ranking metric.")
    gi, si = header.index("gene"), header.index("stat")
    n = 0
    with Path(out_path).open("w") as f:
        for line in lines[1:]:
            cols = line.split("\t")
            if len(cols) <= max(gi, si):
                continue
            stat = cols[si]
            if stat in ("", "NA", "NaN"):
                continue
            try:
                float(stat)
            except ValueError:
                continue
            f.write(f"{cols[gi]}\t{stat}\n")
            n += 1
    return n


def invert_to_gmt(gene2set: dict[str, set[str]], meta: dict[str, tuple[str, str]],
                  out_path: Path) -> int:
    """gen→set haritasını GMT'ye ters çevir: `set_id\\tname\\tgene1\\tgene2…`. Set sayısını döndürür."""
    set2genes: dict[str, set[str]] = {}
    for gene, sets in gene2set.items():
        for sid in sets:
            set2genes.setdefault(sid, set()).add(gene)
    n = 0
    with Path(out_path).open("w") as f:
        for sid in sorted(set2genes):
            genes = set2genes[sid]
            if not genes:
                continue
            name = meta.get(sid, ("", sid))[1] or sid
            f.write(f"{sid}\t{name}\t" + "\t".join(sorted(genes)) + "\n")
            n += 1
    return n


def run_gsea_r(rnk: Path, gmt: Path, gene_map: Path, out_dir: Path, collection: str,
               min_size: int, max_size: int, title: str, env: str = "rnaforge-de") -> str:
    """gsea.R'i çalıştır (fgsea + NES dot-plot). stdout/stderr döndür, hatada gürültülü yüksel."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = ["conda", "run", "-n", env, "Rscript", str(_SCRIPT),
           str(rnk), str(gmt), str(gene_map), str(out_dir), collection,
           str(min_size), str(max_size), title]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"gsea.R failed (exit {r.returncode}):\n{r.stderr}")
    return (r.stdout or "") + (r.stderr or "")
