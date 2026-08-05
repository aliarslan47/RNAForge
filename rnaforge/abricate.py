"""m13 — abricate (AMR/virülans) parser + koordinat eşleme + DE overlay. Saf Python + abricate runner.

abricate genome.fa'yı CARD/VFDB'ye tarar; hit'ler koordinatla GFF genlerine (locus_tag) eşlenir; sonra
DESeq2 DE durumu overlay edilir. Eşleme gen adıyla değil KOORDİNATLA (abricate adı sembolle eşleşmeyebilir).
"""
from __future__ import annotations

import csv
import subprocess
from pathlib import Path


def run_abricate(genome_fa: Path, db: str, out_tsv: Path, env: str = "rnaforge-amr") -> str:
    """abricate --db <db> <genome> -> out_tsv. stdout/stderr döndür, hatada gürültülü yüksel."""
    out_tsv = Path(out_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["conda", "run", "-n", env, "abricate", "--db", db, "--quiet", str(genome_fa)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"abricate (--db {db}) failed (exit {r.returncode}):\n{r.stderr}")
    out_tsv.write_text(r.stdout)
    return r.stderr or ""


def parse_abricate(tsv: Path, min_identity: float = 80.0, min_coverage: float = 80.0) -> list[dict]:
    """abricate 1.4 TSV -> hit'ler. %identity/%coverage filtreli. Başlık `#FILE…` (DictReader)."""
    tsv = Path(tsv)
    if not tsv.exists():
        return []
    hits = []
    with tsv.open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            try:
                pct_id = float(row.get("%IDENTITY", ""))
                pct_cov = float(row.get("%COVERAGE", ""))
                start, end = int(row["START"]), int(row["END"])
            except (ValueError, KeyError, TypeError):
                continue
            if pct_id < min_identity or pct_cov < min_coverage:
                continue
            hits.append({
                "contig": row.get("SEQUENCE", ""), "start": start, "end": end,
                "gene": row.get("GENE", ""), "pct_id": pct_id, "pct_cov": pct_cov,
                "db": row.get("DATABASE", ""), "product": row.get("PRODUCT", ""),
                "resistance": row.get("RESISTANCE", ""),
            })
    return hits


def gene_coords(gff: Path) -> list[dict]:
    """GFF `gene` feature'larından {contig, start, end, locus_tag, symbol}."""
    out = []
    for line in Path(gff).read_text().splitlines():
        if not line or line.startswith("#") or "\tgene\t" not in line:
            continue
        cols = line.split("\t")
        if len(cols) < 9:
            continue
        attrs = dict(kv.split("=", 1) for kv in cols[8].split(";") if "=" in kv)
        lt = attrs.get("locus_tag")
        if not lt:
            continue
        out.append({"contig": cols[0], "start": int(cols[3]), "end": int(cols[4]),
                    "locus_tag": lt, "symbol": attrs.get("gene", "")})
    return out


def map_hits_to_genes(hits: list[dict], genes: list[dict]) -> tuple[list[dict], int]:
    """Her hit'i aynı contig'te EN ÇOK örtüşen gene ata. Eşleşmeyen sayısı ayrı döner.
    Aynı locus_tag'e birden çok hit → en yüksek %identity tutulur (dedup)."""
    by_contig: dict[str, list[dict]] = {}
    for g in genes:
        by_contig.setdefault(g["contig"], []).append(g)

    best: dict[str, dict] = {}
    n_unmapped = 0
    for h in hits:
        cand = by_contig.get(h["contig"], [])
        best_gene, best_ov = None, 0
        for g in cand:
            ov = min(h["end"], g["end"]) - max(h["start"], g["start"])
            if ov > best_ov:
                best_ov, best_gene = ov, g
        if best_gene is None:
            n_unmapped += 1
            continue
        lt = best_gene["locus_tag"]
        rec = dict(h); rec["locus_tag"] = lt; rec["symbol"] = best_gene["symbol"]
        if lt not in best or rec["pct_id"] > best[lt]["pct_id"]:
            best[lt] = rec
    return list(best.values()), n_unmapped


def _read_de(deseq_tsv: Path) -> dict[str, tuple[float | None, float | None]]:
    lines = Path(deseq_tsv).read_text().splitlines() if Path(deseq_tsv).exists() else []
    if not lines:
        return {}
    h = lines[0].split("\t")
    gi, li, pi = h.index("gene"), h.index("log2FoldChange"), h.index("padj")
    out = {}

    def num(x):
        try:
            return float(x)
        except ValueError:
            return None
    for line in lines[1:]:
        c = line.split("\t")
        if len(c) > max(gi, li, pi):
            out[c[gi]] = (num(c[li]), num(c[pi]))
    return out


def overlay_de(mapped: list[dict], deseq_tsv: Path, fdr: float, lfc: float) -> list[dict]:
    """Eşlenmiş hit'lere log2FC/padj/de_status ekle. Sıralama: up/down önce, sonra padj."""
    de = _read_de(deseq_tsv)
    order = {"up": 0, "down": 1, "ns": 2, "untested": 3}
    out = []
    for m in mapped:
        l2fc, padj = de.get(m["locus_tag"], (None, None))
        if padj is None or l2fc is None:
            status = "untested"
        elif padj < fdr and l2fc >= lfc:
            status = "up"
        elif padj < fdr and l2fc <= -lfc:
            status = "down"
        else:
            status = "ns"
        rec = dict(m); rec["log2fc"] = l2fc; rec["padj"] = padj; rec["de_status"] = status
        out.append(rec)
    out.sort(key=lambda r: (order[r["de_status"]],
                            r["padj"] if r["padj"] is not None else 1.0))
    return out
