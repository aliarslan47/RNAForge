"""m14 — Operon tahmini (intergenik-mesafe sezgiseli) + operon-düzeyi DE koordinasyonu. Saf Python.

Aynı yönde bitişik + gap ≤ max_gap genler aynı operon (Moreno-Hagelsieb & Collado-Vides 2002).
Operonlar TAHMİN (deneysel değil). DB/araç yok, organizma-agnostik.
"""
from __future__ import annotations

from pathlib import Path


def _genes_from_gff(gff: Path) -> list[dict]:
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
                    "strand": cols[6], "locus_tag": lt, "symbol": attrs.get("gene", "")})
    return out


def predict_operons(gff: Path, max_gap: int = 50) -> list[dict]:
    """Aynı contig + aynı strand + `start − prev_end − 1 ≤ max_gap` ardışık genleri operona topla.
    Strand/contig değişimi veya büyük gap → operon sınırı. Tek-genli operonlar da üretilir."""
    genes = sorted(_genes_from_gff(gff), key=lambda g: (g["contig"], g["start"]))
    operons: list[dict] = []
    cur: list[dict] | None = None

    def flush():
        if cur:
            operons.append({
                "operon_id": f"operon_{len(operons) + 1}",
                "contig": cur[0]["contig"], "strand": cur[0]["strand"],
                "locus_tags": [g["locus_tag"] for g in cur],
                "symbols": [g["symbol"] or g["locus_tag"] for g in cur],
                "size": len(cur),
            })

    for g in genes:
        if cur is None:
            cur = [g]
            continue
        prev = cur[-1]
        gap = g["start"] - prev["end"] - 1
        if g["contig"] == prev["contig"] and g["strand"] == prev["strand"] and gap <= max_gap:
            cur.append(g)
        else:
            flush(); cur = [g]
    flush()
    return operons


def _read_de(deseq_tsv: Path) -> dict[str, tuple[float | None, float | None]]:
    lines = Path(deseq_tsv).read_text().splitlines() if Path(deseq_tsv).exists() else []
    if not lines:
        return {}
    h = lines[0].split("\t")
    gi, li, pi = h.index("gene"), h.index("log2FoldChange"), h.index("padj")

    def num(x):
        try:
            return float(x)
        except ValueError:
            return None
    out = {}
    for line in lines[1:]:
        c = line.split("\t")
        if len(c) > max(gi, li, pi):
            out[c[gi]] = (num(c[li]), num(c[pi]))
    return out


def aggregate_operon_de(operons: list[dict], deseq_tsv: Path, fdr: float, lfc: float) -> list[dict]:
    """Operon başına DE metrikleri. coordinated = ≥2 gen & ≥2 DEG & hepsi aynı yön."""
    de = _read_de(deseq_tsv)
    out = []
    for op in operons:
        l2fcs, n_up, n_down = [], 0, 0
        for lt in op["locus_tags"]:
            l2fc, padj = de.get(lt, (None, None))
            if l2fc is None or padj is None:
                continue
            l2fcs.append(l2fc)
            if padj < fdr and l2fc >= lfc:
                n_up += 1
            elif padj < fdr and l2fc <= -lfc:
                n_down += 1
        n_deg = n_up + n_down
        mean_l2fc = sum(l2fcs) / len(l2fcs) if l2fcs else None
        coordinated = op["size"] >= 2 and n_deg >= 2 and (n_up == n_deg or n_down == n_deg)
        rec = dict(op)
        rec.update({"n_tested": len(l2fcs), "n_deg": n_deg, "n_up": n_up, "n_down": n_down,
                    "mean_log2fc": mean_l2fc, "coordinated": coordinated})
        out.append(rec)
    out.sort(key=lambda r: (not r["coordinated"], -r["n_deg"]))
    return out
