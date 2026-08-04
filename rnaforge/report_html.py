"""m08 — HTML report builder. Pure Python + stdlib; assembles a single self-contained report.html
from m06/m07 output contracts. No new data gate; verdict carries over from the confidence card."""
from __future__ import annotations
import base64
import csv
import html
import json
from pathlib import Path

N_SECTIONS = 9


def _num(v):
    if v is None:
        return None
    s = str(v).strip()
    if s == "" or s.upper() == "NA":
        return None
    try:
        return float(s)
    except ValueError:
        return None


def parse_de_results(path: Path) -> list[dict]:
    rows: list[dict] = []
    with Path(path).open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            row = dict(r)
            for k in ("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"):
                if k in row:
                    row[k] = _num(row[k])
            rows.append(row)
    return rows


def load_gene_map(path: Path) -> dict[str, str]:
    path = Path(path)
    if not path.exists():
        return {}
    out: dict[str, str] = {}
    with path.open() as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader, None)  # header
        for row in reader:
            if len(row) >= 2 and row[0]:
                out[row[0]] = row[1]
    return out


def top_degs(de: list[dict], gene_map: dict, fdr: float, lfc: float, n: int = 50) -> list[dict]:
    sig = [r for r in de
           if r.get("padj") is not None and r["padj"] < fdr
           and r.get("log2FoldChange") is not None and abs(r["log2FoldChange"]) >= lfc]
    sig.sort(key=lambda r: r["padj"])
    out = []
    for r in sig[:n]:
        gid = r["gene"]
        out.append({
            "gene": gene_map.get(gid, gid),
            "log2fc": r["log2FoldChange"],
            "padj": r["padj"],
            "base_mean": r.get("baseMean"),
            "direction": "Up" if r["log2FoldChange"] > 0 else "Down",
        })
    return out


def embed_png(path: Path) -> str:
    data = Path(path).read_bytes()
    return "data:image/png;base64," + base64.b64encode(data).decode("ascii")


def _read_json(path: Path) -> dict:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"m08 report input missing: {path}")
    return json.loads(path.read_text())


def load_report_inputs(run_dir: Path) -> dict:
    run_dir = Path(run_dir)
    stats = run_dir / "statistics"
    de_tsv = run_dir / "differential_expression" / "deseq2_results.tsv"
    if not de_tsv.exists():
        raise FileNotFoundError(f"m08 report input missing: {de_tsv}")
    figures_dir = run_dir / "figures"
    return {
        "raw": _read_json(stats / "raw_statistics.json"),
        "qc": _read_json(stats / "qc_statistics.json"),
        "trimming": _read_json(stats / "trimming_statistics.json"),
        "alignment": _read_json(stats / "alignment_statistics.json"),
        "count": _read_json(stats / "count_statistics.json"),
        "de": _read_json(stats / "de_statistics.json"),
        "confidence": _read_json(run_dir / "quality" / "confidence_card.json"),
        "figures": _read_json(figures_dir / "manifest.json"),
        "de_results": parse_de_results(de_tsv),
        "gene_map": load_gene_map(figures_dir / "gene_map.tsv"),
        "figures_dir": figures_dir,
    }


LABELS: dict[str, dict[str, str]] = {
    "tr": {
        "confidence": "Güvence Kartı", "dataset": "Veri Kümesi ve Örnekler",
        "quality": "Kalite ve İşleme", "de": "Diferansiyel Ekspresyon",
        "figures": "Figürler", "table": "En Güçlü DEG'ler", "methods": "Yöntemler",
        "references": "Kaynaklar", "verdict": "Karar", "no_degs": "Anlamlı DEG bulunamadı.",
        "sample": "Örnek", "condition": "Koşul", "batch": "Batch", "paired": "Eşleşmiş",
        "read_len": "Ort. okuma uzunluğu", "quality_col": "Ort. kalite",
        "alignment_rate": "Hizalama oranı", "assignment_rate": "Atama oranı",
        "organism": "Organizma", "platform": "Platform", "design": "Tasarım",
        "gate": "Kapı", "status": "Durum", "measured": "Ölçülen", "threshold": "Eşik",
        "profile": "Profil", "contrast": "Kontrast", "n_genes": "Gen sayısı",
        "n_sig": "Anlamlı gen", "gene": "Gen", "log2fc": "log2FC", "padj": "padj",
        "base_mean": "baseMean", "direction": "Yön", "up": "Artan", "down": "Azalan",
        "full_table_note": "Tam tablo: differential_expression/deseq2_results.tsv",
        "min_length": "Min uzunluk", "aggressive": "Agresif kalite trimming",
        "summary": "Özet",
    },
    "en": {
        "confidence": "Confidence Card", "dataset": "Dataset and Samples",
        "quality": "Quality and Processing", "de": "Differential Expression",
        "figures": "Figures", "table": "Top DEGs", "methods": "Methods",
        "references": "References", "verdict": "Verdict", "no_degs": "No significant DEGs found.",
        "sample": "Sample", "condition": "Condition", "batch": "Batch", "paired": "Paired",
        "read_len": "Mean read length", "quality_col": "Mean quality",
        "alignment_rate": "Alignment rate", "assignment_rate": "Assignment rate",
        "organism": "Organism", "platform": "Platform", "design": "Design",
        "gate": "Gate", "status": "Status", "measured": "Measured", "threshold": "Threshold",
        "profile": "Profile", "contrast": "Contrast", "n_genes": "Genes",
        "n_sig": "Significant genes", "gene": "Gene", "log2fc": "log2FC", "padj": "padj",
        "base_mean": "baseMean", "direction": "Direction", "up": "Up", "down": "Down",
        "full_table_note": "Full table: differential_expression/deseq2_results.tsv",
        "min_length": "Min length", "aggressive": "Aggressive quality trimming",
        "summary": "Summary",
    },
}


def _esc(x) -> str:
    if x is None:
        return "—"
    return html.escape(str(x))


def _table(headers: list[str], rows: list[list]) -> str:
    head = "".join(f"<th>{_esc(h)}</th>" for h in headers)
    body = "".join("<tr>" + "".join(f"<td>{_esc(c)}</td>" for c in row) + "</tr>" for row in rows)
    return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def _pct(v) -> str:
    return "—" if v is None else f"{float(v) * 100:.1f}%"


def section_confidence(conf: dict, L: dict) -> str:
    verdict = conf.get("verdict", "UNKNOWN")
    counts = conf.get("counts", {})
    prof = conf.get("profile", {})
    gate_rows = [[g.get("name"), g.get("status"), g.get("measured"), g.get("threshold")]
                 for g in conf.get("gates", [])]
    gate_tbl = _table([L["gate"], L["status"], L["measured"], L["threshold"]], gate_rows)
    overrides = prof.get("overrides") or {}
    ov = "" if not overrides else f"<p>overrides: {_esc(overrides)}</p>"
    return (
        f'<section id="confidence"><h2>{_esc(L["confidence"])}</h2>'
        f'<div class="banner verdict-{verdict.lower()}"><strong>{_esc(L["verdict"])}: {_esc(verdict)}</strong>'
        f' &nbsp; PASS={_esc(counts.get("PASS", 0))} WARN={_esc(counts.get("WARN", 0))} '
        f'FAIL={_esc(counts.get("FAIL", 0))}</div>'
        f'<p>{_esc(L["profile"])}: {_esc(prof.get("name"))}</p>{ov}{gate_tbl}</section>'
    )


def section_dataset(raw: dict, L: dict) -> str:
    conds = ", ".join(f"{k}={v}" for k, v in (raw.get("conditions") or {}).items())
    meta = (f'<p>{_esc(L["organism"])}: {_esc(raw.get("organism"))} · '
            f'{_esc(L["platform"])}: {_esc(raw.get("platform"))} · '
            f'{_esc(L["design"])}: {_esc(raw.get("design"))} · {_esc(conds)}</p>')
    rows = [[s.get("sample_id"), s.get("condition"), s.get("batch"), s.get("paired"),
             s.get("mean_read_length"), s.get("mean_quality")] for s in raw.get("samples", [])]
    tbl = _table([L["sample"], L["condition"], L["batch"], L["paired"], L["read_len"], L["quality_col"]], rows)
    return f'<section id="dataset"><h2>{_esc(L["dataset"])}</h2>{meta}{tbl}</section>'


def section_quality(align: dict, count: dict, trimming_cfg: dict, L: dict) -> str:
    trim = (f'<p>{_esc(L["min_length"])}: {_esc(trimming_cfg.get("min_length"))} · '
            f'{_esc(L["aggressive"])}: {_esc(trimming_cfg.get("aggressive"))}</p>')
    asamp = align.get("samples", {})
    csamp = count.get("samples", {})
    rows = [[sid, _pct(asamp.get(sid, {}).get("alignment_rate")),
             _pct(csamp.get(sid, {}).get("assignment_rate"))] for sid in asamp]
    tbl = _table([L["sample"], L["alignment_rate"], L["assignment_rate"]], rows)
    return f'<section id="quality"><h2>{_esc(L["quality"])}</h2>{trim}{tbl}</section>'
