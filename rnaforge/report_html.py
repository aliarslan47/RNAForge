"""m08 — HTML report builder. Pure Python + stdlib; assembles a single self-contained report.html
from m06/m07 output contracts. No new data gate; verdict carries over from the confidence card."""
from __future__ import annotations
import base64
import csv
import html
import json
from datetime import datetime
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


def top_degs_by_direction(de: list[dict], gene_map: dict, fdr: float, lfc: float,
                          direction: str, n: int = 25) -> list[dict]:
    rows = [r for r in top_degs(de, gene_map, fdr, lfc, n=10**9) if r["direction"] == direction]
    return rows[:n]


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
        "up_table": "En Güçlü 25 Artan (Up)", "down_table": "En Güçlü 25 Azalan (Down)",
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
        "up_table": "Top 25 Up-regulated", "down_table": "Top 25 Down-regulated",
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


_VERDICT_CLASSES = {"TRUSTWORTHY", "SUSPECT", "INVALID", "UNKNOWN"}


def section_confidence(conf: dict, L: dict) -> str:
    verdict = conf.get("verdict", "UNKNOWN")
    # verdict feeds a CSS class; keep it to the known enum so it can never break out of the attribute.
    vclass = verdict.lower() if verdict in _VERDICT_CLASSES else "unknown"
    counts = conf.get("counts", {})
    prof = conf.get("profile", {})
    gate_rows = [[g.get("name"), g.get("status"), g.get("measured"), g.get("threshold")]
                 for g in conf.get("gates", [])]
    gate_tbl = _table([L["gate"], L["status"], L["measured"], L["threshold"]], gate_rows)
    overrides = prof.get("overrides") or {}
    ov = "" if not overrides else f"<p>overrides: {_esc(overrides)}</p>"
    return (
        f'<section id="confidence"><h2>{_esc(L["confidence"])}</h2>'
        f'<div class="banner verdict-{vclass}"><strong>{_esc(L["verdict"])}: {_esc(verdict)}</strong>'
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


def section_de(de: dict, L: dict) -> str:
    n_sig = de.get("n_significant", 0)
    summary = (f'<p class="summary">{_esc(L["summary"])}: {_esc(n_sig)} / '
               f'{_esc(de.get("n_genes"))} — {_esc(de.get("contrast"))} '
               f'(FDR<{_esc(de.get("fdr_threshold"))}, |log2FC|>={_esc(de.get("log2fc_threshold"))}).</p>')
    rows = [
        [L["contrast"], de.get("contrast")],
        [L["n_genes"], de.get("n_genes")],
        [L["n_sig"], n_sig],
        ["min replicate corr.", de.get("min_replicate_correlation")],
    ]
    tbl = _table([" ", " "], rows)
    return f'<section id="de"><h2>{_esc(L["de"])}</h2>{summary}{tbl}</section>'


def section_figures(figures_manifest: dict, figures_dir: Path, L: dict) -> str:
    figures_dir = Path(figures_dir)
    blocks = []
    for fig in figures_manifest.get("figures", []):
        png = figures_dir / fig["png"]
        if not png.exists():
            raise FileNotFoundError(f"m08: figure PNG missing for report: {png}")
        blocks.append(f'<figure><img src="{embed_png(png)}" alt="{_esc(fig.get("title"))}"/>'
                      f'<figcaption>{_esc(fig.get("title"))}</figcaption></figure>')
    return f'<section id="figures"><h2>{_esc(L["figures"])}</h2>{"".join(blocks)}</section>'


def _deg_table(rows: list[dict], L: dict) -> str:
    body = [[r["gene"],
             f'{r["log2fc"]:.2f}' if r["log2fc"] is not None else "—",
             f'{r["padj"]:.2e}' if r["padj"] is not None else "—",
             f'{r["base_mean"]:.1f}' if r["base_mean"] is not None else "—"] for r in rows]
    return _table([L["gene"], L["log2fc"], L["padj"], L["base_mean"]], body)


def section_table(de_results: list, gene_map: dict, fdr: float, lfc: float, L: dict) -> str:
    up = top_degs_by_direction(de_results, gene_map, fdr, lfc, "Up", n=25)
    down = top_degs_by_direction(de_results, gene_map, fdr, lfc, "Down", n=25)
    if not up and not down:
        return f'<section id="table"><h2>{_esc(L["table"])}</h2><p>{_esc(L["no_degs"])}</p></section>'
    up_html = (f'<h3>{_esc(L["up_table"])}</h3>{_deg_table(up, L)}' if up
               else f'<h3>{_esc(L["up_table"])}</h3><p>{_esc(L["no_degs"])}</p>')
    down_html = (f'<h3>{_esc(L["down_table"])}</h3>{_deg_table(down, L)}' if down
                 else f'<h3>{_esc(L["down_table"])}</h3><p>{_esc(L["no_degs"])}</p>')
    return (f'<section id="table"><h2>{_esc(L["table"])}</h2>{up_html}{down_html}'
            f'<p class="note">{_esc(L["full_table_note"])}</p></section>')


def section_methods(config, L: dict) -> str:
    t = config.trimming
    q = config.quantification
    d = config.de
    text = (
        f'FastQC + fastp (min_length={t.min_length}, aggressive_quality={t.aggressive_quality}; '
        f'Williams et al. 2016). bowtie2 alignment; featureCounts '
        f'(feature_type={q.feature_type}, attribute={q.attribute}). '
        f'DESeq2 (design={d.design}, FDR<{d.fdr_threshold}, |log2FC|>={d.log2fc_threshold}). '
        f'Figures: ggplot2. RNAForge pipeline.'
    )
    return f'<section id="methods"><h2>{_esc(L["methods"])}</h2><p>{_esc(text)}</p></section>'


_REFERENCES = [
    "Langmead B, Salzberg SL. Fast gapped-read alignment with Bowtie 2. Nat Methods. 2012.",
    "Chen S et al. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics. 2018.",
    "Liao Y et al. featureCounts. Bioinformatics. 2014.",
    "Love MI et al. DESeq2. Genome Biol. 2014.",
    "Wickham H. ggplot2: Elegant Graphics for Data Analysis. Springer. 2016.",
    "Williams CR et al. Trimming of sequence reads alters RNA-Seq gene expression estimates. "
    "BMC Bioinformatics. 2016.",
]


def section_references(L: dict) -> str:
    items = "".join(f"<li>{_esc(r)}</li>" for r in _REFERENCES)
    return f'<section id="references"><h2>{_esc(L["references"])}</h2><ol>{items}</ol></section>'


_CSS = """
body{font-family:system-ui,-apple-system,Segoe UI,Roboto,sans-serif;max-width:960px;margin:2rem auto;
padding:0 1rem;color:#1a1a1a;line-height:1.5}
h1{font-size:1.7rem} h2{margin-top:2rem;border-bottom:2px solid #eee;padding-bottom:.3rem}
table{border-collapse:collapse;width:100%;margin:.6rem 0;font-size:.9rem}
th,td{border:1px solid #ddd;padding:.35rem .5rem;text-align:left}
th{background:#f6f6f6} .banner{padding:.8rem 1rem;border-radius:6px;margin:.6rem 0;font-size:1.05rem}
.verdict-trustworthy{background:#e6f4ea;border:1px solid #34a853}
.verdict-suspect{background:#fef7e0;border:1px solid #f9ab00}
.verdict-invalid{background:#fce8e6;border:1px solid #d93025}
.verdict-unknown{background:#f1f3f4;border:1px solid #9aa0a6}
figure{margin:1rem 0;text-align:center} img{max-width:100%;height:auto}
figcaption{color:#555;font-size:.9rem;margin-top:.3rem}
.note{color:#666;font-size:.85rem} .summary{font-size:1.05rem;font-weight:600}
@media print{body{max-width:none} h2{page-break-after:avoid}}
"""


def render_report(inputs: dict, config, version: str, run_id: str = "") -> str:
    lang = config.report.language
    L = LABELS.get(lang, LABELS["tr"])
    raw = inputs["raw"]
    trimming_cfg = {"min_length": config.trimming.min_length,
                    "aggressive": config.trimming.aggressive_quality}
    generated = datetime.now().isoformat(timespec="seconds")
    run_part = f'{_esc(run_id)} · ' if run_id else ""
    header = (f'<h1>RNAForge — {_esc(raw.get("organism"))}</h1>'
              f'<p class="note">{run_part}{_esc(generated)} · v{_esc(version)}</p>')
    body = "".join([
        header,
        section_confidence(inputs["confidence"], L),
        section_dataset(raw, L),
        section_quality(inputs["alignment"], inputs["count"], trimming_cfg, L),
        section_de(inputs["de"], L),
        section_figures(inputs["figures"], inputs["figures_dir"], L),
        section_table(inputs["de_results"], inputs["gene_map"],
                      config.de.fdr_threshold, config.de.log2fc_threshold, L),
        section_methods(config, L),
        section_references(L),
    ])
    return (f'<!doctype html><html lang="{_esc(lang)}"><head><meta charset="utf-8">'
            f'<meta name="viewport" content="width=device-width, initial-scale=1">'
            f'<title>RNAForge report</title><style>{_CSS}</style></head>'
            f'<body>{body}</body></html>')
