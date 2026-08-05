"""m08 — HTML report builder. Pure Python + stdlib; assembles a single self-contained report.html
from m06/m07 output contracts. No new data gate; verdict carries over from the confidence card."""
from __future__ import annotations
import base64
import csv
import html
import json
from datetime import datetime
from pathlib import Path

N_SECTIONS = 10


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


def parse_enrichment_tsv(path: Path) -> list[dict]:
    """m09 enrichment_{up,down}.tsv -> tipli satırlar. Yoksa/başlık-only -> boş liste."""
    path = Path(path)
    if not path.exists():
        return []
    rows: list[dict] = []
    with path.open() as f:
        for r in csv.DictReader(f, delimiter="\t"):
            row = dict(r)
            for k in ("study_count", "study_n", "bg_count", "bg_n"):
                row[k] = int(row[k]) if row.get(k) not in (None, "") else None
            for k in ("expected", "fold_enrichment", "p_value", "p_adj"):
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
            "gene_id": gid,
            "gene": gene_map.get(gid, gid),
            "log2fc": r["log2FoldChange"],
            "padj": r["padj"],
            "base_mean": r.get("baseMean"),
            "direction": "Up" if r["log2FoldChange"] > 0 else "Down",
        })
    return out


def parse_coldata(path: Path) -> list[tuple[str, str]]:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"m08 report input missing: {path}")
    out: list[tuple[str, str]] = []
    with path.open() as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader, None)  # header: sample, condition[, batch]
        for row in reader:
            if len(row) >= 2 and row[0]:
                out.append((row[0], row[1]))
    return out


def parse_normalized_counts(path: Path) -> dict[str, dict[str, float]]:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"m08 report input missing: {path}")
    out: dict[str, dict[str, float]] = {}
    with path.open() as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)
        samples = header[1:] if header else []
        for row in reader:
            if not row or not row[0]:
                continue
            out[row[0]] = {s: _num(v) for s, v in zip(samples, row[1:])}
    return out


def condition_layout(coldata: list[tuple[str, str]]) -> tuple[list[str], dict[str, list[str]]]:
    """Koşulları ilk-görülme sırasında döndür + koşul -> örnek listesi."""
    order: list[str] = []
    samples: dict[str, list[str]] = {}
    for sample, cond in coldata:
        if cond not in samples:
            samples[cond] = []
            order.append(cond)
        samples[cond].append(sample)
    return order, samples


def cond_mean(gene_id: str, samples: list[str], norm_counts: dict) -> float | None:
    row = norm_counts.get(gene_id)
    if not row:
        return None
    vals = [row[s] for s in samples if row.get(s) is not None]
    return sum(vals) / len(vals) if vals else None


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
    de_dir = run_dir / "differential_expression"
    de_tsv = de_dir / "deseq2_results.tsv"
    if not de_tsv.exists():
        raise FileNotFoundError(f"m08 report input missing: {de_tsv}")
    figures_dir = run_dir / "figures"
    enrich_dir = run_dir / "enrichment"
    enrich_manifest = enrich_dir / "manifest.json"
    return {
        "norm_counts": parse_normalized_counts(de_dir / "normalized_counts.tsv"),
        "coldata": parse_coldata(de_dir / "coldata.tsv"),
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
        # m09 zenginleştirme — opsiyonel: çalıştırılmadıysa None (rapor dürüstçe not düşer, kırılmaz).
        "enrichment_up": parse_enrichment_tsv(enrich_dir / "enrichment_up.tsv")
        if (enrich_dir / "enrichment_up.tsv").exists() else None,
        "enrichment_down": parse_enrichment_tsv(enrich_dir / "enrichment_down.tsv")
        if (enrich_dir / "enrichment_down.tsv").exists() else None,
        "enrichment_manifest": json.loads(enrich_manifest.read_text())
        if enrich_manifest.exists() else None,
        "enrichment_dir": enrich_dir,
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
        "mean_suffix": "ort.",
        "enrichment": "Fonksiyonel Zenginleştirme (GO)",
        "go_term": "GO terimi", "go_id": "GO id", "namespace": "Kategori",
        "fold": "Kat-zenginleşme", "study_bg": "Set / arka plan",
        "enrichment_up": "Artan genlerde zenginleşen GO terimleri",
        "enrichment_down": "Azalan genlerde zenginleşen GO terimleri",
        "no_enrichment": "Bu yönde anlamlı zenginleşen GO terimi bulunamadı.",
        "enrichment_not_run": "GO zenginleştirme bu koşuda çalıştırılmadı "
                              "(rnaforge enrich ile aynı --run-id üzerinde üretilir).",
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
        "mean_suffix": "mean",
        "enrichment": "Functional Enrichment (GO)",
        "go_term": "GO term", "go_id": "GO id", "namespace": "Category",
        "fold": "Fold enrichment", "study_bg": "Study / background",
        "enrichment_up": "GO terms enriched among up-regulated genes",
        "enrichment_down": "GO terms enriched among down-regulated genes",
        "no_enrichment": "No significantly enriched GO terms in this direction.",
        "enrichment_not_run": "GO enrichment was not run for this run "
                              "(produced by rnaforge enrich on the same --run-id).",
    },
}

FIGURE_CAPTIONS: dict[str, dict[str, str]] = {
    "tr": {
        "pca": "En değişken 500 genin ana bileşen izdüşümü. Aynı koşulun replikaları kümelenmeli; koşullar ayrışmalı.",
        "sample_correlation": "Örnekler arası Pearson korelasyonu (log2 normalize). Yüksek blok = tutarlı replikalar; sapan örnek burada görünür.",
        "expression_dist": "Örnek başına log2 normalize ekspresyon dağılımı. Kutular benzer olmalı; normalizasyonun dengeli olduğunu gösterir.",
        "dispersion": "DESeq2 dispersiyon tahmini: gen-bazlı (gri) vs uyum (mavi) vs son (turuncu). Nokta bulutu uyum eğrisine çökmeli.",
        "pval_histogram": "Ham p-değerlerinin dağılımı. Düz + 0'a yakın tepe sağlıklıdır; anormal biçim model/veri sorununa işaret eder.",
        "volcano": "log2 kat değişimi vs -log10 padj. Sağ üst = anlamlı artan, sol üst = anlamlı azalan genler.",
        "ma": "Ortalama ekspresyon vs log2 kat değişimi. Anlamlı genler renkli; düşük sayımda dağılım genişler.",
        "heatmap": "En güçlü 40 DEG'in örnek-başı z-skoru. Koşullar arası zıt renk blokları beklenir.",
    },
    "en": {
        "pca": "Principal-component projection of the 500 most variable genes. Replicates should cluster; conditions should separate.",
        "sample_correlation": "Between-sample Pearson correlation (log2 normalized). High blocks = consistent replicates; an outlier stands out here.",
        "expression_dist": "Per-sample log2 normalized expression distribution. Boxes should be similar, indicating balanced normalization.",
        "dispersion": "DESeq2 dispersion estimates: gene-wise (grey) vs fit (blue) vs final (orange). The cloud should shrink toward the fit.",
        "pval_histogram": "Distribution of raw p-values. Flat with a peak near 0 is healthy; an odd shape signals a model/data issue.",
        "volcano": "log2 fold change vs -log10 padj. Top-right = significant up, top-left = significant down genes.",
        "ma": "Mean expression vs log2 fold change. Significant genes coloured; spread widens at low counts.",
        "heatmap": "Per-sample z-scores of the top 40 DEGs. Contrasting colour blocks between conditions are expected.",
        "enrichment_up": "GO terms enriched among up-regulated genes; x = fold enrichment, point size = gene count, colour = padj.",
        "enrichment_down": "GO terms enriched among down-regulated genes; x = fold enrichment, point size = gene count, colour = padj.",
    },
}
FIGURE_CAPTIONS["tr"].update({
    "enrichment_up": "Artan genlerde zenginleşen GO terimleri; x = kat-zenginleşme, nokta boyutu = gen sayısı, renk = padj.",
    "enrichment_down": "Azalan genlerde zenginleşen GO terimleri; x = kat-zenginleşme, nokta boyutu = gen sayısı, renk = padj.",
})

SECTION_INTRO: dict[str, dict[str, str]] = {
    "tr": {
        "confidence": "Bu koşulun kalite kapılarının özeti. FAIL varsa sonuç geçersizdir; WARN varsa sonuç şüpheli damgalıdır.",
        "dataset": "Analiz edilen organizma, platform, deney tasarımı ve örnekler.",
        "quality": "Okuma işleme, hizalama ve gene atama oranları — verinin analize uygunluğu.",
        "de": "Koşullar arası diferansiyel ekspresyon özeti (DESeq2).",
        "figures": "Kalite, model tanısı ve sonuç görselleri. Her figürün altında nasıl okunacağı açıklanmıştır.",
        "table": "İstatistiksel eşiği geçen en güçlü artan ve azalan genler.",
        "enrichment": "Artan ve azalan DEG'lerde aşırı temsil edilen GO terimleri (hipergeometrik ORA, BH-FDR). "
                      "Anlamlı terimler (padj<0.05) hangi biyolojik süreçlerin değiştiğini özetler.",
        "methods": "Kullanılan araçlar ve parametreler.",
        "references": "Yöntemlerin dayandığı yayınlar.",
    },
    "en": {
        "confidence": "Summary of this run's quality gates. A FAIL invalidates the result; a WARN stamps it as suspect.",
        "dataset": "Organism, platform, experimental design and samples analysed.",
        "quality": "Read processing, alignment and gene-assignment rates — the data's fitness for analysis.",
        "de": "Summary of differential expression between conditions (DESeq2).",
        "figures": "Quality, model-diagnostic and result figures. Each figure includes how to read it.",
        "table": "The strongest up- and down-regulated genes passing the statistical threshold.",
        "enrichment": "GO terms over-represented among up- and down-regulated DEGs (hypergeometric ORA, BH-FDR). "
                      "Significant terms (padj<0.05) summarise which biological processes changed.",
        "methods": "Tools and parameters used.",
        "references": "Publications the methods are based on.",
    },
}


def _intro(section_id: str, L: dict) -> str:
    # L is a language dict; find which language by identity to pick the intro set.
    lang = "en" if L is LABELS["en"] else "tr"
    text = SECTION_INTRO[lang].get(section_id, "")
    return f'<p class="intro">{_esc(text)}</p>' if text else ""


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
        f'<section id="confidence"><h2>{_esc(L["confidence"])}</h2>{_intro("confidence", L)}'
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
    return f'<section id="dataset"><h2>{_esc(L["dataset"])}</h2>{_intro("dataset", L)}{meta}{tbl}</section>'


def section_quality(align: dict, count: dict, trimming_cfg: dict, L: dict) -> str:
    trim = (f'<p>{_esc(L["min_length"])}: {_esc(trimming_cfg.get("min_length"))} · '
            f'{_esc(L["aggressive"])}: {_esc(trimming_cfg.get("aggressive"))}</p>')
    asamp = align.get("samples", {})
    csamp = count.get("samples", {})
    rows = [[sid, _pct(asamp.get(sid, {}).get("alignment_rate")),
             _pct(csamp.get(sid, {}).get("assignment_rate"))] for sid in asamp]
    tbl = _table([L["sample"], L["alignment_rate"], L["assignment_rate"]], rows)
    return f'<section id="quality"><h2>{_esc(L["quality"])}</h2>{_intro("quality", L)}{trim}{tbl}</section>'


def section_de(de: dict, L: dict) -> str:
    n_sig = de.get("n_significant", 0)
    summary = (f'<p class="summary">{_esc(L["summary"])}: {_esc(n_sig)} / '
               f'{_esc(de.get("n_genes"))} — {_esc(de.get("contrast"))} '
               f'(FDR<{_esc(de.get("fdr_threshold"))}, |log2FC|>={_esc(de.get("log2fc_threshold"))}).</p>')
    rows = [
        [L["contrast"], de.get("contrast")],
        [L["n_genes"], de.get("n_genes")],
        [L["n_sig"], n_sig],
        [L["up"], de.get("n_up")],
        [L["down"], de.get("n_down")],
        ["min replicate corr.", de.get("min_replicate_correlation")],
    ]
    tbl = _table([" ", " "], rows)
    return f'<section id="de"><h2>{_esc(L["de"])}</h2>{_intro("de", L)}{summary}{tbl}</section>'


def section_figures(figures_manifest: dict, figures_dir: Path, L: dict, lang: str = "tr") -> str:
    figures_dir = Path(figures_dir)
    caps = FIGURE_CAPTIONS.get(lang, FIGURE_CAPTIONS["tr"])
    blocks = []
    for fig in figures_manifest.get("figures", []):
        png = figures_dir / fig["png"]
        if not png.exists():
            raise FileNotFoundError(f"m08: figure PNG missing for report: {png}")
        cap = caps.get(fig.get("id"), "")
        cap_html = (f'<figcaption><strong>{_esc(fig.get("title"))}</strong> — {_esc(cap)}</figcaption>'
                    if cap else f'<figcaption>{_esc(fig.get("title"))}</figcaption>')
        blocks.append(f'<figure><img src="{embed_png(png)}" alt="{_esc(fig.get("title"))}"/>{cap_html}</figure>')
    return f'<section id="figures"><h2>{_esc(L["figures"])}</h2>{_intro("figures", L)}{"".join(blocks)}</section>'


def _deg_table(rows: list[dict], L: dict, cond_ctx: dict | None = None) -> str:
    cond_order = cond_ctx["order"] if cond_ctx else []
    headers = [L["gene"], L["log2fc"], L["padj"]]
    headers += [f'{c} {L["mean_suffix"]}' for c in cond_order]
    headers += [L["base_mean"]]
    body = []
    for r in rows:
        row = [r["gene"],
               f'{r["log2fc"]:.2f}' if r["log2fc"] is not None else "—",
               f'{r["padj"]:.2e}' if r["padj"] is not None else "—"]
        for c in cond_order:
            m = cond_mean(r["gene_id"], cond_ctx["samples"][c], cond_ctx["norm_counts"])
            row.append(f'{m:.1f}' if m is not None else "—")
        row.append(f'{r["base_mean"]:.1f}' if r["base_mean"] is not None else "—")
        body.append(row)
    return _table(headers, body)


def section_table(de_results: list, gene_map: dict, fdr: float, lfc: float, L: dict,
                  cond_ctx: dict | None = None) -> str:
    up = top_degs_by_direction(de_results, gene_map, fdr, lfc, "Up", n=25)
    down = top_degs_by_direction(de_results, gene_map, fdr, lfc, "Down", n=25)
    if not up and not down:
        return f'<section id="table"><h2>{_esc(L["table"])}</h2>{_intro("table", L)}<p>{_esc(L["no_degs"])}</p></section>'
    up_html = (f'<h3>{_esc(L["up_table"])}</h3>{_deg_table(up, L, cond_ctx)}' if up
               else f'<h3>{_esc(L["up_table"])}</h3><p>{_esc(L["no_degs"])}</p>')
    down_html = (f'<h3>{_esc(L["down_table"])}</h3>{_deg_table(down, L, cond_ctx)}' if down
                 else f'<h3>{_esc(L["down_table"])}</h3><p>{_esc(L["no_degs"])}</p>')
    return (f'<section id="table"><h2>{_esc(L["table"])}</h2>{_intro("table", L)}{up_html}{down_html}'
            f'<p class="note">{_esc(L["full_table_note"])}</p></section>')


_NS_ORDER = ["BP", "MF", "CC"]


def _enrichment_dir_table(rows: list[dict], L: dict, top_n: int = 10) -> str:
    """Anlamlı (padj<0.05) terimleri namespace başına top_n ile tablola. Boşsa not."""
    sig = [r for r in rows if (r.get("p_adj") is not None and r["p_adj"] < 0.05)]
    if not sig:
        return f'<p>{_esc(L["no_enrichment"])}</p>'
    headers = [L["go_id"], L["go_term"], L["namespace"], L["study_bg"], L["fold"], L["padj"]]
    body: list[list] = []
    for ns in _NS_ORDER:
        grp = [r for r in sig if r.get("namespace") == ns][:top_n]
        for r in grp:
            body.append([
                r["go_id"], r["term"], r["namespace"],
                f'{r["study_count"]}/{r["bg_count"]}',
                f'{r["fold_enrichment"]:.1f}' if r["fold_enrichment"] is not None else "—",
                f'{r["p_adj"]:.2e}' if r["p_adj"] is not None else "—",
            ])
    return _table(headers, body)


def section_enrichment(inputs: dict, L: dict, lang: str = "tr") -> str:
    up, down = inputs.get("enrichment_up"), inputs.get("enrichment_down")
    if up is None and down is None:      # m09 çalıştırılmadı — dürüst not, kırılmaz
        return (f'<section id="enrichment"><h2>{_esc(L["enrichment"])}</h2>'
                f'<p class="note">{_esc(L["enrichment_not_run"])}</p></section>')
    caps = FIGURE_CAPTIONS.get(lang, FIGURE_CAPTIONS["tr"])
    figs_dir = Path(inputs.get("enrichment_dir", "."))
    fig_by_id = {f["id"]: f for f in (inputs.get("enrichment_manifest") or {}).get("figures", [])}
    blocks = []
    for direction, rows, fid in (("up", up or [], "enrichment_up"),
                                 ("down", down or [], "enrichment_down")):
        blocks.append(f'<h3>{_esc(L["enrichment_" + direction])}</h3>')
        blocks.append(_enrichment_dir_table(rows, L))
        fig = fig_by_id.get(fid)
        if fig and (figs_dir / fig["png"]).exists():
            cap = caps.get(fid, "")
            blocks.append(
                f'<figure><img src="{embed_png(figs_dir / fig["png"])}" '
                f'alt="{_esc(fig.get("title"))}"/>'
                f'<figcaption>{_esc(cap)}</figcaption></figure>')
    return (f'<section id="enrichment"><h2>{_esc(L["enrichment"])}</h2>'
            f'{_intro("enrichment", L)}{"".join(blocks)}</section>')


# Yöntem anlatısı — DESeq2 (Love ve ark. 2014) ve standart bulk RNA-seq pratiğinden;
# config parametreleriyle doldurulur. Çift dilli. {aggr} = agresif-trimming ifadesi.
_METHODS_TEXT: dict[str, str] = {
    "tr": (
        "Ham okumaların kalitesi FastQC ile değerlendirildi (taban kalitesi, adaptör ve GC içeriği). "
        "Adaptör dizileri ve kısa okumalar fastp ile kırpıldı (asgari okuma uzunluğu {min_len} nt); "
        "gen ekspresyon tahminlerini yanlıladığı için agresif kalite kırpması {aggr}. Kırpılan okumalar "
        "referans genoma Bowtie2 ile (uçtan-uca) hizalandı. Hizalanan okumalar featureCounts ile "
        "özniteliklere atanarak gen×örnek sayım matrisi oluşturuldu (öznitelik tipi {feature_type}, "
        "kimlik özniteliği {attribute}). Diferansiyel ekspresyon DESeq2 ile hesaplandı: kütüphane "
        "büyüklüğü medyan-oran yöntemiyle normalize edildi (boyut faktörleri); gen-bazlı dispersiyonlar "
        "kestirilip ortalama-dispersiyon eğilimine doğru empirical-Bayes ile büzüldü; koşul etkisi negatif "
        "binom genelleştirilmiş doğrusal modelle test edildi (Wald testi) ve p-değerleri çoklu-test için "
        "Benjamini–Hochberg (FDR) ile düzeltildi. Tasarım formülü {design}. Bir gen, düzeltilmiş p (padj) "
        "< {fdr} ve |log2 kat değişimi| ≥ {lfc} ise anlamlı diferansiyel eksprese sayıldı. Tüm görseller "
        "ggplot2 ile üretildi."
    ),
    "en": (
        "Raw read quality was assessed with FastQC (per-base quality, adapter and GC content). Adapter "
        "sequences and short reads were trimmed with fastp (minimum read length {min_len} nt); aggressive "
        "quality trimming was {aggr} because it biases gene-expression estimates. Trimmed reads were "
        "aligned to the reference genome with Bowtie2 (end-to-end). Aligned reads were assigned to "
        "features with featureCounts to build a gene×sample count matrix (feature type {feature_type}, "
        "identifier attribute {attribute}). Differential expression was computed with DESeq2: library "
        "sizes were normalized by the median-of-ratios method (size factors); gene-wise dispersions were "
        "estimated and shrunk toward the mean–dispersion trend by empirical Bayes; the condition effect "
        "was tested with a negative-binomial generalized linear model (Wald test) and p-values were "
        "adjusted for multiple testing by Benjamini–Hochberg (FDR). Design formula {design}. A gene was "
        "called significantly differentially expressed when adjusted p (padj) < {fdr} and |log2 fold "
        "change| ≥ {lfc}. All figures were produced with ggplot2."
    ),
}
_METHODS_AGGR: dict[str, dict[bool, str]] = {
    "tr": {False: "kapatıldı (Williams ve ark. 2016)", True: "uygulandı"},
    "en": {False: "disabled (Williams et al. 2016)", True: "enabled"},
}


def section_methods(config, L: dict) -> str:
    lang = "en" if L is LABELS["en"] else "tr"
    t = config.trimming
    q = config.quantification
    d = config.de
    text = _METHODS_TEXT[lang].format(
        min_len=t.min_length, aggr=_METHODS_AGGR[lang][bool(t.aggressive_quality)],
        feature_type=q.feature_type, attribute=q.attribute,
        design=d.design, fdr=d.fdr_threshold, lfc=d.log2fc_threshold,
    )
    return f'<section id="methods"><h2>{_esc(L["methods"])}</h2>{_intro("methods", L)}<p>{_esc(text)}</p></section>'


# (atıf, DOI/URL). DOI'ler doğrulandı (doi.org). FastQC bir araçtır; DOI'si yoktur → proje URL'si.
_REFERENCES: list[tuple[str, str]] = [
    ("Andrews S. FastQC: a quality control tool for high throughput sequence data. 2010.",
     "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/"),
    ("Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. "
     "Bioinformatics. 2018;34(17):i884–i890.", "https://doi.org/10.1093/bioinformatics/bty560"),
    ("Langmead B, Salzberg SL. Fast gapped-read alignment with Bowtie 2. "
     "Nat Methods. 2012;9(4):357–359.", "https://doi.org/10.1038/nmeth.1923"),
    ("Liao Y, Smyth GK, Shi W. featureCounts: an efficient general purpose program for assigning "
     "sequence reads to genomic features. Bioinformatics. 2014;30(7):923–930.",
     "https://doi.org/10.1093/bioinformatics/btt656"),
    ("Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data "
     "with DESeq2. Genome Biol. 2014;15:550.", "https://doi.org/10.1186/s13059-014-0550-8"),
    ("Williams CR, Baccarella A, Parrish JZ, Kim CC. Trimming of sequence reads alters RNA-Seq gene "
     "expression estimates. BMC Bioinformatics. 2016;17:103.",
     "https://doi.org/10.1186/s12859-016-0956-2"),
    ("Wickham H. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York; 2016.",
     "https://doi.org/10.1007/978-3-319-24277-4"),
]


def _ref_link_label(url: str) -> str:
    return url.split("doi.org/", 1)[1] if "doi.org/" in url else url


def section_references(L: dict) -> str:
    items = "".join(
        f'<li>{_esc(cite)} '
        f'<a href="{_esc(url)}" target="_blank" rel="noopener">'
        f'{"doi:" if "doi.org/" in url else ""}{_esc(_ref_link_label(url))}</a></li>'
        for cite, url in _REFERENCES
    )
    return f'<section id="references"><h2>{_esc(L["references"])}</h2>{_intro("references", L)}<ol>{items}</ol></section>'


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
.intro{color:#444;font-size:.95rem;margin:.2rem 0 .6rem}
@media print{body{max-width:none} h2{page-break-after:avoid}}
"""


def render_report(inputs: dict, config, version: str, run_id: str = "") -> str:
    lang = config.report.language
    L = LABELS.get(lang, LABELS["tr"])
    raw = inputs["raw"]
    trimming_cfg = {"min_length": config.trimming.min_length,
                    "aggressive": config.trimming.aggressive_quality}
    cond_order, cond_samples = condition_layout(inputs.get("coldata", []))
    cond_ctx = {"norm_counts": inputs.get("norm_counts", {}),
                "order": cond_order, "samples": cond_samples}
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
        section_figures(inputs["figures"], inputs["figures_dir"], L, lang),
        section_table(inputs["de_results"], inputs["gene_map"],
                      config.de.fdr_threshold, config.de.log2fc_threshold, L, cond_ctx),
        section_enrichment(inputs, L, lang),
        section_methods(config, L),
        section_references(L),
    ])
    return (f'<!doctype html><html lang="{_esc(lang)}"><head><meta charset="utf-8">'
            f'<meta name="viewport" content="width=device-width, initial-scale=1">'
            f'<title>RNAForge report</title><style>{_CSS}</style></head>'
            f'<body>{body}</body></html>')
