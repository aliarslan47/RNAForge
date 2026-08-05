"""m08 — HTML report builder. Pure Python + stdlib; assembles a single self-contained report.html
from m06/m07 output contracts. No new data gate; verdict carries over from the confidence card."""
from __future__ import annotations
import base64
import csv
import html
import json
from datetime import datetime
from pathlib import Path

N_SECTIONS = 11


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


def parse_gsea_tsv(path: Path) -> list[dict]:
    """m11 gsea_<coll>.tsv -> tipli satırlar. Yoksa/başlık-only -> boş liste."""
    path = Path(path)
    if not path.exists():
        return []
    rows: list[dict] = []
    with path.open() as f:
        for r in csv.DictReader(f, delimiter="\t"):
            row = dict(r)
            row["size"] = int(row["size"]) if row.get("size") not in (None, "") else None
            for k in ("ES", "NES", "pval", "padj"):
                row[k] = _num(row.get(k))
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
    kegg_dir = run_dir / "kegg"
    kegg_manifest = kegg_dir / "manifest.json"
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
        # m10 KEGG — opsiyonel: çalıştırılmadıysa None.
        "kegg_up": parse_enrichment_tsv(kegg_dir / "kegg_up.tsv")
        if (kegg_dir / "kegg_up.tsv").exists() else None,
        "kegg_down": parse_enrichment_tsv(kegg_dir / "kegg_down.tsv")
        if (kegg_dir / "kegg_down.tsv").exists() else None,
        "kegg_manifest": json.loads(kegg_manifest.read_text())
        if kegg_manifest.exists() else None,
        "kegg_dir": kegg_dir,
        # m11 GSEA — opsiyonel: çalıştırılmadıysa None.
        "gsea_go": parse_gsea_tsv(run_dir / "gsea" / "gsea_go.tsv")
        if (run_dir / "gsea" / "gsea_go.tsv").exists() else None,
        "gsea_kegg": parse_gsea_tsv(run_dir / "gsea" / "gsea_kegg.tsv")
        if (run_dir / "gsea" / "gsea_kegg.tsv").exists() else None,
        "gsea_manifest": json.loads((run_dir / "gsea" / "manifest.json").read_text())
        if (run_dir / "gsea" / "manifest.json").exists() else None,
        "gsea_dir": run_dir / "gsea",
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
        "go_heading": "Gene Ontology (GO)", "kegg_heading": "KEGG Yolakları",
        "kegg_up": "Artan genlerde zenginleşen KEGG yolakları",
        "kegg_down": "Azalan genlerde zenginleşen KEGG yolakları",
        "no_enrichment": "Bu yönde anlamlı zenginleşen terim bulunamadı.",
        "enrichment_not_run": "GO zenginleştirme bu koşuda çalıştırılmadı "
                              "(rnaforge enrich ile aynı --run-id üzerinde üretilir).",
        "enrichment_legend": (
            "Sütunlar — <b>GO id</b>: terim kimliği · <b>GO terimi</b>: terim adı · "
            "<b>Kategori</b>: GO alanı (<b>BP</b> = Biyolojik Süreç, <b>MF</b> = Moleküler İşlev, "
            "<b>CC</b> = Hücresel Bileşen) · <b>Set / arka plan</b>: terimdeki DEG sayısı / "
            "arka plandaki (test edilen, anotasyonlu) gen sayısı · <b>Kat-zenginleşme</b>: "
            "gözlenen oranın beklenene oranı (>1 = zenginleşme) · <b>padj</b>: Benjamini–Hochberg ile "
            "çoklu-test düzeltilmiş p-değeri. Yalnız anlamlı terimler (padj &lt; 0,05), kategori başına "
            "en güçlü ilk 10 gösterilir; tam liste enrichment/ ve kegg/ TSV dosyalarındadır. "
            "KEGG bölümünde Kategori sütunu 'KEGG' (yolak) anlamına gelir."),
        "gsea": "Gen Seti Zenginleştirme (GSEA)", "nes": "NES", "set_size": "Set boyutu",
        "leading_edge": "Öncü genler",
        "gsea_pos": "Artan tarafta zenginleşen (NES > 0)", "gsea_neg": "Azalan tarafta zenginleşen (NES < 0)",
        "gsea_not_run": "GSEA bu koşuda çalıştırılmadı (rnaforge gsea ile aynı --run-id üzerinde üretilir).",
        "gsea_legend": (
            "GSEA tüm genleri Wald istatistiğine göre sıralar ve gen setlerinin sıralamanın hangi ucunda "
            "yoğunlaştığını ölçer. <b>NES</b> = normalize zenginleşme skoru: pozitif → set artan (yüksek) "
            "genlerde, negatif → azalan genlerde zenginleşmiş. <b>Öncü genler</b> = skora en çok katkı veren "
            "çekirdek genler. Yalnız padj &lt; 0,05 gösterilir; tam liste gsea/ TSV dosyalarındadır."),
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
        "go_heading": "Gene Ontology (GO)", "kegg_heading": "KEGG Pathways",
        "kegg_up": "KEGG pathways enriched among up-regulated genes",
        "kegg_down": "KEGG pathways enriched among down-regulated genes",
        "no_enrichment": "No significantly enriched terms in this direction.",
        "enrichment_not_run": "GO enrichment was not run for this run "
                              "(produced by rnaforge enrich on the same --run-id).",
        "enrichment_legend": (
            "Columns — <b>GO id</b>: term identifier · <b>GO term</b>: term name · "
            "<b>Category</b>: GO domain (<b>BP</b> = Biological Process, <b>MF</b> = Molecular Function, "
            "<b>CC</b> = Cellular Component) · <b>Study / background</b>: DEGs in the term / genes in the "
            "background (tested, annotated) · <b>Fold enrichment</b>: observed-to-expected ratio "
            "(>1 = enriched) · <b>padj</b>: Benjamini–Hochberg multiple-testing adjusted p-value. "
            "Only significant terms (padj &lt; 0.05), top 10 per category, are shown; the full lists are in "
            "the enrichment/ and kegg/ TSV files. In the KEGG section the Category column reads 'KEGG' (pathway)."),
        "gsea": "Gene Set Enrichment (GSEA)", "nes": "NES", "set_size": "Set size",
        "leading_edge": "Leading edge",
        "gsea_pos": "Enriched on the up side (NES > 0)", "gsea_neg": "Enriched on the down side (NES < 0)",
        "gsea_not_run": "GSEA was not run for this run (produced by rnaforge gsea on the same --run-id).",
        "gsea_legend": (
            "GSEA ranks all genes by the Wald statistic and measures where a gene set concentrates along "
            "that ranking. <b>NES</b> = normalized enrichment score: positive → the set is enriched among "
            "up-regulated (high) genes, negative → among down-regulated genes. <b>Leading edge</b> = the "
            "core genes driving the score. Only padj &lt; 0.05 shown; full lists in the gsea/ TSV files."),
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
    "kegg_up": "Artan genlerde zenginleşen KEGG yolakları; x = kat-zenginleşme, nokta boyutu = gen sayısı, renk = padj.",
    "kegg_down": "Azalan genlerde zenginleşen KEGG yolakları; x = kat-zenginleşme, nokta boyutu = gen sayısı, renk = padj.",
})
FIGURE_CAPTIONS["en"].update({
    "kegg_up": "KEGG pathways enriched among up-regulated genes; x = fold enrichment, point size = gene count, colour = padj.",
    "kegg_down": "KEGG pathways enriched among down-regulated genes; x = fold enrichment, point size = gene count, colour = padj.",
})
FIGURE_CAPTIONS["tr"].update({
    "gsea_go": "GSEA (GO) işaretli NES; sağ = artan tarafta, sol = azalan tarafta zenginleşen setler. Renk = padj, boyut = set boyutu.",
    "gsea_kegg": "GSEA (KEGG) işaretli NES; sağ = artan tarafta, sol = azalan tarafta zenginleşen yolaklar. Renk = padj, boyut = set boyutu.",
})
FIGURE_CAPTIONS["en"].update({
    "gsea_go": "GSEA (GO) signed NES; right = enriched on the up side, left = on the down side. Colour = padj, size = set size.",
    "gsea_kegg": "GSEA (KEGG) signed NES; right = enriched on the up side, left = on the down side. Colour = padj, size = set size.",
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
        "gsea": "Tüm genlerin sıralı listesinde koordineli değişen gen setleri (GSEA, fgsea). "
                "ORA'nın kaçırabildiği zayıf ama tutarlı sinyalleri yakalar.",
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
        "gsea": "Gene sets that change coordinately across the full ranked gene list (GSEA, fgsea). "
                "Captures weak but consistent signals that ORA can miss.",
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
    """Anlamlı (padj<0.05) terimleri namespace başına top_n ile tablola. Boşsa not.
    Namespace GO (BP/MF/CC) veya KEGG olabilir — bilinen sıra önce, gerisi alfabetik."""
    sig = [r for r in rows if (r.get("p_adj") is not None and r["p_adj"] < 0.05)]
    if not sig:
        return f'<p>{_esc(L["no_enrichment"])}</p>'
    present = list(dict.fromkeys(r.get("namespace") for r in sig))
    ordered = [ns for ns in _NS_ORDER if ns in present] + sorted(
        ns for ns in present if ns not in _NS_ORDER)
    headers = [L["go_id"], L["go_term"], L["namespace"], L["study_bg"], L["fold"], L["padj"]]
    body: list[list] = []
    for ns in ordered:
        for r in [x for x in sig if x.get("namespace") == ns][:top_n]:
            body.append([
                r["go_id"], r["term"], r["namespace"],
                f'{r["study_count"]}/{r["bg_count"]}',
                f'{r["fold_enrichment"]:.1f}' if r["fold_enrichment"] is not None else "—",
                f'{r["p_adj"]:.2e}' if r["p_adj"] is not None else "—",
            ])
    return _table(headers, body)


def _enrichment_collection(up, down, manifest, figs_dir: Path, L: dict, lang: str,
                           up_label: str, down_label: str, cap_ids: tuple[str, str]) -> str:
    """Bir zenginleştirme koleksiyonunu (GO veya KEGG) artan/azalan tablo + figürlerle basar."""
    caps = FIGURE_CAPTIONS.get(lang, FIGURE_CAPTIONS["tr"])
    fig_by_id = {f["id"]: f for f in (manifest or {}).get("figures", [])}
    blocks = []
    for rows, lab, fid in ((up or [], up_label, cap_ids[0]),
                           (down or [], down_label, cap_ids[1])):
        blocks.append(f'<h4>{_esc(L[lab])}</h4>')
        blocks.append(_enrichment_dir_table(rows, L))
        fig = fig_by_id.get(fid)
        if fig and (Path(figs_dir) / fig["png"]).exists():
            cap = caps.get(fid, "")
            blocks.append(
                f'<figure><img src="{embed_png(Path(figs_dir) / fig["png"])}" '
                f'alt="{_esc(fig.get("title"))}"/><figcaption>{_esc(cap)}</figcaption></figure>')
    return "".join(blocks)


def section_enrichment(inputs: dict, L: dict, lang: str = "tr") -> str:
    go_up, go_down = inputs.get("enrichment_up"), inputs.get("enrichment_down")
    kegg_up, kegg_down = inputs.get("kegg_up"), inputs.get("kegg_down")
    go_ran = go_up is not None or go_down is not None
    kegg_ran = kegg_up is not None or kegg_down is not None
    if not go_ran and not kegg_ran:      # m09/m10 çalıştırılmadı — dürüst not, kırılmaz
        return (f'<section id="enrichment"><h2>{_esc(L["enrichment"])}</h2>'
                f'<p class="note">{_esc(L["enrichment_not_run"])}</p></section>')
    blocks = []
    if go_ran:
        blocks.append(f'<h3>{_esc(L["go_heading"])}</h3>')
        blocks.append(_enrichment_collection(
            go_up, go_down, inputs.get("enrichment_manifest"), inputs.get("enrichment_dir", "."),
            L, lang, "enrichment_up", "enrichment_down", ("enrichment_up", "enrichment_down")))
    if kegg_ran:
        blocks.append(f'<h3>{_esc(L["kegg_heading"])}</h3>')
        blocks.append(_enrichment_collection(
            kegg_up, kegg_down, inputs.get("kegg_manifest"), inputs.get("kegg_dir", "."),
            L, lang, "kegg_up", "kegg_down", ("kegg_up", "kegg_down")))
    # Tabloların altına sütun + kısaltma (BP/MF/CC, KEGG) açıklaması. Sabit, kontrollü HTML.
    blocks.append(f'<p class="note">{L["enrichment_legend"]}</p>')
    return (f'<section id="enrichment"><h2>{_esc(L["enrichment"])}</h2>'
            f'{_intro("enrichment", L)}{"".join(blocks)}</section>')


def _gsea_table(rows: list[dict], L: dict, top_n: int = 10) -> str:
    """İşaretli NES tablosu: en güçlü pozitif + en güçlü negatif (padj<0.05). Boşsa not."""
    sig = [r for r in rows if (r.get("padj") is not None and r["padj"] < 0.05)]
    if not sig:
        return f'<p>{_esc(L["no_enrichment"])}</p>'
    pos = sorted([r for r in sig if (r.get("NES") or 0) > 0],
                 key=lambda r: -r["NES"])[:top_n]
    neg = sorted([r for r in sig if (r.get("NES") or 0) < 0],
                 key=lambda r: r["NES"])[:top_n]
    headers = [L["go_term"], L["nes"], L["padj"], L["set_size"], L["leading_edge"]]
    body: list[list] = []
    for r in pos + neg:
        le = (r.get("leading_edge") or "").split(";")
        le_txt = ", ".join(le[:8]) + ("…" if len(le) > 8 else "")
        body.append([
            r.get("name") or r.get("pathway_id"),
            f'{r["NES"]:.2f}' if r.get("NES") is not None else "—",
            f'{r["padj"]:.2e}' if r.get("padj") is not None else "—",
            r.get("size"), le_txt,
        ])
    return _table(headers, body)


def section_gsea(inputs: dict, L: dict, lang: str = "tr") -> str:
    go, kegg = inputs.get("gsea_go"), inputs.get("gsea_kegg")
    if go is None and kegg is None:      # m11 çalıştırılmadı — dürüst not
        return (f'<section id="gsea"><h2>{_esc(L["gsea"])}</h2>'
                f'<p class="note">{_esc(L["gsea_not_run"])}</p></section>')
    caps = FIGURE_CAPTIONS.get(lang, FIGURE_CAPTIONS["tr"])
    figs_dir = Path(inputs.get("gsea_dir", "."))
    fig_by_id = {f["id"]: f for f in (inputs.get("gsea_manifest") or {}).get("figures", [])}
    blocks = []
    for coll, rows, heading, fid in (("go", go, "go_heading", "gsea_go"),
                                     ("kegg", kegg, "kegg_heading", "gsea_kegg")):
        if rows is None:
            continue
        blocks.append(f'<h3>{_esc(L[heading])}</h3>')
        blocks.append(_gsea_table(rows, L))
        fig = fig_by_id.get(fid)
        if fig and (figs_dir / fig["png"]).exists():
            cap = caps.get(fid, "")
            blocks.append(
                f'<figure><img src="{embed_png(figs_dir / fig["png"])}" '
                f'alt="{_esc(fig.get("title"))}"/><figcaption>{_esc(cap)}</figcaption></figure>')
    blocks.append(f'<p class="note">{L["gsea_legend"]}</p>')
    return (f'<section id="gsea"><h2>{_esc(L["gsea"])}</h2>'
            f'{_intro("gsea", L)}{"".join(blocks)}</section>')


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

# GO ORA yöntem anlatısı — yalnız m09 çalıştıysa eklenir (koşmadıysa yazmak yanıltıcı olur).
_ENRICHMENT_METHODS: dict[str, str] = {
    "tr": (
        "Diferansiyel eksprese genlerin fonksiyonel yorumu için Gene Ontology (GO) aşırı-temsil "
        "analizi (ORA) uygulandı. Artan ve azalan gen setleri ayrı ayrı, tek-yönlü hipergeometrik "
        "testle değerlendirildi (arka plan: test edilen, GO anotasyonu bulunan tüm genler; en az "
        "{min_term} gene anotlı terimler dikkate alındı). GO anotasyonları birincil olarak referans "
        "genom GFF'inden alındı ve EBI-GOA anotasyonlarıyla (yalnız anotasyonsuz genler, benzersiz gen "
        "sembolü eşleşmesiyle) tamamlandı; terimler Gene Ontology yapısı (go-basic) boyunca ata "
        "terimlere yayıldı (is_a ve part_of ilişkileri). p-değerleri her GO kategorisi (biyolojik "
        "süreç, moleküler işlev, hücresel bileşen) içinde Benjamini–Hochberg yöntemiyle düzeltildi; "
        "düzeltilmiş p < 0,05 olan terimler zenginleşmiş sayıldı."
    ),
    "en": (
        "For functional interpretation of the differentially expressed genes, Gene Ontology (GO) "
        "over-representation analysis (ORA) was performed. Up- and down-regulated gene sets were "
        "tested separately with a one-sided hypergeometric test (background: all tested genes carrying "
        "a GO annotation; terms annotated to at least {min_term} genes were considered). GO annotations "
        "were taken primarily from the reference genome GFF and supplemented with EBI-GOA annotations "
        "(only for otherwise unannotated genes, matched by unique gene symbol); terms were propagated to "
        "ancestor terms along the Gene Ontology structure (go-basic; is_a and part_of relationships). "
        "P-values were adjusted within each GO category (biological process, molecular function, "
        "cellular component) by the Benjamini–Hochberg method; terms with adjusted p < 0.05 were "
        "considered enriched."
    ),
}


_KEGG_METHODS: dict[str, str] = {
    "tr": (
        "Aynı over-representation çerçevesi KEGG yolakları için de uygulandı: gen→yolak eşlemesi KEGG "
        "veritabanından (organizma {org}) alındı ve gen sembolü ile referans genoma eşlendi; çok geniş "
        "genel/özet haritalar dışlandı. Zenginleşme hipergeometrik testle hesaplanıp Benjamini–Hochberg "
        "ile düzeltildi (padj < 0,05)."
    ),
    "en": (
        "The same over-representation framework was applied to KEGG pathways: gene-to-pathway mapping "
        "was taken from the KEGG database (organism {org}) and matched to the reference genome by gene "
        "symbol; broad global/overview maps were excluded. Enrichment was computed with the hypergeometric "
        "test and adjusted by Benjamini–Hochberg (padj < 0.05)."
    ),
}


_GSEA_METHODS: dict[str, str] = {
    "tr": (
        "Gen seti zenginleştirme analizi (GSEA) genlerin DESeq2 Wald istatistiğine göre sıralanmış tam "
        "listesi üzerinde fgsea (multilevel) ile yapıldı; GO ve KEGG gen setleri kullanıldı (set boyutu "
        "{min_size}–{max_size} genle sınırlandı). Normalize zenginleşme skorunun (NES) işareti yönü verir "
        "(pozitif = artan tarafta) ve p-değerleri Benjamini–Hochberg ile düzeltildi."
    ),
    "en": (
        "Gene set enrichment analysis (GSEA) was run over the full gene list ranked by the DESeq2 Wald "
        "statistic with fgsea (multilevel), using GO and KEGG gene sets (set size restricted to "
        "{min_size}–{max_size} genes). The sign of the normalized enrichment score (NES) gives direction "
        "(positive = up side) and p-values were adjusted by Benjamini–Hochberg."
    ),
}


def section_methods(config, L: dict, enrichment_ran: bool = False, kegg_ran: bool = False,
                    gsea_ran: bool = False) -> str:
    lang = "en" if L is LABELS["en"] else "tr"
    t = config.trimming
    q = config.quantification
    d = config.de
    text = _METHODS_TEXT[lang].format(
        min_len=t.min_length, aggr=_METHODS_AGGR[lang][bool(t.aggressive_quality)],
        feature_type=q.feature_type, attribute=q.attribute,
        design=d.design, fdr=d.fdr_threshold, lfc=d.log2fc_threshold,
    )
    paras = f'<p>{_esc(text)}</p>'
    if enrichment_ran:
        go = _ENRICHMENT_METHODS[lang].format(min_term=config.enrichment.min_term_size)
        paras += f'<p>{_esc(go)}</p>'
    if kegg_ran:
        kegg = _KEGG_METHODS[lang].format(org=config.enrichment.kegg_organism or "—")
        paras += f'<p>{_esc(kegg)}</p>'
    if gsea_ran:
        g = _GSEA_METHODS[lang].format(min_size=config.enrichment.gsea_min_size,
                                       max_size=config.enrichment.gsea_max_size)
        paras += f'<p>{_esc(g)}</p>'
    return f'<section id="methods"><h2>{_esc(L["methods"])}</h2>{_intro("methods", L)}{paras}</section>'


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

# GO ORA kaynakları — yalnız m09 çalıştıysa eklenir (DOI'ler doğrulandı).
_ENRICHMENT_REFERENCES: list[tuple[str, str]] = [
    ("Ashburner M, Ball CA, Blake JA, et al. Gene ontology: tool for the unification of biology. "
     "Nat Genet. 2000;25(1):25–29.", "https://doi.org/10.1038/75556"),
    ("The Gene Ontology Consortium. The Gene Ontology knowledgebase in 2023. "
     "Genetics. 2023;224(1):iyad031.", "https://doi.org/10.1093/genetics/iyad031"),
    ("Huntley RP, Sawford T, Mutowo-Meullenet P, et al. The GOA database: Gene Ontology annotation "
     "updates for 2015. Nucleic Acids Res. 2015;43(D1):D1057–D1063.",
     "https://doi.org/10.1093/nar/gku1113"),
    ("Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach "
     "to multiple testing. J R Stat Soc Series B. 1995;57(1):289–300.",
     "https://doi.org/10.1111/j.2517-6161.1995.tb02031.x"),
]

# KEGG kaynağı — yalnız m10 çalıştıysa (DOI doğrulandı).
_KEGG_REFERENCES: list[tuple[str, str]] = [
    ("Kanehisa M, Goto S. KEGG: Kyoto Encyclopedia of Genes and Genomes. "
     "Nucleic Acids Res. 2000;28(1):27–30.", "https://doi.org/10.1093/nar/28.1.27"),
]

# GSEA kaynakları — yalnız m11 çalıştıysa (DOI doğrulandı).
_GSEA_REFERENCES: list[tuple[str, str]] = [
    ("Subramanian A, Tamayo P, Mootha VK, et al. Gene set enrichment analysis: a knowledge-based "
     "approach for interpreting genome-wide expression profiles. Proc Natl Acad Sci USA. "
     "2005;102(43):15545–15550.", "https://doi.org/10.1073/pnas.0506580102"),
    ("Korotkevich G, Sukhov V, Budin N, et al. Fast gene set enrichment analysis. bioRxiv. 2021.",
     "https://doi.org/10.1101/060012"),
]


def _ref_link_label(url: str) -> str:
    return url.split("doi.org/", 1)[1] if "doi.org/" in url else url


def section_references(L: dict, enrichment_ran: bool = False, kegg_ran: bool = False,
                       gsea_ran: bool = False) -> str:
    refs = (_REFERENCES + (_ENRICHMENT_REFERENCES if enrichment_ran else [])
            + (_KEGG_REFERENCES if kegg_ran else [])
            + (_GSEA_REFERENCES if gsea_ran else []))
    items = "".join(
        f'<li>{_esc(cite)} '
        f'<a href="{_esc(url)}" target="_blank" rel="noopener">'
        f'{"doi:" if "doi.org/" in url else ""}{_esc(_ref_link_label(url))}</a></li>'
        for cite, url in refs
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
    enrichment_ran = (inputs.get("enrichment_up") is not None
                      or inputs.get("enrichment_down") is not None)
    kegg_ran = (inputs.get("kegg_up") is not None
                or inputs.get("kegg_down") is not None)
    gsea_ran = (inputs.get("gsea_go") is not None
                or inputs.get("gsea_kegg") is not None)
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
        section_gsea(inputs, L, lang),
        section_methods(config, L, enrichment_ran, kegg_ran, gsea_ran),
        section_references(L, enrichment_ran, kegg_ran, gsea_ran),
    ])
    return (f'<!doctype html><html lang="{_esc(lang)}"><head><meta charset="utf-8">'
            f'<meta name="viewport" content="width=device-width, initial-scale=1">'
            f'<title>RNAForge report</title><style>{_CSS}</style></head>'
            f'<body>{body}</body></html>')
