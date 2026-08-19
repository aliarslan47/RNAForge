"""m08 — HTML report builder. Pure Python + stdlib; assembles a single self-contained report.html
from m06/m07 output contracts. No new data gate; verdict carries over from the confidence card."""
from __future__ import annotations
import base64
import csv
import html
import json
import re
from datetime import datetime
from pathlib import Path

N_SECTIONS = 17


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


def parse_reduced_tsv(path: Path) -> list[dict]:
    """m12 reduced_*.tsv -> tipli satırlar (go_id, namespace, term, padj, n_collapsed, members)."""
    path = Path(path)
    if not path.exists():
        return []
    rows: list[dict] = []
    with path.open() as f:
        for r in csv.DictReader(f, delimiter="\t"):
            row = dict(r)
            row["padj"] = _num(row.get("padj"))
            row["n_collapsed"] = int(row["n_collapsed"]) if row.get("n_collapsed") else 1
            rows.append(row)
    return rows


def parse_amr_tsv(path: Path) -> list[dict]:
    """m13 amr_genes.tsv / virulence_genes.tsv -> tipli satırlar."""
    path = Path(path)
    if not path.exists():
        return []
    rows: list[dict] = []
    with path.open() as f:
        for r in csv.DictReader(f, delimiter="\t"):
            row = dict(r)
            for k in ("pct_identity", "pct_coverage", "log2fc", "padj"):
                row[k] = _num(row.get(k))
            rows.append(row)
    return rows


def parse_operon_tsv(path: Path) -> list[dict]:
    """m14 operons.tsv -> tipli satırlar."""
    path = Path(path)
    if not path.exists():
        return []
    rows: list[dict] = []
    with path.open() as f:
        for r in csv.DictReader(f, delimiter="\t"):
            row = dict(r)
            for k in ("size", "n_tested", "n_deg", "n_up", "n_down"):
                row[k] = int(row[k]) if row.get(k) not in (None, "") else 0
            row["mean_log2fc"] = _num(row.get("mean_log2fc"))
            row["coordinated"] = row.get("coordinated") == "yes"
            rows.append(row)
    return rows


def parse_community_tsv(path: Path) -> list[dict]:
    """m15 communities.tsv -> tipli satırlar."""
    path = Path(path)
    if not path.exists():
        return []
    rows: list[dict] = []
    with path.open() as f:
        for r in csv.DictReader(f, delimiter="\t"):
            row = dict(r)
            for k in ("size", "n_up", "n_down"):
                row[k] = int(row[k]) if row.get(k) not in (None, "") else 0
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
        "figure_errors": collect_figure_errors(run_dir),
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
        # m12 REVIGO — opsiyonel: çalıştırılmadıysa None.
        "reduced_ora_up": parse_reduced_tsv(run_dir / "semantic" / "reduced_ora_up.tsv")
        if (run_dir / "semantic" / "reduced_ora_up.tsv").exists() else None,
        "reduced_ora_down": parse_reduced_tsv(run_dir / "semantic" / "reduced_ora_down.tsv")
        if (run_dir / "semantic" / "reduced_ora_down.tsv").exists() else None,
        "reduced_gsea_go": parse_reduced_tsv(run_dir / "semantic" / "reduced_gsea_go.tsv")
        if (run_dir / "semantic" / "reduced_gsea_go.tsv").exists() else None,
        "semantic_manifest": json.loads((run_dir / "semantic" / "manifest.json").read_text())
        if (run_dir / "semantic" / "manifest.json").exists() else None,
        "semantic_dir": run_dir / "semantic",
        # m13 AMR/virülans — opsiyonel: çalıştırılmadıysa None.
        "amr_genes": parse_amr_tsv(run_dir / "amr" / "amr_genes.tsv")
        if (run_dir / "amr" / "amr_genes.tsv").exists() else None,
        "virulence_genes": parse_amr_tsv(run_dir / "amr" / "virulence_genes.tsv")
        if (run_dir / "amr" / "virulence_genes.tsv").exists() else None,
        "amr_stats": _read_json(stats / "amr_statistics.json")
        if (stats / "amr_statistics.json").exists() else None,
        # m14 operon — opsiyonel: çalıştırılmadıysa None.
        "operons": parse_operon_tsv(run_dir / "operon" / "operons.tsv")
        if (run_dir / "operon" / "operons.tsv").exists() else None,
        "operon_stats": _read_json(stats / "operon_statistics.json")
        if (stats / "operon_statistics.json").exists() else None,
        "operon_manifest": json.loads((run_dir / "operon" / "manifest.json").read_text())
        if (run_dir / "operon" / "manifest.json").exists() else None,
        "operon_dir": run_dir / "operon",
        # m15 PPI — opsiyonel: çalıştırılmadıysa None.
        "communities": parse_community_tsv(run_dir / "ppi" / "communities.tsv")
        if (run_dir / "ppi" / "communities.tsv").exists() else None,
        "ppi_stats": _read_json(stats / "ppi_statistics.json")
        if (stats / "ppi_statistics.json").exists() else None,
        # m16 seqqc — opsiyonel: çalıştırılmadıysa None.
        "seqqc": _read_json(stats / "seqqc_statistics.json")
        if (stats / "seqqc_statistics.json").exists() else None,
        # m17 hizalama-sonrası QC (insert-size/read-dist/coverage) — opsiyonel.
        "alignqc": _read_json(stats / "alignqc_statistics.json")
        if (stats / "alignqc_statistics.json").exists() else None,
        # m18 MultiQC toplu görünüm — opsiyonel.
        "multiqc": _read_json(stats / "multiqc_statistics.json")
        if (stats / "multiqc_statistics.json").exists() else None,
        "ppi_manifest": json.loads((run_dir / "ppi" / "manifest.json").read_text())
        if (run_dir / "ppi" / "manifest.json").exists() else None,
        "ppi_dir": run_dir / "ppi",
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
        "rrna_pct": "rRNA %", "rrna_mean": "Ortalama rRNA", "strandedness": "Strandedness (çıkarılan)",
        "dedup_pct": "Benzersiz % (dedup)", "rd_group": "Genomik bölge", "rd_pct": "Okuma %",
        "cap_per_base": "Şekil — Per-base baz kompozisyonu (A/T/G/C, örnek ortalaması)",
        "cap_insert_size": "Şekil — Insert-size (fragment uzunluğu) dağılımı",
        "cap_coverage": "Şekil — Kontig başına ortalama okuma derinliği (coverage)",
        "cap_read_dist": "Okumaların genomik özelliklere dağılımı (RSeQC read distribution)",
        "multiqc_note": "Tüm araç çıktılarının toplu MultiQC görünümü:",
        "organism": "Organizma", "platform": "Platform", "read_type": "Okuma tipi", "design": "Tasarım",
        "gate": "Kapı", "status": "Durum", "measured": "Ölçülen", "threshold": "Eşik",
        "profile": "Profil", "contrast": "Kontrast", "n_genes": "Gen sayısı",
        "n_sig": "Anlamlı gen", "gene": "Gen", "log2fc": "log2FC", "padj": "padj",
        "base_mean": "baseMean", "direction": "Yön", "up": "Artan", "down": "Azalan",
        "full_table_note": "Tam tablo: differential_expression/deseq2_results.tsv",
        "min_length": "Min uzunluk", "aggressive": "Agresif kalite trimming",
        "summary": "Özet",
        "cap_gates": "Kalite kapıları", "cap_samples": "Örnekler",
        "cap_quality": "Hizalama, atama ve rRNA oranları", "cap_de": "Diferansiyel ekspresyon özeti",
        "cap_operon": "Koordineli operonlar", "cap_ppi": "Protein etkileşim modülleri",
        "software": "Yazılım ve Veritabanları", "sw_tool": "Yazılım", "sw_version": "Sürüm",
        "sw_purpose": "Amaç", "db_name": "Veritabanı", "db_version": "Sürüm / kaynak", "db_purpose": "Amaç",
        "cap_software": "Kullanılan yazılımlar", "cap_database": "Kullanılan veritabanları",
        "expr_note": "Gen ekspresyon değerleri TPM ve FPKM olarak quantification/tpm.tsv ve fpkm.tsv'de üretildi.",
        "up_table": "En Güçlü 25 Artan (Up)", "down_table": "En Güçlü 25 Azalan (Down)",
        "mean_suffix": "ort.",
        "enrichment": "Fonksiyonel Zenginleştirme (GO)",
        "go_section": "GO Zenginleştirme (ORA)", "kegg_section": "KEGG Yolak Zenginleştirme (ORA)",
        "kegg_not_run": "KEGG zenginleştirme bu koşuda çalıştırılmadı (rnaforge kegg ile üretilir).",
        "kegg_legend": (
            "Artan/azalan DEG'lerde aşırı temsil edilen KEGG yolakları (hipergeometrik ORA + BH-FDR). "
            "<b>Set / arka plan</b>: yolaktaki DEG / arka plandaki gen; <b>Kat-zenginleşme</b>: gözlenen/beklenen; "
            "<b>padj</b>: BH-düzeltilmiş p. Yalnız anlamlı (padj &lt; 0,05), en güçlü ilk 10; tam liste kegg/ TSV'de."),
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
        "semantic": "Anlamsal İndirgeme (REVIGO)", "n_collapsed": "Temsil ettiği terim",
        "sem_ora_up": "Artan GO terimleri (temsilciler)", "sem_ora_down": "Azalan GO terimleri (temsilciler)",
        "sem_gsea_go": "GSEA GO terimleri (temsilciler)",
        "sem_summary": "{n} anlamlı terim → {m} temsilci",
        "semantic_not_run": "Anlamsal indirgeme bu koşuda çalıştırılmadı (rnaforge semantic ile üretilir).",
        "semantic_legend": (
            "Fazlalık GO terimleri (parent/child, benzer süreçler) Lin semantik benzerliğiyle kümelenip "
            "her kümeden en iyi padj'li <b>temsilci</b> tutulur (REVIGO fikri). <b>Temsil ettiği terim</b> = "
            "o temsilcinin altında toplanan terim sayısı. Namespace (BP/MF/CC) ayrı işlenir; tam üye "
            "listesi semantic/ TSV dosyalarındadır."),
        "amr": "Direnç ve Virülans (AMR / VFDB)", "amr_genes": "Direnç genleri (AMR)",
        "vir_genes": "Virülans genleri (VFDB)", "de_status": "DE durumu",
        "identity": "%kimlik", "amr_label": "Sınıf / faktör",
        "card_col": "CARD", "amrfinder_col": "AMRFinderPlus",
        "amr_not_run": "AMR/virülans taraması bu koşuda çalıştırılmadı (rnaforge amr ile üretilir).",
        "amr_more": "… ve {n} DE-olmayan gen daha (tam liste amr/ TSV dosyalarında).",
        "amr_legend": (
            "Direnç genleri iki bağımsız araçla <b>yan yana</b> gösterilir: <b>CARD</b> (abricate; geniş, "
            "intrinsik efflux dahil) ve <b>AMRFinderPlus</b> (NCBI; küratörlü, yüksek-güven). Her ikisinin "
            "bulduğu gen konkordans işareti; yalnız birinin bulduğu araç farkını gösterir. Virülans <b>VFDB</b> "
            "(abricate). Genler koordinatla locus_tag'e eşlenip <b>DE durumu</b> eklendi. Veritabanları "
            "araçlarla paket gelir (abricate-get_db / amrfinder -u ile yenilenir)."),
        "operon": "Operon Analizi", "operon_genes": "Genler", "operon_size": "Boyut",
        "operon_dir": "Yön", "operon_not_run": "Operon analizi bu koşuda çalıştırılmadı "
                                               "(rnaforge operon ile üretilir).",
        "operon_summary": "{n} operon tahmin edildi · {m} çok-genli · {k} koordineli DE",
        "operon_legend": (
            "Operonlar <b>tahmin</b>: aynı yönde bitişik ve intergenik boşluğu ≤ {gap} bp olan genler aynı "
            "transkripsiyon birimi sayıldı (Moreno-Hagelsieb ve Collado-Vides 2002). Deneysel doğrulanmamıştır. "
            "<b>Koordineli</b> = ≥2 geni ve ≥2 DEG'i aynı yönde değişen operon (birlikte-düzenlenen yanıt). "
            "Tam liste operon/operons.tsv'de."),
        "ppi": "Protein Etkileşim Modülleri (STRING)", "ppi_module": "Modül", "ppi_size": "Boyut",
        "ppi_not_run": "PPI/community analizi bu koşuda çalıştırılmadı (rnaforge ppi ile üretilir).",
        "ppi_summary": "{net}/{deg} DEG ağda · {edges} kenar · {mods} modül",
        "ppi_legend": (
            "Diferansiyel eksprese genler STRING protein-etkileşim kenarlarıyla (combined score ≥ {score}) "
            "bağlandı ve ağ Louvain yöntemiyle modüllere ayrıldı. STRING etkileşimleri kanıt-skorlu "
            "tahminlerdir, hepsi deneysel doğrulanmış değildir. Tam liste ppi/communities.tsv'de."),
    },
    "en": {
        "confidence": "Confidence Card", "dataset": "Dataset and Samples",
        "quality": "Quality and Processing", "de": "Differential Expression",
        "figures": "Figures", "table": "Top DEGs", "methods": "Methods",
        "references": "References", "verdict": "Verdict", "no_degs": "No significant DEGs found.",
        "sample": "Sample", "condition": "Condition", "batch": "Batch", "paired": "Paired",
        "read_len": "Mean read length", "quality_col": "Mean quality",
        "alignment_rate": "Alignment rate", "assignment_rate": "Assignment rate",
        "rrna_pct": "rRNA %", "rrna_mean": "Mean rRNA", "strandedness": "Strandedness (inferred)",
        "dedup_pct": "Unique % (dedup)", "rd_group": "Genomic region", "rd_pct": "Read %",
        "cap_per_base": "Figure — Per-base sequence composition (A/T/G/C, sample mean)",
        "cap_insert_size": "Figure — Insert-size (fragment length) distribution",
        "cap_coverage": "Figure — Mean read depth (coverage) per contig",
        "cap_read_dist": "Read distribution across genomic features (RSeQC)",
        "multiqc_note": "Aggregate MultiQC view of all tool outputs:",
        "organism": "Organism", "platform": "Platform", "read_type": "Read type", "design": "Design",
        "gate": "Gate", "status": "Status", "measured": "Measured", "threshold": "Threshold",
        "profile": "Profile", "contrast": "Contrast", "n_genes": "Genes",
        "n_sig": "Significant genes", "gene": "Gene", "log2fc": "log2FC", "padj": "padj",
        "base_mean": "baseMean", "direction": "Direction", "up": "Up", "down": "Down",
        "full_table_note": "Full table: differential_expression/deseq2_results.tsv",
        "min_length": "Min length", "aggressive": "Aggressive quality trimming",
        "summary": "Summary",
        "cap_gates": "Quality gates", "cap_samples": "Samples",
        "cap_quality": "Alignment, assignment and rRNA rates", "cap_de": "Differential expression summary",
        "cap_operon": "Coordinated operons", "cap_ppi": "Protein interaction modules",
        "software": "Software and Databases", "sw_tool": "Software", "sw_version": "Version",
        "sw_purpose": "Purpose", "db_name": "Database", "db_version": "Version / source", "db_purpose": "Purpose",
        "cap_software": "Software used", "cap_database": "Databases used",
        "expr_note": "Gene expression values were produced as TPM and FPKM in quantification/tpm.tsv and fpkm.tsv.",
        "up_table": "Top 25 Up-regulated", "down_table": "Top 25 Down-regulated",
        "mean_suffix": "mean",
        "enrichment": "Functional Enrichment (GO)",
        "go_section": "GO Enrichment (ORA)", "kegg_section": "KEGG Pathway Enrichment (ORA)",
        "kegg_not_run": "KEGG enrichment was not run for this run (produced by rnaforge kegg).",
        "kegg_legend": (
            "KEGG pathways over-represented among up/down DEGs (hypergeometric ORA + BH-FDR). "
            "<b>Study / background</b>: DEGs in the pathway / genes in the background; <b>Fold enrichment</b>: "
            "observed/expected; <b>padj</b>: BH-adjusted p. Only significant (padj &lt; 0.05), top 10; full list in kegg/ TSV."),
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
        "amr": "Resistance and Virulence (AMR / VFDB)", "amr_genes": "Resistance genes (AMR)",
        "vir_genes": "Virulence genes (VFDB)", "de_status": "DE status",
        "identity": "%identity", "amr_label": "Class / factor",
        "card_col": "CARD", "amrfinder_col": "AMRFinderPlus",
        "amr_not_run": "AMR/virulence scan was not run for this run (produced by rnaforge amr).",
        "amr_more": "… and {n} more non-DE genes (full list in the amr/ TSV files).",
        "amr_legend": (
            "Resistance genes are shown <b>side by side</b> from two independent tools: <b>CARD</b> "
            "(abricate; broad, includes intrinsic efflux) and <b>AMRFinderPlus</b> (NCBI; curated, "
            "high-confidence). A gene found by both is a concordance mark; a gene found by only one shows the "
            "tool difference. Virulence uses <b>VFDB</b> (abricate). Hits were mapped to locus tags by "
            "coordinate and annotated with their <b>DE status</b>. Databases ship with the tools "
            "(refresh with abricate-get_db / amrfinder -u)."),
        "operon": "Operon Analysis", "operon_genes": "Genes", "operon_size": "Size",
        "operon_dir": "Direction", "operon_not_run": "Operon analysis was not run for this run "
                                                     "(produced by rnaforge operon).",
        "operon_summary": "{n} operons predicted · {m} multi-gene · {k} coordinated DE",
        "operon_legend": (
            "Operons are <b>predicted</b>: genes on the same strand with an intergenic gap ≤ {gap} bp were "
            "taken as one transcription unit (Moreno-Hagelsieb and Collado-Vides 2002). Not experimentally "
            "validated. <b>Coordinated</b> = an operon with ≥2 genes and ≥2 DEGs changing in the same direction "
            "(a co-regulated response). Full list in operon/operons.tsv."),
        "ppi": "Protein Interaction Modules (STRING)", "ppi_module": "Module", "ppi_size": "Size",
        "ppi_not_run": "PPI/community analysis was not run for this run (produced by rnaforge ppi).",
        "ppi_summary": "{net}/{deg} DEGs in the network · {edges} edges · {mods} modules",
        "ppi_legend": (
            "The differentially expressed genes were connected by STRING protein-interaction edges (combined "
            "score ≥ {score}) and the network was partitioned into modules by the Louvain method. STRING "
            "interactions are evidence-scored predictions, not all experimentally validated. Full list in "
            "ppi/communities.tsv."),
        "semantic": "Semantic Reduction (REVIGO)", "n_collapsed": "Terms represented",
        "sem_ora_up": "Up GO terms (representatives)", "sem_ora_down": "Down GO terms (representatives)",
        "sem_gsea_go": "GSEA GO terms (representatives)",
        "sem_summary": "{n} significant terms → {m} representatives",
        "semantic_not_run": "Semantic reduction was not run for this run (produced by rnaforge semantic).",
        "semantic_legend": (
            "Redundant GO terms (parent/child, similar processes) are clustered by Lin semantic similarity "
            "and the best-padj <b>representative</b> is kept per cluster (the REVIGO idea). <b>Terms "
            "represented</b> = how many terms collapsed into that representative. Namespaces (BP/MF/CC) are "
            "processed separately; full member lists are in the semantic/ TSV files."),
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
        "kegg": "Artan ve azalan DEG'lerde aşırı temsil edilen KEGG yolakları (hipergeometrik ORA, BH-FDR). "
                "Değişen metabolik/sinyal yolaklarını özetler.",
        "gsea": "Tüm genlerin sıralı listesinde koordineli değişen gen setleri (GSEA, fgsea). "
                "ORA'nın kaçırabildiği zayıf ama tutarlı sinyalleri yakalar.",
        "semantic": "Zenginleşen GO terimlerinin fazlalığı, semantik benzerlikle temsilcilere indirgenmiş "
                    "öz görünümü. Uzun listeleri okunur kılar.",
        "amr": "Suşun taşıdığı antibiyotik direnç (CARD) ve virülans (VFDB) genleri; her biri için "
               "diferansiyel ekspresyon durumu — tedavi altında indüklenen direnç/virülans yanıtını gösterir.",
        "operon": "Tahmin edilen operonlar ve içlerindeki genlerin koordineli değişimi — birlikte "
                  "düzenlenen transkripsiyon birimlerini gen-düzeyinin üstünde gösterir.",
        "ppi": "DEG'lerin STRING protein etkileşim ağındaki modülleri — birlikte işleyen fonksiyonel "
               "gen gruplarını (kompleksler/yolaklar) ağ topolojisinden çıkarır.",
        "software": "Bu koşuda kullanılan yazılımların sürümleri ve referans veritabanları (tekrarlanabilirlik için).",
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
        "kegg": "KEGG pathways over-represented among up- and down-regulated DEGs (hypergeometric ORA, BH-FDR). "
                "Summarises which metabolic/signalling pathways changed.",
        "gsea": "Gene sets that change coordinately across the full ranked gene list (GSEA, fgsea). "
                "Captures weak but consistent signals that ORA can miss.",
        "semantic": "A concise view of the enriched GO terms, reduced to representatives by semantic "
                    "similarity. Makes long lists readable.",
        "amr": "Antibiotic-resistance (CARD) and virulence (VFDB) genes carried by the strain, each with "
               "its differential-expression status — showing the resistance/virulence response under treatment.",
        "operon": "Predicted operons and the coordinated change of their genes — revealing co-regulated "
                  "transcription units above the gene level.",
        "ppi": "Modules of the DEGs in the STRING protein-interaction network — revealing co-functioning "
               "gene groups (complexes/pathways) from network topology.",
        "software": "Versions of the software used in this run and the reference databases (for reproducibility).",
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


def _table(headers: list[str], rows: list[list], caption: str = "") -> str:
    head = "".join(f"<th>{_esc(h)}</th>" for h in headers)
    body = "".join("<tr>" + "".join(f"<td>{_esc(c)}</td>" for c in row) + "</tr>" for row in rows)
    cap = f"<caption>{_esc(caption)}</caption>" if caption else ""
    return f"<table>{cap}<thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


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
    gate_tbl = _table([L["gate"], L["status"], L["measured"], L["threshold"]], gate_rows, L["cap_gates"])
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
    # read_type rozeti (kısa/uzun okuma) — hangi araç zincirinin koştuğunu gösterir.
    rt = raw.get("read_type")
    lang_tr = L is LABELS["tr"]
    rt_gloss = {"short": ("kısa" if lang_tr else "short"),
                "long": ("uzun" if lang_tr else "long")}
    rt_disp = (f'{_esc(L["read_type"])}: {_esc(rt_gloss.get(rt, rt))} · ' if rt else "")
    meta = (f'<p>{_esc(L["organism"])}: {_esc(raw.get("organism"))} · '
            f'{_esc(L["platform"])}: {_esc(raw.get("platform"))} · '
            f'{rt_disp}'
            f'{_esc(L["design"])}: {_esc(raw.get("design"))} · {_esc(conds)}</p>')
    rows = [[s.get("sample_id"), s.get("condition"), s.get("batch"), s.get("paired"),
             s.get("mean_read_length"), s.get("mean_quality")] for s in raw.get("samples", [])]
    tbl = _table([L["sample"], L["condition"], L["batch"], L["paired"], L["read_len"], L["quality_col"]],
                 rows, L["cap_samples"])
    return f'<section id="dataset"><h2>{_esc(L["dataset"])}</h2>{_intro("dataset", L)}{meta}{tbl}</section>'


def section_quality(align: dict, count: dict, trimming_cfg: dict, L: dict,
                    seqqc: dict | None = None, qc: dict | None = None,
                    figures_dir: Path | None = None, alignqc: dict | None = None,
                    multiqc: dict | None = None) -> str:
    trim = (f'<p>{_esc(L["min_length"])}: {_esc(trimming_cfg.get("min_length"))} · '
            f'{_esc(L["aggressive"])}: {_esc(trimming_cfg.get("aggressive"))}</p>')
    asamp = align.get("samples", {})
    csamp = count.get("samples", {})
    seqqc = seqqc or {}
    qc = qc or {}
    rrna = (seqqc.get("rrna_per_sample") or {})
    dedup = (qc.get("deduplication") or {})
    rows = [[sid, _pct(asamp.get(sid, {}).get("alignment_rate")),
             _pct(csamp.get(sid, {}).get("assignment_rate")),
             _pct(rrna[sid]) if sid in rrna else "—",
             (f'{dedup[sid]:.1f}%' if isinstance(dedup.get(sid), (int, float)) else "—")]
            for sid in asamp]
    tbl = _table([L["sample"], L["alignment_rate"], L["assignment_rate"], L["rrna_pct"], L["dedup_pct"]],
                 rows, L["cap_quality"])
    # rRNA% + strandedness özet satırı (m16 çalıştıysa)
    seq_line = ""
    if seqqc:
        match = ("uyumlu" if seqqc.get("strandedness_match") else "UYUŞMUYOR") if L is LABELS["tr"] \
            else ("match" if seqqc.get("strandedness_match") else "MISMATCH")
        seq_line = (f'<p>{_esc(L["rrna_mean"])}: {_pct(seqqc.get("mean_rrna_fraction"))} · '
                    f'{_esc(L["strandedness"])}: {_esc(seqqc.get("inferred_strandedness"))} '
                    f'(≟ {_esc(seqqc.get("declared_strandedness"))} — {_esc(match)})</p>')
    figs = _quality_figures(qc, figures_dir, alignqc, L)
    rdist = _read_distribution_table(alignqc, L)
    mqc = _multiqc_note(multiqc, L)
    return (f'<section id="quality"><h2>{_esc(L["quality"])}</h2>'
            f'{_intro("quality", L)}{trim}{seq_line}{tbl}{rdist}{figs}{mqc}</section>')


def _fig_block(png: Path, title: str) -> str:
    return (f'<figure><img src="{embed_png(png)}" alt="{_esc(title)}"/>'
            f'<figcaption>{_esc(title)}</figcaption></figure>')


def _quality_figures(qc: dict, figures_dir: Path | None,
                     alignqc: dict | None, L: dict) -> str:
    """Kalite bölümü QC figürleri: per-base kompozisyon (F1), insert-size (F2),
    coverage (F4). Yalnız üretilmiş olanlar gömülür."""
    if figures_dir is None:
        return ""
    figures_dir = Path(figures_dir)
    blocks = []
    alignqc = alignqc or {}
    for name, title in (
        (qc.get("per_base_composition_figure"), L["cap_per_base"]),
        (alignqc.get("insert_size_figure"), L["cap_insert_size"]),
        (alignqc.get("coverage_figure"), L["cap_coverage"]),
    ):
        if name:
            png = figures_dir / name
            if png.exists():
                blocks.append(_fig_block(png, title))
    return "".join(blocks)


def _read_distribution_table(alignqc: dict | None, L: dict) -> str:
    """RSeQC read-distribution (F3): okumaların genomik özelliklere yüzde dağılımı."""
    alignqc = alignqc or {}
    rd = alignqc.get("read_distribution")
    if not rd:
        return ""
    rows = [[grp, f'{pct:.1f}%'] for grp, pct in rd.items()]
    return _table([L["rd_group"], L["rd_pct"]], rows, L["cap_read_dist"])


def _multiqc_note(multiqc: dict | None, L: dict) -> str:
    """MultiQC toplu görünüm (F5) — göreli link + kısa not."""
    multiqc = multiqc or {}
    rel = multiqc.get("report_relpath")
    if not rel:
        return ""
    return (f'<p class="note">{_esc(L["multiqc_note"])} '
            f'<a href="{_esc(rel)}">{_esc(rel)}</a></p>')


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
    tbl = _table([" ", " "], rows, L["cap_de"])
    expr = f'<p class="note">{_esc(L["expr_note"])}</p>'
    return f'<section id="de"><h2>{_esc(L["de"])}</h2>{_intro("de", L)}{summary}{tbl}{expr}</section>'


def collect_figure_errors(run_dir: Path) -> list[str]:
    """statistics/*.json içindeki `figure_errors` (best-effort figürlerin başarısızlıkları)
    tek bir insan-okunur listede toplar. Bu figürler TANISALDIR; verdict'i etkilemez ama
    üretilememeleri müşteri çıktısında GÖRÜNMELİ (yalnız log'da kalmasın)."""
    stats_dir = Path(run_dir) / "statistics"
    if not stats_dir.exists():
        return []
    out: list[str] = []
    for path in sorted(stats_dir.glob("*.json")):
        try:
            data = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        errors = data.get("figure_errors") if isinstance(data, dict) else None
        if isinstance(errors, dict):
            source = path.stem.replace("_statistics", "")
            for name, msg in errors.items():
                out.append(f"{source}/{name}: {msg}")
    return out


def section_figures(figures_manifest: dict, figures_dir: Path, L: dict, lang: str = "tr",
                    figure_errors: list[str] | None = None) -> str:
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
    note = ""
    if figure_errors:
        if lang == "en":
            lead = (f"Note: {len(figure_errors)} diagnostic figure(s) could not be generated "
                    "(this does not affect the results):")
        else:
            lead = (f"Not: {len(figure_errors)} tanısal figür üretilemedi "
                    "(sonucu etkilemez):")
        items = "".join(f"<li>{_esc(e)}</li>" for e in figure_errors)
        note = f'<div class="note"><p>{_esc(lead)}</p><ul>{items}</ul></div>'
    return (f'<section id="figures"><h2>{_esc(L["figures"])}</h2>{_intro("figures", L)}'
            f'{note}{"".join(blocks)}</section>')


def _deg_table(rows: list[dict], L: dict, cond_ctx: dict | None = None, caption: str = "") -> str:
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
    return _table(headers, body, caption)


def section_table(de_results: list, gene_map: dict, fdr: float, lfc: float, L: dict,
                  cond_ctx: dict | None = None) -> str:
    up = top_degs_by_direction(de_results, gene_map, fdr, lfc, "Up", n=25)
    down = top_degs_by_direction(de_results, gene_map, fdr, lfc, "Down", n=25)
    if not up and not down:
        return f'<section id="table"><h2>{_esc(L["table"])}</h2>{_intro("table", L)}<p>{_esc(L["no_degs"])}</p></section>'
    up_html = (f'<h3>{_esc(L["up_table"])}</h3>{_deg_table(up, L, cond_ctx, L["up_table"])}' if up
               else f'<h3>{_esc(L["up_table"])}</h3><p>{_esc(L["no_degs"])}</p>')
    down_html = (f'<h3>{_esc(L["down_table"])}</h3>{_deg_table(down, L, cond_ctx, L["down_table"])}' if down
                 else f'<h3>{_esc(L["down_table"])}</h3><p>{_esc(L["no_degs"])}</p>')
    return (f'<section id="table"><h2>{_esc(L["table"])}</h2>{_intro("table", L)}{up_html}{down_html}'
            f'<p class="note">{_esc(L["full_table_note"])}</p></section>')


_NS_ORDER = ["BP", "MF", "CC"]


def _enrichment_dir_table(rows: list[dict], L: dict, top_n: int = 10, caption: str = "") -> str:
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
    return _table(headers, body, caption)


def _enrichment_collection(up, down, manifest, figs_dir: Path, L: dict, lang: str,
                           up_label: str, down_label: str, cap_ids: tuple[str, str]) -> str:
    """Bir zenginleştirme koleksiyonunu (GO veya KEGG) artan/azalan tablo + figürlerle basar."""
    caps = FIGURE_CAPTIONS.get(lang, FIGURE_CAPTIONS["tr"])
    fig_by_id = {f["id"]: f for f in (manifest or {}).get("figures", [])}
    blocks = []
    for rows, lab, fid in ((up or [], up_label, cap_ids[0]),
                           (down or [], down_label, cap_ids[1])):
        blocks.append(f'<h4>{_esc(L[lab])}</h4>')
        blocks.append(_enrichment_dir_table(rows, L, caption=L[lab]))
        fig = fig_by_id.get(fid)
        if fig and (Path(figs_dir) / fig["png"]).exists():
            cap = caps.get(fid, "")
            blocks.append(
                f'<figure><img src="{embed_png(Path(figs_dir) / fig["png"])}" '
                f'alt="{_esc(fig.get("title"))}"/><figcaption>{_esc(cap)}</figcaption></figure>')
    return "".join(blocks)


def section_go(inputs: dict, L: dict, lang: str = "tr") -> str:
    """GO zenginleştirme (ORA) — ayrı bölüm."""
    up, down = inputs.get("enrichment_up"), inputs.get("enrichment_down")
    if up is None and down is None:
        return (f'<section id="go"><h2>{_esc(L["go_section"])}</h2>'
                f'<p class="note">{_esc(L["enrichment_not_run"])}</p></section>')
    block = _enrichment_collection(
        up, down, inputs.get("enrichment_manifest"), inputs.get("enrichment_dir", "."),
        L, lang, "enrichment_up", "enrichment_down", ("enrichment_up", "enrichment_down"))
    return (f'<section id="go"><h2>{_esc(L["go_section"])}</h2>'
            f'{_intro("enrichment", L)}{block}<p class="note">{L["enrichment_legend"]}</p></section>')


def section_kegg(inputs: dict, L: dict, lang: str = "tr") -> str:
    """KEGG yolak zenginleştirme (ORA) — ayrı bölüm."""
    up, down = inputs.get("kegg_up"), inputs.get("kegg_down")
    if up is None and down is None:
        return (f'<section id="kegg"><h2>{_esc(L["kegg_section"])}</h2>'
                f'<p class="note">{_esc(L["kegg_not_run"])}</p></section>')
    block = _enrichment_collection(
        up, down, inputs.get("kegg_manifest"), inputs.get("kegg_dir", "."),
        L, lang, "kegg_up", "kegg_down", ("kegg_up", "kegg_down"))
    return (f'<section id="kegg"><h2>{_esc(L["kegg_section"])}</h2>'
            f'{_intro("kegg", L)}{block}<p class="note">{L["kegg_legend"]}</p></section>')


def _gsea_table(rows: list[dict], L: dict, top_n: int = 10, caption: str = "") -> str:
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
    return _table(headers, body, caption)


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
        blocks.append(_gsea_table(rows, L, caption=f'GSEA — {L[heading]}'))
        fig = fig_by_id.get(fid)
        if fig and (figs_dir / fig["png"]).exists():
            cap = caps.get(fid, "")
            blocks.append(
                f'<figure><img src="{embed_png(figs_dir / fig["png"])}" '
                f'alt="{_esc(fig.get("title"))}"/><figcaption>{_esc(cap)}</figcaption></figure>')
    blocks.append(f'<p class="note">{L["gsea_legend"]}</p>')
    return (f'<section id="gsea"><h2>{_esc(L["gsea"])}</h2>'
            f'{_intro("gsea", L)}{"".join(blocks)}</section>')


def _reduced_table(rows: list[dict], L: dict, caption: str = "") -> str:
    """Temsilci GO terimleri: term, namespace, padj, temsil ettiği terim sayısı."""
    if not rows:
        return f'<p>{_esc(L["no_enrichment"])}</p>'
    n_terms = sum(r.get("n_collapsed", 1) for r in rows)
    summary = f'<p class="note">{_esc(L["sem_summary"].format(n=n_terms, m=len(rows)))}</p>'
    headers = [L["go_term"], L["namespace"], L["padj"], L["n_collapsed"]]
    body = []
    for r in sorted(rows, key=lambda x: -(x.get("n_collapsed") or 1)):
        body.append([
            r.get("term") or r.get("go_id"), r.get("namespace"),
            f'{r["padj"]:.2e}' if r.get("padj") is not None else "—",
            r.get("n_collapsed"),
        ])
    return summary + _table(headers, body, caption)


def section_semantic(inputs: dict, L: dict, lang: str = "tr") -> str:
    sources = [("sem_ora_up", inputs.get("reduced_ora_up")),
               ("sem_ora_down", inputs.get("reduced_ora_down")),
               ("sem_gsea_go", inputs.get("reduced_gsea_go"))]
    if all(rows is None for _, rows in sources):
        return (f'<section id="semantic"><h2>{_esc(L["semantic"])}</h2>'
                f'<p class="note">{_esc(L["semantic_not_run"])}</p></section>')
    blocks = []
    for label, rows in sources:
        if rows is None:
            continue
        blocks.append(f'<h3>{_esc(L[label])}</h3>')
        blocks.append(_reduced_table(rows, L, caption=L[label]))
    blocks.append(_embed_first_figure(inputs.get("semantic_manifest"), inputs.get("semantic_dir", ".")))
    blocks.append(f'<p class="note">{L["semantic_legend"]}</p>')
    return (f'<section id="semantic"><h2>{_esc(L["semantic"])}</h2>'
            f'{_intro("semantic", L)}{"".join(blocks)}</section>')


def _amr_table(rows: list[dict], L: dict, cap: int = 40, caption: str = "") -> str:
    """AMR (CARD↔AMRFinderPlus yan yana) veya virülans (tek araç) tablosu. up/down önce; uzunsa kırp."""
    if not rows:
        return f'<p>{_esc(L["no_degs"])}</p>'
    de = [r for r in rows if r.get("de_status") in ("up", "down")]
    shown = rows[:cap]
    side_by_side = "card" in rows[0] or "amrfinder" in rows[0]   # AMR tablosu iki araç sütunlu
    if side_by_side:
        headers = [L["gene"], L["card_col"], L["amrfinder_col"], L["identity"],
                   L["log2fc"], L["padj"], L["de_status"]]
        body = [[
            r.get("gene"), r.get("card") or "—", r.get("amrfinder") or "—",
            f'{r["pct_identity"]:.1f}' if r.get("pct_identity") is not None else "—",
            f'{r["log2fc"]:.2f}' if r.get("log2fc") is not None else "—",
            f'{r["padj"]:.2e}' if r.get("padj") is not None else "—",
            r.get("de_status"),
        ] for r in shown]
    else:
        headers = [L["gene"], L["amr_label"], L["identity"], L["log2fc"], L["padj"], L["de_status"]]
        body = [[
            r.get("gene"), r.get("label"),
            f'{r["pct_identity"]:.1f}' if r.get("pct_identity") is not None else "—",
            f'{r["log2fc"]:.2f}' if r.get("log2fc") is not None else "—",
            f'{r["padj"]:.2e}' if r.get("padj") is not None else "—",
            r.get("de_status"),
        ] for r in shown]
    tbl = _table(headers, body, caption)
    note = ""
    if len(rows) > cap:
        note = f'<p class="note">{_esc(L["amr_more"].format(n=len(rows) - cap))}</p>'
    summary = f'<p class="note">{len(rows)} gen · {len(de)} DE</p>'
    return summary + tbl + note


def section_amr(inputs: dict, L: dict, lang: str = "tr") -> str:
    amr, vir = inputs.get("amr_genes"), inputs.get("virulence_genes")
    if amr is None and vir is None:
        return (f'<section id="amr"><h2>{_esc(L["amr"])}</h2>'
                f'<p class="note">{_esc(L["amr_not_run"])}</p></section>')
    blocks = []
    if amr is not None:
        blocks.append(f'<h3>{_esc(L["amr_genes"])}</h3>{_amr_table(amr, L, caption=L["amr_genes"])}')
    if vir is not None:
        blocks.append(f'<h3>{_esc(L["vir_genes"])}</h3>{_amr_table(vir, L, caption=L["vir_genes"])}')
    blocks.append(f'<p class="note">{L["amr_legend"]}</p>')
    return (f'<section id="amr"><h2>{_esc(L["amr"])}</h2>'
            f'{_intro("amr", L)}{"".join(blocks)}</section>')


def _embed_first_figure(manifest: dict | None, figs_dir) -> str:
    """Manifest'teki figürleri (varsa) gömer. Yoksa/dosya yoksa boş döner (tolerant)."""
    if not manifest:
        return ""
    figs_dir = Path(figs_dir)
    blocks = []
    for fig in manifest.get("figures", []):
        png = figs_dir / fig.get("png", "")
        if png.exists():
            blocks.append(f'<figure><img src="{embed_png(png)}" alt="{_esc(fig.get("title"))}"/>'
                          f'<figcaption><strong>{_esc(fig.get("title"))}</strong></figcaption></figure>')
    return "".join(blocks)


def section_operon(inputs: dict, L: dict, lang: str = "tr", cap: int = 30) -> str:
    operons = inputs.get("operons")
    if operons is None:
        return (f'<section id="operon"><h2>{_esc(L["operon"])}</h2>'
                f'<p class="note">{_esc(L["operon_not_run"])}</p></section>')
    stats = inputs.get("operon_stats") or {}
    gap = stats.get("max_gap", 50)
    summary = f'<p class="note">{_esc(L["operon_summary"].format(n=stats.get("n_operons", len(operons)), m=stats.get("n_multi_gene", 0), k=stats.get("n_coordinated", 0)))}</p>'
    coord = [o for o in operons if o.get("coordinated")][:cap]
    if not coord:
        body_html = f'<p>{_esc(L["no_degs"])}</p>'
    else:
        headers = [L["operon_genes"], L["operon_size"], L["n_sig"], L["operon_dir"], L["log2fc"]]
        rows = []
        for o in coord:
            direction = L["up"] if o.get("n_up", 0) >= o.get("n_down", 0) else L["down"]
            rows.append([
                ", ".join(o.get("genes", "").split(";")), o.get("size"), o.get("n_deg"),
                direction, f'{o["mean_log2fc"]:.2f}' if o.get("mean_log2fc") is not None else "—",
            ])
        body_html = _table(headers, rows, L["cap_operon"])
    fig_html = _embed_first_figure(inputs.get("operon_manifest"), inputs.get("operon_dir", "."))
    legend = f'<p class="note">{L["operon_legend"].format(gap=gap)}</p>'
    return (f'<section id="operon"><h2>{_esc(L["operon"])}</h2>'
            f'{_intro("operon", L)}{summary}{body_html}{fig_html}{legend}</section>')


def section_ppi(inputs: dict, L: dict, lang: str = "tr", cap: int = 20) -> str:
    comms = inputs.get("communities")
    if comms is None:
        return (f'<section id="ppi"><h2>{_esc(L["ppi"])}</h2>'
                f'<p class="note">{_esc(L["ppi_not_run"])}</p></section>')
    stats = inputs.get("ppi_stats") or {}
    summary = f'<p class="note">{_esc(L["ppi_summary"].format(net=stats.get("n_deg_in_network", 0), deg=stats.get("n_deg", 0), edges=stats.get("n_edges", 0), mods=stats.get("n_communities", len(comms))))}</p>'
    if not comms:
        body_html = f'<p>{_esc(L["no_degs"])}</p>'
    else:
        headers = [L["ppi_module"], L["ppi_size"], L["operon_dir"], L["operon_genes"]]
        rows = []
        for c in comms[:cap]:
            direction = L["up"] if c.get("n_up", 0) >= c.get("n_down", 0) else L["down"]
            rows.append([c.get("community_id"), c.get("size"), direction,
                         ", ".join(c.get("genes", "").split(";"))])
        body_html = _table(headers, rows, L["cap_ppi"])
    fig_html = _embed_first_figure(inputs.get("ppi_manifest"), inputs.get("ppi_dir", "."))
    score = stats.get("min_score", 700)
    legend = f'<p class="note">{L["ppi_legend"].format(score=score)}</p>'
    return (f'<section id="ppi"><h2>{_esc(L["ppi"])}</h2>'
            f'{_intro("ppi", L)}{summary}{fig_html}{body_html}{legend}</section>')


# Kullanılan yazılım sürümleri (envs/*.yml ile eşleşir). (araç, sürüm, amaç_tr, amaç_en, koşul).
_SOFTWARE: list[tuple] = [
    ("Python", "3.11", "Orkestrasyon", "Orchestration", None),
    # Kısa-okuma araç zinciri (cond "short")
    ("FastQC", "0.12.1", "Ham okuma QC", "Raw read QC", "short"),
    ("fastp", "1.3.6", "Kırpma", "Trimming", "short"),
    ("Bowtie2", "2.5.5", "Hizalama", "Alignment", "short"),
    # Uzun-okuma araç zinciri (cond "long")
    ("NanoPlot", "1.47.1", "Uzun-okuma QC", "Long-read QC", "long"),
    ("Pychopper", "2.7.10", "Tam-boy cDNA yönlendirme/kırpma", "Full-length cDNA orient/trim", "long"),
    ("chopper", "0.13.0", "Uzun-okuma uzunluk/kalite filtresi", "Long-read length/quality filter", "long"),
    ("minimap2", "2.31", "Uzun-okuma hizalama (-ax map-ont/map-hifi)", "Long-read alignment (-ax map-ont/map-hifi)", "long"),
    ("Subread/featureCounts", "2.1.1", "Sayım (TPM/FPKM; uzun-okumada -L)", "Quantification (TPM/FPKM; -L for long-read)", None),
    ("DESeq2 (R)", "1.50.2", "Diferansiyel ekspresyon", "Differential expression", None),
    ("ggplot2 (R)", "4.0", "Figürler", "Figures", None),
    ("fgsea (R)", "1.36.2", "GSEA", "GSEA", "gsea"),
    ("abricate", "1.4.0", "AMR/virülans (CARD/VFDB)", "AMR/virulence (CARD/VFDB)", "amr"),
    ("AMRFinderPlus", "4.2.7", "AMR (küratörlü)", "AMR (curated)", "amr"),
    ("SortMeRNA", "7.0.0", "rRNA%", "rRNA%", "seqqc"),
    ("RSeQC", "5.0.5", "Strandedness", "Strandedness", "seqqc"),
    ("networkx", "3.6.1", "PPI community", "PPI community", "ppi"),
]


def section_software(config, L: dict, inputs: dict, flags: dict) -> str:
    lang = "en" if L is LABELS["en"] else "tr"
    # F2.2: runtime'da yakalanan gerçek sürümler (varsa) curated fallback'i ezer.
    runtime = inputs.get("software_versions") or {}
    sw_rows = [[t, runtime.get(t, v), (ptr if lang == "tr" else pen)]
               for t, v, ptr, pen, cond in _SOFTWARE
               if cond is None or flags.get(cond)]
    sw_tbl = _table([L["sw_tool"], L["sw_version"], L["sw_purpose"]], sw_rows, L["cap_software"])
    # Veritabanları — yalnız koşan analizler
    e = config.enrichment
    db_rows = []
    if flags.get("enrichment"):
        db_rows.append(["Gene Ontology (go-basic)", "obo", "GO anotasyonu" if lang == "tr" else "GO annotation"])
        db_rows.append(["EBI-GOA", "GAF", "GO anotasyonu (tamamlayıcı)" if lang == "tr" else "GO annotation (supplement)"])
    if flags.get("kegg"):
        db_rows.append(["KEGG", f"REST ({e.kegg_organism or '—'})", "Yolak" if lang == "tr" else "Pathways"])
    if flags.get("gsea") or flags.get("enrichment"):
        pass
    if flags.get("amr"):
        db_rows.append(["CARD", "abricate 1.4.0", "Direnç" if lang == "tr" else "Resistance"])
        db_rows.append(["VFDB", "abricate 1.4.0", "Virülans" if lang == "tr" else "Virulence"])
        db_rows.append(["AMRFinderPlus DB", "2026-05-15", "Direnç (küratörlü)" if lang == "tr" else "Resistance (curated)"])
    if flags.get("ppi"):
        db_rows.append(["STRING", f"v12.0 ({config.ppi.taxid or '—'})", "PPI ağı" if lang == "tr" else "PPI network"])
    if inputs.get("seqqc"):
        db_rows.append(["rRNA referansı", "referans genom" if lang == "tr" else "reference genome",
                        "rRNA%" if lang == "tr" else "rRNA%"])
    db_html = _table([L["db_name"], L["db_version"], L["db_purpose"]], db_rows, L["cap_database"]) if db_rows else ""
    return (f'<section id="software"><h2>{_esc(L["software"])}</h2>'
            f'{_intro("software", L)}{sw_tbl}{db_html}</section>')


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

# Uzun-okuma (ONT/PacBio) yöntem anlatısı — kısa-okuma araç zinciri yerine NanoPlot/
# Pychopper+chopper/minimap2/featureCounts -L. DESeq2 ve sonrası kuyruğu kısa ile aynı.
# {aggr} bilinçli olarak kullanılmaz (uzun-okumada agresif-kırpma kavramı yok); .format
# fazladan kwarg'ı yok sayar.
_METHODS_TEXT_LONG: dict[str, str] = {
    "tr": (
        "Ham uzun okumaların kalitesi NanoPlot ile değerlendirildi (okuma uzunluğu dağılımı, N50, "
        "okuma kalitesi). cDNA kütüphaneleri için tam-boy transkriptler Pychopper ile yönlendirilip "
        "kesildi; okumalar chopper ile uzunluk/kalite açısından filtrelendi (asgari uzunluk {min_len} nt). "
        "Filtrelenen uzun okumalar referans genoma minimap2 ile hizalandı (uzun-okuma ön ayarı "
        "-ax map-ont/map-hifi). Hizalanan okumalar featureCounts uzun-okuma modunda (-L) özniteliklere "
        "atanarak gen×örnek sayım matrisi oluşturuldu (öznitelik tipi {feature_type}, kimlik özniteliği "
        "{attribute}). Diferansiyel ekspresyon DESeq2 ile hesaplandı: kütüphane büyüklüğü medyan-oran "
        "yöntemiyle normalize edildi (boyut faktörleri); gen-bazlı dispersiyonlar kestirilip ortalama-"
        "dispersiyon eğilimine doğru empirical-Bayes ile büzüldü; koşul etkisi negatif binom genelleştirilmiş "
        "doğrusal modelle test edildi (Wald testi) ve p-değerleri Benjamini–Hochberg (FDR) ile düzeltildi. "
        "Tasarım formülü {design}. Bir gen, düzeltilmiş p (padj) < {fdr} ve |log2 kat değişimi| ≥ {lfc} ise "
        "anlamlı diferansiyel eksprese sayıldı. Uzun-okuma kalite eşikleri bilinçli olarak permissive'dir "
        "(prokaryote_long profili, rapora damgalanır). Tüm görseller ggplot2 ile üretildi."
    ),
    "en": (
        "Raw long-read quality was assessed with NanoPlot (read-length distribution, N50, read quality). "
        "For cDNA libraries, full-length transcripts were oriented and trimmed with Pychopper, and reads were "
        "length/quality filtered with chopper (minimum length {min_len} nt). Filtered long reads were aligned "
        "to the reference genome with minimap2 (long-read preset -ax map-ont/map-hifi). Aligned reads were "
        "assigned to features with featureCounts in long-read mode (-L) to build a gene×sample count matrix "
        "(feature type {feature_type}, identifier attribute {attribute}). Differential expression was computed "
        "with DESeq2: library sizes were normalized by the median-of-ratios method (size factors); gene-wise "
        "dispersions were estimated and shrunk toward the mean–dispersion trend by empirical Bayes; the "
        "condition effect was tested with a negative-binomial generalized linear model (Wald test) and p-values "
        "were adjusted by Benjamini–Hochberg (FDR). Design formula {design}. A gene was called significantly "
        "differentially expressed when adjusted p (padj) < {fdr} and |log2 fold change| ≥ {lfc}. Long-read "
        "quality thresholds are deliberately permissive (prokaryote_long profile, stamped in the report). All "
        "figures were produced with ggplot2."
    ),
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


_SEMANTIC_METHODS: dict[str, str] = {
    "tr": (
        "Zenginleşen GO terimlerinin fazlalığı, terim çiftleri arasındaki Lin semantik benzerliğiyle "
        "(bilgi içeriği arka plan anotasyonundan türetildi) azaltıldı: benzerliği {thr} eşiğini aşan "
        "terimler, her namespace içinde en iyi düzeltilmiş p-değerli temsilci altında toplandı (REVIGO yaklaşımı)."
    ),
    "en": (
        "Redundancy among enriched GO terms was reduced using Lin semantic similarity between term pairs "
        "(information content derived from the background annotation): terms exceeding a similarity of "
        "{thr} were collapsed, within each namespace, under the representative with the best adjusted "
        "p-value (the REVIGO approach)."
    ),
}


_AMR_METHODS: dict[str, str] = {
    "tr": (
        "Suşun referans genomu abricate ({amr_db} ve {vir_db} veritabanları) ile tarandı; en az %{min_id} "
        "kimlik ve %{min_cov} kapsama sağlayan gen isabetleri, genom koordinatlarıyla anotasyon "
        "genlerine (locus_tag) eşlendi ve diferansiyel ekspresyon durumları eklendi. Veritabanları "
        "abricate ile paket halinde gelir."
    ),
    "en": (
        "The strain reference genome was scanned with abricate (databases {amr_db} and {vir_db}); gene "
        "hits passing at least {min_id}% identity and {min_cov}% coverage were mapped to annotation genes "
        "(locus tags) by genome coordinates and annotated with their differential-expression status. "
        "The databases are bundled with abricate."
    ),
}


_OPERON_METHODS: dict[str, str] = {
    "tr": (
        "Operon yapısı, referans anotasyonundaki gen koordinatlarından tahmin edildi: aynı yönde bitişik "
        "ve intergenik boşluğu en fazla {gap} bç olan genler tek bir operon (transkripsiyon birimi) sayıldı "
        "(Moreno-Hagelsieb ve Collado-Vides 2002). Her operon için genlerinin diferansiyel ekspresyon yönü "
        "değerlendirilerek koordineli değişen (birlikte-düzenlenen) operonlar belirlendi."
    ),
    "en": (
        "Operon structure was predicted from gene coordinates in the reference annotation: genes on the same "
        "strand with an intergenic gap of at most {gap} bp were taken as one operon (transcription unit) "
        "(Moreno-Hagelsieb and Collado-Vides 2002). For each operon the differential-expression direction of "
        "its genes was assessed to identify coordinately (co-regulated) changing operons."
    ),
}


_PPI_METHODS: dict[str, str] = {
    "tr": (
        "Diferansiyel eksprese genler arasındaki protein-protein etkileşimleri STRING veritabanından "
        "(combined score ≥ {score}) alınarak bir ağ kuruldu ve ağ, modülarite temelli Louvain yöntemiyle "
        "(networkx) topluluklara (modüllere) ayrıldı. STRING etkileşimleri kanıt-skorlu tahminlerdir."
    ),
    "en": (
        "Protein–protein interactions among the differentially expressed genes were taken from the STRING "
        "database (combined score ≥ {score}) to build a network, which was partitioned into communities "
        "(modules) by the modularity-based Louvain method (networkx). STRING interactions are "
        "evidence-scored predictions."
    ),
}


def section_methods(config, L: dict, enrichment_ran: bool = False, kegg_ran: bool = False,
                    gsea_ran: bool = False, semantic_ran: bool = False, amr_ran: bool = False,
                    operon_ran: bool = False, ppi_ran: bool = False,
                    read_type: str = "short") -> str:
    lang = "en" if L is LABELS["en"] else "tr"
    t = config.trimming
    q = config.quantification
    d = config.de
    # read_type'a göre araç-zinciri anlatısı (uzun: NanoPlot/Pychopper/minimap2/-L).
    methods_src = _METHODS_TEXT_LONG if read_type == "long" else _METHODS_TEXT
    text = methods_src[lang].format(
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
    if semantic_ran:
        s = _SEMANTIC_METHODS[lang].format(thr=config.enrichment.revigo_similarity)
        paras += f'<p>{_esc(s)}</p>'
    if amr_ran:
        am = _AMR_METHODS[lang].format(amr_db=config.amr.amr_db, vir_db=config.amr.virulence_db,
                                       min_id=config.amr.min_identity, min_cov=config.amr.min_coverage)
        if config.amr.amrfinder_organism:
            am += (" AMR genleri ayrıca NCBI AMRFinderPlus ile bağımsız olarak tespit edilip CARD "
                   "sonuçlarıyla yan yana raporlandı (konkordans)." if lang == "tr" else
                   " AMR genes were additionally detected independently with NCBI AMRFinderPlus and "
                   "reported side by side with the CARD results (concordance).")
        paras += f'<p>{_esc(am)}</p>'
    if operon_ran:
        op = _OPERON_METHODS[lang].format(gap=config.operon.max_gap)
        paras += f'<p>{_esc(op)}</p>'
    if ppi_ran:
        pp = _PPI_METHODS[lang].format(score=config.ppi.min_score)
        paras += f'<p>{_esc(pp)}</p>'
    return f'<section id="methods"><h2>{_esc(L["methods"])}</h2>{_intro("methods", L)}{paras}</section>'


# (atıf, DOI/URL). DOI'ler doğrulandı (doi.org). FastQC bir araçtır; DOI'si yoktur → proje URL'si.
# Kısa-okuma araç zincirine özgü atıflar (yalnız read_type=short raporunda).
_REFERENCES_SHORT: list[tuple[str, str]] = [
    ("Andrews S. FastQC: a quality control tool for high throughput sequence data. 2010.",
     "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/"),
    ("Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. "
     "Bioinformatics. 2018;34(17):i884–i890.", "https://doi.org/10.1093/bioinformatics/bty560"),
    ("Langmead B, Salzberg SL. Fast gapped-read alignment with Bowtie 2. "
     "Nat Methods. 2012;9(4):357–359.", "https://doi.org/10.1038/nmeth.1923"),
    ("Williams CR, Baccarella A, Parrish JZ, Kim CC. Trimming of sequence reads alters RNA-Seq gene "
     "expression estimates. BMC Bioinformatics. 2016;17:103.",
     "https://doi.org/10.1186/s12859-016-0956-2"),
]
# Uzun-okuma (ONT/PacBio) araç zincirine özgü atıflar (yalnız read_type=long raporunda).
_REFERENCES_LONG: list[tuple[str, str]] = [
    ("De Coster W, Rademakers R. NanoPack2: population-scale evaluation of long-read sequencing data. "
     "Bioinformatics. 2023;39(5):btad311.", "https://doi.org/10.1093/bioinformatics/btad311"),
    ("Li H. Minimap2: pairwise alignment for nucleotide sequences. "
     "Bioinformatics. 2018;34(18):3094–3100.", "https://doi.org/10.1093/bioinformatics/bty191"),
    ("Oxford Nanopore Technologies. Pychopper: identify, orient and trim full-length cDNA reads. 2023.",
     "https://github.com/epi2me-labs/pychopper"),
]
# Okuma-tipinden bağımsız (her raporda) atıflar.
_REFERENCES: list[tuple[str, str]] = [
    ("Liao Y, Smyth GK, Shi W. featureCounts: an efficient general purpose program for assigning "
     "sequence reads to genomic features. Bioinformatics. 2014;30(7):923–930.",
     "https://doi.org/10.1093/bioinformatics/btt656"),
    ("Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data "
     "with DESeq2. Genome Biol. 2014;15:550.", "https://doi.org/10.1186/s13059-014-0550-8"),
    ("Wickham H. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York; 2016.",
     "https://doi.org/10.1007/978-3-319-24277-4"),
    # Genel RNA-seq analiz/raporlama metodolojisi (kullanıcı sağladı).
    ("Dawadi P, Pokharel B, Shrestha A, et al. From bench to bytes: a practical guide to RNA sequencing "
     "data analysis. Front Genet. 2025;16:1697922.", "https://doi.org/10.3389/fgene.2025.1697922"),
    ("Deshpande D, Chhugani K, Chang Y, et al. RNA-seq data science: from raw data to effective "
     "interpretation. Front Genet. 2023;14:997383.", "https://doi.org/10.3389/fgene.2023.997383"),
    ("Pola-Sánchez E, Hernández-Martínez KM, Pérez-Estrada R, et al. RNA-Seq data analysis: a practical "
     "guide for model and non-model organisms. Curr Protoc. 2024;4(4):e1054.",
     "https://doi.org/10.1002/cpz1.1054"),
    ("Claussen H. Bulk RNAseq Methodology and Analysis Guide. Emory Integrated Computational Core "
     "(EICC); 2023.", "https://www.cores.emory.edu/eicc/"),
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

# Semantic reduction kaynakları — yalnız m12 çalıştıysa (DOI doğrulandı).
_SEMANTIC_REFERENCES: list[tuple[str, str]] = [
    ("Lin D. An information-theoretic definition of similarity. Proc 15th Int Conf Machine Learning "
     "(ICML). 1998:296–304.", "https://dl.acm.org/doi/10.5555/645527.657297"),
    ("Supek F, Bošnjak M, Škunca N, Šmuc T. REVIGO summarizes and visualizes long lists of gene "
     "ontology terms. PLoS One. 2011;6(7):e21800.", "https://doi.org/10.1371/journal.pone.0021800"),
]

# AMR/virülans kaynakları — yalnız m13 çalıştıysa (DOI doğrulandı).
_AMR_REFERENCES: list[tuple[str, str]] = [
    ("Seemann T. Abricate: mass screening of contigs for antimicrobial and virulence genes. "
     "GitHub.", "https://github.com/tseemann/abricate"),
    ("Alcock BP, Raphenya AR, Lau TTY, et al. CARD 2020: antibiotic resistome surveillance with the "
     "comprehensive antibiotic resistance database. Nucleic Acids Res. 2020;48(D1):D517–D525.",
     "https://doi.org/10.1093/nar/gkz935"),
    ("Liu B, Zheng D, Jin Q, Chen L, Yang J. VFDB 2019: a comparative pathogenomic platform with an "
     "interactive web interface. Nucleic Acids Res. 2019;47(D1):D687–D692.",
     "https://doi.org/10.1093/nar/gky1080"),
    ("Feldgarden M, Brover V, Gonzalez-Escalona N, et al. AMRFinderPlus and the Reference Gene Catalog "
     "facilitate examination of the genomic links among antimicrobial resistance, stress response, and "
     "virulence. Sci Rep. 2021;11:12728.", "https://doi.org/10.1038/s41598-021-91456-0"),
]

# Operon kaynağı — yalnız m14 çalıştıysa (DOI doğrulandı).
_OPERON_REFERENCES: list[tuple[str, str]] = [
    ("Moreno-Hagelsieb G, Collado-Vides J. A powerful non-homology method for the prediction of operons "
     "in prokaryotes. Bioinformatics. 2002;18(Suppl 1):S329–S336.",
     "https://doi.org/10.1093/bioinformatics/18.suppl_1.S329"),
]

# PPI/community kaynakları — yalnız m15 çalıştıysa (DOI doğrulandı).
_PPI_REFERENCES: list[tuple[str, str]] = [
    ("Szklarczyk D, Kirsch R, Koutrouli M, et al. The STRING database in 2023: protein–protein "
     "association networks and functional enrichment analyses. Nucleic Acids Res. 2023;51(D1):D638–D646.",
     "https://doi.org/10.1093/nar/gkac1000"),
    ("Blondel VD, Guillaume JL, Lambiotte R, Lefebvre E. Fast unfolding of communities in large networks. "
     "J Stat Mech. 2008;2008(10):P10008.", "https://doi.org/10.1088/1742-5468/2008/10/P10008"),
]


def _ref_link_label(url: str) -> str:
    return url.split("doi.org/", 1)[1] if "doi.org/" in url else url


def section_references(L: dict, enrichment_ran: bool = False, kegg_ran: bool = False,
                       gsea_ran: bool = False, semantic_ran: bool = False,
                       amr_ran: bool = False, operon_ran: bool = False,
                       ppi_ran: bool = False, read_type: str = "short") -> str:
    # Okuma-tipine özgü araç atıfları önce (kullanılmayan aracı atıflamaz — dürüstlük).
    platform_refs = _REFERENCES_LONG if read_type == "long" else _REFERENCES_SHORT
    refs = (platform_refs + _REFERENCES + (_ENRICHMENT_REFERENCES if enrichment_ran else [])
            + (_KEGG_REFERENCES if kegg_ran else [])
            + (_GSEA_REFERENCES if gsea_ran else [])
            + (_SEMANTIC_REFERENCES if semantic_ran else [])
            + (_AMR_REFERENCES if amr_ran else [])
            + (_OPERON_REFERENCES if operon_ran else [])
            + (_PPI_REFERENCES if ppi_ran else []))
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
th{background:#f6f6f6}
caption{caption-side:top;text-align:left;font-weight:600;font-size:.9rem;color:#333;margin-bottom:.25rem}
.banner{padding:.8rem 1rem;border-radius:6px;margin:.6rem 0;font-size:1.05rem}
.verdict-trustworthy{background:#e6f4ea;border:1px solid #34a853}
.verdict-suspect{background:#fef7e0;border:1px solid #f9ab00}
.verdict-invalid{background:#fce8e6;border:1px solid #d93025}
.verdict-unknown{background:#f1f3f4;border:1px solid #9aa0a6}
figure{margin:1rem 0;text-align:center} img{max-width:100%;height:auto}
figcaption{color:#555;font-size:.9rem;margin-top:.3rem}
.note{color:#666;font-size:.85rem} .summary{font-size:1.05rem;font-weight:600}
.intro{color:#444;font-size:.95rem;margin:.2rem 0 .6rem}
section{padding-top:.5rem}
/* Her analiz bir sayfada: yazdırma/PDF'te her bölüm yeni sayfada başlar. */
@media print{
  body{max-width:none}
  section{break-before:page;page-break-before:always}
  section#confidence{break-before:auto;page-break-before:avoid}
  h2{page-break-after:avoid} figure{page-break-inside:avoid} img{max-width:100%}
}
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
    semantic_ran = any(inputs.get(k) is not None for k in
                       ("reduced_ora_up", "reduced_ora_down", "reduced_gsea_go"))
    amr_ran = (inputs.get("amr_genes") is not None
               or inputs.get("virulence_genes") is not None)
    operon_ran = inputs.get("operons") is not None
    ppi_ran = inputs.get("communities") is not None
    # read_type: raw_statistics otoriter; alignment stats yedek; varsayılan short.
    read_type = raw.get("read_type") or (inputs.get("alignment") or {}).get("read_type") or "short"
    generated = datetime.now().isoformat(timespec="seconds")
    run_part = f'{_esc(run_id)} · ' if run_id else ""
    header = (f'<h1>RNAForge — {_esc(raw.get("organism"))}</h1>'
              f'<p class="note">{run_part}{_esc(generated)} · v{_esc(version)}</p>')
    # Bölüm sırası kullanıcının standart RNA-seq listesine göre: QC/işleme → DESeq2 → figürler
    # (PCA…heatmap) → Top DEG → GO → KEGG → GSEA → REVIGO → (ek: AMR/operon/PPI) → yöntem/kaynak.
    body = "".join([
        header,
        section_confidence(inputs["confidence"], L),
        section_dataset(raw, L),
        section_quality(inputs["alignment"], inputs["count"], trimming_cfg, L, inputs.get("seqqc"),
                        inputs.get("qc"), inputs.get("figures_dir"),
                        inputs.get("alignqc"), inputs.get("multiqc")),
        section_de(inputs["de"], L),
        section_figures(inputs["figures"], inputs["figures_dir"], L, lang,
                        inputs.get("figure_errors")),
        section_table(inputs["de_results"], inputs["gene_map"],
                      config.de.fdr_threshold, config.de.log2fc_threshold, L, cond_ctx),
        section_go(inputs, L, lang),
        section_kegg(inputs, L, lang),
        section_gsea(inputs, L, lang),
        section_semantic(inputs, L, lang),
        section_amr(inputs, L, lang),
        section_operon(inputs, L, lang),
        section_ppi(inputs, L, lang),
        section_software(config, L, inputs, {
            "enrichment": enrichment_ran, "kegg": kegg_ran, "gsea": gsea_ran,
            "amr": amr_ran, "operon": operon_ran, "ppi": ppi_ran,
            "seqqc": inputs.get("seqqc") is not None,
            "short": read_type == "short", "long": read_type == "long"}),
        section_methods(config, L, enrichment_ran, kegg_ran, gsea_ran, semantic_ran, amr_ran, operon_ran, ppi_ran,
                        read_type=read_type),
        section_references(L, enrichment_ran, kegg_ran, gsea_ran, semantic_ran, amr_ran, operon_ran, ppi_ran,
                           read_type=read_type),
    ])
    # Figürleri belge sırasına göre numaralandır: "Şekil N." / "Figure N."
    fig_word = "Şekil" if lang == "tr" else "Figure"
    _fig = {"n": 0}

    def _numbered(_m):
        _fig["n"] += 1
        return f'<figcaption>{fig_word} {_fig["n"]}. '
    body = re.sub(r'<figcaption>', _numbered, body)
    # Tabloları belge sırasına göre numaralandır: "Tablo N. <ad>" / "Table N. <name>"
    tbl_word = "Tablo" if lang == "tr" else "Table"
    _tbl = {"n": 0}

    def _numbered_tbl(_m):
        _tbl["n"] += 1
        return f'<caption>{tbl_word} {_tbl["n"]}. '
    body = re.sub(r'<caption>', _numbered_tbl, body)
    return (f'<!doctype html><html lang="{_esc(lang)}"><head><meta charset="utf-8">'
            f'<meta name="viewport" content="width=device-width, initial-scale=1">'
            f'<title>RNAForge report</title><style>{_CSS}</style></head>'
            f'<body>{body}</body></html>')
