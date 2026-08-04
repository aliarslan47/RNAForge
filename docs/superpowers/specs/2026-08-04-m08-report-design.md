# m08 — HTML Rapor · Tasarım Spec'i

**Tarih:** 2026-08-04 · **Dal:** `feat/m08-report` · **Referans:** `PLAN.md` v1.3 §9, §11
**Ön koşul:** m07 (figures) `main`'de (merge `c46baf5`). Zincir: m01→m03→m04→m05→m06→m07→**m08**.

## 1. Amaç

m06/m07 çıktılarını tek, **self-contained**, çift dilli (`tr`|`en`) bir `report.html`'de birleştirmek —
her koşuda OTOMATİK. Rapor koşunun GERÇEK verisinden üretilir; **uydurma biyolojik yorum yoktur**
(Discussion/Conclusion otomatik yazılmaz — bkz. `feedback_uydurma_kavram`). MVP: yalnız HTML (PDF Faz 2+).

## 2. İlkeler / Kısıtlar

- **Yeni veri-kapısı YOK.** Rapor görsel/özet katmanıdır; biyolojiyi geçersiz kılmaz. Verdict
  güvence kartından (m06 gate'leri) **değişmeden** taşınır (m07 gibi).
- **Yüksek sesle hata:** zorunlu bir girdi (aşağıdaki sözleşmeler) eksik/bozuksa → `FileNotFoundError`
  (veya net `ValueError`), exit 1. Sessiz yarım rapor YASAK (`feedback_gurultulu_hata`).
- **Kenar durumları çökmez:** 0 anlamlı DEG → tablo "anlamlı gen yok" notu; <50 DEG → mevcut kadarı.
  (m07 heatmap dersinin tekrarı: geçerli biyoloji raporu düşürmez.)
- **Self-contained:** figür PNG'leri base64 olarak HTML'e gömülür → tek dosya taşınabilir (~2–3 MB).
- **Dil:** kod/log EN; rapor metni config `report.language` (`tr`|`en`). Veri değerleri dil-nötr.
- **Bağımlılık yok:** saf Python + stdlib (`base64`, `html`, `json`, `string`). jinja2/R YOK.
- `python -m rnaforge.cli` ÇALIŞMAZ; entry point `rnaforge`. Test:
  `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest -q`.

## 3. Mimari (m06/m07 desenini izler)

- **`rnaforge/report_html.py`** — saf yardımcılar (R runner YOK):
  - `load_report_inputs(run_dir) -> dict` — tüm istatistik JSON'ları + güvence kartı + figür manifesti +
    `deseq2_results.tsv` + `gene_map.tsv`'yi okur; zorunlu eksikse yüksek sesle hata.
  - `top_degs(de_results, gene_map, n=50) -> list[dict]` — padj'ye göre sıralı ilk n anlamlı gen;
    `{gene, log2fc, padj, base_mean, direction}`; gen adı locus_tag→gene (yoksa locus_tag).
  - `embed_png(path) -> str` — PNG → `data:image/png;base64,…`.
  - Bölüm kurucular: `section_confidence`, `section_dataset`, `section_quality`, `section_de`,
    `section_figures`, `section_table`, `section_methods`, `section_references` — her biri HTML string döner.
  - `render_report(inputs, language, version) -> str` — tam self-contained HTML (inline CSS) monte eder.
  - `LABELS: dict[str, dict[str, str]]` — `{"tr": {...}, "en": {...}}` bölüm başlıkları + sabit metinler.
- **`rnaforge/modules/m08_report.py`** — `run_report(config, metadata_path, run_dir, force=False) -> dict`:
  ön koşul m07 (`m07_figures` done değilse `ValueError`); resume/heartbeat; `report/report.html` +
  `statistics/report_statistics.json` yazar; `state.mark_done("m08_report", ...)`. `MODULE_NAME="m08_report"`.
  Döner: `{"report": "report/report.html", "language", "n_sections", "resumed"?}`. YENİ gate YOK.
- **`rnaforge/cli.py`** — `report` subparser + `_cmd_report` + dispatch (m07 `figures` deseni birebir);
  "report OK: <path>" + run dir + güvence verdict basar.

## 4. Tüketilen sözleşmeler (hepsi zaten üretiliyor)

| Girdi | Bölüm |
|---|---|
| `statistics/raw_statistics.json` | Dataset & Örnek Bilgisi |
| `statistics/qc_statistics.json`, `trimming_statistics.json` | Kalite & İşleme |
| `statistics/alignment_statistics.json` | hizalama oranı (örnek-başı) |
| `statistics/count_statistics.json` | atama oranı, n_genes |
| `statistics/de_statistics.json` | DE özeti (n_significant, kontrast, eşikler) |
| `quality/confidence_card.json` | Güvence banner (verdict + kapılar) |
| `figures/manifest.json` (+ PNG'ler) | Figürler (base64 gömülü) |
| `differential_expression/deseq2_results.tsv` | Top 50 DEG tablosu |
| `figures/gene_map.tsv` | tablo gen adları (locus_tag→gene) |

## 5. Rapor bölümleri (sıra)

1. **Başlık** — organizma · run id · üretim tarihi · pipeline sürümü.
2. **Güvence Kartı (belirgin, renk-kodlu banner)** — verdict (TRUSTWORTHY=yeşil / SUSPECT=amber /
   INVALID=kırmızı / UNKNOWN=gri) + PASS/WARN/FAIL sayıları + kapı tablosu (ad, durum, ölçüm, eşik) +
   profil adı + ezilen eşikler. En üstte; kalite felsefesinin merkezi.
3. **Dataset & Örnek Bilgisi** — organizma/platform/tasarım/koşullar + örnek tablosu
   (id, koşul, batch, eşleşme, ort. okuma uzunluğu, ort. kalite).
4. **Kalite & İşleme** — QC/trimming özeti (min_length, agresif=kapalı) + örnek-başı hizalama oranı +
   atama oranı (kapı durumu etiketiyle).
5. **DE Sonuçları** — kontrast, FDR/log2FC eşikleri, n_genes, **n_significant (kaç UP / kaç DOWN)**,
   min replika korelasyonu (WARN ise damga) + tek cümle **sayısal** otomatik özet.
6. **Figürler** — PCA · Volcano · Heatmap · MA gömülü `<img>` + TR/EN başlık/altyazı.
7. **Top 50 DEG tablosu** — gen adı · log2FC · padj · baseMean · yön (Up/Down). Statik (JS yok).
   Tam tablo için `differential_expression/deseq2_results.tsv` notu. 0 DEG → "anlamlı gen yok" notu.
8. **Methods** — config'ten parametrize metin (araç+sürüm, nazik trimming gerekçesi, bowtie2,
   featureCounts feature_type/attribute, DESeq2 tasarımı, eşikler).
9. **References** — sabit atıf listesi (bowtie2, fastp, featureCounts, DESeq2, ggplot2, Williams 2016).

## 6. Stil

Gömülü inline `<style>` (tek dosya CSS), modern/temiz/okunur, yazdırmaya uygun (`@media print`),
açık tema. Verdict banner renk-kodlu. Tüm kullanıcı metni `html.escape` ile kaçışlı (enjeksiyon yok).

## 7. Çıktı sözleşmesi

- `runs/<ts>_<id>/report/report.html` — tek self-contained dosya.
- `runs/<ts>_<id>/statistics/report_statistics.json` — `{"report","language","n_sections"}` (resume için).
- `runs/<ts>_<id>/logs/report.log`.

## 8. Test stratejisi

- **Birim (`tests/test_report_html.py`):** `top_degs` (padj sıralama, gen-adı eşleme, yön sınıflaması,
  <50 ve 0-DEG kenar); `embed_png` (data-URI prefix); her bölüm kurucu (girdi dict → HTML beklenen değeri
  içeriyor, `html.escape` uygulanıyor); `render_report` (verdict + 4 `<img data:image/png;base64` + tüm
  bölüm başlıkları var; dil `tr` vs `en` başlıkları değiştiriyor).
- **Orkestrasyon + entegrasyon (`tests/test_m08_report.py`):** ön koşul m07 yoksa `ValueError`; sentetik
  run dizininden tam rapor → tek dosya, verdict + tüm bölümler + 4 gömülü img; resume `True`.
- Env-gated GEREKMEZ (saf Python; R/dış araç yok). Tüm testler `rnaforge-core`'da koşar.

## 9. Kapsam DIŞI (Faz 2+)

PDF, MultiQC gömme, workflow diyagramı, interaktif/JS tablo, GO/KEGG, dashboard. MVP: statik HTML.

## 10. Kabul kriteri

`rnaforge report --config … --metadata … --run-id …` → `report OK: …/report/report.html`;
tek dosya tarayıcıda açılıyor, 4 figür gömülü görünüyor, verdict m06/m07'deki ile AYNI; 0-DEG koşusunda
çökmeden "anlamlı gen yok" raporu üretiliyor; `report.language: en` config'i İngilizce başlık veriyor.
