# Rapor & Figür Zenginleştirme · Tasarım Spec'i

**Tarih:** 2026-08-04 · **Dal:** `feat/report-enrichment` · **Referans:** `PLAN.md` v1.3 §9
**Ön koşul:** m06/m07/m08 `main`'de (prokaryot MVP tamam). Bu iş m06+m07+m08'i zenginleştirir.

## 1. Amaç

Güncel DE raporlama konvansiyonlarına (nf-core/differentialabundance, DEGreport, EnhancedVolcano,
DESeq2 vignette) hizalanmak: **4 yeni standart figür**, **up/down ayrı tabloları** ve **açıklayıcı
içerik** (figür-altı + bölüm-başı) eklemek. Yalnız koşunun GERÇEK verisi; **uydurma biyolojik yorum
YOK** (`feedback_uydurma_kavram`). Değişmeyen ilkeler: yeni veri-kapısı yok, verdict m06'dan taşınır,
self-contained tek HTML, `tr|en`, saf Python+R (yeni bağımlılık yok), TDD.

## 2. Kapsam — üç modül

### 2.1 m06 — `rnaforge/scripts/deseq2.R` (dispersiyon verisi)
- `dds`'ten ek çıktı: `dispersions.tsv` sütunlar `gene_id`, `baseMean`, `dispGeneEst`, `dispFit`,
  `dispFinal`. `mcols(dds)$dispFit` fallback (gene-est) yolunda olmayabilir → yoksa `NA` yaz.
- Karar/kapı YOK; yalnız ek dosya. Diğer çıktılar (deseq2_results, normalized_counts, de_metrics)
  değişmez. Ön koşul zinciri ve `run_deseq2` Python imzası değişmez (aynı `out_dir`'e bir dosya daha).
- `de_statistics.json`'a `n_up`/`n_down` eklenir (m06 `modules/m06_de.py`): sırasıyla
  `padj<fdr & log2FC>=lfc` ve `padj<fdr & log2FC<=-lfc` sayıları. (Rapor up/down özetinde kullanılır.)

### 2.2 m07 — `rnaforge/scripts/figures.R` + `rnaforge/figures.py` (4 yeni figür)
- `FIGURE_SPECS` 4→8, **anlatı sırasıyla** (id, basename, title):
  1. `("pca","01_pca","PCA")`
  2. `("sample_correlation","02_sample_correlation","Örnek korelasyonu")`
  3. `("expression_dist","03_expression_dist","Ekspresyon dağılımı")`
  4. `("dispersion","04_dispersion","Dispersiyon")`
  5. `("pval_histogram","05_pval_histogram","p-değeri dağılımı")`
  6. `("volcano","06_volcano","Volcano")`
  7. `("ma","07_ma","MA plot")`
  8. `("heatmap","08_heatmap","Heatmap")`
  (Mevcut 4'ün basename'leri yeniden numaralandırılır — `runs/` gitignore'lu, dış etki yok.)
- Yeni figürler (`figures.R`):
  - **sample_correlation:** `cor(log2(nc+1))` Pearson matrisi + `hclust` sıralı tile heatmap
    (renk 0.9–1.0 vurgulu). <2 örnek → boş-durum paneli (m07 dersi: çökme yok).
  - **expression_dist:** örnek-başı `log2(nc+1)` boxplot (koşula göre renkli).
  - **dispersion:** `dispersions.tsv` log-log — `baseMean` vs `dispGeneEst` (nokta) + `dispFinal`
    (çizgi/nokta); DESeq2 plotDispEsts tarzı. Dosya yoksa yüksek sesle hata (m06 artık üretir).
  - **pval_histogram:** `de$pvalue` histogramı (0–1, 40 kova); `pvalue` NA'ları atılır.
- `run_figures_r` yeni argüman alır: `dispersions.tsv` yolu (`de_dir/dispersions.tsv`). figures.R arg
  sırası buna göre güncellenir (runner ile birebir).
- Kenar: 0 anlamlı DEG mevcut heatmap/volcano korumasını (Task m07) korur; yeni figürler DEG'e
  bağlı değil (QC/diagnostik) → doğal çalışır.

### 2.3 m08 — `rnaforge/report_html.py` (tablolar + açıklama)
- **Up/Down ayrı tabloları:** `section_table` → `top_degs`'i yön filtresiyle kullanan iki alt-tablo.
  Yeni yardımcı `top_degs_by_direction(de, gene_map, fdr, lfc, direction, n=25)`; bölüm tek
  `<h2>En Güçlü DEG'ler</h2>` altında **ARTAN (Up) top-25** + **AZALAN (Down) top-25** iki `<table>`.
  Her iki yön de boşsa `no_degs` notu; biri boşsa o alt-tabloda "yok" notu.
- **Figür-altı caption:** `FIGURE_CAPTIONS: dict[str,dict[str,str]]` (`{"tr":{id:metin},"en":{...}}`),
  her figür için "ne gösterir + nasıl okunur" (sabit, çift dilli). `section_figures` her figürün
  altına manifest `title` + caption basar. Manifest'te olmayan id için caption boş geçilir (kırılmaz).
- **Bölüm-başı açıklama:** `SECTION_INTRO: dict[str,dict[str,str]]` (bölüm id → çift dilli 1-2 cümle
  sabit metin). Her `section_*` başlığından hemen sonra `<p class="intro">`. Yalnız sabit/sayısal;
  uydurma biyoloji yok.
- `N_SECTIONS` değişmez (9); yalnız içerik zenginleşir.

## 3. Tüketilen/üretilen sözleşmeler (değişiklikler)

| Dosya | Değişiklik |
|---|---|
| `differential_expression/dispersions.tsv` | **YENİ** (m06 üretir, m07 tüketir) |
| `statistics/de_statistics.json` | `n_up`, `n_down` **eklenir** (m06) |
| `figures/manifest.json` | 4→8 figür (m07) |
| `report/report.html` | up/down 2 tablo + caption + bölüm-intro (m08) |

## 4. Test stratejisi

- **m06:** `test_deseq2.py`/`test_m06_de.py` — `dispersions.tsv` üretiliyor + beklenen sütunlar;
  `de_statistics.json` `n_up`+`n_down` = `n_significant` (yön ayrımı toplamı tutarlı).
- **m07:** `test_figures.py` — `FIGURE_SPECS` 8 eleman + id/sıra; `run_figures_r` yeni argüman.
  `test_m07_figures.py` env-gated entegrasyon — 8 figür PNG+SVG üretiliyor (>1KB), 0-DEG dahil.
- **m08:** `test_report_html.py` — `top_degs_by_direction` (yön filtresi + sıra); iki alt-tablo
  (Up/Down başlıkları + gen adları); `FIGURE_CAPTIONS`/`SECTION_INTRO` tr≠en, `section_figures`
  caption basıyor, `render_report` intro paragrafları içeriyor.
- Tümü mevcut env'lerde (`rnaforge-core` + env-gated `rnaforge-de`).

## 5. Kapsam DIŞI (Faz 2+)

GO/KEGG enrichment, pathway (Manhattan/gprofiler2), 3D PCA, sample-distance dendrogram ayrı figür,
interaktif/JS grafik, PDF. Otomatik biyolojik yorum (Discussion/Conclusion) — hâlâ YOK.

## 6. Kabul kriteri

Gerçek GSE300731'de (m06 `--force` → m07 → m08 zinciri): `report.html` **8 gömülü figür** + her figür
altında açıklama + her bölüm başında intro + **ARTAN 25 / AZALAN 25 iki tablo**; verdict SUSPECT
(değişmez); `dispersions.tsv` ve `de_statistics.json` `n_up`/`n_down` mevcut; 0-DEG sentetik koşusunda
hiçbir figür çökmez. `report.language: en` İngilizce caption/intro verir.
