# m09 — GO Fonksiyonel Zenginleştirme (ORA) · Tasarım Spec'i

**Tarih:** 2026-08-05 · **Dal:** `feat/m09-go-enrichment` · **Referans:** `PLAN.md` v1.3
**Ön koşul:** Prokaryot MVP tamam (`main`): m01→m08 uçtan uca, zengin HTML rapor.
**Yol:** A yalın yol (Ali onayladı 2026-08-04). Ağır DB yolu B (eggNOG) BİLİNÇLİ sonraya.

## 1. Amaç

DESeq2 sonuçlarından **artan** ve **azalan** DEG'ler için ayrı ayrı **GO over-representation
analizi (ORA)** yürütmek; sonuç tablosu + figür üretmek ve **m08 HTML rapora yeni bölüm** eklemek.
Yorumlama değeri yüksek; disk yükü ~50 MB. Değişmeyen ilkeler:

- **YENİ veri-kapısı YOK.** Zenginleştirme yorumlayıcıdır; kötü biyolojiyi geçersiz kılmaz.
  Verdict m06/m07'den **değişmeden taşınır** (m06/m07 gibi).
- **Yalancı sonuç yok** ([[feedback_dogruluk_kontrol]]): annotation birleştirmesi belirsiz eşleşmeyi
  **tahmin etmez, atar**; her GO kaydı kaynak-damgalıdır.
- **Sessiz hata yok** ([[feedback_gurultulu_hata]]): eksik referans (obo/GAF) → yüksek sesle hata + talimat.
- **Uydurma yorum yok** ([[feedback_uydurma_kavram]]): rapor yalnız koşunun gerçek sayılarını gösterir.
- Saf Python (ORA) + R/ggplot2 (figür, m07 deseni); yeni ağır bağımlılık yok (stdlib `math.comb`).
- TDD; modül deseni m07'yle bire bir (saf yardımcı + runner + orkestratör + CLI).

## 2. Annotation birleştirme (`rnaforge/go_annotation.py` — saf Python)

Gen (locus_tag) → GO eşlemesini üç kaynaktan kurar:

### 2.1 GFF otorite (birincil)
`config.reference.annotation_gff`'in CDS satırlarından her `locus_tag` için:
- `Ontology_term=GO:...` → GO id listesi.
- `go_process=` / `go_function=` / `go_component=` → her GO id'nin **namespace**'i
  (BP/MF/CC) ve **okunabilir adı** (`ad|id||kanıt` formatından ayrıştırılır).

GSE300731 referansında kapsama: 4416 CDS'in **2278'i** GO'lu (doğrulandı).

### 2.2 GAF doldurma (tamamlayıcı)
EBI-GOA GAF (E. coli), **yalnızca GFF'te hiç GO'su olmayan** genlere eklenir:
- Birleştirme anahtarı **tam + benzersiz gen sembolü** (GFF `gene=` ↔ GAF `db_object_symbol`, sütun 3).
- Belirsiz eşleşme (sembol GFF'te ≥2 locus_tag'e karşılık geliyorsa **veya** GAF'ta ≥2 farklı
  UniProt id'ye) → **ATILIR** (tahmin yok).
- GAF'tan gelen GO'nun namespace'i sütun 9 (`P`/`F`/`C` → BP/MF/CC); adı obo'dan alınır.
- Her GO kaydı `source ∈ {GFF, GOA}` damgalı tutulur.

### 2.3 Propagation (`go-basic.obo`)
`go-basic.obo` DAG'ı ayrıştırılır (`[Term]` blokları: `id`, `name`, `namespace`, `is_a`,
`relationship: part_of`). Her genin doğrudan GO seti, `is_a` + `part_of` kenarları izlenerek tüm
**ata-terimlere** yayılır (geçişli kapanış, döngü-korumalı). `obsolete: true` terimler atlanır.

**Çıktı:** bellek-içi `gene2go: dict[str, set[str]]` (yayılmış) + `go_meta: dict[str, (namespace, name)]`
+ denetim izi `runs/.../enrichment/gene2go.tsv` (gene, go_id, namespace, name, source, direct|propagated).

## 3. ORA motoru (`rnaforge/enrichment.py` — saf Python, stdlib)

- **Setler:** UP (`padj < de.fdr_threshold & log2FoldChange ≥ +de.log2fc_threshold`) ve DOWN (`≤ −lfc`),
  `deseq2_results.tsv`'den. Eşikler m06 config'inden (tek kaynak, drift yok).
- **Arka plan (background):** DE matrisinde **GO'su olan** tüm genler (test edilen evren).
- **Test:** her (set × namespace) için, her GO terimi için tek-yönlü **hipergeometrik** over-representation
  p-değeri, `math.comb` ile:
  p = Σ_{i=k}^{min(K,n)} C(K,i)·C(N−K,n−i) / C(N,n)
  (N=arka plan, K=terimdeki arka plan geni, n=set büyüklüğü, k=terimdeki set geni).
- **Gürültü filtresi:** arka planda `< enrichment.min_term_size` (varsayılan 3) gene sahip terimler atlanır.
- **Çoklu-test:** **BH FDR**, her (set × namespace) grubu içinde ayrı.
- **Fold-enrichment:** (k/n) / (K/N); expected = n·K/N.

**Çıktı:** `enrichment/enrichment_up.tsv`, `enrichment/enrichment_down.tsv`. Sütunlar:
`go_id, namespace, term, study_count(k), study_n(n), bg_count(K), bg_n(N), expected, fold_enrichment,
p_value, p_adj, genes(sembol;…)`. padj artan sıralı. Set boşsa başlıklı boş dosya (çökme yok).

## 4. Figür (`rnaforge/scripts/enrichment.R`, `rnaforge-de` env — m07 deseni)

- UP ve DOWN için **top-N (varsayılan 15) zenginleşmiş terim dot plot**: x=fold_enrichment,
  boyut=study_count, renk=p_adj; namespace'e (BP/MF/CC) göre facet. `p_adj < 0.05` süzülür,
  padj'e göre sıralı. PNG 300dpi + SVG → `enrichment/`.
- **Boş durum** (anlamlı terim yok / set boş) → **çökme yok**, dürüst boş-durum paneli (m07 kritik-bug dersi).
- `manifest.json` (m07 `build_manifest` deseni yeniden kullanılır): `enrichment_up`, `enrichment_down`.

## 5. Orkestrasyon (`rnaforge/modules/m09_enrichment.py`)

`run_enrichment(config, metadata_path, run_dir, force=False) -> dict` (m07 imza deseni):
- Ön koşul **m06** (`state.is_done("m06_de")`); değilse `ValueError` (net mesaj, aynı `--run-id`).
- Resume: `state.is_done("m09_enrichment")` + stats varsa `resumed=True` döndür.
- Sıra: annotation kur → ORA (up/down TSV) → figür → manifest → `statistics/enrichment_statistics.json`
  (`n_terms_up`, `n_terms_down`, `n_sig_up`, `n_sig_down`, `background_size`, `n_annotated`, kaynak sayıları).
- **Gate YOK.** Güvence kartı verdict m06/m07'den taşınır (m09 dokunmaz).
- Eksik `obo`/`gaf` dosyası → yüksek sesle hata (`FileNotFoundError`/`ValueError`) + nereden alınacağı talimatı.
- `logs/enrichment.log`'a R stdout/stderr + özet.

**Zincir:** m06 → m07 → **m09** → m08. m09 ön koşulu **m06** (m07'ye bağlı değil; figürleri bağımsız).
m08 ön koşulu **m07 + m09** olur (rapor GO bölümünü gömer).

## 6. CLI (`rnaforge enrich`)

Yeni subcommand `rnaforge enrich` (`--config`, `--metadata`, `--run-id`, `--force`), m01 ön koşullu
(aynı `--run-id` → aynı `run_dir`). `_cmd_enrich` → `run_enrichment`; güvence kartı basılır (verdict taşınır).
`main()` dispatch'ine eklenir.

## 7. Config & referanslar

- Yeni üst-seviye anahtar **`enrichment`** → `KNOWN_TOP_LEVEL_KEYS`'e eklenir (bilinmeyen anahtar
  reddi korunur, sessiz yutma yok). Yeni `EnrichmentConfig` (frozen dataclass):
  - `min_term_size: int = 3`
  - `top_n: int = 15` (figür)
  - `obo: Path | None` (varsayılan `references/go/go-basic.obo`)
  - `gaf: Path | None` (opsiyonel; verilmezse GAF doldurma atlanır — GFF+propagation yeterli, sessiz DEĞİL: log'a "GAF yok, yalnız GFF" yazılır)
- **Referans dosyaları (ikisi de `.gitignore`'lu):**
  - `references/go/go-basic.obo` (~35-40 MB) — kaynak `http://purl.obolibrary.org/obo/go/go-basic.obo`
  - `references/ecoli_bw25113/ecoli.gaf` — kaynak EBI-GOA (E. coli); indirme adımında **doğrulanır**
    (körlemesine güvenme — dosya formatı/sütunları kontrol edilir).
- **İndirme = ayrı hazırlık adımı** (`docs`'ta belgelenir; opsiyonel yardımcı script). Modül dosyaları
  yalnız **okur**; yoksa gürültülü hata.

## 8. m08 rapor bölümü (`rnaforge/report_html.py`)

- Yeni bölüm **"Fonksiyonel Zenginleştirme (GO)"** (DE Sonuçları'ndan sonra, Methods'tan önce).
- `load_report_inputs`'a `enrichment_up.tsv`/`enrichment_down.tsv` + `enrichment_statistics.json`
  + enrichment figür manifesti eklenir (yoksa bölüm dürüstçe "zenginleştirme çalıştırılmadı" notu — kırılmaz).
- İçerik: UP ve DOWN için namespace başına **top-N anlamlı terim tablosu** (go_id, term, namespace,
  study/bg, fold-enrichment, padj) + gömülü dot-plot figürleri. Çift dilli caption + sabit/sayısal intro
  (uydurma biyoloji yok). m08 **organizma-agnostik** kalır; verdict değişmez; yeni kapı yok.

## 9. Tüketilen/üretilen sözleşmeler

**Tüketilen:** `differential_expression/deseq2_results.tsv` (m06) · `config.reference.annotation_gff` ·
`references/go/go-basic.obo` · opsiyonel `references/.../ecoli.gaf` · `config.de.{fdr,log2fc}_threshold`.
**Üretilen:** `enrichment/{gene2go.tsv, enrichment_up.tsv, enrichment_down.tsv, *.png, *.svg, manifest.json}` ·
`statistics/enrichment_statistics.json`. m08 bunları tüketir.

## 10. Doğrulama

- **Birim testler** (`rnaforge-core`, monkeypatch):
  - GFF GO parser (Ontology_term + go_process/function/component → id/namespace/name).
  - GAF doldurma: tam-benzersiz eşleşme eklenir; belirsiz eşleşme atılır; kaynak damgası.
  - obo parser + propagation (is_a/part_of geçişli kapanış, döngü/obsolete koruması).
  - hipergeometrik p (bilinen küçük örnekte elle hesapla) + fold-enrichment + expected.
  - BH FDR (bilinen p vektöründe).
  - min_term_size filtresi; boş set → başlıklı boş TSV, çökme yok.
  - orkestrasyon: m06 ön koşul yoksa ValueError; resume; gate yok (verdict taşınır).
  - config: `enrichment` anahtarı parse; bilinmeyen alt-anahtar reddi.
- **Entegrasyon** (env-gated skip): gerçek `enrichment.R` bir küçük TSV'yi PNG'ye render eder.
- **GSE300731 canlı smoke:** obo+GAF indir → `rnaforge enrich` → beklenti: DOWN'da asit-direnç
  (gad*/hde*) ve UP'ta zarf-stres/kapsül ilgili GO terimleri anlamlı zenginleşsin (biyolojik akıl kontrolü).
  m08 raporunda GO bölümü + figürler gömülü; verdict SUSPECT (değişmez).

## 11. İş akışı

brainstorm ✅ → **bu spec** → writing-plans → TDD (~8-9 task) → requesting-code-review → GSE300731 smoke
→ `feat/m09-go-enrichment` → `main` merge + push → DURUM.md + bellek güncelle.
