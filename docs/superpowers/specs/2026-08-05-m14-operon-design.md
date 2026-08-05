# m14 — Operon Analizi · Tasarım Spec'i

**Tarih:** 2026-08-05 · **Dal:** `feat/m14-operon` · **Referans:** `PLAN.md` v1.3
**Ön koşul:** Dalga 2 #A (m13 AMR) `main`'de. Bu iş **Dalga 2 #B**.
**İlke:** intergenik-mesafe sezgiseli (saf Python, GFF'ten); organizma-agnostik; DB/araç yok. Gate YOK.

## 1. Amaç

Operon-düzeyi DE koordinasyonu: aynı yönde bitişik + kısa intergenik boşluklu genler aynı operon
(birlikte transkribe). Tedaviye **birlikte yanıt veren transkripsiyon birimlerini** ortaya çıkarır
(gen-düzeyi DE'yi tamamlar). Değişmeyen ilkeler: gate YOK (verdict m06'dan taşınır), uydurma yok
(operonlar **tahmin**, deneysel değil — rapora dürüstçe yazılır), TDD, saf Python.

**Neden intergenik-mesafe:** operon yapısının en güçlü tek belirleyicisi (Moreno-Hagelsieb & Collado-Vides
2002). Harici araçlar (Operon-mapper = web servisi, veri dışarı; Rockhopper = ağır) yerel/tekrarlanabilir
felsefeye uymaz. Ali kararı: yalın yerel sezgisel.

## 2. Operon tahmini (`rnaforge/operon.py`)

- `predict_operons(gff, max_gap=50) -> list[dict]`: GFF `gene` feature'ları `(contig, start)` sıralı;
  ardışık genler **aynı strand** ve `gap = next.start − cur.end − 1 ≤ max_gap` ise aynı operonda toplanır;
  strand değişimi / `gap > max_gap` / contig değişimi → operon sınırı. Her operon:
  `{operon_id, contig, strand, locus_tags:[…], symbols:[…], size}`. Tek-genli operonlar da üretilir.
  (Örtüşen genler gap<0 → aynı operon.) thrABC gibi bilinen operonları yakalar.

## 3. DE koordinasyonu (`rnaforge/operon.py`)

- `aggregate_operon_de(operons, deseq_tsv, fdr, lfc) -> list[dict]`: operon başına deseq2'den
  `log2FC/padj`; metrikler: `size`, `n_tested`, `n_deg` (padj<fdr & |log2FC|≥lfc), `n_up`, `n_down`,
  `mean_log2fc` (test edilenlerin), `coordinated` (≥2 gen & ≥2 DEG & hepsi aynı yön → birlikte-düzenlenen).
  Koordineli + `n_deg` azalan sıralı.

## 4. Modül (`rnaforge/modules/m14_operon.py`) + CLI

`run_operon(config, metadata_path, run_dir, force=False) -> dict` (m13 deseni):
- Ön koşul **m06** (DE); GFF gerekli. **Gate YOK.** Resume/heartbeat.
- Çıktı `operon/operons.tsv` (operon_id, contig, strand, genes, size, n_tested, n_deg, n_up, n_down,
  mean_log2fc, coordinated). `statistics/operon_statistics.json` (n_operons, n_multi_gene, n_coordinated,
  max_gap). Yeni `rnaforge operon` subcommand. Zincir m06 tabanlı.

## 5. Rapor (`report_html.py`)

Yeni **"Operon Analizi"** bölümü (AMR'dan sonra, Methods'tan önce). Top **koordineli** operonlar tablosu
(üyeler, boyut, n_DEG, yön, mean log2FC) + "N operon tahmin edildi, M çok-genli, K koordineli DE" özeti.
**Dürüst not:** operonlar intergenik-mesafeyle **tahmin** edildi (deneysel doğrulanmadı). Çift dilli,
tolerant. Yöntem/Kaynak: Moreno-Hagelsieb & Collado-Vides 2002. verdict değişmez.

## 6. Config

Yeni top-level `operon` (`OperonConfig`): `max_gap: int = 50` (bp). `KNOWN_TOP_LEVEL_KEYS`'e.

## 7. Doğrulama

- **Birim** (`rnaforge-core`): `predict_operons` (aynı-yön birleşme, gap sınırı, strand/contig kırılımı,
  tek-gen), `aggregate_operon_de` (koordineli/karışık-yön/tek-gen; n_deg; mean), config, rapor bölümü tolerant.
- **GSE300731 canlı smoke:** thrABC operonu tahmin edilmeli; enterobaktin (ent*) + efflux (acrAB/emrAB)
  operonları **koordineli DE** çıkmalı (m13 AMR sonuçlarıyla uyumlu — birlikte indüklenen birimler).
  Rapor bölümü; verdict SUSPECT değişmez.

## 8. İş akışı

spec → writing-plans → TDD (~6 task) → GSE300731 smoke → `feat/m14-operon` → `main` merge + push →
DURUM + bellek. Sonra Dalga 2 #C (PPI + community). Bkz. [[rnaforge-project]].
