# m11 — GSEA (Gene Set Enrichment Analysis) · Tasarım Spec'i

**Tarih:** 2026-08-05 · **Dal:** `feat/m11-gsea` · **Referans:** `PLAN.md` v1.3
**Ön koşul:** m09 (GO ORA) + m10 (KEGG ORA) `main`'de. Downstream Dalga 1 #2.
**İlke:** m09/m10 gen-seti kurucularını yeniden kullan; motor fgsea (altın standart); organizma-agnostik.

## 1. Amaç

DESeq2 sonuçlarından **GSEA** — tüm genlerin **sıralı (ranked) listesi** üzerinde gen-seti zenginleştirme.
ORA'dan farklı: anlamlı DEG alt kümesini değil, **tüm ranked listeyi** kullanır → koordineli zayıf sinyalleri
yakalar. Değişmeyen ilkeler (m09/m10 ile aynı): **gate YOK** (verdict m06'dan taşınır), yalancı sonuç yok,
sessiz hata yok, uydurma yorum yok, TDD, organizma-agnostik.

**Neden fgsea:** GSEA istatistiği (normalize ES, permütasyon/multilevel p) ince; elle yazmak ince-hata riski
([[feedback_dogruluk_kontrol]]). fgsea (Bioconductor) altın standart + hızlı → `rnaforge-de` env'ine kurulur
(DESeq2/ggplot2 zaten orada). DESeq2 kararıyla tutarlı (izole conda + altın standart, Ali onaylı).

## 2. Girdi hazırlığı (`rnaforge/gsea.py` — saf Python)

- **Ranked liste:** `deseq2_results.tsv` → `gene(locus_tag) → stat` (Wald istatistiği; yön+anlamlılığı
  birlikte kodlar, `stat = log2FC/lfcSE`). NA `stat` atılır. `gsea/ranked.rnk` (`gene\tstat`).
- **Gen setleri (GMT):** m09/m10 kurucuları **YENİDEN KULLANILIR** — `go_annotation.build_gene2go`
  (GO, propagate'li) + `kegg_annotation.build_gene2pathway` (KEGG). Gen→set **ters çevrilir** → set→genler
  (locus_tag). GMT satırı: `<set_id>\t<set_name>\t<gene1>\t<gene2>…`. GO ve KEGG **ayrı koleksiyon**
  (`gsea/go.gmt`, `gsea/kegg.gmt`). Boyut filtresi fgsea'da (min/max).
- Fonksiyonlar: `write_rnk(deseq_tsv, out)`, `invert_to_gmt(gene2set, meta, out)`,
  `run_gsea_r(rnk, gmt, gene_map, out_dir, collection, min_size, max_size, title, env)`.

## 3. Motor (`rnaforge/scripts/gsea.R`, `rnaforge-de` env)

- `fgsea::fgsea(pathways, stats, minSize, maxSize)` (multilevel; modern varsayılan). `stats` = adlandırılmış
  vektör (locus_tag→stat), azalan sıralı. `pathways` = GMT'den okunan liste.
- Çıktı `gsea/gsea_<collection>.tsv`: `pathway_id, name, namespace, size, ES, NES, pval, padj, leading_edge`
  (öncü genler locus_tag→sembol eşlenip `;` ile). **NES işaretli** (+ artan tarafta, − azalan tarafta) —
  ayrı up/down koşu YOK.
- **Figür:** koleksiyon başına **NES dot-plot** (ggplot): en güçlü ±NES terimler; x=NES, renk=padj,
  boyut=set boyutu; sıfır çizgisi. Boş-durum panelli. `gsea/gsea_<collection>.png/svg`.

## 4. Orkestrasyon (`rnaforge/modules/m11_gsea.py`) + CLI

`run_gsea(config, metadata_path, run_dir, force=False) -> dict` (m09/m10 deseni):
- Ön koşul **m06**; resume/heartbeat; **GATE YOK**.
- Koleksiyonlar: GO (obo varsa) + KEGG (kegg_organism varsa). **En az biri** gerekli; hiçbiri yoksa net hata.
- `deseq2_results.tsv`'de `stat` yoksa gürültülü hata. Referans eksikse (obo/kegg) o koleksiyon
  atlanır ama **sessiz değil** (log + stats'ta işaretli); hiç koleksiyon kalmazsa hata.
- Çıktı: `gsea/{ranked.rnk, go.gmt, kegg.gmt, gsea_go.tsv, gsea_kegg.tsv, *.png/svg, manifest.json}` +
  `statistics/gsea_statistics.json` (koleksiyon başına n_sets, n_sig_pos, n_sig_neg).
- Zincir: m06→m07→**m09/m10/m11**→m08 (hepsi m06'ya bağlı, bağımsız). Yeni `rnaforge gsea` subcommand.

## 5. Config

`Enrichment`'e `gsea_min_size:int=15`, `gsea_max_size:int=500` (fgsea eşikleri = veri/config ilkesi).
obo/gaf/kegg_organism/kegg_dir/min_term_size yeniden kullanılır. `envs/rnaforge-de.yml`'e
`bioconductor-fgsea` eklenir.

## 6. Rapor (`report_html.py`)

Yeni **"Gen Seti Zenginleştirme (GSEA)"** bölümü (Fonksiyonel Zenginleştirme'den sonra, Methods'tan önce).
`load_report_inputs`'a `gsea_go/gsea_kegg` tsv + manifest (yoksa None). Her koleksiyon için: en güçlü
**pozitif** (artan) ve **negatif** (azalan) NES terim tablosu (term, namespace, NES, padj, size, **öncü genler**)
+ gömülü NES dot-plot. Çift dilli, tolerant. Yöntemler'e GSEA paragrafı + Kaynaklar'a **Subramanian ve ark.
2005** + **Korotkevich ve ark. 2021 (fgsea)** (yalnız `gsea_ran`). m08 agnostik kalır; verdict değişmez.

## 7. Tüketilen/üretilen sözleşmeler

**Tüketilen:** `deseq2_results.tsv` (m06, `stat`) · `annotation_gff` · obo/gaf · KEGG dosyaları ·
`config.enrichment.{gsea_min_size,gsea_max_size,...}`. **Üretilen:** `gsea/*` + `statistics/gsea_statistics.json`.

## 8. Doğrulama

- **Birim** (`rnaforge-core`): `write_rnk` (NA atma, stat), `invert_to_gmt` (ters çevirme + GMT format),
  orkestrasyon (m06 ön koşul, stat yok → hata, koleksiyon yok → hata, resume, gate yok), config alanları,
  rapor GSEA bölümü tolerant.
- **Entegrasyon** (env-gated): gsea.R küçük rnk+gmt'de fgsea NES üretir.
- **GSE300731 canlı smoke:** `rnaforge gsea`. Beklenti — GSEA, ORA'nın zarf-stres/respirasyon temalarını
  işaretli NES ile doğrular (UP: kolanik asit/peptidoglikan pozitif NES; DOWN: oksidatif fosforilasyon
  negatif NES) + koordineli zayıf setleri yakalayabilir. Rapor GSEA bölümü; verdict SUSPECT değişmez.

## 9. İş akışı

spec → writing-plans → fgsea kur + CLI doğrula → TDD (~7-8 task) → GSE300731 smoke → `feat/m11-gsea` →
`main` merge + push → DURUM.md + bellek. Sonra Dalga 1 #3: semantic similarity + REVIGO. Bkz. [[rnaforge-project]].
