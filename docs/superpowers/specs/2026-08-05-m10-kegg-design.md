# m10 — KEGG Pathway ORA · Tasarım Spec'i

**Tarih:** 2026-08-05 · **Dal:** `feat/m10-kegg` · **Referans:** `PLAN.md` v1.3
**Ön koşul:** m09 GO enrichment `main`'de. Bu iş Dalga 1'in ilk parçası (KEGG → GSEA → semantic).
**İlke:** m09 motorunu (`enrichment.py`) DEĞİŞTİRMEDEN yeniden kullan; organizma-agnostik (ökaryota taşınır).

## 1. Amaç

DESeq2 DEG'leri için **KEGG pathway over-representation analizi (ORA)**. Artan/azalan ayrı; sonuç
tablosu + figür + m08 rapora KEGG alt-bölümü. Değişmeyen ilkeler (m09 ile aynı):
- **YENİ veri-kapısı YOK**; verdict m06/m07'den taşınır.
- **Yalancı anotasyon yok** ([[feedback_dogruluk_kontrol]]): sembol join belirsizse ATILIR (tahmin yok).
- **Sessiz hata yok** ([[feedback_gurultulu_hata]]): eksik KEGG dosyası → gürültülü hata + talimat.
- **Uydurma yorum yok**: rapor yalnız gerçek sayılar.
- Saf Python + mevcut R; yeni ağır bağımlılık yok. TDD.
- **Motor yeniden kullanımı:** `deg_sets`/`run_ora`/`bh_fdr`/`write_ora_tsv` (`enrichment.py`) DEĞİŞMEZ.

## 2. Annotation (`rnaforge/kegg_annotation.py` — saf Python)

KEGG REST 3 dosyasından gen(locus_tag)→pathway kurar (prep'te indirilir, gitignore'lu):
- `link/pathway/<org>` → `eco:b0114 \t path:eco00010` (b-number → pathway).
- `list/pathway/<org>` → `eco00010 \t Glycolysis - <organism>` (pathway → ad; " - <organism>" eki kırpılır).
- `list/<org>` → `eco:b0002 \t CDS \t ... \t thrA; fused ...` (b-number → **sembol**; sütun 4'ün ilk
  `;`/`,` öncesi jetonu).

**Join (m09 deseni):** KEGG b-number → sembol → bizim `locus_tag`. Ters harita `gene_symbol`'den
(`go_annotation._symbol_to_locus` YENİDEN KULLANILIR); **çoklu locus'lu sembol ATILIR**. `gene_symbol`
GFF'ten (`go_annotation.parse_gff_go` üçüncü dönüşü).

**Global/overview map filtresi:** ORA'yı bozmamak için KEGG "global and overview" haritaları hariç
tutulur (numerik id sabit seti, organizma-agnostik): `{01100,01110,01120,01200,01210,01212,01220,
01230,01232,01240,01250}`.

**Çıktı:** `gene2pathway: dict[locus_tag, set[pathway_id]]` + `pathway_meta: dict[id, ("KEGG", ad)]`.
Namespace **tek grup "KEGG"** (tek BH grubu — pathway ORA standardı; BRITE kategorileri SONRAYA).

## 3. Motor (yeniden kullanım — YENİ KOD YOK)

`enrichment.py`'deki `deg_sets` (m06 eşikleri), `run_ora(gene_set, background, gene2pathway,
pathway_meta, gene_symbol, min_term_size)`, `bh_fdr`, `write_ora_tsv` **aynen** kullanılır (jenerik:
gen→set + set→(namespace,ad) alır). Arka plan = KEGG-anotasyonlu test edilen genler. `all_tested_genes`
yeniden kullanılır. Çıktı `kegg/kegg_{up,down}.tsv` (m09 ile aynı sütun sözleşmesi).

## 4. Figür (`enrichment.R` parametrize — m09 R yolu bozulmaz)

`run_enrichment_r`'a opsiyonel `title_prefix` + `basename_prefix` argümanı eklenir; `enrichment.R`
sonda iki opsiyonel arg okur (yoksa m09 varsayılanı: "GO zenginleştirme"/"enrichment"). m10 çağrısı
`title_prefix="KEGG Pathway"`, `basename_prefix="kegg"` → `kegg_up.png/svg`, `kegg_down.png/svg`.
`build_enrichment_manifest` parametrik basename ile yeniden kullanılır. Boş-durum paneli korunur.

## 5. Orkestrasyon (`rnaforge/modules/m10_kegg.py`) + CLI

`run_kegg(config, metadata_path, run_dir, force=False) -> dict` (m09 imza deseni):
- Ön koşul **m06**; resume; heartbeat; **GATE YOK** (verdict taşınır).
- `config.enrichment.kegg_organism` yoksa → net ValueError (KEGG kodu zorunlu, ör. "eco").
- KEGG dosyaları (`kegg_dir`) eksikse → FileNotFoundError + REST indirme talimatı (sessiz skip yok).
- Sıra: `build_gene2pathway` → `deg_sets` → `run_ora`(up/down) → `write_ora_tsv` → `run_enrichment_r`
  (KEGG prefix) → manifest → `statistics/kegg_statistics.json` (n_terms/n_sig up-down, background, n_annotated).
- **Zincir:** m06→m07→**m09/m10**→m08 (m10, m09'dan bağımsız; ön koşul yalnız m06).
- CLI: yeni `rnaforge kegg` subcommand (config/metadata/run-id/force) + `_cmd_kegg` (güvence kartı, verdict taşınır).

## 6. Config & referans

`Enrichment` dataclass'a eklenir (KNOWN_TOP_LEVEL_KEYS'te `enrichment` zaten var):
- `kegg_organism: str | None` (ör. "eco"/"hsa"; yoksa m10 çalışmaz — net hata).
- `kegg_dir: Path | None` (varsayılan `references/kegg/<org>/`).
- `min_term_size`/`top_n` yeniden kullanılır.
Referans dosyaları (gitignore'lu, `references/kegg/<org>/`): `pathway_links.tsv`, `pathway_names.tsv`,
`gene_list.tsv` — prep'te KEGG REST'ten indirilir. KEGG akademik kullanım; dosyalar redistribute edilmez.

## 7. Rapor (`report_html.py`)

"Fonksiyonel Zenginleştirme" bölümü **GO + KEGG alt-bölümleri** gösterir (her biri tolerant; koşmadıysa
dürüst not). `load_report_inputs`'a `kegg_{up,down}.tsv` + kegg manifest eklenir (yoksa None).
`section_enrichment` GO'dan sonra KEGG alt-bölümünü (tablo + gömülü figür) basar. Yöntemler'e KEGG cümlesi
+ Kaynaklar'a **Kanehisa & Goto 2000 (KEGG)** (yalnız `kegg_ran`). Tablo açıklaması KEGG'i de kapsar
(pathway = KEGG yolağı). m08 organizma-agnostik kalır; verdict değişmez.

## 8. Tüketilen/üretilen sözleşmeler

**Tüketilen:** `deseq2_results.tsv` (m06) · `config.reference.annotation_gff` · `references/kegg/<org>/*` ·
`config.de.{fdr,log2fc}_threshold` · `config.enrichment.{kegg_organism,kegg_dir,min_term_size,top_n}`.
**Üretilen:** `kegg/{gene2pathway.tsv, kegg_up.tsv, kegg_down.tsv, kegg_*.png/svg, manifest.json}` ·
`statistics/kegg_statistics.json`. m08 tüketir.

## 9. Doğrulama

- **Birim** (`rnaforge-core`): KEGG parser (link/list/gene_list → gene2pathway; " - org" eki kırpma;
  path: prefix normalize); sembol join belirsizlik atma; global-map filtresi; motor yeniden kullanımı
  (run_ora KEGG haritasıyla); orkestrasyon (m06 ön koşul, kegg_organism yok → hata, dosya yok → gürültülü,
  resume, gate yok); rapor KEGG alt-bölümü + tolerant.
- **Entegrasyon** (env-gated): enrichment.R KEGG prefiksiyle PNG üretir.
- **GSE300731 canlı smoke:** `rnaforge kegg --run-id GSE300731`. Beklenti — DOWN: oksidatif fosforilasyon/
  sitrat döngüsü/respirasyon; UP: two-component system / katyonik AMP direnci / biyofilm gibi antibiyotik
  yanıt pathway'leri. Rapor GO+KEGG birlikte; verdict SUSPECT değişmez.

## 10. İş akışı

spec → writing-plans → TDD (~7-8 task) → GSE300731 smoke → `feat/m10-kegg` → `main` merge + push →
DURUM.md + bellek. Sonra Dalga 1 #2: GSEA. Bkz. [[rnaforge-project]].
