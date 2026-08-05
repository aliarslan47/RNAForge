# m12 — Semantic Reduction (REVIGO-benzeri) · Tasarım Spec'i

**Tarih:** 2026-08-05 · **Dal:** `feat/m12-semantic` · **Referans:** `PLAN.md` v1.3
**Ön koşul:** m09 (GO ORA) + m11 (GSEA) `main`'de. Downstream Dalga 1 #3 (son).
**İlke:** obo + m09 kurucusunu yeniden kullan; saf Python stdlib (numpy YOK); organizma-agnostik. Gate YOK.

## 1. Amaç

GO zenginleştirmesi çok **fazlalık** terim üretir (parent/child, benzer süreçler; GSE300731: 58+51 ORA GO).
Semantik benzerlikle kümeleyip her kümeden **en iyi skorlu temsilciyi** tutarak öz, okunur liste üretmek
(REVIGO fikri). Kaynaklar: m09 ORA GO (up/down) + m11 GSEA GO. Değişmeyen ilkeler: gate YOK (verdict m06'dan
taşınır), uydurma yok, sessiz hata yok, saf Python, TDD.

## 2. IC + Lin benzerliği (`rnaforge/semantic.py`)

- **IC(t)** = −log( count(t) / N ); count(t) = arka planda (propagate'li `build_gene2go`) t'ye anotlı gen
  sayısı, N = toplam anotlı gen. m09 `build_gene2go` yeniden kullanılır (obo+gff+gaf).
- **Lin(a,b)** = 2·IC(MICA)/(IC(a)+IC(b)); MICA = en yüksek IC'li ortak ata (`go_annotation._ancestors`
  yeniden kullanılır). [0,1]. a==b→1; payda 0 (kök)→0.
- Benzerlik **YALNIZ aynı namespace içinde** anlamlı (BP/MF/CC ayrı — farklı namespace ortak ata yok, kök hariç).

## 3. REVIGO-benzeri indirgeme

`reduce_terms(terms, obo, ic, threshold) -> list[dict]`:
- Girdi: `[{go_id, namespace, term, padj}]` (yalnız anlamlı, padj<0.05). Namespace başına grupla.
- Greedy: padj artan sırala (en iyi önce). Her terim mevcut temsilcilere **max Lin ≥ threshold** ise o kümeye
  (fazlalık), değilse yeni temsilci. `threshold = config.enrichment.revigo_similarity` (varsayılan **0.7**).
- Çıktı: temsilci satırı + `n_collapsed` (temsil ettiği terim sayısı) + `members` (go_id;…).

## 4. Modül (`rnaforge/modules/m12_semantic.py`) + CLI

- Ön koşul: **m09 VEYA m11** çalışmış + `config.enrichment.obo`; hiç GO kaynağı/çıktısı yoksa net hata.
- Kaynak çıktıları: `enrichment/enrichment_{up,down}.tsv` (namespace TSV'de) + `gsea/gsea_go.tsv` (namespace
  `go_meta`'dan — build_gene2go zaten üretir). Her birini indirger.
- Çıktı: `semantic/reduced_{ora_up,ora_down,gsea_go}.tsv` (sütunlar: go_id, namespace, term, padj, n_collapsed,
  members) + `statistics/semantic_statistics.json` (kaynak başına n_terms→n_representatives). **Gate YOK.**
- **Figür YOK** (MDS scatter numpy/R ister → bilinçli sonraya); değer öz tabloda. Mevcut tam figürler kalır.
- Yeni `rnaforge semantic` subcommand. Zincir: m09/m11 sonrası; m06 tabanlı.

## 5. Config

`Enrichment`'e `revigo_similarity: float = 0.7`. obo/gaf yeniden kullanılır.

## 6. Rapor (`report_html.py`)

Yeni **"Anlamsal İndirgeme (REVIGO)"** bölümü (GSEA'dan sonra, Methods'tan önce). `load_report_inputs`'a
`reduced_{ora_up,ora_down,gsea_go}` (yoksa None). Kaynak başına temsilci tablosu (term, namespace, padj,
"temsil ettiği terim sayısı") + "N terim → M temsilci" özeti. Çift dilli, tolerant. Yöntem/Kaynak: Lin 1998 +
Supek ve ark. 2011 (REVIGO) — yalnız `semantic_ran`. verdict değişmez.

## 7. Doğrulama

- **Birim** (`rnaforge-core`): `compute_ic` (bilinen sayımda), `lin_similarity` (küçük DAG'da elle MICA/Lin),
  `reduce_terms` (benzer terimler tek temsilcide; namespace ayrımı; boş liste; eşik davranışı), config alanı,
  rapor bölümü tolerant.
- **GSE300731 canlı smoke:** `rnaforge semantic`. Beklenti — 58 ORA-up terimi birkaç temsilciye insin
  (ör. "polysaccharide biosynthetic/metabolic process" ailesi tek temsilcide); rapor öz tablo; verdict SUSPECT değişmez.

## 8. İş akışı

spec → writing-plans → TDD (~6 task) → GSE300731 smoke → `feat/m12-semantic` → `main` merge + push →
DURUM.md + bellek. **Dalga 1 BİTER.** Sonra Dalga 2 (AMR/virulence…) veya ökaryot. Bkz. [[rnaforge-project]].
