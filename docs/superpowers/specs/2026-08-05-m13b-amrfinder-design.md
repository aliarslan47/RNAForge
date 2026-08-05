# m13b — AMRFinderPlus ikinci araç (yan-yana konkordans) · Tasarım Spec'i

**Tarih:** 2026-08-05 · **Dal:** `feat/m13-amrfinder` · **Referans:** m13 AMR modülü
**Amaç:** m13'e ikinci bağımsız AMR aracı (AMRFinderPlus) ekleyip AMR tablosunda abricate/CARD ile
**yan yana** göstermek (araç konkordansı). Ali isteği + reviewer geri bildirimi.

## 1. Neden

CARD (abricate) geniş taban atar (intrinsik efflux dahil → BW25113'te 43 gen); **AMRFinderPlus** NCBI'ın
küratörlü, yüksek-güven aracıdır (K-12'de ~5: blaEC β-laktamaz, acrF/mdtM efflux). İki aracı yan yana
göstermek **konkordans + güven** katmanı ekler: hangi geni hangi araç çağırıyor. Kapsam m13'e sınırlı;
gate YOK; verdict m06'dan taşınır.

## 2. Araç koşusu (`rnaforge/abricate.py`)

- `run_amrfinder(genome_fa, out_tsv, organism, env="ali-amrfinder") -> str`: `conda run -n <env> amrfinder
  --nucleotide <genome> --organism <organism> --plus --quiet`. `--plus` (K-12'de core=0; küratörlü ek AMR/
  stres). Hata gürültülü.
- `parse_amrfinder(tsv, min_identity, min_coverage) -> list[dict]`: DictReader (tab); kolonlar "Contig id",
  "Start", "Stop", "Element symbol", "Type", "Class", "% Identity/Coverage of reference". Type ∈ {AMR, STRESS}
  tut (direnç ilişkili). Hit **abricate ile aynı dict şekli** → `{contig, start, end, gene, pct_id, pct_cov,
  db="amrfinderplus", product=Element name, resistance=Class}` (map_hits_to_genes/overlay_de yeniden kullanılır).

## 3. Birleştirme (m13_amr.py)

- CARD: abricate → mapped_card (locus_tag → hit). AMRFinderPlus (organism verildiyse): parse → map → mapped_afp.
- **Birleşim:** locus_tag birliği; her gen için `card` = CARD Class/resistance (yoksa "—"), `amrfinder` =
  AMRFinderPlus Class (yoksa "—"), `pct_identity` = ikisinden max, DE overlay birleşik gen setine tek sefer.
- `amr/amr_genes.tsv` sütunları: `gene, locus_tag, card, amrfinder, pct_identity, log2fc, padj, de_status`.
  de_status'a göre sıralı. Virülans (VFDB) tablosu DEĞİŞMEZ (tek araç).
- amrfinder_organism yoksa: eski davranış (yalnız CARD); `amrfinder` sütunu tümü "—" (veya sütun atlanır).
  Basitlik: sütun hep var, AMRFinderPlus koşmadıysa hep "—".
- stats: `n_amr_card`, `n_amr_amrfinder`, `n_amr_both` (konkordans), `amrfinder_organism`.

## 4. Config

`AMR`'e: `amrfinder_organism: str | None` (verilirse AMRFinderPlus koşar; ör. "Escherichia"),
`amrfinder_env: str = "ali-amrfinder"`. Diğer alanlar korunur.

## 5. Rapor

AMR tablosu sütunları: **Gen · CARD · AMRFinderPlus · %kimlik · log2FC · padj · DE** (yan yana). Yöntem
metnine AMRFinderPlus cümlesi + Kaynak (Feldgarden ve ark. 2021, AMRFinderPlus). DB tarih notu her iki araç
için. Virülans bölümü aynı.

## 6. Doğrulama

Birim (amrfinder parser: kolon/Type filtre/dict şekli; merge: birlik/yan-yana/konkordans; sütunlar);
GSE300731 canlı: CARD ~43, AMRFinderPlus ~5 (blaEC/acrF/mdtM); yan-yana tabloda ikisi görünür; efflux
genleri her ikisinde, blaEC yalnız AMRFinderPlus'ta. verdict SUSPECT değişmez.

## 7. İş akışı

spec → TDD (~5 task: config, parse_amrfinder+runner, merge, rapor, smoke) → GSE300731 → merge → DURUM/bellek.
