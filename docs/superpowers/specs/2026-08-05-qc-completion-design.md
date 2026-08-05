# QC Tamamlama — Tasarım (5 düşük-öncelik eksik)

**Tarih:** 2026-08-05
**Amaç:** RNAForge'un kalite/tanı katmanındaki 5 kalan eksiği kapatmak. Hepsi
*diagnostik* (verdict'i yalnız WARN ile etkiler, asla FAIL üretmez — m02/m16 ilkesi).

## Kapsam (sırayla)

1. **Per-base baz kompozisyonu + duplikasyon** — mevcut FastQC zip çıktısından
   (yeni araç yok). Per-base A/T/G/C içeriği + `Total Deduplicated Percentage`.
2. **Insert-size dağılımı** — `samtools stats` (BAM). Yalnız paired-end anlamlı;
   single-end'de zarifçe atlanır (skip + not, fail değil).
3. **RSeQC read-distribution** — `read_distribution.py` (BAM + BED). Okumaların
   genomik özelliklere (CDS exons, UTR, intron, intergenic) dağılımı.
4. **Coverage görselleştirme** — `samtools depth` → kontig başına ortalama
   coverage figürü (matplotlib).
5. **MultiQC toplu görünüm** — `multiqc` ile FastQC/fastp/featureCounts/samtools
   stats/RSeQC çıktılarını tek toplu görünümde birleştirir.

## Mimari — modül yerleşimi

- **F1 → m02 genişletir** (FastQC orada yaşıyor). `rnaforge/fastqc.py`'ye saf
  parser'lar; m02 stats JSON'a `per_base_composition` + `deduplication` alanları.
- **F2,F3,F4 → yeni modül `m17_alignqc`** (hizalama-sonrası BAM tanısı — üçü de
  sıralı BAM okur). `rnaforge/alignqc.py` (saf parser + runner'lar), figürler
  `figures/` altına, stats `statistics/alignqc_statistics.json`.
- **F5 → yeni modül `m18_multiqc`** (kapstone toplayıcı; en son koşar).
  `rnaforge/multiqc.py`. Çıktı `multiqc/` altında kendi HTML'i + general-stats
  TSV; rapora general-stats tablosu + göreli link gömülür.

Her modül mevcut deseni izler: `run_x(config, metadata_path, run_dir, force=False)`
+ RunState resume/heartbeat + önkoşul kontrolü (m04) + stats JSON + CLI subcommand
+ unit testler (saf parser) + rapor bölümü.

## Kapı politikası

- F1: yeni **duplikasyon WARN kapısı** — `Total Deduplicated Percentage` çok
  düşükse (yüksek duplikasyon) WARN. m02 sözleşmesi gereği asla FAIL.
- F2/F3/F4: **kapı yok** — saf diagnostik figür/tablo (raporda gösterilir).
- F5: kapı yok — toplayıcı görünüm.

## Rapor

`report_html.py` kalite bölümü (section_quality) genişler: per-base kompozisyon
figürü, duplikasyon satırı, insert-size figürü/satırı, read-distribution tablosu,
coverage figürü, MultiQC linki. Figürler "Şekil N", tablolar "Tablo N" olarak
mevcut post-pass ile numaralanır. `N_SECTIONS` etkilenmez (kalite bölümü içine).

## Config

Yeni top-level anahtar yok. Diagnostikler her zaman açık; girdi yoksa (SE insert-
size gibi) zarifçe atlanır. Duplikasyon WARN eşiği profile eklenir
(`dedup_fraction`, prokaryot varsayılan makul bir alt sınır).

## Test / doğrulama

- Saf parser'lar için hızlı, deterministik unit testler (örnek çıktı string'leri).
- Canlı smoke: mevcut tamamlanmış run `runs/20260805_160103_GSE300731_final`
  üzerinde her yeni modülü koştur (BAM/zip zaten var); PASS/FAIL gözle.
- Tüm test paketi yeşil kalmalı.

## Sıra ve commit

Tek dal `feat/qc-completion`. Her feature ayrı commit (F1→F5), sonra main'e merge,
rapor yeniden üretilir, DURUM.md + bellek güncellenir.
