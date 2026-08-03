# DURUM — RNAForge

> Bu dosya "nerede kaldık" anlık görüntüsüdür. Tüm karar detayı Claude belleğindedir
> (`rnaforge-project` memory). Claude bunu anlamlı her durakta ve `/clear` öncesi günceller.

**Konum:** `/home/ali/rnaforge-pipeline/` (git deposu)
**GitHub:** `github.com/aliarslan47/RNAForge` — **PRIVATE**, remote `origin` (SSH)
**Referans doküman:** `PLAN.md` **v1.3** (tek referans — Kural 1)
**Son güncelleme:** 2026-08-03

## Şu an nerede kaldık
- **GERÇEK YAYIN-VERİSİ DOĞRULAMASI BİTTİ (2026-08-03) — açık konu KAPANDI.** Veri: **GSE300731**
  (Nature Microbiology 2025, "enterololin" dar-spektrum antibiyotik, Brown Lab). *E. coli* BW25113
  ΔbamBΔtolC, **5× MIC enterololin vs kontrol, 4h, 3'er replika** (6 örnek, PRJNA1281986). Referans
  BW25113 GCF_000750555.1. Girdiler gitignore'lu `raw/GSE300731/` + `references/ecoli_bw25113/`.
  Config: prok GFF → `feature_type=CDS`, `attribute=locus_tag`; `de.reference=control`.
- **Uçtan uca canlı koştu** (`runs/20260803_143036_GSE300731/`): hizalama %98.9–99.3, atama %85–86,
  replika kor. min 0.98, **DE 1634 anlamlı/4398 gen**. Tek WARN: ctrl_rep3 ham-QC GC → dürüstçe
  SUSPECT damgalandı (kapı sistemi çalışıyor). Kalan tüm kapılar PASS.
- **KONKORDANS (Katman A/B) — güçlü:** makalenin kendi kallisto abundance'ları (GEO suppl) → aynı
  DESeq2 → bizim ham-FASTQ (bowtie2+featureCounts) LFC'leriyle **Pearson r=0.972 / Spearman 0.958**
  (3592 ortak gen), **makale DEG recall %92.6**, **yön uyumu %99.9**. Biyolojik: top UP = zarf-stres
  imzası (pspA/pspC/spy + Rcs/kapsül rcsA/ugd/gmd/wzc) = makalenin doğruladığı **LolCDE** (lipoprotein
  taşıma) hedefiyle birebir; top DOWN = asit-direnç (gad*/hde*) baskılanması.
- **Analiz/figür script'leri scratchpad'de** (repo dışı): `run_all.sh`, `figures.R`, paper_de.
  Konkordans + figürler tek seferlik prototip — **ASIL İŞ: bunları pipeline'a m07/m08 olarak kur.**
- **YENİ İSTEK (kullanıcı, 2026-08-03):** figürler+rapor pipeline'ın kalıcı parçası olsun, **her koşuda
  OTOMATİK** çıksın; **güncel/modern, güzel, YÜKSEK ÇÖZÜNÜRLÜK** görseller. → m07 (figürler) + m08 (rapor)
  şimdi sıradaki iş; tasarım kalitesi çıtası yüksek. Plotlama: **rnaforge-de'de ggplot2 KURULU** (ggrepel
  YOK — ya kur ya kaçın). matplotlib hiçbir env'de yok. Konkordans figürü STANDART DEĞİL (referans DEG
  tablosu gerektirir; müşteri koşusunda olmaz) → yalnız doğrulama aracı, m07'ye girmez.
- **m07 SPEC + PLAN YAZILDI ve commit'lendi (2026-08-03), dal `feat/m07-figures`.** Kod HENÜZ YOK.
  - Spec: `docs/superpowers/specs/2026-08-03-m07-figures-design.md` · Plan: `docs/superpowers/plans/2026-08-03-m07-figures.md` (6 task, TDD).
  - Onaylı kararlar: figürler **PCA·Volcano·Heatmap·MA** (koşunun KENDİ verisinden); **statik yüksek çöz — PNG 300dpi + SVG**;
    R/ggplot2 `rnaforge-de` env (ggrepel+svglite EKLENECEK); m07 **gate yok/FAIL yok** (m06 gibi); ön koşul m06; çıktı `runs/.../figures/` + `manifest.json` (m08 tüketir).
  - **Konkordans/makale figürü m07'ye GİRMEZ** (Ali netleştirdi) — yalnız tek seferlik doğrulama aracıydı.
  - **DEVAM: Task 1'den başla** (`gene_name_map`, saf TDD) → executing-plans/inline. `conda run -n rnaforge-core --cwd <repo> python -m pytest -q`.
- (Önceki) **m06 DESeq2 `main`'de**, 181 test. PROKARYOT MVP DE zinciri tam.
- **AÇIK KONU (ARTIK KAPALI):** ~~GERÇEK yayımlanmış veri seti ile doğrulama henüz YOK.~~
  Testler + canlı smoke SENTETİK (Kural 8: gerçek/müşteri verisi repo'da yok — doğru). Ama
  Katman A/B doğruluğu için gerçek bakteri veri seti seçilmeli (bkz. "Açık konu" bölümü,
  demo veri kriterleri). SIRADAKİ GERÇEK İŞ olabilir.
- Önceki: m01, m02, m03, m04(prok), m05(prok) `main`'de.
- **Test komutu:** `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest -q`
  (repo dışından çağırırsan `tests.conftest` importu kırılır — yanlış alarm verir.)
- **m05 detay:** featureCounts TÜM BAM'lere tek çağrı → native matris; sütun→sample_id KONUMLA
  (BAM adıyla değil). Veri kapısı `assignment_rate` (ÜÇÜNCÜ FAIL kapısı). featureCounts params
  config-driven: `quantification.feature_type`(exon)/`attribute`(gene_id) — prok GFF3 için ez
  (CDS/locus_tag). `quantification` KNOWN_TOP_LEVEL_KEYS'e eklendi. Yanlış feature_type → yüksek
  sesle hata, exit 1 (sessiz boş matris YOK). Ön koşul m04; zincir m01→m03→m04→m05.
- **Çalıştırma NOTU:** `python -m rnaforge.cli` ÇALIŞMAZ (main-guard yok); entry point
  `rnaforge` kullan. Referans yolları config'te göreliyse CWD'ye bağlı — smoke'ta cd gerekti.

## Tamamlanan (kalite kapıları planı)
| Task | İş | Durum |
|---|---|---|
| 1 | `gates.py` sözleşmesi (PASS/WARN/FAIL, gates.json) | ✅ re-review temiz |
| 2 | Profiller veri olarak (`profiles/*.yml`) + eşik ezme | ✅ re-review temiz |
| 3 | `subject` sütunu + `looks_paired()` detektörü | ✅ review temiz |
| 4 | `validate_design` → `GateResult` döndürüyor (B kararı) | ✅ re-review temiz |
| 5 | m01 kapıları yazıp zorluyor | ✅ review temiz (mutasyon testli) |
| 6 | Güvence kartı + CLI verdict | ✅ controller + final review (m01 gates.json'ı enforce'tan ÖNCE yazıyor → FAIL'de kart INVALID, doğrulandı) |
| 7 | Dokümantasyon (README/PLAN v1.4) | ⏸ **ERTELENDİ** (m02-m03 sonrası) |

**2026-07-30 tamamlanan (bu oturum):**
- **config sertleştirme:** bilinmeyen üst-seviye anahtar (`design:` → `de.design`, yazım
  hataları) artık `ConfigError` ile reddediliyor (`KNOWN_TOP_LEVEL_KEYS`). Sessiz yutma bitti.
- **test temizliği:** `test_m01_validate.py`'de ikinci `_illumina` gölgeleme kaldırıldı,
  tek kanonik helper (200-read), gates importu tepeye taşındı.
- **Task 6 final review:** temiz. m02'ye geçiş serbest.

Reviewer'ların yakaladığı 3 gerçek bug (hepsi düzeltildi):
1. `gates.json` bozulunca önceki modüllerin kaydı **sessizce siliniyordu** (atomik yazma yok).
2. `~subject + condition` doygun tasarımda **üç kapı da yeşil yanıyordu** → DESeq2 kriptik
   çökme + "TRUSTWORTHY" damgası. (Spec'in deliğiydi, implementer hatası değil.)
3. Güvence kartı **FAIL'de hiç yazılmıyordu** — en çok gerektiği anda yok.

## SIRADAKİ İŞ
1. ~~config.py sessiz hata fix'i~~ ✅ BİTTİ (2026-07-30)
2. ~~Final whole-branch review + `_illumina` gölgeleme~~ ✅ BİTTİ (2026-07-30)
3. ~~`feat/kalite-kapilari` → `main` merge~~ ✅ BİTTİ (2026-07-30)
4. ~~m02 = FastQC~~ ✅ BİTTİ (2026-07-30) — `main`'de (merge `db4b501`).
5. ~~m03 = fastp (nazik trimming)~~ ✅ BİTTİ — `main`'de (`25b3ca2`).
6. ~~m04 = prokaryot hizalama (bowtie2)~~ ✅ BİTTİ — `main`'de (`7595d66`).
7. ~~m05 = prokaryot count matrisi (featureCounts)~~ ✅ BİTTİ — `main`'de (`caba117`).
8. ~~m06 = DESeq2~~ ✅ BİTTİ (2026-07-30) — spec+plan `.../2026-07-30-m06-deseq2*`; env `rnaforge-de`
   (bioconductor-deseq2 1.50.2), betik `rnaforge/scripts/deseq2.R` (dispersiyon fallback'li).
9. **`feat/m06-deseq2` → `main` merge + push. ← ŞİMDİ BURADAYIZ.**
10. **SIRADAKİ — iki yön:**
    a) **GERÇEK VERİ doğrulaması (Katman A/B):** yayımlanmış bakteri RNA-seq (count matrisi + ham
       FASTQ) ile pipeline'ı doğrula — demo veri seti seçimi (kriterler "Açık konu"da). Kullanıcı
       gerçek veri istedi; muhtemelen öncelik.
    b) **m07 figürler** (PCA/Volcano/Heatmap; `deseq2_results.tsv`+`normalized_counts.tsv` tüketir)
       → **m08 HTML rapor** → MVP tamam. VE/VEYA ökaryot yolu (m04-euk salmon + m05-euk tximport).
    - m06 = ilk ORTAK (organizma-agnostik) adım; asla FAIL üretmez (replicate_correlation WARN).

## Kalite kapıları — Ali ile onaylanan kararlar (2026-07-20)
Gerekçe: *"Yalancı sonuç asla istemem, müşteri güvenceli alsın"* + *"Sorun varsa sorun var
densin, patladıysa bilelim."* Pipeline doğru olsa bile kötü girdi MAKUL GÖRÜNEN sahte sonuç
üretir (%8 hizalama da count matrisi + p-değeri + rapor üretir).
- **Hedef kitle: D (karma)** — bakteriyel derin, altyapı genele açık.
- **İkili politika:** FAIL = sonuç GEÇERSİZ (durur, biyolojik çıktı YOK) · WARN = ŞÜPHELİ
  (üretilir + damgalanır).
- **Eşikler veri:** `profiles/{prokaryote,eukaryote}.yml`, `organism_type` seçer.
  Ökaryot BİLİNÇLİ gevşek + "geniş toleranslı" damgası (ökaryot doğrulaması YOK; uydurma
  eşik kapı sistemini itibarsızlaştırır). Ezilen eşik rapora YAZILIR — sessiz gevşetme yok.
- **FAIL çıktısı:** teşhis raporu (hangi kapı, ölçüm, eşik, sorumlu örnek, ne yapılmalı).
  "Damgalı ama üretilmiş sonuç" REDDEDİLDİ (grafik kopyalanır, damga kaybolur).
- **Eşleşmiş tasarım (öncesi/sonrası):** metadata'ya `subject` sütunu. Tespit VAR, karar YOK —
  eşleşmiş görünüp design kullanmıyorsa DUR ve sor; `paired: false` ile bilerek geçilir.
- **Kapsam:** çerçeve + tasarım kapıları (m01, metadata'dan çıkar) ŞİMDİ; veri kapıları
  kendi modülüyle (m04 yazılırken `alignment_rate` da yazılır).

## Ortam
- 54 GB RAM · 16 çekirdek · R sistemde kurulu (`/usr/bin/Rscript`)
- `rnaforge-core` (python 3.11 + pyyaml + pytest, `pip install -e .`)
- **Araç env'leri KURULDU (2026-07-20):** `rnaforge-qc` (FastQC 0.12.1, fastp 1.3.6) ·
  `rnaforge-quant-prok` (bowtie2 2.5.5, samtools, featureCounts 2.1.1) ·
  `rnaforge-quant-euk` (salmon 2.3.4)
- **DİKKAT:** salmon **2.3.4** geldi, PLAN 1.x varsayıyordu. m04 yazılmadan ÖNCE
  index/CLI davranışı doğrulanmalı — körlemesine güvenme.

## Onaylanan kararlar
1. **Orkestrasyon:** Python paketi, `ali-wgs-pipeline` deseni (config.yaml, runs/, conda env'ler).
2. **Yönlendirme:** `organism_type` (prokaryote|eukaryote) ZORUNLU, varsayılanı yok.
   - prokaryote → bowtie2 + featureCounts · eukaryote → Salmon + tximport/tx2gene
   - Ayrım YALNIZCA `m04`/`m05`'te; ikisi de **gen × örnek count matrisi** sözleşmesinde
     buluşur. `m06`/`m07`/`m08` organizma tipini bilmez.
3. **Platform:** FASTQ'dan tespit (Illumina/ONT/PacBio). MVP **Illumina-only**;
   ONT/PacBio tespit edilir + net hatayla REDDEDİLİR. Kütüphane kimyası
   (stranded/rRNA-polyA) FASTQ'dan tespit EDİLEMEZ → config'ten.
4. **Trimming:** fastp NAZİK — adapter + zorunlu min-uzunluk filtresi;
   agresif kalite trimming KAPALI. Gerekçe: Williams et al. 2016 (BMC Bioinformatics)
   agresif trimming genlerin >%10'unun ekspresyon tahminini bozuyor.
5. **DE motoru:** R/Bioconductor DESeq2 birincil + **pydeseq2 çapraz kontrol**
   (yalnız doğrulamada, her run'da DEĞİL).
6. **Doğrulama iki katmanlı:** Katman A (yayımlanmış count matrisi → DE → birebir yakın
   beklenir) · Katman B (ham FASTQ → tüm pipeline → yayınla KONKORDANS; birebir değil).
7. **Dosyalama:** deney-bazlı `runs/<ts>_<id>/{raw_qc,quantification,...}` (PLAN §14).
   WGS'deki örnek-bazlı `<id>/analiz/` düzeninden BİLİNÇLİ ayrılış — DE/PCA/volcano
   tek örneğe değil deneye ait.
8. **Dil:** kod/log EN · `README.md` (EN) + `README.tr.md` (TR) · PLAN/DURUM/spec TR ·
   HTML rapor dili config'ten (`tr`|`en`).

## MVP kapsamı (PLAN §2.1)
VAR: doğrulama+platform tespiti → FastQC → fastp → quant (2 yol) → count matrisi →
DESeq2 → PCA/Volcano/Heatmap → HTML rapor.
YOK (Faz 2+): ONT, MultiQC, GO/KEGG, GSEA, workflow diyagramı, dashboard, PDF.

## Açık konu
- **Demo veri seti SEÇİLMEDİ.** Kriterler PLAN §3.2'de yazılı: bakteriyel, public
  Illumina FASTQ, üst segment dergi (Frontiers/PLOS One/Sci Rep HARİÇ), 2024-2026
  (tercih 2025-26), GÜÇLÜ sinyal, temiz tasarım, yayında açık DEG tablosu.
- Literatür taraması yapıldı, henüz tüm kriterleri sağlayan net kazanan yok.
  En güçlü aday sinyal açısından: S. coelicolor Scr1 overekspresyonu (1308 DEG) ama
  dergisi (Microbial Cell Factories) üst segment değil.
- MVP kodunu BLOKE ETMEZ (pipeline organizma-agnostik) — Katman B doğrulamasından
  önce seçilmeli.

## Ortam
- 54 GB RAM · 16 çekirdek · 654 GB boş · R sistemde kurulu (`/usr/bin/Rscript`)
- Conda env **`rnaforge-core` KURULU** (python 3.11 + pyyaml + pytest; `pip install -e .` yapıldı).
  Test: `conda run -n rnaforge-core python -m pytest -q`
- Plan B araç env'leri YOK: **fastqc, fastp, salmon, bowtie2 KURULU DEĞİL**
  (featureCounts ve samtools sistemde var).

## Kritik kurallar
- Müşteri verisi/PII asla commit edilmez (`runs/`, `raw/`, `references/` .gitignore'da).
- PLAN.md yeniden yazılmaz, yalnız revize + sürüm yükseltilir (Kural 3).
- Tespit etmek ≠ desteklemek: desteklenmeyen girdi sessizce işlenmez (Kural 7).
