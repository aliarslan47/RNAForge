# DURUM — RNAForge

> Bu dosya "nerede kaldık" anlık görüntüsüdür. Tüm karar detayı Claude belleğindedir
> (`rnaforge-project` memory). Claude bunu anlamlı her durakta ve `/clear` öncesi günceller.

**Konum:** `/home/ali/rnaforge-pipeline/` (git deposu)
**GitHub:** `github.com/aliarslan47/RNAForge` — **PRIVATE**, remote `origin` (SSH), yalnız `main`
**Referans doküman:** `PLAN.md` **v1.3** (tek referans — Kural 1)
**Son güncelleme:** 2026-07-20

## Şu an nerede kaldık
- **Plan A BİTTİ ve `main`'e MERGE EDİLDİ** (`e3e68ec`). 52/52 test geçiyor.
  Final review'un 6 Important'ının **6'sı da düzeltildi** ve canlı doğrulandı (`654add4`).
- **HENÜZ PUSH EDİLMEDİ** — `main` yerelde origin'in önünde. Ali onayı bekliyor.
- `rnaforge validate --config X --metadata Y` çalışıyor. Resume artık GERÇEKTEN çalışıyor:
  aynı `--run-id` var olan dizini yeniden kullanır, biten modülü atlar ve **atladığını yazar**,
  `--force` ezer. Log satır satır flush edilir → koşu patlasa da neden diskte kalır.
- **Sıradaki iş: kalite kapıları çerçevesi.** Plan hazır ve onaylı:
  `docs/superpowers/plans/2026-07-20-kalite-kapilari-cerceve.md` (7 task, TDD adımlı).
  Tasarım: `docs/superpowers/specs/2026-07-20-kalite-kapilari-design.md`.
- **İlerleme ledger'ı: `.superpowers/sdd/progress.md`** (git-ignored — devam ederken ÖNCE oku).

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
