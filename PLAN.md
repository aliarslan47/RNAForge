# PLAN.md

**Version:** v1.4\
**Status:** Active\
**Purpose:** Single Source of Truth

> **Changelog (v1.3 → v1.4)**
> - **Uzun-okuma (ONT/PacBio) kolu TAMAMLANDI** (7 adım, hepsi TDD+merge). ONT/PacBio artık
>   REDDEDİLMİYOR, **yönlendiriliyor**: yeni **read_type (short|long)** boyutu m02'den itibaren
>   akışı dallandırır. Kısa: FastQC→fastp→Bowtie2. Uzun: NanoPlot→Pychopper+chopper→minimap2
>   (`-ax map-ont`/`map-hifi`)→featureCounts `-L`. İki boyut (`organism_type × read_type`) aynı
>   gen×örnek count matrisinde buluşur → m06+ değişmez. Yeni `library.chemistry` config
>   (cdna|direct_rna; ONT-long için zorunlu, FASTQ'tan tespit edilemez). Yeni **`prokaryote_long`**
>   kalite profili (bilinçli permissive + damgalı; alignment 0.50 FAIL, survival/assignment WARN).
>   Rapor read_type-farkında (rozet + long araç/yöntem/atıflar). Bölüm 2.1, 4, 5, 11.
> - **Downstream fonksiyonel analiz katmanı** (m09–m15) eklendi: GO ORA, KEGG ORA, GSEA (fgsea),
>   REVIGO-benzeri semantik indirgeme, AMR/virülans (CARD+AMRFinderPlus/VFDB), operon, STRING PPI+
>   community. Hepsi organizma/okuma-tipi agnostik, yeni FAIL kapısı YOK.
> - **QC tamamlama** (m16 seqqc rRNA%/strandedness · m17 alignqc insert-size/coverage/read-dist ·
>   m18 MultiQC) + TPM/FPKM + yazılım/veritabanı sürüm tabloları. Tümü tanısal.
> - **m00 basecall (ham sinyal → FASTQ):** FAST5/POD5 girdisi artık desteklenir — `rnaforge basecall`
>   dorado (GPU, `hac` model) ile FASTQ üretir, m01 çözülmüş metadata'yı devralır; FASTQ girdisinde
>   atlanır. GPU zorunlu (RTX 4050 ile doğrulandı: POD5→FASTQ canlı). Yeni `rnaforge-basecall` env
>   (pod5) + ayrı dorado binary. Bölüm 4.
> - **Doğrulama:** kısa-okuma GSE300731 uçtan uca (konkordans r=0.972); uzun-okuma microbepore
>   canlı smoke (araç zinciri doğrulandı; DE-sinyalli bakteri ONT seti henüz seçilmedi).
>
> **Changelog (v1.2 → v1.3)**
> - Trimming varsayılanı testi netleştirildi: `min_length >= 1` iddiası **vacuous**'tu
>   (config zaten `< 1`'de `ConfigError` atıyor, assertion asla düşemezdi). Williams 2016
>   gerekçesini gerçekten koruyan iddia **`min_length == 36`**'dır (Bölüm 4.2).
> - Kapatma dayanıklılığı iddiası koda hizalandı: resume yalnız kayıt değil, **çalışan**
>   davranıştır — aynı `run_id` var olan run dizinini yeniden kullanır, biten modül
>   atlanır (`--force` ile ezilebilir), heartbeat 10 sn kadansını gerçekten uygular
>   (Bölüm 15).
>
> **Changelog (v1.1 → v1.2)**
> - Organizma tipi girişe alındı: prokaryot/ökaryot seçimi akışı dallandırır
>   (Bölüm 2.1, 4, 5, 11). Bakterilerde intron yok → tx2gene ~1:1; Salmon+tximport
>   yerine genom hizalama + featureCounts yolu.
> - Platform tespiti eklendi: Illumina/ONT/PacBio FASTQ'dan tespit edilir.
>   MVP yalnızca Illumina destekler; ONT/PacBio tespit edilip reddedilir (Bölüm 2.1, 4.1, 11).
> - Trimming politikası eklendi (literatür temelli): agresif kalite trimming
>   ekspresyon tahminlerini bozar → nazik varsayılan + zorunlu min-uzunluk filtresi
>   (Bölüm 4.1, 11).
> - Doğrulama iki katmanlı tanımlandı: Katman A (count matrisi → birebir yakın),
>   Katman B (FASTQ → konkordans). Birebir yayın tekrarı hedef değildir (Bölüm 3).
> - Çıktı yapısının adım-bazlı olduğu gerekçesiyle netleştirildi (Bölüm 14).
> - Çift dilli (TR/EN) dokümantasyon ve rapor dili config'e bağlandı (Bölüm 11, 15).
> - Detaylı MVP tasarımı: `docs/superpowers/specs/2026-07-16-rnaforge-mvp-design.md`

> **Changelog (v1.0 → v1.1)**
> - MVP-First build önceliği eklendi (Bölüm 2.1).
> - Development Policy'deki doğrulama verisi kısıtı gevşetildi (yayıncı markası
>   yerine "known published ground-truth") (Bölüm 3).
> - Workflow düzeltildi: MultiQC, quantification'dan **sonraya** alındı;
>   tximport/tx2gene adımı eklendi (Bölüm 4).
> - Method Notes bölümü eklendi: Salmon = pseudo-alignment, ORA vs GSEA,
>   DESeq2 design formülü (Bölüm 4.1).
> - Configuration, transkriptom/tx2gene ve design formülü/batch içerecek
>   şekilde güncellendi (Bölüm 11).
> - Yeni bölüm: Data Security & Retention (Bölüm 16).

------------------------------------------------------------------------

# 1. Project Overview

Bu proje; kullanıcı tarafından sağlanan Bulk RNA-seq verilerini uçtan
uca analiz eden, yeniden üretilebilir (reproducible), modüler ve
otomatik raporlama yapabilen profesyonel bir analiz pipeline'ı
geliştirmeyi amaçlamaktadır.

Pipeline'ın temel çıktısı; bilimsel standartlara uygun analiz sonuçları,
yayın kalitesinde görseller, kapsamlı bir rapor ve interaktif bir
dashboard olacaktır.

------------------------------------------------------------------------

# 2. Objectives

-   Modüler Bulk RNA-seq pipeline geliştirmek.
-   Analiz sürecini tamamen otomatik hale getirmek.
-   Publication-quality rapor üretmek.
-   İnteraktif dashboard oluşturmak.
-   Tekrar üretilebilir analiz sağlamak.

## 2.1 Build Priority (MVP-First)

Bu belge kapsamlıdır; tehlike, ilk müşteriden önce haftalarca cilalama
tuzağına düşmektir. Bu yüzden geliştirme **çalışan ince bir dilimle** başlar,
süsler sonra eklenir.

-   **MVP (Faz 1 hedefi):** FASTQ → doğrulama + **platform tespiti** → QC →
    nazik trimming → **kantifikasyon (prokaryot | ökaryot yolu)** → gen-seviyesi
    count matrisi → DESeq2 → 3 temel figür (PCA, Volcano, Heatmap) → basit HTML
    rapor. Demo veride uçtan uca **koşmalı.**
-   **MVP sınırları (TARİHSEL):** MVP yalnız **Illumina** kısa-okumaydı; ONT/PacBio o aşamada
    reddediliyordu. **v1.4'te bu değişti:** uzun-okuma (ONT/PacBio) kolu tamamlandı — artık
    reddedilmiyor, read_type ile yönlendiriliyor (yukarıdaki changelog). Yalnız **tanımlanamayan**
    platform reddedilir.
-   **Faz 2+ (ÇOĞU TAMAMLANDI, v1.4):** ✅ ONT yolu (minimap2), ✅ MultiQC, ✅ GSEA, ✅ KEGG/GO (ORA),
    ✅ genişletilmiş figür/tablo seti + downstream (AMR/operon/PPI/REVIGO). **KALAN:** ökaryot yolu
    canlı doğrulaması (m04-euk salmon + m05-euk tximport), otomatik workflow diyagramı, interaktif
    dashboard, PDF rapor.
-   Kural: "Çalışan MVP" > "mükemmel plan". Her modül önce en basit çalışan
    haliyle bitirilir, sonra zenginleştirilir.

### Organizma tipi (girişte seçilir)

Pipeline organizma-agnostiktir; müşteri işi herhangi bir organizma olabilir.
Organizma tipi **config'ten gelir ve zorunludur** (varsayılanı yoktur — yanlış
varsayım sessiz hataya yol açar):

-   `prokaryote`: genom hizalama (bowtie2) + featureCounts
-   `eukaryote`: Salmon + tximport/tx2gene

Bu ayrım yalnızca kantifikasyon ve count modüllerine hapsedilir; her iki yol da
aynı **gen × örnek count matrisi** sözleşmesinde buluşur. DE, figürler ve rapor
matrisin hangi yoldan geldiğini bilmez.

------------------------------------------------------------------------

# 3. Development Policy

Pipeline yalnızca bir kez doğrulama amacıyla public bir veri seti
üzerinde test edilecektir.

Doğrulama veri seti;

-   GEO, SRA veya ArrayExpress kaynaklı,
-   **açık erişimli (public),**
-   **ve sonuçları yayımlanmış/bilinen (known published ground-truth)** —
    yani kendi çıktımızı karşılaştırıp doğrulayabileceğimiz bir referans sonucu olan

olmalıdır.

> Not: Önemli olan, beklenen sonucun bilinir olması ve pipeline çıktısıyla
> karşılaştırılabilmesidir. Buna karşılık veri seti **üst segment bir dergiden**
> gelmelidir (Nature Communications, Nature Microbiology, Science, Cell, PNAS,
> EMBO, mBio, NAR vb.); Frontiers, PLOS One ve Scientific Reports dahil değildir.
> Tarih aralığı 2024–2026 (tercih 2025–2026).

## 3.1 Doğrulama Katmanları

Bir yayının sayılarını **birebir** yeniden üretmek genelde mümkün değildir: yayınlar
farklı araç, farklı anotasyon sürümü ve farklı filtreler kullanır. Gerçekçi ve
savunulabilir hedef **konkordanstır**. Bu ayrım yapılmazsa var olmayan bir hata aranır.

-   **Katman A — sayısal doğruluk:** Yayımlanmış count matrisi → DE modülü →
    yayımlanmış DE sonucuyla karşılaştır. Kantifikasyon değişkenliği devre dışı;
    yalnızca DE modülü sınanır → **birebir yakın eşleşme beklenir.** Tutmuyorsa
    DE modülü bozuktur. **pydeseq2 çapraz kontrolü burada çalışır:** R/DESeq2
    birincil sonuç, pydeseq2 bağımsız ikinci uygulama; anlamlı sapma bir alarmdır.
-   **Katman B — uçtan uca konkordans:** Ham FASTQ → tüm pipeline → yayının DEG
    tablosu. Hedef: anlamlı genlerin örtüşmesi ve etki yönü uyumu. Birebir log2FC
    eşleşmesi **beklenmez.**

pydeseq2 çapraz kontrolü yalnızca doğrulamada çalışır, her müşteri run'ında değil.

## 3.2 Demo Veri Seti Seçim Kriterleri

1.  Bakteriyel (prokaryot yolunu sınar)
2.  Ham FASTQ public (SRA/GEO/ENA), Illumina
3.  Üst segment dergi (yukarıdaki liste)
4.  2024–2026 (tercih 2025–2026)
5.  **Güçlü ve net sinyal** — knockout/overekspresyon veya güçlü muamele; yüzlerce
    DEG. Gerekçe: zayıf sinyalli sette doğru kurulmuş bir pipeline bile yayından
    sapar ve "sapma bizden mi veriden mi" ayırt edilemez.
6.  Temiz tasarım: az değişken, replikalı
7.  Yayında karşılaştırılabilir açık sonuç (DEG tablosu / supplementary)

Bu veri seti yalnızca geliştirme ve test amacıyla kullanılacaktır.

Gerçek kullanımda tüm analizler yalnızca müşterinin sağladığı veri
üzerinden gerçekleştirilecek ve oluşturulan tüm çıktılar yalnızca o
çalışmaya ait olacaktır.

------------------------------------------------------------------------

# 4. Pipeline Workflow

``` text
Input FASTQ
↓
Metadata Validation + Platform & read_type Detection (Illumina→short | ONT/PacBio→long)
│   └─ tanımlanamayan platform → REJECT (yalnız unknown reddedilir; v1.4)
↓
Quality Control (short: FastQC · long: NanoPlot)
↓
Read Trimming (short: fastp nazik · long: Pychopper+chopper / chopper)
↓
┌─────────────── organism_type ───────────────┐
│                                             │
prokaryote                                eukaryote
│                                             │
Alignment (short: bowtie2 · long: minimap2) Quantification (Salmon, transcriptome index)
│                                             │
featureCounts (GFF/GTF; long: -L)         tximport + tx2gene (transcript → gene)
│                                             │
└─────────────────┬───────────────────────────┘
↓
Count Matrix (gen × örnek — ortak sözleşme)
↓
Differential Expression (DESeq2, esnek design formülü)
↓
Functional Enrichment (GO/KEGG ORA, GSEA, REVIGO, AMR/operon/PPI)  [✅ v1.4]
↓
Visualization (PCA, Volcano, Heatmap, MA, dispersion, …)
↓
MultiQC (tüm adımların QC'sini toplar)             [✅ v1.4]
↓
Interactive Dashboard                              [Faz 2+]
↓
Scientific Report (MVP: HTML; PDF Faz 2+)
```

Akış `organism_type`'a göre dallanır ve **count matrisi sözleşmesinde birleşir.**
Dallanma yalnızca kantifikasyon/count modüllerindedir; sonraki tüm adımlar ortaktır.

Pipeline tamamlandıktan sonra bu akış ayrıca SVG, PNG ve PDF
formatlarında otomatik bir Workflow Diagram olarak oluşturulmalıdır.

## 4.1 Method Notes (kritik teknik notlar)

-   **MultiQC yeri:** MultiQC bir toplayıcıdır; FastQC + fastp + Salmon
    loglarını **tek raporda birleştirebilmesi için quantification'dan sonra**
    çalışır. (v1.0'da yanlış konumdaydı.)
-   **Salmon = alignment DEĞİL:** Salmon *pseudo-alignment* ile kantifikasyon
    yapar ve **genom değil transkriptom** indeksi kullanır. Gen seviyesine çıkmak
    için **tximport + tx2gene** eşlemesi zorunludur.
-   **STAR opsiyonu:** STAR + featureCounts alternatif akıştır ama insan genomunda
    ~30GB RAM ister → WSL2/laptop için ağır. Bulut opsiyonu olarak tutulur;
    varsayılan Salmon'dur.
-   **DESeq2 design:** Sabit `~condition` yerine esnek formül; batch/confounder
    varsa `~batch + condition`. Design, config'ten gelmelidir.
-   **ORA vs GSEA ayrımı:** GO/KEGG **ORA** anlamlı DEG *alt kümesini* kullanır;
    **GSEA** *tüm genlerin sıralı (ranked) listesini* kullanır. İkisi karıştırılmaz.
-   **Prokaryot ≠ ökaryot:** Bakterilerde intron yoktur ve gen başına tek transkript
    vardır → tx2gene ~1:1, tximport bir formaliteye döner. Prokaryotta doğru yol
    genom hizalama + featureCounts'tur. Ayrıca rRNA deplesyonu polyA seçiminin
    yerini alır; GO/KEGG anotasyonu farklı kaynaklardan gelir (org.Hs.eg.db yok →
    eggNOG/KEGG gerekir, Faz 2+).
-   **Platform tespiti ve sınırı:** FASTQ başlık formatı + read uzunluk dağılımından
    platform **güvenilir tespit edilir** (Illumina / ONT / PacBio). Buna karşılık
    kütüphane kimyası (rRNA-deplesyon mu polyA mı, stranded mı) FASTQ'da **yoktur,
    tespit edilemez** → config'ten gelmelidir. Strandedness hizalama sonrası
    çıkarsanabilir (Salmon `-l A`) ama tahmin edilmez; config ile çelişirse uyarılır.

## 4.2 Trimming Politikası (literatür temelli)

**Agresif kalite trimming RNA-seq ekspresyon tahminlerini bozar.** Williams et al.
2016 (BMC Bioinformatics, doi:10.1186/s12859-016-0956-2) agresif trimming ile
genlerin **%10'undan fazlasının** ekspresyon tahmininin anlamlı şekilde değiştiğini
gösterir; üç veri setinde ve farklı DE pipeline'larında tekrarlanmış, sonuçlar
microarray'e karşı doğrulanmıştır. Sebep: kısalan read'lerin yanlış hizalanması
(spurious mapping). Sapmanın büyük kısmı **minimum uzunluk filtresiyle** ortadan
kalkar. Yazarların sonucu: *no or modest trimming results in the most biologically
accurate gene expression estimates.*

Ayrıca Salmon gibi araçlar soft-clipping yaptığı için adapter trimming'in faydası
yalnızca gerçek adapter kontaminasyonu varsa ortaya çıkar.

**Karar:** fastp kullanılır, ama nazik varsayılanlarla:

-   Adapter tespiti/temizliği: **açık** (fastp otomatik)
-   Agresif kalite trimming (sliding-window `--cut_right` vb.): **kapalı**
-   Minimum uzunluk filtresi: **açık ve zorunlu**
-   Kalite filtresi: fastp varsayılanı (nazik)

Eşikler config'ten ayarlanabilir; güvenli varsayılan daima nazik taraftadır.
Gerekçe doğrulama hedefini de korur: agresif trimming ile koşulsaydı sonuçlar
referans yayından sapar ve bu sapma pipeline hatası sanılırdı.

------------------------------------------------------------------------

# 5. Pipeline Modules

-   `m01` Input Validation + **Platform Detection**
-   `m02` Quality Control (FastQC)
-   `m03` Read Trimming (fastp — nazik, bkz. Bölüm 4.2)
-   `m04` Quantification **ROUTER** — prokaryote: bowtie2 | eukaryote: Salmon
-   `m05` Count Matrix — prokaryote: featureCounts | eukaryote: tximport + tx2gene
-   `m06` Differential Expression (DESeq2, esnek design formülü, batch desteği)
-   `m07` Visualization (MVP: PCA, Volcano, Heatmap)
-   `m08` Reporting (MVP: HTML)
-   Functional Enrichment (ORA: GO/KEGG — GSEA ayrı)   [Faz 2+]
-   MultiQC Aggregation                                 [Faz 2+]
-   Dashboard                                           [Faz 2+]

`m04`/`m05` dışındaki hiçbir modül organizma tipini bilmez.

------------------------------------------------------------------------

# 6. Statistics

## Raw Statistics

-   Sample sayısı
-   Grup dağılımı
-   FASTQ bilgileri
-   Read sayısı
-   Read uzunluğu
-   GC oranı
-   Kalite skorları

## Final Statistics

-   Mapping Rate
-   Assigned Reads
-   Gene Count
-   Significant DEGs
-   GO Summary
-   KEGG Summary
-   GSEA Summary

Başlangıç ve sonuç istatistikleri karşılaştırmalı olarak sunulmalıdır.

------------------------------------------------------------------------

# 7. Figures

-   Figure 1. Analysis Workflow
-   Figure 2. QC Summary
-   Figure 3. MultiQC
-   Figure 4. PCA
-   Figure 5. Heatmap
-   Figure 6. Volcano Plot
-   Figure 7. MA Plot
-   Figure 8. GO Enrichment
-   Figure 9. KEGG Enrichment
-   Figure 10. GSEA
-   Figure 11. Dashboard Overview

Her figür numaralandırılmalı ve açıklayıcı bir caption içermelidir.

> MVP kapsamı (bkz. 2.1): önce Figure 4 (PCA), 5 (Heatmap), 6 (Volcano).
> Geri kalanı Faz 2+.

------------------------------------------------------------------------

# 8. Tables

-   Table 1. Dataset Summary
-   Table 2. Sample Information
-   Table 3. Raw Statistics
-   Table 4. QC Summary
-   Table 5. Differential Expression Results
-   Table 6. GO Results
-   Table 7. KEGG Results
-   Table 8. GSEA Results
-   Table 9. Final Statistics

------------------------------------------------------------------------

# 9. Scientific Report

PDF ve HTML formatında otomatik rapor oluşturulmalıdır.
(MVP: önce HTML; PDF Faz 2+.)

Bölümler:

-   Dataset Information
-   Sample Information
-   Raw Statistics
-   Workflow
-   Methods
-   Results
-   Figures
-   Tables
-   Final Statistics
-   Discussion
-   Conclusion
-   References

------------------------------------------------------------------------

# 10. Dashboard

-   Gene Search
-   Sample Search
-   PCA Viewer
-   Volcano Explorer
-   Heatmap Viewer
-   GO Browser
-   KEGG Browser
-   GSEA Browser
-   Download Center

------------------------------------------------------------------------

# 11. Configuration

-   Organism
-   **Organism Type** — `prokaryote` | `eukaryote` — **ZORUNLU, varsayılanı yok**
-   **Platform** — `auto` (varsayılan) | `illumina`. `auto` tespit eder; açık değer
    verilirse tespit yine çalışır ve çelişkide hata verir.
-   Reference (organism_type'a göre zorunlu olur):
    -   prokaryote: **Genome FASTA** + **Annotation GFF/GTF**
    -   eukaryote: **Transcriptome FASTA** (Salmon index) + **tx2gene** mapping
-   **Library** — `strandedness` (unstranded | stranded | reverse),
    `selection` (rrna_depletion | polya). FASTQ'dan tespit edilemez (Bölüm 4.1).
-   **Trimming** — `min_length` (zorunlu filtre), `aggressive_quality: false`
    (literatür temelli varsayılan, bkz. Bölüm 4.2)
-   **Design Formula** (ör. `~condition` veya `~batch + condition`)
-   **Batch / Covariates** (opsiyonel)
-   FDR Threshold
-   Log2FC Threshold
-   **Report Language** — `tr` | `en`
-   Thread
-   Memory

Opsiyonel/Faz 2: STAR (bulut), ONT yolu (minimap2).

------------------------------------------------------------------------

# 12. Logging

-   validation.log
-   qc.log
-   quantification.log
-   deseq2.log
-   enrichment.log
-   report.log

------------------------------------------------------------------------

# 13. Error Handling

Kontrol edilmesi gereken durumlar:

-   Eksik FASTQ
-   Hatalı metadata
-   Eşleşmeyen sample
-   Eksik referans transkriptom / tx2gene
-   Yetersiz disk
-   Yetersiz RAM
-   Araç hataları

------------------------------------------------------------------------

# 14. Output Structure

``` text
runs/<timestamp>_<run_id>/
├── raw_qc/                    # örnek-bazlı: her FASTQ'nun FastQC çıktısı
├── quantification/            # örnek-bazlı quant + counts.tsv (birleşik matris)
├── differential_expression/   # deney-bazlı: DESeq2 sonuç tablosu
├── enrichment/                # [Faz 2+]
├── figures/                   # deney-bazlı: PCA, Volcano, Heatmap (≥300 DPI)
├── tables/
├── statistics/                # raw + final istatistikler (Bölüm 6)
├── workflow/                  # [Faz 2+]
├── report/                    # MVP: HTML
├── dashboard/                 # [Faz 2+]
└── logs/
```

**Bölme adım-bazlıdır, örnek-bazlı değil.** Gerekçe: DE, PCA ve volcano tek bir
örneğe değil *deneye* aittir — örnek-bazlı klasörde bu çıktıların konacağı yer
yoktur. Örnek kırılımı adım klasörlerinin *içinde* korunur (`raw_qc/`,
`quantification/`). Bu, örnek-bazlı WGS düzeninden bilinçli bir ayrılıştır:
WGS'de her örnek bağımsız analiz edilir, RNA-seq DE ise doğası gereği
örnekler-arası bir analizdir.

`runs/` .gitignore'dadır — müşteri verisi asla commit edilmez (Bölüm 16).

------------------------------------------------------------------------

# 15. Quality Standards

-   Reproducible
-   Docker & Conda destekli
-   Modüler mimari
-   Publication-quality görseller (≥300 DPI)
-   Otomatik workflow diyagramı
-   Otomatik rapor
-   Otomatik istatistik özeti
-   Numaralandırılmış tablo ve şekiller
-   GitHub'a uygun proje yapısı (**private** depo)
-   **Çift dilli (TR/EN) dokümantasyon:** `README.md` (EN) + `README.tr.md` (TR),
    içerik eşdeğer. Kod, değişken adları ve log mesajları İngilizce; PLAN/DURUM ve
    spec'ler Türkçe (çalışma dili). HTML rapor dili config'ten (`tr` | `en`).
-   **Kapatma dayanıklılığı:** uzun işler kapanmaya dayanıklı — 10 sn heartbeat +
    her modül bitiminde kalıcı durum kaydı → yeniden başlatmada resume.
    Somut sözleşme: aynı `--run-id` **var olan** run dizinini yeniden kullanır (yeni
    zaman damgalı dizin açmaz, aksi halde `state.json` erişilemez kalır ve resume
    sessizce çöker); tamamlanmış modül atlanır ve bu **ekrana yazılır** (`--force`
    bilerek yeniden koşturur); modül log'ları satır satır flush edilir, böylece koşu
    hata verip düştüğünde nedeni diskte kalır.

------------------------------------------------------------------------

# 16. Data Security & Retention

Müşteri verisi işleneceği için gizlilik ve saklama politikası açıktır
(bu aynı zamanda bir satış argümanıdır):

-   Müşteri verisi yalnızca güvenli kanaldan transfer edilir.
-   Analiz çıktıları yalnızca ilgili çalışmaya aittir; başka işte kullanılmaz.
-   İş teslim edilip onaylandıktan sonra ham veri, kararlaştırılan süre içinde
    güvenli şekilde silinir (retention süresi sözleşmede belirtilir).
-   KVKK / GDPR farkındalığı; talep halinde NDA imzalanır.
-   Ödeme ve iletişim yalnızca resmi platform üzerinden yürütülür.

------------------------------------------------------------------------

# 17. Claude Code Rules

1.  Bu belge tek referans dokümandır.
2.  Yeni gereksinimler mevcut belge üzerinde revize edilerek
    işlenecektir.
3.  PLAN.md yeniden oluşturulmayacak, yalnızca sürümü güncellenecektir.
4.  Proje yapısı, rapor standardı ve çıktı formatı korunacaktır.
5.  Geliştirme MVP-First ilkesine uyar (bkz. Bölüm 2.1): önce çalışan ince dilim,
    sonra zenginleştirme.
6.  Organizma tipi ayrımı `m04`/`m05`'e hapsedilir; diğer modüllere sızdırılmaz.
    Yeni bir yol (ONT, STAR) eklenirse yalnızca bu iki modül değişir.
7.  Tespit etmek ≠ desteklemek. Desteklenmeyen girdi sessizce işlenmez; net bir
    hatayla reddedilir.
8.  Müşteri verisi ve PII asla commit edilmez, teste veya dokümana yazılmaz.
