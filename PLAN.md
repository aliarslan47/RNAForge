# PLAN.md

**Version:** v1.1\
**Status:** Active\
**Purpose:** Single Source of Truth

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

-   **MVP (Faz 1 hedefi):** FASTQ → QC → Salmon → DESeq2 → 3-4 temel figür
    (PCA, Volcano, Heatmap) → basit HTML rapor. Demo veride uçtan uca **koşmalı.**
-   **Faz 2+:** GSEA, KEGG, 11 figür/9 tablo tam seti, otomatik workflow diyagramı,
    9 modüllü dashboard, PDF rapor.
-   Kural: "Çalışan MVP" > "mükemmel plan". Her modül önce en basit çalışan
    haliyle bitirilir, sonra zenginleştirilir.

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

> Not: Yayıncı markası (Nature/Springer vb.) bir kriter değildir. Önemli olan,
> beklenen sonucun bilinir olması ve pipeline çıktısıyla karşılaştırılabilmesidir.
> Kriteri gereksiz daraltmak zaman kaybettirir.

Bu veri seti yalnızca geliştirme ve test amacıyla kullanılacaktır.

Gerçek kullanımda tüm analizler yalnızca müşterinin sağladığı veri
üzerinden gerçekleştirilecek ve oluşturulan tüm çıktılar yalnızca o
çalışmaya ait olacaktır.

------------------------------------------------------------------------

# 4. Pipeline Workflow

``` text
Input FASTQ
↓
Metadata Validation
↓
Quality Control (FastQC)
↓
Read Trimming (fastp)
↓
Quantification (Salmon, transcriptome index)
↓
tximport + tx2gene (transcript → gene level)
↓
Count Matrix
↓
Differential Expression (DESeq2)
↓
Functional Enrichment (GO/KEGG = ORA, GSEA)
↓
Visualization
↓
MultiQC (tüm adımların QC'sini toplar)
↓
Interactive Dashboard
↓
Scientific Report
```

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

------------------------------------------------------------------------

# 5. Pipeline Modules

-   Input Validation
-   Quality Control
-   Read Trimming
-   Quantification (Salmon)
-   Transcript-to-Gene (tximport + tx2gene)
-   Differential Expression (esnek design formülü, batch desteği)
-   Functional Enrichment (ORA: GO/KEGG — GSEA ayrı)
-   Visualization
-   MultiQC Aggregation
-   Reporting
-   Dashboard

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
-   Reference **Transcriptome** (Salmon index)
-   Annotation / **tx2gene** mapping
-   Aligner/Quantifier (default: Salmon; opsiyonel: STAR — bulut)
-   **Design Formula** (ör. `~condition` veya `~batch + condition`)
-   **Batch / Covariates** (opsiyonel)
-   FDR Threshold
-   Log2FC Threshold
-   Thread
-   Memory

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
results/
├── raw_qc/
├── quantification/
├── differential_expression/
├── enrichment/
├── figures/
├── tables/
├── statistics/
├── workflow/
├── report/
├── dashboard/
└── logs/
```

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
-   GitHub'a uygun proje yapısı

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
