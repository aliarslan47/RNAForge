# Uzun-Okuma (ONT/PacBio) Yolu — Tasarım

**Tarih:** 2026-08-05 (kararlar onaylandı; **uygulama 2026-08-06'da başlar**)
**Durum:** Tasarım onaylı, spec/plan + kod YARIN. Ali: "kaydet yarın başla".

## Amaç
Pipeline veriyi **kısa vs uzun okuma** olarak otomatik algılayıp yönlendirsin; her tür
kendi QC/ön-işleme/hizalama/sayım araçlarını kullansın, DESeq2 ve sonrası (m06–m18)
**değişmeden** ortak count matrisinde buluşsun. Şu an uzun okuma tespit edilip
REDDEDİLİYOR (`platform.py require_supported`) — bunu **yönlendirmeye** çeviriyoruz.

## Ali'nin onayladığı kararlar (2026-08-05)
- **Platform/kimya:** ONT cDNA + ONT direct-RNA + PacBio HiFi/Iso-Seq (üçü de).
- **Nicelik seviyesi:** **gen-seviyesi** (minimap2 → featureCounts → mevcut DESeq2 downstream).
  İzoform-seviyesi HAYIR — bakteride splicing yok, getiri düşük.
- **Sıra:** **prokaryot uzun-okuma ÖNCE**; ökaryot short-read (salmon, bekleyen) ve
  euk-uzun sonra. [[reminder_rnaforge_eukaryote]]

## Literatür temeli (2026-08-05 tarandı)
- Uzun-okuma literatürü ağırlıkla **ökaryot izoform/splicing** (IsoQuant Nat Biotech 2023;
  TALON; LRGASP Nat Methods 2024; SG-NEx Nat Methods 2025). **Bakteride uzun-okumanın değeri
  izoform değil**: tam-boy transkript, **operon sınırları**, TSS/TTS (NAR 2024 bakteri
  direct-RNA; *E. coli*/*S. aureus* operon sınırlarını genişletti).
- **LRGASP kilit bulgu:** iyi anotasyonlu genomda **referans-tabanlı** araçlar en iyi;
  uzun+doğru okuma keşif için, derinlik nicelik için önemli.
- **Nicelik:** NanoCount ONT'de Illumina ile **gen-seviyesi r=0.927**; featureCounts long-read
  `-L` moduyla gen-seviyesi sayım yapar (bizim seçtiğimiz yol).
- **Araçlar (yerleşik):** hizalama **minimap2** (`-ax map-ont`/`map-hifi`; bakteri splice yok);
  QC **NanoPlot**; cDNA ön-işleme **Pychopper** (tam-boy yönlendir/kes) + **chopper/NanoFilt**
  (uzunluk/kalite filtre); direct-RNA'da Pychopper yok.

## Mimari — yeni boyut: read_type (short|long)
Bugün ayrım tek boyutlu (organizma → m04/m05). Ekleniyor: **read_type**, m02'den itibaren.
`organism_type × read_type` → hepsi **gen×örnek count matrisi** sözleşmesinde buluşur → m06–m18 aynen.

| Aşama | short (mevcut) | long (yeni, prok-önce) |
|-------|----------------|------------------------|
| m02 QC | FastQC | NanoPlot (uzunluk/N50/kalite) |
| m03 ön-işleme | fastp | cDNA→Pychopper+chopper; direct-RNA→chopper; HiFi→filtre |
| m04 hizalama | bowtie2 (prok) | minimap2 (`-ax map-ont`/`map-hifi`) |
| m05 sayım | featureCounts | featureCounts `-L` (long-read modu) |

## Tespit (sağlamlaştırma)
`platform.py` zaten illumina/ont/pacbio_hifi ayırıyor → **red yerine route**. Uzunluk eşiğine
**uzunluk-dağılımı** (N50/medyan, %>500bp payı) eklenir (ONT kısa molekül de üretebilir).
cDNA vs direct-RNA FASTQ'dan tespit EDİLEMEZ → yeni config `library.chemistry`
(`cdna|direct_rna`; HiFi platformdan çıkar).

## Config + kapılar
- Yeni `library.chemistry` alanı (Library dataclass; `library` zaten KNOWN top-level).
- Long-read profili: farklı eşikler (ONT kalite ~Q10-15, Q30 değil) — read_type'a göre seçilir.
  Kapı ilkesi aynı: FAIL=geçersiz (durur), WARN=şüpheli (damgalı). [[rnaforge-project]]

## Yeni ortam
`rnaforge-longread`: minimap2 + NanoPlot + pychopper + chopper + samtools.

## Uygulama sırası (parça parça, her biri TDD + merge)
1. **Tespit+yönlendirme** (red→route) + `library.chemistry` config
2. m02-long: NanoPlot
3. m03-long: Pychopper (cDNA) + chopper filtre
4. m04-long: minimap2
5. m05-long: featureCounts `-L`
6. Long-read profili/kapıları
7. Rapor read_type rozeti + uçtan-uca canlı smoke

## Test + doğrulama
- Saf parser'lar TDD (NanoPlot/pychopper/minimap2 log/çıktıları).
- **Canlı smoke: gerçek bakteri ONT RNA-seq** veri seti seçilecek (ör. microbepore *E. coli*
  ONT dRNA/cDNA, SRA). Kısa-okuma GSE300731 gibi referans koşu.

## Riskler / açık noktalar
- Gerçek ONT doğrulama veri seti henüz SEÇİLMEDİ (yarın ilk iş adaylardan biri).
- featureCounts `-L` gen-seviyesi DE için yeterli mi → smoke'ta NanoCount ile çapraz kontrol opsiyonu.
- read_type × organism matrisi modül isimlendirmesi (m04-prok-long vs dispatcher) — spec'te netleşecek.
