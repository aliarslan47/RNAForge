# RNAForge MVP — Tasarım Dokümanı (Design Spec)

**Tarih:** 2026-07-16
**Durum:** Onay bekliyor
**Kaynak plan:** `PLAN.md` v1.1, Bölüm 2.1 (MVP-First)
**Kapsam:** Bulk RNA-seq pipeline, MVP ince dilimi

> Bu belge PLAN.md'nin yerini almaz. PLAN.md tek referans dokümandır (Kural 1).
> Bu spec, PLAN.md Bölüm 2.1'deki MVP kapsamının nasıl inşa edileceğini tanımlar.
> Bu spec'teki kararlar PLAN.md v1.2 revizyonunu gerektirir (bkz. Bölüm 9).

---

## 1. Amaç ve MVP sınırı

PLAN.md Bölüm 2.1 MVP'yi şöyle tanımlar: FASTQ → QC → kantifikasyon → DESeq2 →
3 temel figür (PCA, Volcano, Heatmap) → basit HTML rapor; demo veride uçtan uca koşmalı.

**MVP'de VAR:**
- Metadata + FASTQ doğrulama, platform tespiti
- QC (FastQC), nazik trimming (fastp)
- Kantifikasyon (prokaryot ve ökaryot yolları)
- Gen-seviyesi count matrisi
- Diferansiyel ekspresyon (DESeq2, esnek design formülü)
- PCA, Volcano, Heatmap (≥300 DPI)
- HTML rapor
- İki katmanlı doğrulama (Bölüm 6)

**MVP'de YOK (Faz 2+):**
- GSEA, GO/KEGG (ORA)
- 11 figür / 9 tablo tam seti
- Otomatik workflow diyagramı (SVG/PNG/PDF)
- 9 modüllü interaktif dashboard
- PDF rapor
- ONT / PacBio desteği
- MultiQC agregasyonu

---

## 2. Mimari

`ali-wgs-pipeline` deseni izlenir: Python paketi, YAML config, zaman damgalı run
dizinleri, araç grubu başına ayrı conda ortamı.

```
rnaforge-pipeline/
├── PLAN.md                  # tek referans doküman (v1.2)
├── README.md                # EN
├── README.tr.md             # TR
├── DURUM.md                 # "nerede kaldık" anlık görüntüsü
├── config/
│   └── config.yaml
├── rnaforge/
│   ├── __init__.py
│   ├── cli.py               # komut satırı girişi
│   ├── config.py            # config yükleme + şema doğrulama
│   ├── run.py               # orkestrasyon, resume, heartbeat
│   ├── contracts.py         # modüller arası veri sözleşmeleri
│   └── modules/
│       ├── m01_validate.py  # metadata + FASTQ doğrulama, platform tespiti
│       ├── m02_qc.py        # FastQC
│       ├── m03_trim.py      # fastp (nazik)
│       ├── m04_quant.py     # YÖNLENDİRİCİ: prokaryot | ökaryot
│       ├── m05_counts.py    # gen-seviyesi count matrisi
│       ├── m06_de.py        # DESeq2 (Rscript köprüsü)
│       ├── m07_figures.py   # PCA, Volcano, Heatmap
│       └── m08_report.py    # HTML rapor
├── r/
│   ├── tximport_deseq2.R    # ökaryot: tximport + DESeq2
│   └── deseq2_counts.R      # prokaryot: hazır count matrisi + DESeq2
├── envs/                    # conda ortam tanımları
├── tests/
├── docs/
└── runs/                    # .gitignore — müşteri verisi asla commit edilmez
```

### 2.1 Yönlendirme (routing) — tasarımın kilit fikri

Prokaryot ve ökaryot RNA-seq temelde farklıdır: bakterilerde intron yoktur, gen başına
tek transkript vardır (tx2gene ~1:1), rRNA deplesyonu polyA seçiminin yerini alır.

Bu fark **yalnızca `m04` ve `m05`'e hapsedilir**:

| Adım | Prokaryot | Ökaryot |
|---|---|---|
| Kantifikasyon (`m04`) | bowtie2 → genom BAM | Salmon → transkript quant |
| Gen seviyesi (`m05`) | featureCounts (GFF/GTF) | tximport + tx2gene |
| Referans | genom FASTA + GFF | transkriptom FASTA + tx2gene |

**Her iki yol da aynı sözleşmede buluşur: gen × örnek count matrisi.**

Sonuç: `m06` (DE), `m07` (figürler), `m08` (rapor) matrisin hangi yoldan geldiğini
bilmez ve umursamaz. Organizma tipi ayrımı pipeline'ın geri kalanına sızmaz.
Faz 2'de ONT yolu (minimap2) eklenirse yine sadece `m04`/`m05` değişir.

Yönlendirme `config.organism_type` (`prokaryote` | `eukaryote`) ile yapılır. Bu alan
zorunludur; varsayılanı yoktur — yanlış varsayım sessiz hataya yol açar.

### 2.2 Modül sözleşmeleri

Her modül şu sözleşmeye uyar: girdi olarak run dizini + config alır, çıktısını
diske yazar, bir durum kaydı bırakır. Modüller birbirini doğrudan çağırmaz;
`run.py` sırayı yönetir. Böylece her modül tek başına test edilebilir.

Kritik sözleşme — **count matrisi** (`m05` çıktısı, `m06` girdisi):
`quantification/counts.tsv`, satırlar gen ID, sütunlar örnek ID, değerler tamsayı
raw count. Örnek ID'leri metadata ile birebir eşleşir. Bu dosya `m04`'ün hangi
yolu izlediğine dair hiçbir iz taşımaz.

---

## 3. Platform tespiti ve sınırları

`m01` FASTQ başlık formatı + read uzunluk dağılımından platformu tespit eder.

**Güvenilir tespit edilir:**
- Illumina: sabit/dar read uzunluğu (~50–150 bp), tanınan başlık formatı
- ONT: geniş uzunluk dağılımı, uzun read'ler
- PacBio: uzun read'ler, PacBio başlık deseni

**Güvenilir tespit EDİLEMEZ** (FASTQ'da bu bilgi yoktur — config'ten gelmelidir):
- Kütüphane kimyası: rRNA deplesyonu mu polyA seçimi mi
- Strandedness (stranded / unstranded / reverse)

Strandedness hizalama sonrası çıkarsanabilir (Salmon `-l A`), ama tahmin edilmez;
config'ten alınır, tespit ile çelişirse uyarı verilir.

**MVP davranışı:** Yalnızca Illumina desteklenir. ONT/PacBio tespit edilir ve
**net bir hata mesajıyla reddedilir** — sessizce yanlış araçla koşulmaz.
Bu bilinçli bir sınırdır: tespit etmek, desteklemek değildir.

`config.platform` = `auto` (varsayılan) | `illumina`. `auto` tespit eder; açık değer
verilirse tespit yine çalışır ve çelişki varsa hata verir.

---

## 4. Trimming politikası (literatür temelli)

**Bulgu:** Agresif kalite trimming RNA-seq ekspresyon tahminlerini bozar.
Williams et al. 2016 (BMC Bioinformatics, doi:10.1186/s12859-016-0956-2) agresif
trimming ile genlerin %10'undan fazlasının ekspresyon tahmininin anlamlı şekilde
değiştiğini gösterir; üç veri setinde ve farklı DE pipeline'larında tekrarlanmış,
sonuçlar microarray'e karşı doğrulanmıştır. Sebep: kısalan read'lerin yanlış
hizalanması (spurious mapping). Sapmanın büyük kısmı **minimum uzunluk filtresiyle**
ortadan kalkar. Yazarların sonucu: *no or modest trimming results in the most
biologically accurate gene expression estimates.*

**Karar:** fastp kullanılır (PLAN §4), ama nazik varsayılanlarla:

| Ayar | Varsayılan | Gerekçe |
|---|---|---|
| Adapter tespiti/temizliği | açık | fastp otomatik tespit; kontaminasyon varsa gerekir |
| Agresif kalite trimming (sliding-window) | **kapalı** | literatür: doğruluğu bozar |
| Minimum uzunluk filtresi | **açık, zorunlu** | sapmayı engelleyen asıl mekanizma |
| Kalite filtresi | fastp varsayılanı (nazik) | agresif değil |

Eşikler config'ten ayarlanabilir; güvenli varsayılan daima nazik taraftadır.

Bu politika doğrulama hedefini de korur: agresif trimming ile koşulsaydı sonuçlar
referans yayından sapar ve bu sapma pipeline hatası sanılırdı.

---

## 5. Konfigürasyon

PLAN.md Bölüm 11 temel alınır, aşağıdaki alanlarla:

```yaml
organism: "Escherichia coli"
organism_type: "prokaryote"        # prokaryote | eukaryote — ZORUNLU, varsayılan yok
platform: "auto"                   # auto | illumina

reference:
  # prokaryote:
  genome_fasta: "..."
  annotation_gff: "..."
  # eukaryote:
  transcriptome_fasta: "..."
  tx2gene: "..."

library:
  strandedness: "unstranded"       # unstranded | stranded | reverse
  selection: "rrna_depletion"      # rrna_depletion | polya

trimming:
  min_length: 36
  aggressive_quality: false        # literatür temelli varsayılan (Bölüm 4)

de:
  design: "~condition"             # esnek formül; ör. "~batch + condition"
  fdr_threshold: 0.05
  log2fc_threshold: 1.0

report:
  language: "tr"                   # tr | en

resources:
  threads: 8
  memory_gb: 32
```

`organism_type`'a göre `reference` altındaki ilgili alanlar zorunlu olur; ilgisiz
alanlar yok sayılır. Doğrulama `m01`'den önce, config yüklenirken yapılır.

---

## 6. Doğrulama stratejisi (iki katmanlı)

Bir yayının sayılarını **birebir** yeniden üretmek genelde mümkün değildir: yayınlar
farklı araç, farklı anotasyon sürümü ve farklı filtreler kullanır. Gerçekçi ve
savunulabilir hedef konkordanstır. Bu ayrım yapılmazsa var olmayan bir hata aranır.

**Katman A — sayısal doğruluk (birebir yakın beklenir):**
Yayımlanmış count matrisi → `m06` → yayımlanmış DE sonuçlarıyla karşılaştır.
Kantifikasyon değişkenliği devre dışı; sadece DE modülü sınanır. Tutmuyorsa
`m06` bozuktur. pydeseq2 çapraz kontrolü burada çalışır: R/DESeq2 birincil sonuç,
pydeseq2 bağımsız ikinci uygulama; ikisi arasında anlamlı sapma bir alarmdır.

**Katman B — uçtan uca konkordans:**
Ham FASTQ → tüm pipeline → yayının DEG tablosu. Hedef: anlamlı genlerin örtüşmesi
ve etki yönü uyumu. Birebir log2FC eşleşmesi beklenmez.

pydeseq2 çapraz kontrolü **yalnızca doğrulamada** çalışır, her müşteri run'ında değil
(MVP-First: doğrulama aracı, üretim yükü değil).

### 6.1 Demo veri seti seçim kriterleri

Veri seti henüz seçilmedi. Seçim, Katman B doğrulaması başlamadan önce yapılır;
MVP kodunu bloke etmez (pipeline organizma-agnostiktir).

Zorunlu kriterler:
1. Bakteriyel (prokaryot yolunu sınar)
2. Ham FASTQ public (SRA/GEO/ENA)
3. Üst segment dergi: Nature Communications, Nature Microbiology, Science, Cell,
   PNAS, EMBO, mBio, Nucleic Acids Research vb. — Frontiers, PLOS One,
   Scientific Reports dahil DEĞİL
4. Yayın tarihi 2024–2026 (tercih 2025–2026)
5. **Güçlü ve net sinyal** — ör. knockout/overekspresyon veya güçlü muamele;
   yüzlerce-binlerce DEG. Gerekçe: zayıf sinyalli sette doğru kurulmuş bir pipeline
   bile yayından sapar ve "sapma bizden mi veriden mi" ayırt edilemez.
6. Temiz tasarım: az sayıda değişken, tercihen eşleştirilmiş/replikalı
7. Yayında karşılaştırılabilir açık sonuç (DEG tablosu / supplementary)

Illumina olmalıdır (MVP sınırı, Bölüm 3).

---

## 7. Çıktı yapısı

PLAN.md Bölüm 14 izlenir. Bölme **adım-bazlıdır**, örnek-bazlı değil — çünkü DE,
PCA ve volcano tek bir örneğe değil *deneye* aittir; örnek-bazlı klasörde bu
çıktıların yeri yoktur. Örnek kırılımı adım klasörlerinin *içinde* korunur.

```
runs/20260716_143022_<run_id>/
├── raw_qc/                    # örnek-bazlı: her FASTQ'nun FastQC çıktısı
├── quantification/            # örnek-bazlı quant + counts.tsv (birleşik matris)
├── differential_expression/   # deney-bazlı: DESeq2 sonuç tablosu
├── figures/                   # deney-bazlı: PCA, Volcano, Heatmap (≥300 DPI)
├── tables/
├── statistics/                # raw + final istatistikler (PLAN §6)
├── report/                    # HTML rapor
└── logs/                      # validation, qc, quantification, deseq2, report
```

Not: Bu, `ali-wgs-pipeline`'daki örnek-bazlı (`<id>/analiz/`) düzenden bilinçli
bir ayrılıştır. WGS'de her örnek bağımsız analiz edilir; RNA-seq DE ise doğası
gereği örnekler-arası bir analizdir.

---

## 8. Dayanıklılık, hata yönetimi, test

**Kapatma dayanıklılığı:** Uzun işler kapanmaya dayanıklı olmalıdır. `run.py`
10 saniyelik heartbeat yazar ve her modül bittiğinde kalıcı durum kaydı bırakır.
Yeniden başlatıldığında tamamlanmış modüller atlanır (resume).

**Hata yönetimi** (PLAN §13): eksik FASTQ, hatalı metadata, eşleşmeyen örnek,
eksik referans, yetersiz disk/RAM, araç hataları. Hepsi `m01`'de veya ilgili
modülün başında yakalanır ve *net mesajla* durur — sessiz devam yok.

**Test (TDD):** Her modül için önce test. Küçük sentetik FASTQ fixture'ları ile
hızlı birim testleri; gerçek veri yalnızca doğrulamada. Modül sözleşmeleri
(özellikle count matrisi) test edilir.

---

## 9. PLAN.md v1.2 revizyonu (gerekli)

Bu spec'teki kararlar PLAN.md v1.1'i aşar. Kural 2/3 gereği belge yeniden
yazılmaz, revize edilip sürümü yükseltilir:

- **Bölüm 2.1:** MVP kapsamına platform tespiti + prokaryot/ökaryot yönlendirmesi
- **Bölüm 4:** Akış tek hat değil; `m04`/`m05`'te dallanır (prok/ökaryot)
- **Bölüm 4.1:** Trimming politikası (literatür temelli, nazik) eklenir
- **Bölüm 5:** Modül listesi yönlendirmeyi yansıtacak şekilde güncellenir
- **Bölüm 11:** `organism_type`, `platform`, `library`, `trimming`, `report.language`
- **Bölüm 15:** Çift dilli (TR/EN) dokümantasyon

---

## 10. Dil politikası (TR/EN)

- Kod, değişken adları, log mesajları: İngilizce
- `README.md` İngilizce, `README.tr.md` Türkçe — içerik eşdeğer tutulur
- `PLAN.md`, `DURUM.md`, spec'ler: Türkçe (çalışma dili)
- HTML rapor: `config.report.language` ile `tr` | `en` — müşteriye göre seçilir

---

## 11. Gizlilik

- Depo **private**. Müşteri verisi asla commit edilmez (`runs/`, `raw/`,
  `references/` .gitignore'da).
- Gerçek hasta/müşteri ismi veya PII hiçbir dosyaya, teste veya dokümana yazılmaz.
- Demo veri public olduğu için doğrulama çıktıları paylaşılabilir.
- PLAN.md Bölüm 16 (Data Security & Retention) bağlayıcıdır.
