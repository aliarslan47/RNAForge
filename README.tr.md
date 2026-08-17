# RNAForge

Tekrarlanabilir, modüler **bulk RNA-seq** analiz pipeline'ı — ham FASTQ'tan tek, kendi kendine yeten
HTML rapora; diferansiyel ekspresyonun üstünde tam bir fonksiyonel-analiz katmanıyla.

İngilizce sürüm: [README.md](README.md) · Referans doküman: [PLAN.md](PLAN.md) (v1.4)

## Ne yapar

Ham okumaları biyolojiye taşıyan aşamalı bir pipeline:

```
validate → qc → trim → quant → counts → de → figures → report
                                          └→ enrich · kegg · gsea · semantic · amr · operon · ppi
```

- **Çekirdek**: girdi/tasarım doğrulama, QC, nazik kırpma, hizalama/kantifikasyon, DESeq2 diferansiyel
  ekspresyon, yayın kalitesinde figürler ve çift dilli (`tr`/`en`) kendi kendine yeten HTML rapor.
  QC → trim → hizalama araç zinciri **okuma tipine** göre otomatik seçilir (aşağıya bak);
  kantifikasyon `organism_type` ile yönlendirilir (prokaryot: Bowtie2/minimap2 + featureCounts ·
  ökaryot: Salmon + tximport). İki boyut da aynı gen × örnek sayım matrisinde buluşur.
- **Okuma tipleri (kısa / uzun)**: Illumina okumaları kısa-okuma zincirini çalıştırır (FastQC → fastp →
  Bowtie2); ONT/PacBio okumaları uzun-okuma zincirini (NanoPlot → Pychopper+chopper → minimap2 →
  featureCounts `-L`). Okuma tipi m01'de otomatik tespit edilir; m05'ten itibaren (DESeq2 ve tüm
  fonksiyonel analiz) pipeline okuma-tipinden bağımsızdır.
- **Fonksiyonel analiz** (hepsi opsiyonel, organizma-agnostik, hiçbiri yeni FAIL kapısı üretmez):
  GO aşırı-temsil (ORA), KEGG yolak ORA, GSEA (fgsea), REVIGO-benzeri semantik indirgeme,
  AMR + virülans overlay (abricate/CARD/VFDB), operon tahmini + koordinasyon ve STRING
  protein-etkileşim modülleri (Louvain community detection).

### Kalite kapıları (sonuçlar neden güvenilir)

Doğru bir pipeline bile kötü girdiden makul görünen ama **sahte** bir sonuç üretir. RNAForge kapıları
ikili politikayla zorlar:

- **FAIL** → sonuç **geçersiz**: koşu durur, biyolojik çıktı üretilmez (exit 1).
- **WARN** → sonuç **şüpheli**: üretilir ama öyle damgalanır.

Eşikler veridir (`profiles/{prokaryote,eukaryote,prokaryote_long}.yml`); ezilen eşik rapora yazılır
(sessiz gevşetme yok). Uzun-okuma koşuları `prokaryote_long` profilini kullanır; eşikleri bilinçli
olarak permissive ve damgalıdır (ONT kalitesi ~Q10–15, Q30 değil): yalnız katastrofik hizalama hatası
(yanlış referans) FAIL verir, ONT'de doğal olan düşük survival/assignment WARN olur. Her koşu ayrıca
bir güvence kartı yazar (`UNKNOWN`/`INVALID`/`SUSPECT`/`TRUSTWORTHY`).

## Pipeline modülleri

| Aşama | Subcommand | Ne yapar |
|---|---|---|
| m01 | `validate` | Config + metadata + FASTQ doğrulama, platform + okuma-tipi tespiti, tasarım kapıları |
| m02 | `qc` | kısa: FastQC · uzun: NanoPlot (teşhis; koşuyu asla durdurmaz) |
| m03 | `trim` | kısa: fastp (nazik) · uzun: Pychopper+chopper (cDNA) / chopper (direct-RNA) |
| m04 | `quant` | Hizalama (prok kısa: Bowtie2 · uzun: minimap2 · ökaryot: Salmon) |
| m05 | `counts` | featureCounts (uzun-okumada `-L`) → gen × örnek sayım matrisi |
| m06 | `de` | DESeq2 diferansiyel ekspresyon |
| m07 | `figures` | PCA, volcano, MA, heatmap, dispersiyon, … (PNG 300dpi + SVG) |
| m08 | `report` | Tek, kendi kendine yeten çift dilli HTML rapor (okuma-tipi farkında) |
| m09 | `enrich` | GO aşırı-temsil (ORA), hipergeometrik + BH |
| m10 | `kegg` | KEGG yolak ORA |
| m11 | `gsea` | Sıralı gen listesinde GSEA (fgsea) |
| m12 | `semantic` | GO terimlerinin REVIGO-benzeri semantik indirgemesi (+ MDS haritası) |
| m13 | `amr` | AMR (CARD + AMRFinderPlus) + virülans (VFDB) genlerinin DE'ye overlay'i |
| m14 | `operon` | Operon tahmini (intergenik mesafe) + DE koordinasyonu |
| m15 | `ppi` | STRING PPI alt-ağı + Louvain community modülleri |
| m16 | `seqqc` | rRNA% (SortMeRNA) + strandedness (RSeQC) — WARN kapıları |
| m17 | `alignqc` | insert-size + coverage + read-distribution (samtools/RSeQC) |
| m18 | `multiqc` | koşu genelinde toplu MultiQC görünümü (en son) |

Downstream analizler (m09–m15) organizma- ve okuma-tipi-agnostiktir ve bir koşuyu asla geçersiz kılmaz —
verdict kalite kapılarından değişmeden taşınır. QC ekleri (m16–m18) tanısaldır.

## Kurulum

```bash
conda env create -f envs/rnaforge-core.yml     # orkestrasyon (Python) + networkx/numpy/scipy
conda activate rnaforge-core
pip install -e .
```

Araç ortamları (bir kez kurulur, modüller tarafından kullanılır):

```bash
conda env create -f envs/rnaforge-qc.yml           # FastQC, fastp
conda env create -f envs/rnaforge-quant-prok.yml   # Bowtie2, samtools, featureCounts
conda env create -f envs/rnaforge-quant-euk.yml    # Salmon
conda env create -f envs/rnaforge-longread.yml     # minimap2, NanoPlot, Pychopper, chopper, samtools
conda env create -f envs/rnaforge-de.yml           # R: DESeq2, ggplot2, fgsea
conda env create -f envs/rnaforge-amr.yml          # abricate (CARD/VFDB)
conda env create -f envs/rnaforge-seqqc.yml        # SortMeRNA, RSeQC, MultiQC (m16/m18)
```

## Kullanım

```bash
# çekirdek zincir (baştan sona aynı --run-id)
rnaforge validate --config config/config.yaml --metadata samples.tsv --run-id demo
rnaforge qc       --config config/config.yaml --metadata samples.tsv --run-id demo
rnaforge trim     --config config/config.yaml --metadata samples.tsv --run-id demo
rnaforge quant    --config config/config.yaml --metadata samples.tsv --run-id demo
rnaforge counts   --config config/config.yaml --metadata samples.tsv --run-id demo
rnaforge de       --config config/config.yaml --metadata samples.tsv --run-id demo
rnaforge figures  --config config/config.yaml --metadata samples.tsv --run-id demo

# opsiyonel QC / tanısal (m04 gerekir; tanısal figür/tablo üretir, asla FAIL vermez)
rnaforge seqqc    --config config/config.yaml --metadata samples.tsv --run-id demo  # rRNA% + strandedness (m16)
rnaforge alignqc  --config config/config.yaml --metadata samples.tsv --run-id demo  # insert-size + coverage + read-distribution (m17)
rnaforge multiqc  --config config/config.yaml --metadata samples.tsv --run-id demo  # toplu MultiQC görünümü (m18, en son)

# opsiyonel fonksiyonel analizler (herhangi alt-küme; her biri referans verisi ister — aşağıya bak)
rnaforge enrich   --config config/config.yaml --metadata samples.tsv --run-id demo
rnaforge kegg     --config config/config.yaml --metadata samples.tsv --run-id demo
rnaforge gsea     --config config/config.yaml --metadata samples.tsv --run-id demo
rnaforge semantic --config config/config.yaml --metadata samples.tsv --run-id demo
rnaforge amr      --config config/config.yaml --metadata samples.tsv --run-id demo
rnaforge operon   --config config/config.yaml --metadata samples.tsv --run-id demo
rnaforge ppi      --config config/config.yaml --metadata samples.tsv --run-id demo

# raporu en son üret — hangi analizler koştuysa onları gömer
rnaforge report   --config config/config.yaml --metadata samples.tsv --run-id demo
```

> Not: `python -m rnaforge.cli` ÇALIŞMAZ (main-guard yok); kurulu `rnaforge` entry point'ini kullan.

### Metadata biçimi (TSV)

| Sütun | Zorunlu | Açıklama |
|---|---|---|
| `sample_id` | evet | Benzersiz örnek kimliği |
| `condition` | evet | Deney grubu; ≥2 seviye ve her birinde ≥2 replika gerekir |
| `fastq_1` | evet | R1 yolu (veya single-end okumalar) |
| `fastq_2` | hayır | Paired-end için R2 yolu |
| `subject` | hayır | Eşleşmiş/subject id; tespit edilir, eşleşmiş görünüyorsa bilinçli ele alınmalı |
| `batch` | hayır | Batch/kovaryat; tasarım formülü `batch` kullanıyorsa zorunlu |

## Referans verisi (tek sefer hazırlık, git-ignore'lu)

Fonksiyonel analizler `references/` altındaki yerel dosyaları okur (asla commit edilmez). Organizmanız
için bir kez indirin (*E. coli* K-12 örnekleri):

```bash
# GO ontolojisi (m09/m12) + organizma GO anotasyonu (EBI-GOA)
curl -L -o references/go/go-basic.obo http://purl.obolibrary.org/obo/go/go-basic.obo
curl -L https://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/18.E_coli_MG1655.goa \
     -o references/ecoli_bw25113/ecoli.gaf

# KEGG (m10) — organizma-başı REST dosyaları
curl -s https://rest.kegg.jp/link/pathway/eco > references/kegg/eco/pathway_links.tsv
curl -s https://rest.kegg.jp/list/pathway/eco > references/kegg/eco/pathway_names.tsv
curl -s https://rest.kegg.jp/list/eco        > references/kegg/eco/gene_list.tsv

# STRING (m15) — taxon-başı ağ
curl -s https://stringdb-downloads.org/download/protein.info.v12.0/511145.protein.info.v12.0.txt.gz \
     -o references/string/511145/protein.info.txt.gz
curl -s https://stringdb-downloads.org/download/protein.links.v12.0/511145.protein.links.v12.0.txt.gz \
     -o references/string/511145/protein.links.txt.gz
```

AMR/virülans (m13) abricate'in paketli CARD/VFDB veritabanlarını kullanır (ayrı indirme yok).

## Temel tasarım kararları

- **`organism_type` zorunlu ve varsayılanı yok** (`prokaryote` | `eukaryote`). Yalnız kantifikasyonu
  (m04/m05) yönlendirir; iki yol da aynı gen × örnek sayım matrisinde buluşur, böylece tüm downstream
  adımlar (m06–m15) organizma-agnostiktir.
- **İki yönlendirme boyutu: `organism_type` × okuma tipi.** Okuma tipi (kısa/uzun) m01'de FASTQ'tan
  otomatik tespit edilir (Illumina → kısa; ONT/PacBio → uzun) ve m02–m05'i yönlendirir; `organism_type`
  m04/m05'i yönlendirir. Tanımlanamayan platformlar yine net hatayla reddedilir (asla yanlış yola sessizce
  sokulmaz).
- **ONT uzun okumalar için `library.chemistry` zorunlu** (`cdna` | `direct_rna`) — FASTQ'tan tespit
  edilemez ve m03 uzun-okuma ön-işlemesini seçer (cDNA → Pychopper+chopper; direct-RNA → yalnız chopper).
  PacBio HiFi bunu gerektirmez.
- **Kırpma bilinçli olarak nazik.** Agresif kalite kırpması ekspresyon tahminlerini bozar
  ([Williams ve ark. 2016](https://doi.org/10.1186/s12859-016-0956-2)); bozulmayı asgari-uzunluk
  filtresi engeller.
- **Uydurma sonuç yok.** Belirsiz annotation birleştirmeleri (gen sembolüyle) tahmin edilmez, atılır;
  tahmin edilen yapılar (operon, STRING etkileşimleri) rapora tahmin olarak damgalanır.

## Geliştirme

```bash
conda run -n rnaforge-core --cwd "$(pwd)" python -m pytest -q
```

Depo kökünden çalıştırın (test suit'i `tests.conftest`'i import eder).

## Gizlilik

Müşteri verisi asla commit edilmez. `runs/`, `raw/` ve `references/` git-ignore'ludur.
