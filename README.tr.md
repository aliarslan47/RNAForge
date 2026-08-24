# RNAForge

Tekrarlanabilir, modüler **bulk RNA-seq** pipeline'ı — ham FASTQ'tan tek, kendine yeten bir HTML rapora; diferansiyel ifadenin üzerine kurulu bir fonksiyonel analiz katmanıyla.

[![Pipeline DAG](https://img.shields.io/badge/pipeline-DAG-0d6b8f)](https://aliarslan47.github.io/RNAForge/pipeline_architecture.html)
[![organism](https://img.shields.io/badge/organism-prokaryote%20%C2%B7%20eukaryote%20%C2%B7%20metatranscriptome-2f8f5b)](https://aliarslan47.github.io/RNAForge/pipeline_architecture.html)
[![reads](https://img.shields.io/badge/reads-short%20%C2%B7%20long-c07211)](https://aliarslan47.github.io/RNAForge/pipeline_architecture.html)

**Türkçe** · [English](README.md)

## Nedir?

RNAForge, Forge ailesinin bulk RNA-seq üyesidir — BacForge (bakteri) ve VirusForge (virüs/faj) ile aynı mimari, ancak ayrı ve izole bir kurulum. Ham okumaları tek komutla biyolojiye taşır ve çift dilli (`tr`/`en`), kendine yeten bir HTML raporla sonlanır.

## Ne yapar?

Aşamalı bir pipeline — `validate → qc → trim → quant → counts → de → figures → report` — üzerine opsiyonel fonksiyonel analiz katmanı (GO / KEGG / GSEA / semantik / AMR / operon / PPI).

İki yönlendirme ekseni aynı gen × örnek sayım matrisinde birleşir; böylece DE'den sonrası tüm adımlar bağımsızdır:

- **Organizma** (`organism_type`): prokaryot (Bowtie2/minimap2 + featureCounts) · ökaryot (Salmon + tximport) · metatranskriptom (rRNA arındırma + Kraken2/Bracken + gen-kataloğu).
- **Okuma tipi** (otomatik algılanır): kısa (Illumina: FastQC → fastp → Bowtie2) · uzun (ONT/PacBio: NanoPlot → Pychopper+chopper → minimap2).

Tasarımı gereği güvenilir: iki kademeli kalite kapıları (**FAIL** run'ı durdurur, **WARN** şüpheli sonucu damgalar), eşikler veridir (`profiles/*.yml`) ve her run bir güven kartı yazar. Uydurma sonuç yok.

Etkileşimli çift dilli düğüm grafiği: **[render edilmiş diyagram](https://aliarslan47.github.io/RNAForge/pipeline_architecture.html)**.

## Kurulum

```bash
bash install.sh
conda run -n rnaforge-core rnaforge doctor   # gerekli tüm env'lerin varlığını doğrula
```

İdempotenttir; sürümü sabitlenmiş dokuz conda ortamı yaratır (`envs/*.yml`). `dorado` (ONT ham-sinyal basecalling, m00) conda dışında kurulan yalnızca-GPU bir ikili dosyadır — sadece FAST5/POD5 girdisinde gerekir.

## Kullanım

```bash
# tüm pipeline (FAIL'de durur, resume edilebilir)
rnaforge run --config config/config.yaml --metadata samples.tsv --run-id demo

# rapordan önce opsiyonel aşamalar ekle
rnaforge run ... --include enrich,kegg,gsea,seqqc,alignqc,multiqc

# ya da çekirdek zincirin bir dilimini çalıştır
rnaforge run ... --from trim --to counts
```

Her aşama aynı `--run-id` ile elle de sürülebilir. Kurulu `rnaforge` giriş noktasını kullan (`python -m` değil).

## Modüller

| Kod | Alt komut | Ne yapar |
|---|---|---|
| m00 | `basecall` | ONT ham sinyal (FAST5/POD5) → dorado ile FASTQ (GPU); opsiyonel |
| m01 | `validate` | Config/metadata/FASTQ doğrulama, platform + okuma-tipi algılama |
| m02 | `qc` | kısa: FastQC · uzun: NanoPlot |
| m03 | `trim` | kısa: fastp (yumuşak) · uzun: Pychopper+chopper / chopper |
| m04 | `quant` | Hizalama (prok: Bowtie2/minimap2 · ökar: Salmon) |
| m05 | `counts` | featureCounts → gen × örnek sayım matrisi |
| m06 | `de` | DESeq2 diferansiyel ifade |
| m07 | `figures` | PCA, volkan, MA, ısı haritası, dispersiyon (PNG 300dpi + SVG) |
| m08 | `report` | Tek, kendine yeten çift dilli HTML rapor |
| m09 | `enrich` | GO aşırı-temsil analizi (ORA) |
| m10 | `kegg` | KEGG yolak ORA |
| m11 | `gsea` | Sıralanmış gen listesinde GSEA (fgsea) |
| m12 | `semantic` | GO terimlerinin REVIGO-benzeri semantik indirgemesi |
| m13 | `amr` | AMR (CARD + AMRFinderPlus) + virülans (VFDB) örtüşmesi |
| m14 | `operon` | Operon tahmini + DE koordinasyonu |
| m15 | `ppi` | STRING PPI alt-ağı + Louvain modülleri |
| m16 | `seqqc` | rRNA% (SortMeRNA) + iplik yönü (RSeQC) |
| m17 | `alignqc` | insert boyu + kapsam + okuma dağılımı |
| m18 | `multiqc` | run genelinde toplu MultiQC görünümü (en son) |

Metatranskriptom run'ları `trim` ile `quant` arasına `rrna-deplete` ve `taxonomy` aşamalarını otomatik ekler. Aşağı-akış analizleri (m09–m18) organizma ve okuma-tipinden bağımsızdır ve bir run'ı asla geçersiz kılmaz. Tam tasarım, metadata formatı ve referans-veri hazırlığı `PLAN.md` ve `docs/` içindedir.

---

Forge ailesi: **RNAForge** (bulk RNA-seq) · [BacForge](https://github.com/aliarslan47/BacForge) (bakteri) · [VirusForge](https://github.com/aliarslan47/VirusForge) (virüs/faj) · [PipelineForge](https://github.com/aliarslan47/PipelineForge) (DAG üreticisi). Müşteri verisi asla commit'lenmez (`runs/`, `raw/`, `references/` git-ignore'da).
