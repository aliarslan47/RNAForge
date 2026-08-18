# Ökaryot kısa-okuma yolu (Salmon + tximport) — Tasarım

**Tarih:** 2026-08-18
**Durum:** Onaylı (brainstorming), implementasyona hazır.
**Kapsam:** MVP'nin ikinci ana kolu — ökaryot **kısa-okuma** (Illumina) niceleme.
Ökaryot uzun-okuma, izoform-seviyesi DE, GTF→tx2gene türetme KAPSAM DIŞI (YAGNI).

## 1. Amaç ve bağlam

Şu an `organism_type: eukaryote` → m04'te `NotImplementedError`. Prokaryot kolu
(bowtie2 + featureCounts) ökaryotta yanlış: intron/izoform nedeniyle genoma hizalama +
featureCounts uygun değil. Ökaryot için doğru yol **Salmon** (transkriptom niceleme) +
**tximport** (transkript→gen toplama).

**Değişmezlik ilkesi:** Ayrım YALNIZ m04/m05'te. m06+ (DESeq2, figürler, GO/KEGG/GSEA/
REVIGO, rapor) organizma-agnostiktir ve **değişmez**. İki kol ortak `counts.tsv`
sözleşmesinde (`gene\t<sample_id...>`) buluşur.

## 2. Mimari

`m04_quant.run_quant` mevcut router deseni (organism_type):
- `prokaryote` → `_quant_short` (bowtie2, mevcut) / `_quant_long` (minimap2, mevcut)
- `eukaryote` → **`_quant_euk` (YENİ, Salmon)** — eski `NotImplementedError` kaldırılır

`m05_counts.run_counts` router:
- `prokaryote` → `_counts_short`/`_counts_long` (featureCounts, mevcut)
- `eukaryote` → **`_counts_euk` (YENİ, tximport)**

Ökaryot read_type şimdilik yalnız `short` (long kolu kapsam dışı). Ökaryot + ONT
gelirse m04-euk yüksek sesle `NotImplementedError` verir (sessiz yanlış koşma yok).

## 3. m04-euk — Salmon

**Yeni modül `rnaforge/salmon.py`** (fastp/bowtie2/minimap2 deseni: saf parser + runner).
Env: **`rnaforge-quant-euk`** (salmon 2.3.4 zaten kurulu).

### Index (decoy-aware, genom opsiyonel)
- `reference.genome_fasta` verilmişse → **decoy-aware selective alignment**:
  `decoys.txt` = genom kontig adları; `gentrome.fa` = `transcriptome_fasta` + `genome_fasta`;
  `salmon index -t gentrome.fa -d decoys.txt -i <idx> -k 31`.
  Gerekçe: anotasyonsuz/intergenik bölgelerden gelen sahte eşleşmeleri eler (Salmon önerisi,
  en doğru — [[feedback_dogruluk_kontrol]]).
- `genome_fasta` YOKSA → transkriptom-only index (`salmon index -t transcriptome_fasta`)
  + **yüksek sesle log** ("decoy yok, doğruluk için genome_fasta öner").
- Index bir kez kurulur, run_dir altında cache (varsa yeniden kurmaz).

### Quant
- Her örnek: `salmon quant -i <idx> -l A [-r r1 | -1 r1 -2 r2] -p <threads> -o <sample>`.
  `-l A` = kütüphane tipini otomatik tespit (strandedness). Paired-end: fastq_2 varsa `-1/-2`.
- Trimlenmiş okuma girişi: `m03.trimmed_reads(run_dir, sample)` (tek kaynak sözleşmesi).
- Çıktı: `quantification/<sample>/quant.sf` + `aux_info/meta_info.json` + `lib_format_counts.json`.

### Parser + kapı
- `parse_salmon_meta(meta_info.json)` → `percent_mapped` (salmon mapping rate).
- **`mapping_rate`** = `alignment_rate` kapısına bağlanır (eukaryote.yml eşiği 0.50).
  Prokaryot short deseni: mapping_rate < eşik → **FAIL** (koşu durur, katastrofik).
  eukaryote.yml `permissive: true` damgalı (doğrulanmış set gelene dek).

## 4. m05-euk — tximport

**Yeni R betiği `rnaforge/scripts/tximport.R`** (deseq2.R deseni).
Env: **`rnaforge-de`** (R + Bioconductor + DESeq2 zaten var) — `bioconductor-tximport`
eklenir; ikinci R env AÇILMAZ (reuse). `envs/rnaforge-de.yml` güncellenir.

### Girdi/çıktı
- Girdi: her örnek `quant.sf` + `reference.tx2gene` (config; kullanıcı sağlar, TSV: `tx_id\tgene_id`).
- tximport `type="salmon"`, `tx2gene` ile transkript→gen toplama.
- **⚠️ KRİTİK DOĞRULUK KARARI:** `countsFromAbundance="lengthScaledTPM"`.
  Neden: salmon+tximport'ta ham toplanmış sayımlar DESeq2'ye doğrudan verilmemeli
  (transkript-uzunluğu bias). Doğru yol ya (a) length-offset matrisini DESeq2'ye vermek
  (`DESeqDataSetFromTximport`) ya da (b) `lengthScaledTPM`/`scaledTPM` ile uzunluk-düzeltilmiş
  sayım üretmek. (a) m06'ya offset/organizma bilgisi sızdırır → agnostiklik bozulur.
  **(b) seçildi:** `lengthScaledTPM` uzunluk-düzeltilmiş sayım verir; m06'nın mevcut DESeq2'si
  bunu düz sayım gibi okur ve **istatistiksel olarak doğru olur**, m06 DEĞİŞMEZ.
- Çıktı: `quantification/counts.tsv` (`gene\t<sample_id...>`, sütun→sample_id KONUMLA) —
  mevcut `_write_count_outputs` yardımcısıyla. TPM: salmon quant.sf TPM'lerinden gen-toplamı
  (tximport `abundance`); FPKM ökaryotta atlanabilir (salmon TPM verir, uzunluk normalize).
- `assignment_rate` diagnostik (kapı prokaryot-short'a özgü; ökaryotta salmon zaten atar).

## 5. Config / referans değişiklikleri

- `Reference` dataclass'a opsiyonel **`genome_fasta: str | None = None`** (ökaryot decoy;
  prokaryotta zaten `genome_fasta` var — alan paylaşılır, char yok).
  DİKKAT: prokaryot `REQUIRED_REFERENCE`'te `genome_fasta` zaten zorunlu. Ökaryotta opsiyonel.
- `REQUIRED_REFERENCE[eukaryote]` = `(transcriptome_fasta, tx2gene)` DEĞİŞMEZ (genome opsiyonel).
- `quantification` config ökaryotta kullanılmaz (feature_type/attribute prokaryot-özgü);
  ökaryot yolu bunları okumaz.

## 6. Test stratejisi (modül deseni)

- `tests/test_salmon.py`: `parse_salmon_meta` (meta_info.json fixture → percent_mapped),
  `build_salmon_index` decoy var/yok dallanma (komut kurulumu doğrulama, run monkeypatch),
  `run_salmon_quant` tek/paired komut.
- `tests/test_m04_quant.py`: eukaryote artık NotImplementedError vermez; `_quant_euk`
  salmon çağırır (monkeypatch), mapping_rate kapısı (PASS + FAIL) yazılır.
- `tests/test_m05_counts.py`: `_counts_euk` tximport çağırır (R runner monkeypatch),
  counts.tsv sözleşmesi (gene\t<ids>), tx2gene toplama.
- `tests/test_config.py`: eukaryote genome_fasta opsiyonel; verilmezse None, verilirse yol.
- Gerçek-araç entegrasyonu: küçük sentetik transkriptom + 2 örnek (env yoksa skip).

## 7. Doğrulama (implementasyon SONRASI — ayrı bilimsel adım)

Prokaryot Typhi doğrulamasındaki disiplin: yayınlanmış DEG tablolu, güçlü sinyalli,
üst-segment ökaryot **Illumina** RNA-seq seti seçilir (kör indirme yok — bilimsel karar).
Aday kriterleri: model organizma (insan/fare/maya/Arabidopsis), temiz tedavi-kontrol,
≥3 replika, referans transkriptom (GENCODE/Ensembl) + tx2gene erişilebilir. Uçtan uca
+ konkordans (yön + anahtar gen setleri). Bu adım MVP kodunu bloke etmez.

## 8. Karar özeti (kodda görünmeyen "neden"ler)

- **Salmon (bowtie2 değil):** ökaryotta izoform/intron → genom-hizalama + featureCounts
  yanlış; transkriptom-niceleme doğru. Prokaryotta tersi (tx2gene ~1:1, Salmon anlamsız).
- **decoy opsiyonel:** doğruluk-öncelikli varsayılan ama küçük/decoy'suz kurulumları da
  destekle; sessiz düşürme yok, loglanır.
- **lengthScaledTPM:** m06 agnostikliğini korur (bkz. §4). Alternatif offset-yolu reddedildi.
- **rnaforge-de env'inde tximport:** R+Bioconductor reuse; ikinci R env maliyeti yok.
