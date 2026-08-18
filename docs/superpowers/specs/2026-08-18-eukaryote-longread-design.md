# Ökaryot uzun-okuma yolu (gen düzeyi, transkriptom-hizalama) — Tasarım

**Tarih:** 2026-08-18
**Durum:** Onaylı (brainstorming), implementasyona hazır.
**Kapsam:** Ökaryot ONT/PacBio uzun-okuma → **gen düzeyi** count matrisi. İzoform-DE, genom-splice
hizalama, salmon-alignment-mode KAPSAM DIŞI (YAGNI).

## 1. Amaç ve değişmezlik

Ökaryot + long-read için niceleme. Ayrım YALNIZ m04/m05'te; m06+ organizma/okuma-agnostik,
DEĞİŞMEZ. Ortak sözleşme `counts.tsv` (`gene\t<sample_id...>`). Referans = mevcut ökaryot referansı
(`transcriptome_fasta` + `tx2gene`); genom/GTF gerekmez.

## 2. Yönlendirme (router)

m04/m05'te sıra: önce `organism_type`, ökaryot içinde `read_type`:
- **m04** eukaryote: `short` → `_quant_euk` (salmon, mevcut) · `long` → **`_quant_euk_long` (YENİ)**
- **m05** eukaryote: `short` → `_counts_euk` (tximport, mevcut) · `long` → **`_counts_euk_long` (YENİ)**
Prokaryot dalları (short bowtie2 / long minimap2-genom) DEĞİŞMEZ.

## 3. m04-euk-long — minimap2 → transkriptom

- `run_minimap2` (mevcut) yeniden kullanılır; **referans = `transcriptome_fasta`** (transkript=hedef,
  intron yok → splice gerekmez). Preset `resolve_platform` → ont: `map-ont`, pacbio_hifi: `map-hifi`
  (`minimap2_preset`, mevcut). Trimlenmiş okuma `trimmed_reads` (m03-long çıktısı; ONT tek-uçlu).
- Çıktı BAM: `quantification/<sid>/aligned.sorted.bam` (m05 sözleşme yolu, prok-long ile aynı).
- `mapping_rate` = `parse_flagstat_mapped` (primary-mapped/primary, mevcut). **DIAGNOSTİK — FAIL kapısı
  YOK** (tüm long yolları gibi; ONT düşük oranı yanlış FAIL'lemesin). Summary'ye stat.

## 4. m05-euk-long — primer-hizalama sayımı → tx2gene topla

- **Yeni** `minimap2.count_primary_alignments(bam_path) -> dict[transcript_id, int]`:
  `samtools view -F 2308 <bam>` (0x904 = unmapped+secondary+supplementary hariç → okuma başına tek
  primer hizalama) → hedef (sütun 3) sayımı. Env `rnaforge-longread` (samtools).
- **Yeni** `tximport.parse_tx2gene(path) -> dict[tx_id, gene_id]` (TSV `tx\tgene`; mevcut R tximport'la
  aynı dosya, artık Python da okur).
- `_counts_euk_long`: her örnek BAM → transkript sayımı → tx2gene ile gen'e topla (gen-içi izoform
  çoklu-eşleşmesi aynı gene toplanır). Gen evreni = tüm örneklerde gözlenen genlerin birleşimi (eksik=0).
  `counts.tsv` (`gene\t<sid...>`, sütun→sample_id KONUMLA). Boş matris → yüksek sesle hata. Diagnostik
  (kapı yok — salmon/uzun zaten hizalamada eledi).
- TPM/FPKM ökaryot-long'da atlanır (uzun-okuma uzunluk-normalizasyonu farklı; DE için counts yeterli).

## 5. Config / referans

Değişiklik YOK. Ökaryot referansı zaten `transcriptome_fasta` + `tx2gene`. Platform `ont`/`pacbio_hifi`
config'ten (kısa cDNA yanlış-tespitini ezmek için, commit deba674 deseni). `library.chemistry` ökaryot-long
için m01 zorunlu kılar mı? — ONT-long chemistry zaten m01'de zorunlu (mevcut); ökaryotta da geçerli.

## 6. Test (modül deseni)

- `count_primary_alignments`: küçük SAM→BAM fixture (primer + secondary/supplementary/unmapped karışık) →
  yalnız primer sayılır, hedef başına doğru; env yoksa skip (samtools gerçek).
- `parse_tx2gene`: TSV → map.
- m04: eukaryote+long dispatch → `_quant_euk_long` (run_minimap2 monkeypatch), mapping_rate diagnostik
  (FAIL yok), BAM yolu.
- m05: eukaryote+long dispatch → `_counts_euk_long` (count monkeypatch), tx2gene toplama, counts.tsv.
- Regresyon: ökaryot-short (salmon/tximport) + prokaryot yolları değişmedi.

## 7. Doğrulama (implementasyon SONRASI — bilimsel adım)

Gerçek ökaryot ONT/PacBio cDNA seti (model organizma, temiz DE kontrastı, ≥3 replika, yayınlanmış DEG,
transkriptom+tx2gene). Kör indirme yok. Typhi/airway disiplini. MVP kodunu bloke etmez.

## 8. Karar özeti

- **Gen düzeyi + transkriptom hizalama:** m06+ reuse, ökaryot-short ile referans tutarlı; izoform-DE ayrı
  büyük iş (YAGNI).
- **Primer-sayım (salmon-EM değil):** ONT için kısa-okuma EM varsayımlarından kaçınır; gen-düzeyi toplama
  gen-içi çoklu-eşleşmeyi doğal çözer. Basit, tekrarlanabilir, denetlenebilir.
- **Diagnostik (FAIL yok):** tüm uzun-okuma yollarıyla tutarlı (prokaryote_long deseni).
