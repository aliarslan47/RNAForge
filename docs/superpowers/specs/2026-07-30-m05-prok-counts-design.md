# m05 — Count Matrix (prokaryot: featureCounts) — Tasarım

**Tarih:** 2026-07-30
**Durum:** Onaylandı (Ali — "1", otonom oturum)
**İlgili:** PLAN.md v1.3 §5, §14; m04 spec (aynı modül deseni)

## 1. Problem

m04 örnek başına sıralı BAM üretti. m05 bunları anotasyona (GFF/GTF) göre sayıp
**gen × örnek count matrisi** üretir — PLAN §5'teki **ortak sözleşme**. Bu matris pipeline'ın
dönüm noktasıdır: buradan sonra prokaryot/ökaryot ayrımı biter, m06 (DESeq2) ve sonrası hep
bu matrisi tüketir. Bu spec **prokaryot yolunu (featureCounts)** kapsar; ökaryot (tximport +
tx2gene) ayrı spec.

## 2. Onaylanan kararlar

1. **featureCounts tüm BAM'lere tek çağrı.** featureCounts çok-örnekli girdide native olarak
   gen×örnek matrisi üretir (`counts.txt`: Geneid…Length + BAM başına bir sütun). Sütunlar BAM
   sırasındadır → `sample_id`'ye eşlenir. Temiz `quantification/counts.tsv` (Geneid + örnek
   başına sütun, BAM yolları DEĞİL sample_id başlıkları) yazılır = **ortak sözleşme**.
2. **Veri kapısı: `assignment_rate` — FAIL kapısı.** featureCounts summary'den örnek başına
   `Assigned / toplam(tüm Status satırları)` profil eşiğinin (`assignment_rate`, yeni; prok 0.50)
   altındaysa **FAIL → koşu durur** (çok düşük atama = yanlış anotasyon/tür → güvenilmez sayım;
   m04 alignment kapısı deseniyle aynı). Eşik `profiles/*.yml`'ye eklenir.
3. **featureCounts parametreleri config-driven.** Yeni `quantification` config bölümü:
   `feature_type` (varsayılan `exon`), `attribute` (varsayılan `gene_id`). Anotasyon kaynağı
   değişir (NCBI prokaryot GFF3 tipik `CDS`/`locus_tag`); kullanıcı ezer. **`config.py`'nin
   `KNOWN_TOP_LEVEL_KEYS`'ine `quantification` EKLENİR** (aksi halde config sertleştirme reddeder).
4. **Ön koşul: m04 done** (BAM gerekir). Zincir m01→m03→m04→m05, aynı `--run-id`. m04 yapılmamışsa
   net `ValueError`.
5. **Paired-end destekli.** featureCounts `-p --countReadPairs` (paired) yoksa read-bazlı.
   Örnek metadata'da `fastq_2` varsa paired.

### Reddedilen seçenekler
- **Örnek başına ayrı featureCounts + sonra birleştir:** featureCounts zaten çok-örnekli matris
  üretir; ayrı çağrı gereksiz karmaşa + tutarsızlık riski.
- **assignment_rate = WARN:** Çok düşük atama sayımı GEÇERSİZ kılar → FAIL (alignment ile tutarlı).
- **feature_type/attribute sabit kod:** Anotasyon kaynağı değişir; sabit `exon/gene_id` prokaryot
  GFF3'te sıfır sayım verir → sahte "boş matris". Config-driven doğru.

## 3. Arayüz (public sözleşme)

```python
def run_counts(config: Config, metadata_path: Path, run_dir: Path,
               force: bool = False) -> dict: ...
```
- **Ön koşul:** m04 done değilse `ValueError`. (Router değil — m05 prok; euk tximport ayrı spec.)
- **Dönüş:** özet dict: `n_samples`, `n_genes`, örnek başına `assignment_rate`/`assigned`, `gate_counts`, `resumed?`.

```python
# rnaforge/featurecounts.py
def parse_counts(counts_text: str) -> tuple[list[str], dict[str, list[int]]]: ...
    # (gene_ids, {bam_column_name: counts_list})
def parse_summary(summary_text: str) -> dict[str, float]: ...  # {bam_column: assignment_rate}
def run_featurecounts(bams: list[Path], gff: Path, out_dir: Path, feature_type: str,
                      attribute: str, paired: bool = False, threads: int = 4, env=...) -> FeatureCountsResult: ...
```

## 4. Bileşenler

- **`rnaforge/featurecounts.py`** — saf parserlar + runner (m04 deseni):
  - `parse_counts(counts_text) -> (gene_ids, {column: counts})`: `#` yorum + Geneid başlığını
    atlar; Length'ten sonraki sütunlar (BAM'ler) → sütun adı: sayım listesi.
  - `parse_summary(summary_text) -> {column: assignment_rate}`: her BAM sütunu için
    `Assigned / sum(tüm Status)`; toplam 0 ise 0.0.
  - `run_featurecounts(bams, gff, out_dir, feature_type, attribute, paired, threads, env)
    -> FeatureCountsResult`: featureCounts'u çalıştırır (`counts.txt` + `.summary`), parse eder.
  - `FeatureCountsResult` frozen: `gene_ids: list[str]`, `counts: dict[str, list[int]]`
    (BAM sütun adı → sayımlar), `assignment_rates: dict[str, float]` (BAM sütun adı → oran).
  - `FeatureCountsParseError(ValueError)`, `FeatureCountsRunError(RuntimeError)`.
- **`rnaforge/modules/m05_counts.py`** — `run_counts` + `build_count_gates(assignment_rates:
  dict[str,float], profile) -> list[GateResult]` (m04 alignment kapısı deseni). BAM sütun→sample_id
  eşlemesi **KONUMLA (isimle değil):** `run_counts` BAM'leri metadata örnek sırasında verir;
  featureCounts sütunları aynı sırada döner; `parse_counts` insertion-order koruyan dict verir →
  `zip(sample_ids, counts.values())`. BAM adları sample_id'den farklı olabilir; sıra sözleşmedir.
- **`rnaforge/config.py`** — `quantification` bölümü (`feature_type`, `attribute`) + `KNOWN_TOP_LEVEL_KEYS`.
- **`rnaforge/profiles/{prokaryote,eukaryote}.yml`** — `assignment_rate` eşiği.
- **`rnaforge/cli.py`** — `counts` subcommand + `_cmd_counts`.

## 5. Veri akışı ve çıktı yapısı

```
runs/<ts>_<run_id>/
├── quantification/
│   ├── <id>/aligned.sorted.bam        # m04 girdisi (okunur)
│   ├── _featurecounts/counts.txt(.summary)  # ham featureCounts
│   └── counts.tsv                     # TEMIZ gen×örnek matris (Geneid + sample_id sütunlari) = SOZLESME
├── statistics/count_statistics.json   # örnek başına assignment_rate + assigned; n_genes
├── quality/gates.json                 # m05 kapısı EKLENİR
└── logs/counts.log
```

featureCounts komutu: `featureCounts -a <gff> -o _featurecounts/counts.txt -t <feature_type>
-g <attribute> [-p --countReadPairs] -T <threads> <bam1> <bam2> ...`.

## 6. Hata yönetimi (Kural 7)
- **m04 done değil:** `ValueError` ("run m04 (quant) first").
- **featureCounts binary yok / sıfırdan farklı exit / çıktı yok:** `FeatureCountsRunError` (cmd, exit, stderr).
- **counts.txt/summary bozuk/eksik sütun:** `FeatureCountsParseError`.
- **Sıfır gen sayıldı (boş matris):** yüksek sesle — muhtemelen yanlış `feature_type`/`attribute`;
  `run_counts` `n_genes==0` ise net hata (`ValueError`), sessiz boş matris değil.
- **FAIL yazma sırası (m04 deseni):** counts.tsv + count_statistics.json yaz → `write_gate_results`
  → EN SON `raise_if_failed`; `mark_done` yalnız FAIL yoksa.
- **Resume:** m05 done + `count_statistics.json` varsa tekrar koşmaz (`resumed`); `--force`.

## 7. Test stratejisi (TDD)
**Birim (fixture/monkeypatch, env'siz):**
- `parse_counts`: fixture counts.txt → doğru gene_ids + sütun sayımları; `#`/başlık atlanır.
- `parse_summary`: fixture summary → doğru assignment_rate; toplam 0 → 0.0.
- `build_count_gates`: oran < eşik → FAIL + samples/measured; hepsi ≥ → PASS; config ezimi → overridden.
- `config`: `quantification.feature_type/attribute` yüklenir + varsayılan; bilinmeyen anahtar hâlâ reddedilir.
- `run_counts` (featureCounts monkeypatch): counts.tsv (sample_id sütunları) + count_statistics.json
  + gates; öncekiler korunur; m04 done değilse ValueError; n_genes==0 → ValueError; düşük assignment
  → GateFailure; resume featureCounts'u tekrar çağırmaz.
- CLI: `counts` exit 0 (validate→trim→quant→counts), verdict; düşük assignment exit 1.

**Entegrasyon (gerçek featureCounts, `rnaforge-quant-prok` yoksa skip):**
- Sentetik genom + GTF (2 gen, exon/gene_id) + o bölgelerden okumalar → m04 BAM → `run_featurecounts`
  → 2 gen sayılır, assignment yüksek. Parser'ın featureCounts 2.1.1 uyumu.

## 8. Kapsam dışı (sonraki spec'ler)
- **Ökaryot m05 (tximport + tx2gene):** ayrı spec; aynı `counts.tsv` sözleşmesini üretir.
- **m04 ökaryot yolu (Salmon):** ayrı spec.
- m06 DESeq2 (`counts.tsv`'yi tüketir), normalizasyon, TPM/FPKM (featureCounts ham sayım verir).
- MultiQC, çoklu-anotasyon otomatik tespiti (Faz 2+).
