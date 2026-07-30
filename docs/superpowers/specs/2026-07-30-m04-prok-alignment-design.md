# m04 — Quantification ROUTER (prokaryot: bowtie2 hizalama) — Tasarım

**Tarih:** 2026-07-30
**Durum:** Onaylandı (Ali — "devam", otonom oturum)
**İlgili:** PLAN.md v1.3 §4.1, §5, §14; m02/m03 spec'leri (aynı modül deseni)

## 1. Problem

m03 nazik trimlenmiş okumalar üretti. m04 bunları referansa **hizalar** ve kantifikasyona
hazırlar. PLAN §5: m04 bir ROUTER'dır — `organism_type`'a göre dallanır:
- **prokaryote:** genom hizalama (bowtie2) → BAM (m05 featureCounts ile sayar)
- **eukaryote:** Salmon pseudo-alignment → quant.sf (m05 tximport ile gene'e çıkar)

İki yol **gen × örnek count matrisi** sözleşmesinde birleşir (m05). Bu spec **yalnız
prokaryot yolunu (bowtie2 hizalama)** kapsar; ökaryot (Salmon) yolu ayrı sonraki spec.

**Neden prokaryot ilk:** demo veri seti bakteriyel, çekirdek hedef kitle bakteriyel derinlik
([[rnaforge-project]] hedef kitle D). Prokaryotta intron yok → genom hizalama + featureCounts
doğru yoldur (PLAN §4.1: tx2gene ~1:1, Salmon+tximport anlamsız).

## 2. Onaylanan kararlar

1. **Router, tek arm bağlı.** `run_quant` `organism_type`'a bakar; `prokaryote` → bowtie2 yolu.
   `eukaryote` → **net `NotImplementedError`** ("m04 eukaryote (Salmon) path not yet
   implemented; prokaryote only for now"). Sessiz/sahte davranış YOK (Kural 7). Router yapısı
   şimdi kurulur ki euk yolu eklenince yalnız bir dal dolar.
2. **Veri kapısı: `alignment_rate` — FAIL kapısı.** bowtie2 "overall alignment rate" profil
   eşiğinin (`alignment_rate`, prok = 0.70) altındaysa **FAIL → koşu durur** (sonuç GEÇERSİZ:
   düşük hizalama = güvenilmez count matrisi; PLAN §3 "%8 hizalama da makul görünen sahte
   sonuç üretir"). Eşik zaten `profiles/*.yml`'de. m03 survival kapısı deseniyle aynı.
3. **Girdi = m03'ün trimlenmiş okumaları.** m04 ham FASTQ'yu DEĞİL, `trimmed/<id>/`'deki
   trimlenmiş okumaları hizalar. **Ön koşul: m03 done** (m03 zaten m01 ön koşullu → zincir
   m01→m03→m04). m03 yapılmamışsa net `ValueError`. Bu, ilk **modül-çıktısı-tüketen** bağdır.
4. **Trimlenmiş yol sözleşmesi tek kaynaktan.** m03'ün adlandırma kuralı (`<stem>.trimmed.fastq`)
   TEK yerde tanımlanır (`rnaforge/modules/m03_trim.py::trimmed_reads(run_dir, sample)`); hem m03
   yazar hem m04 okur. İkinci bir kopya sessizce sürüklenirdi.
5. **Paired-end destekli.** bowtie2 `-1/-2` (paired) veya `-U` (unpaired); m03 iki trimlenmiş
   dosya ürettiyse paired.

### Reddedilen seçenekler
- **m04 = hizalama + sayım (m05'i içine al):** PLAN modül sınırını bozar (m04=quant, m05=count).
  Ayrı tutmak her modülü bağımsız test edilebilir bırakır; BAM ara ürün, count matrisi m05'in işi.
- **alignment_rate = WARN:** Düşük hizalama sonucu GEÇERSİZ kılar (PLAN §3 tam bu senaryo) → FAIL.
- **Ham okumayı hizala:** m03 trimming'i atlar; adapter kontaminasyonu hizalamayı bozar.
- **eukaryote'u sessiz atla/stub:** Kural 7 ihlali; net NotImplementedError doğrusu.

## 3. Arayüz (public sözleşme)

```python
# rnaforge/modules/m04_quant.py
def run_quant(config: Config, metadata_path: Path, run_dir: Path,
              force: bool = False) -> dict: ...
```
- **Ön koşul:** m03 done değilse `ValueError`. `organism_type == "eukaryote"` → `NotImplementedError`.
- **Dönüş:** özet dict: `n_samples`, örnek başına `alignment_rate`/`bam` yolu, `gate_counts`, `resumed?`.

```python
# rnaforge/bowtie2.py
def parse_bowtie2_summary(stderr_text: str) -> float: ...   # overall alignment rate (0..1)
def build_index(genome_fasta: Path, index_dir: Path, env=...) -> Path: ...  # index prefix
def run_bowtie2(index_prefix: Path, out_dir: Path, fastq_1: Path,
                fastq_2: Path | None = None, threads: int = 4, env=...) -> AlignmentResult: ...
```

## 4. Bileşenler

- **`rnaforge/bowtie2.py`** — saf parser + runnerlar (m02/m03 deseni):
  - `parse_bowtie2_summary(stderr_text) -> float`: "`NN.NN% overall alignment rate`" satırından
    oranı [0,1] float olarak çıkarır. Satır yoksa `Bowtie2ParseError`.
  - `build_index(genome_fasta, index_dir, env) -> Path`: `bowtie2-build` çalıştırır; index
    prefix döner. Idempotent değil — çağıran resume kontrolü yapar.
  - `run_bowtie2(index_prefix, out_dir, fastq_1, fastq_2=None, threads, env) -> AlignmentResult`:
    bowtie2 → stdout SAM'i `samtools sort` ile pipe → `aligned.sorted.bam`; `samtools index` →
    `.bai`; bowtie2 stderr'i `bowtie2.log`'a yazar ve `parse_bowtie2_summary` ile parse eder.
  - `AlignmentResult` frozen dataclass: `bam: Path`, `alignment_rate: float`.
  - `Bowtie2ParseError(ValueError)`, `Bowtie2RunError(RuntimeError)`.
- **`rnaforge/modules/m04_quant.py`** — `run_quant` router + `build_alignment_gates(
  results: dict[str, AlignmentResult], profile: Profile) -> list[GateResult]` (m03 survival
  kapısı deseni: eşik altı örnek → FAIL, `measured`=en düşük, `samples`=suçlular).
- **`rnaforge/modules/m03_trim.py`** — `trimmed_reads(run_dir, sample) -> tuple[Path, Path|None]`
  ekle (adlandırma kuralının tek kaynağı; m03 zaten `<stem>.trimmed.fastq` yazıyor).
- **`rnaforge/cli.py`** — `quant` subcommand + `_cmd_quant`.

## 5. Veri akışı ve çıktı yapısı (PLAN §14)

```
runs/<ts>_<run_id>/
├── quantification/
│   ├── _index/                  # bowtie2 genome index (run başına bir kez)
│   └── <sample_id>/             # aligned.sorted.bam + .bai + bowtie2.log
├── statistics/alignment_statistics.json   # örnek başına alignment_rate + bam yolu
├── quality/gates.json           # m04 kapısı EKLENİR (öncekilere dokunma)
└── logs/quant.log
```

bowtie2 komutu (unpaired): `bowtie2 -x <prefix> -U <trimmed1> -p <threads> -S <tmp.sam>`
(stderr → `bowtie2.log`, alignment rate oradan parse) → `samtools sort -o aligned.sorted.bam
<tmp.sam>` → `samtools index` → tmp.sam sil. Paired: `-1 <t1> -2 <t2>`. Ara SAM kullanılır
(conda-run üzerinden pipe karmaşası yerine); tek `conda run ... bash -c` içinde zincirlenir.

## 6. Hata yönetimi (Kural 7)
- **m03 done değil:** `ValueError` ("run m03 (trim) first").
- **eukaryote:** `NotImplementedError` (net mesaj; sessiz değil).
- **genome_fasta yok / index kurulamadı / bowtie2|samtools sıfırdan farklı exit:** `Bowtie2RunError`
  (env, cmd, exit, stderr). Sahte BAM üretme.
- **bowtie2 özet satırı yok (parse):** `Bowtie2ParseError`.
- **FAIL yazma sırası (m03 deseni):** tüm örnekler hizalanır → BAM + `alignment_statistics.json`
  yazılır → `write_gate_results` → EN SON `raise_if_failed`; `mark_done` yalnız FAIL yoksa.
- **Resume:** m04 done + `alignment_statistics.json` varsa tekrar koşmaz (`resumed`); `--force`.
  Index resume: `quantification/_index/` prefix dosyaları varsa yeniden kurma.

## 7. Test stratejisi (TDD)
**Birim (fixture/monkeypatch, env'siz):**
- `parse_bowtie2_summary`: fixture stderr → doğru oran; "0.00%"/"100.00%"; satır yoksa hata.
- `build_alignment_gates`: oran < eşik → FAIL + `samples`/`measured`; hepsi ≥ → PASS; config ezimi → `overridden`.
- `trimmed_reads`: metadata örneğinden doğru trimlenmiş yol(lar)ı üretir (single + paired).
- `run_quant` (bowtie2 monkeypatch): quantification/<id>/ + alignment_statistics.json + gates;
  öncekiler korunur; m03 done değilse ValueError; eukaryote → NotImplementedError; düşük
  alignment → GateFailure + gates.json FAIL; resume bowtie2'yi tekrar çağırmaz.
- CLI: `quant` exit 0 (validate→trim→quant), verdict basılır; düşük alignment exit 1.

**Entegrasyon (gerçek bowtie2+samtools, `rnaforge-quant-prok` yoksa skip):**
- Sentetik genom + O GENOMDAN türetilmiş okumalar → `build_index` + `run_bowtie2` → BAM üretilir,
  `.bai` var, `alignment_rate > 0.95` (genomdan türetilen okumalar hizalanır). Rastgele (genom-dışı)
  okumalar → düşük alignment (parser + gate'in gerçek araçla uyumu).

## 8. Kapsam dışı (sonraki spec'ler)
- **Ökaryot yolu (Salmon + tximport):** ayrı spec; router'ın euk armı. salmon 2.3.4 doğrulandı (uyumlu).
- **m05 count matrisi (featureCounts + gen×örnek TSV):** ayrı spec; m04 yalnız BAM üretir.
- MultiQC, STAR (bulut), duplicate marking, BAM istatistik raporu (Faz 2+).
