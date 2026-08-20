# Metatranskriptom Referans-Tabanlı Kısa-Okuma (M1) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** RNAForge'a `organism_type=metatranscriptome` üçüncü kolu olarak referans-tabanlı,
kısa-okuma topluluk RNA-seq yolu eklemek (rRNA depletion → taksonomi → katalog-kuant → mevcut DE/enrichment/rapor).

**Architecture:** Yeni `organism_type` değeri m04/m05 router'larına (mevcut prokaryot/ökaryot
dallanma deseni) yeni dallar ekler. Ön uç (m01/m02/m03) ve downstream (m06+) paylaşılır ve
değişmez. Metatranskriptoma-özel: SortMeRNA ile rRNA'yı **çıkaran** işlem adımı, Kraken2/Bracken
taksonomi (diagnostik), gen-kataloğuna Bowtie2 hizalama + featureCounts (mapping/assignment
DİAGNOSTİK, FAIL kapısı yok), permissive+damgalı profil.

**Tech Stack:** Python 3, conda envs, SortMeRNA, Kraken2, Bracken, Bowtie2, subread/featureCounts,
DESeq2 (mevcut), pytest (TDD).

**Spec:** `docs/superpowers/specs/2026-08-20-metatranscriptome-refbased-design.md`

## Global Constraints

- **Aile standardı:** her modül tek `rnaforge <cmd>` subcommand'ı; per-sample atomic state +
  heartbeat + resume (`RunState`, `mark_done`/`is_done`/`mark_item_done`); hata sessiz yutulmaz,
  yüksek sesle yükselir (`feedback_gurultulu_hata`).
- **Kapı felsefesi:** uydurma eşik yok. Metatranskriptomda katalog eksikliği doğal → alignment/
  assignment **DİAGNOSTİK (FAIL YOK)**; yalnız katastrofik durum WARN. Profil `permissive: true`
  ve DAMGALI (`prokaryote_long.yml` deseni).
- **Downstream değişmezliği:** m06 (DESeq2), m07 (figürler), m09–m12 (KEGG/GO/GSEA/REVIGO) **kod
  değişmeden** ortak `counts.tsv` sözleşmesinde buluşur. Bunları değiştiren hiçbir görev yok.
- **Regresyon yok:** prokaryot ve ökaryot (short/long) yolları DEĞİŞMEDEN geçmeli; her görevde
  tam test paketi yeşil kalır (mevcut 586 test).
- **Rapor dürüstlüğü:** kullanılmayan aracı atıflamaz; permissive rapor "geniş toleranslı" damgası taşır.
- **Dosyalama:** çıktılar run_dir altında sözleşme yollarında; raw'a dokunulmaz.

---

### Task 1: Config — `metatranscriptome` organism_type + referans/taxonomy/rrna bölümleri

**Files:**
- Modify: `rnaforge/config.py` (`ORGANISM_TYPES:9`, `REQUIRED_REFERENCE:27`, `Reference:37`,
  `KNOWN_TOP_LEVEL_KEYS:20`, builder `parse_config`)
- Test: `tests/test_config_metatranscriptome.py`

**Interfaces:**
- Produces: `ORGANISM_TYPES` içerir `"metatranscriptome"`;
  `REQUIRED_REFERENCE["metatranscriptome"] == ("gene_catalog_fasta", "catalog_annotation")`;
  `Reference.gene_catalog_fasta: Path|None`, `Reference.catalog_annotation: Path|None`;
  yeni `Taxonomy` dataclass (`kraken2_db: Path|None`, `bracken_read_len: int=100`,
  `bracken_level: str="S"`, `env: str="rnaforge-meta"`); yeni `Rrna` dataclass
  (`db_fasta: Path|None`, `env: str="rnaforge-seqqc"`); `Config.taxonomy: Taxonomy`,
  `Config.rrna: Rrna`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_config_metatranscriptome.py
from pathlib import Path
import pytest
from rnaforge.config import parse_config, ConfigError, ORGANISM_TYPES

def _base(tmp_path):
    cat = tmp_path / "catalog.fa"; cat.write_text(">g1\nACGT\n")
    ann = tmp_path / "catalog.gff"; ann.write_text("g1\t.\tgene\t1\t4\t.\t+\t.\tID=g1\n")
    return cat, ann

def test_metatranscriptome_is_allowed(tmp_path):
    cat, ann = _base(tmp_path)
    cfg = parse_config({
        "organism": "gut community", "organism_type": "metatranscriptome",
        "platform": "illumina",
        "reference": {"gene_catalog_fasta": str(cat), "catalog_annotation": str(ann)},
        "taxonomy": {"kraken2_db": str(tmp_path), "bracken_read_len": 150},
        "rrna": {"db_fasta": str(cat)},
    })
    assert cfg.organism_type == "metatranscriptome"
    assert cfg.reference.gene_catalog_fasta == cat
    assert cfg.taxonomy.bracken_read_len == 150
    assert cfg.taxonomy.bracken_level == "S"
    assert cfg.rrna.db_fasta == cat

def test_metatranscriptome_requires_catalog(tmp_path):
    with pytest.raises(ConfigError, match="gene_catalog_fasta"):
        parse_config({"organism": "x", "organism_type": "metatranscriptome",
                      "platform": "illumina", "reference": {}})

def test_metatranscriptome_in_organism_types():
    assert "metatranscriptome" in ORGANISM_TYPES
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_config_metatranscriptome.py -v`
Expected: FAIL (metatranscriptome not in ORGANISM_TYPES / KeyError REQUIRED_REFERENCE).

- [ ] **Step 3: Write minimal implementation**

`config.py`:
```python
ORGANISM_TYPES = ("prokaryote", "eukaryote", "metatranscriptome")

REQUIRED_REFERENCE = {
    "prokaryote": ("genome_fasta", "annotation_gff"),
    "eukaryote": ("transcriptome_fasta", "tx2gene"),
    "metatranscriptome": ("gene_catalog_fasta", "catalog_annotation"),
}
```
`Reference` dataclass'ına alanlar ekle:
```python
    gene_catalog_fasta: Path | None = None
    catalog_annotation: Path | None = None
```
Yeni dataclass'lar (Basecall yanına):
```python
@dataclass(frozen=True)
class Taxonomy:
    kraken2_db: Path | None = None
    bracken_read_len: int = 100
    bracken_level: str = "S"
    env: str = "rnaforge-meta"

@dataclass(frozen=True)
class Rrna:
    db_fasta: Path | None = None
    env: str = "rnaforge-seqqc"
```
`KNOWN_TOP_LEVEL_KEYS`'e `"taxonomy", "rrna"` ekle. `_build_reference` içinde yeni Path
alanlarını dahil et (mevcut Path-dönüştürme desenini izle). `parse_config`'e ekle:
```python
        taxonomy=_build_taxonomy(_section(raw, "taxonomy")),
        rrna=_build_rrna(_section(raw, "rrna")),
```
`_build_taxonomy`/`_build_rrna` yardımcıları (Path alanları None-güvenli dönüştürür; sayısal/str
alanlar tip-doğrulanır — mevcut `_build_*` desenini izle).

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_config_metatranscriptome.py -v`
Expected: PASS (3 test).

- [ ] **Step 5: Run full suite (regresyon)**

Run: `pytest -q`
Expected: mevcut testler yeşil (586 + 3).

- [ ] **Step 6: Commit**

```bash
git add rnaforge/config.py tests/test_config_metatranscriptome.py
git commit -m "feat(config): metatranscriptome organism_type + taxonomy/rrna bölümleri"
```

---

### Task 2: Kalite profili — `metatranscriptome.yml`

**Files:**
- Create: `rnaforge/profiles/metatranscriptome.yml`
- Test: `tests/test_profile_metatranscriptome.py`

**Interfaces:**
- Consumes: `quality.load_profile("metatranscriptome")`, `quality.profile_name_for`.
- Produces: `permissive: true` profil; eşikler `read_depth`, `base_quality`, `rrna_depletion_rate`,
  `alignment_rate`, `assignment_rate`, `replicate_correlation`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_profile_metatranscriptome.py
from rnaforge.quality import load_profile, profile_name_for

def test_metatranscriptome_profile_loads_permissive():
    p = load_profile("metatranscriptome")
    assert p.permissive is True
    assert p.thresholds["replicate_correlation"] == 0.75
    # alignment düşük eşik (katalog eksikliği doğal): FAIL yerine tolere
    assert p.thresholds["alignment_rate"] <= 0.10

def test_profile_name_for_metatranscriptome_short():
    assert profile_name_for("metatranscriptome", "short") == "metatranscriptome"
```

- [ ] **Step 2: Run to verify it fails**

Run: `pytest tests/test_profile_metatranscriptome.py -v`
Expected: FAIL (no such profile file).

- [ ] **Step 3: Write the profile**

```yaml
# rnaforge/profiles/metatranscriptome.yml
name: metatranscriptome
permissive: true
description: >
  Referans-tabanlı kısa-okuma metatranskriptom (topluluk RNA-seq). Eşikler BİLİNÇLİ olarak
  permissive ve DAMGALIDIR: gen kataloğu doğası gereği eksiktir (topluluğun tümü referansta
  yoktur) → düşük hizalama/atama NORMALDİR, geçersizlik değildir. Yalnız KATASTROFİK (~0)
  hizalama şüpheli (WARN). rRNA depletion verimi düşükse WARN. replicate_correlation korunur.
  Bu profille üretilen rapor "geniş toleranslı" damgası taşır; temsili bir metatranskriptom
  veri setine kalibre edilene dek sıkılaştırılmaz.
thresholds:
  read_depth: 1000000
  base_quality: 20
  rrna_depletion_rate: 0.30
  alignment_rate: 0.05
  assignment_rate: 0.02
  replicate_correlation: 0.75
```

- [ ] **Step 4: Run to verify it passes**

Run: `pytest tests/test_profile_metatranscriptome.py -v`
Expected: PASS.

Not: `quality.load_profile` gate allow-list'i profil anahtarlarından türetiyorsa
`rrna_depletion_rate` yeni bir gate adı olarak tanınmalı; `load_profile`'daki bilinmeyen-gate
kontrolü (`quality.py:72`) profilin kendi anahtarlarını kabul eder, ek kod gerekmez. Eğer
merkezi bir "known gates" kümesi varsa oraya `rrna_depletion_rate` eklenir (grep ile doğrula).

- [ ] **Step 5: Commit**

```bash
git add rnaforge/profiles/metatranscriptome.yml tests/test_profile_metatranscriptome.py
git commit -m "feat(profile): metatranscriptome permissive+damgalı kalite profili"
```

---

### Task 3: rRNA depletion runner — `rrna_deplete.py`

**Files:**
- Create: `rnaforge/rrna_deplete.py`
- Test: `tests/test_rrna_deplete.py`

**Interfaces:**
- Produces: `run_sortmerna_deplete(reads: list[Path], rrna_db: Path, workdir: Path, paired: bool,
  threads: int=8, env: str="rnaforge-seqqc") -> dict` → `{"other": [Path,...] (rRNA'sız FASTQ),
  "depletion_rate": float, "aligned_log": Path}`; `parse_depletion_rate(log_path: Path) -> float`.

- [ ] **Step 1: Write the failing test (parser, saf birim)**

```python
# tests/test_rrna_deplete.py
from pathlib import Path
from rnaforge.rrna_deplete import parse_depletion_rate

def test_parse_depletion_rate(tmp_path):
    log = tmp_path / "aligned.log"
    log.write_text(
        "Total reads = 1000\n"
        "Total reads passing E-value threshold = 120 (12.00)\n")
    # depletion_rate = rRNA fraksiyonu = çıkarılan pay
    assert abs(parse_depletion_rate(log) - 0.12) < 1e-6

def test_parse_depletion_rate_missing(tmp_path):
    log = tmp_path / "x.log"; log.write_text("garbage")
    assert parse_depletion_rate(log) == 0.0
```

- [ ] **Step 2: Run to verify it fails**

Run: `pytest tests/test_rrna_deplete.py -v`
Expected: FAIL (module missing).

- [ ] **Step 3: Write minimal implementation**

```python
# rnaforge/rrna_deplete.py
"""Metatranskriptom rRNA depletion: SortMeRNA ile rRNA okumalarını ÇIKAR (--other).
seqqc.run_sortmerna yalnız ÖLÇER; bu modül rRNA'sız FASTQ üretir (downstream girdisi)."""
import re
import shutil
import subprocess
from pathlib import Path

def parse_depletion_rate(log_path: Path) -> float:
    """aligned.log'tan rRNA fraksiyonu (0-1) = çıkarılan pay. seqqc.parse_sortmerna_log deseni."""
    text = Path(log_path).read_text() if Path(log_path).exists() else ""
    m = re.search(r"passing E-value threshold\s*=\s*\d+\s*\(([\d.]+)%?\)", text)
    if m:
        return float(m.group(1)) / 100.0
    total = re.search(r"Total reads\s*=\s*(\d+)", text)
    passed = re.search(r"passing E-value threshold\s*=\s*(\d+)", text)
    if total and passed and int(total.group(1)) > 0:
        return int(passed.group(1)) / int(total.group(1))
    return 0.0

def run_sortmerna_deplete(reads, rrna_db, workdir, paired, threads=8,
                          env="rnaforge-seqqc") -> dict:
    """SortMeRNA --fastx --other ile rRNA'sız okumaları üret. Hatada gürültülü yüksel."""
    workdir = Path(workdir)
    if workdir.exists():
        shutil.rmtree(workdir)
    workdir.mkdir(parents=True)
    aligned = workdir / "out" / "aligned"
    other = workdir / "out" / "other"
    cmd = ["conda", "run", "-n", env, "sortmerna", "--ref", str(rrna_db),
           "--workdir", str(workdir), "--fastx",
           "--aligned", str(aligned), "--other", str(other),
           "--threads", str(threads), "-v"]
    if paired:
        cmd += ["--paired_in", "--out2"]
    for r in reads:
        cmd += ["--reads", str(r)]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"sortmerna deplete failed (exit {res.returncode}):\n{res.stderr[-2000:]}")
    out_dir = workdir / "out"
    others = sorted(out_dir.glob("other*.f*q*"))
    if not others:
        raise RuntimeError(f"sortmerna produced no --other output in {out_dir}")
    log = out_dir / "aligned.log"
    return {"other": others, "depletion_rate": parse_depletion_rate(log), "aligned_log": log}
```

- [ ] **Step 4: Run to verify it passes**

Run: `pytest tests/test_rrna_deplete.py -v`
Expected: PASS (parser testleri).

- [ ] **Step 5: (Opsiyonel) Gerçek SortMeRNA entegrasyon testi**

Küçük sentetik rRNA + non-rRNA FASTQ ile `run_sortmerna_deplete` çağır; conda/sortmerna yoksa
`pytest.importorskip`/skip. `other` FASTQ'nun oluştuğunu ve depletion_rate'in [0,1] olduğunu doğrula.

- [ ] **Step 6: Commit**

```bash
git add rnaforge/rrna_deplete.py tests/test_rrna_deplete.py
git commit -m "feat(rrna): SortMeRNA --other ile rRNA depletion runner"
```

---

### Task 4: `m_rrna_deplete` modülü + CLI

**Files:**
- Create: `rnaforge/modules/m_rrna_deplete.py`
- Modify: `rnaforge/cli.py` (yeni `rrna-deplete` subcommand), `rnaforge/modules/__init__.py`
- Test: `tests/test_m_rrna_deplete.py`

**Interfaces:**
- Consumes: `rrna_deplete.run_sortmerna_deplete`, `RunState`, `config.rrna`.
- Produces: `run_rrna_deplete(config, metadata_path, run_dir, force=False) -> dict`; çıktı
  `rrna_depleted/<sid>/other_*.fastq.gz`; `statistics/rrna_depletion.json`
  (`{sid: {"depletion_rate": float}}`); state `m_rrna_deplete` + per-item.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_m_rrna_deplete.py — run_sortmerna_deplete monkeypatch'lenir (araç gerektirmez)
import json
from pathlib import Path
from rnaforge.modules.m_rrna_deplete import run_rrna_deplete
# ... küçük config + metadata fixture (mevcut m03 testlerindeki deseni izle);
# monkeypatch run_sortmerna_deplete → sahte other fastq yaz + depletion_rate=0.4 döndür.
def test_rrna_deplete_writes_stats_and_output(meta_cfg, monkeypatch, tmp_path):
    # assert: rrna_depletion.json var, her örnek için depletion_rate; other fastq sözleşme yolunda;
    # state.is_done("m_rrna_deplete"); WARN (0.4>=0.30) verdict'i FAIL yapmaz.
    ...
```

- [ ] **Step 2: Run to verify it fails**

Run: `pytest tests/test_m_rrna_deplete.py -v`
Expected: FAIL (module missing).

- [ ] **Step 3: Write minimal implementation**

`m_rrna_deplete.py` — m03_trim.py'nin per-sample döngü + state + heartbeat desenini izle:
her örnek için `run_sortmerna_deplete` çağır (paired = config.paired/metadata layout);
`other` FASTQ'ları `rrna_depleted/<sid>/` altına gzip'le taşı/yaz; depletion_rate'i topla;
profil `rrna_depletion_rate` WARN kapısı (`build_*_gates` + `warn_only`); `statistics/
rrna_depletion.json` yaz; `state.mark_item_done`/`mark_done`. m01 tamamlanmış olmalı (guard).

- [ ] **Step 4: Wire CLI**

`cli.py`: `rrna-deplete` subparser (`--run-id`, `--force`), `run_rrna_deplete`'i çağırır
(mevcut `seqqc`/`trim` subcommand kayıt desenini izle). `%` içeren help metninde `%%` kullan.

- [ ] **Step 5: Run tests**

Run: `pytest tests/test_m_rrna_deplete.py tests/test_cli_help.py -v`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add rnaforge/modules/m_rrna_deplete.py rnaforge/cli.py rnaforge/modules/__init__.py tests/test_m_rrna_deplete.py
git commit -m "feat(m_rrna_deplete): rRNA depletion modülü + CLI"
```

---

### Task 5: Kraken2 + Bracken runner ve parser'lar — `kraken2.py`

**Files:**
- Create: `rnaforge/kraken2.py`
- Test: `tests/test_kraken2.py`

**Interfaces:**
- Produces: `run_kraken2(reads, db, out_prefix, paired, threads, env) -> Path` (report yolu);
  `parse_kraken2_report(path) -> list[dict]` (`{rank, taxid, name, reads, fraction}`);
  `run_bracken(kraken_report, db, out_path, read_len, level, env) -> Path`;
  `parse_bracken(path) -> dict[str, float]` (taxon adı → fraksiyon).

- [ ] **Step 1: Write the failing test (parser'lar, saf birim)**

```python
# tests/test_kraken2.py
from rnaforge.kraken2 import parse_kraken2_report, parse_bracken

def test_parse_kraken2_report(tmp_path):
    r = tmp_path / "k.report"
    # Kraken2 report: fraction, clade_reads, taxon_reads, rank, taxid, name
    r.write_text(" 50.00\t500\t500\tS\t562\tEscherichia coli\n"
                 " 30.00\t300\t300\tS\t1280\tStaphylococcus aureus\n")
    rows = parse_kraken2_report(r)
    assert rows[0]["name"] == "Escherichia coli"
    assert abs(rows[0]["fraction"] - 0.50) < 1e-6
    assert rows[0]["rank"] == "S"

def test_parse_bracken(tmp_path):
    b = tmp_path / "b.bracken"
    b.write_text("name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\t"
                 "added_reads\tnew_est_reads\tfraction_total_reads\n"
                 "Escherichia coli\t562\tS\t400\t100\t500\t0.55\n")
    d = parse_bracken(b)
    assert abs(d["Escherichia coli"] - 0.55) < 1e-6
```

- [ ] **Step 2: Run to verify it fails**

Run: `pytest tests/test_kraken2.py -v`
Expected: FAIL (module missing).

- [ ] **Step 3: Write minimal implementation**

`kraken2.py`: parser'lar TSV satır-ayrıştırma (yukarıdaki sütun düzeni); `run_kraken2`
`conda run -n <env> kraken2 --db <db> --report <prefix>.report [--paired] r1 r2`; `run_bracken`
`bracken -d <db> -i <report> -o <out> -r <read_len> -l <level>`. Hatada gürültülü yüksel
(`bowtie2.py`/`minimap2.py` runner deseni).

- [ ] **Step 4: Run to verify it passes**

Run: `pytest tests/test_kraken2.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add rnaforge/kraken2.py tests/test_kraken2.py
git commit -m "feat(kraken2): Kraken2+Bracken runner ve parser'lar"
```

---

### Task 6: `m_taxonomy` modülü + CLI

**Files:**
- Create: `rnaforge/modules/m_taxonomy.py`
- Modify: `rnaforge/cli.py` (`taxonomy` subcommand)
- Test: `tests/test_m_taxonomy.py`

**Interfaces:**
- Consumes: `kraken2.run_kraken2/run_bracken/parse_bracken`, `config.taxonomy`.
- Produces: `run_taxonomy(config, metadata_path, run_dir, force=False) -> dict`; çıktı
  `taxonomy/<sid>.bracken` + birleşik `taxonomy/abundance_matrix.tsv` (satır=taxon, sütun=örnek);
  **kapı YOK (diagnostik)**. Girdi: rRNA'sız okumalar (m_rrna_deplete çıktısı; yoksa m03 trimmed).

- [ ] **Step 1: Write the failing test** (run_kraken2/run_bracken monkeypatch'li; abundance matris birleştirme doğru).

- [ ] **Step 2: Run to verify it fails** — `pytest tests/test_m_taxonomy.py -v`.

- [ ] **Step 3: Implement** — per-sample Kraken2→Bracken; `parse_bracken` ile fraksiyonlar;
  tüm örnekleri `abundance_matrix.tsv`'ye birleştir (eksik taxon=0); state + heartbeat; diagnostik.

- [ ] **Step 4: Wire CLI** — `taxonomy` subcommand.

- [ ] **Step 5: Run tests** — `pytest tests/test_m_taxonomy.py tests/test_cli_help.py -v` → PASS.

- [ ] **Step 6: Commit**

```bash
git commit -am "feat(m_taxonomy): Kraken2/Bracken topluluk kompozisyonu modülü + CLI"
```

---

### Task 7: m04 metatranscriptome dalı — `_quant_meta`

**Files:**
- Modify: `rnaforge/modules/m04_quant.py` (`run_quant` router + yeni `_quant_meta`)
- Test: `tests/test_m04_meta.py`

**Interfaces:**
- Consumes: `bowtie2.build_index/run_bowtie2/parse_bowtie2_summary`, `config.reference.gene_catalog_fasta`,
  rRNA'sız okumalar (`rrna_depleted/<sid>/`).
- Produces: BAM `quantification/<sid>/aligned.sorted.bam`; `alignment_statistics.json` içinde
  `"organism_type": "metatranscriptome"`, per-sample `alignment_rate` (DİAGNOSTİK).

- [ ] **Step 1: Write the failing test**

```python
# tests/test_m04_meta.py — küçük katalog fasta + sentetik rRNA'sız fastq; bowtie2 monkeypatch
# (veya gerçek bowtie2 skip'li). assert: BAM sözleşme yolunda; alignment düşük olsa bile
# GateFailure YÜKSELMEZ (permissive metatranscriptome profili, FAIL kapısı yok).
def test_quant_meta_no_fail_gate_on_low_alignment(meta_cfg, monkeypatch):
    ...
```

- [ ] **Step 2: Run to verify it fails** — `pytest tests/test_m04_meta.py -v`.

- [ ] **Step 3: Implement router + `_quant_meta`**

`run_quant` router'a (m04_quant.py:96 civarı, `eukaryote` dalının yanına):
```python
    if config.organism_type == "metatranscriptome":
        summary = _quant_meta(config, metadata_path, run_dir,
                              quant_dir, stats_dir, logs_dir, state, force)
    elif config.organism_type == "eukaryote":
        ...
```
`_quant_meta`: gen kataloğuna Bowtie2 index (`build_index`) + per-sample `run_bowtie2`
(rRNA'sız okumalar); `parse_bowtie2_summary` → alignment_rate. Profil
`metatranscriptome` (permissive); `build_alignment_gates` **warn_only/FAIL yok** — düşük oran
diagnostik. `_quant_short` gövdesini şablon al, referansı `gene_catalog_fasta` yap, girdi
`rrna_depleted/<sid>/` (yoksa trimmed'e düş + yüksek sesle logla). BAM sözleşme yoluna.

- [ ] **Step 4: Run tests + regresyon** — `pytest tests/test_m04_meta.py -q && pytest -q` → PASS,
  prokaryot/ökaryot m04 testleri DEĞİŞMEDEN yeşil.

- [ ] **Step 5: Commit**

```bash
git commit -am "feat(m04): metatranscriptome katalog-hizalama dalı (diagnostik, FAIL yok)"
```

---

### Task 8: m05 metatranscriptome dalı — `_counts_meta`

**Files:**
- Modify: `rnaforge/modules/m05_counts.py` (`run_counts` router + yeni `_counts_meta`)
- Test: `tests/test_m05_meta.py`

**Interfaces:**
- Consumes: `featurecounts.run_featurecounts`, `_write_count_outputs`, `config.reference.catalog_annotation`.
- Produces: `quantification/counts.tsv` (+ tpm/fpkm) — özellik = katalog geni/KO; ortak sözleşme
  (m06 girdisi). `count_statistics.json` içinde `assignment_rate` (DİAGNOSTİK).

- [ ] **Step 1: Write the failing test** — küçük BAM + katalog GFF/SAF; `run_featurecounts`
  monkeypatch veya gerçek subread skip'li; assert counts.tsv üretilir, assignment düşükse
  GateFailure YÜKSELMEZ.

- [ ] **Step 2: Run to verify it fails** — `pytest tests/test_m05_meta.py -v`.

- [ ] **Step 3: Implement router + `_counts_meta`**

`run_counts` router'a (m05_counts.py:79 civarı) `metatranscriptome` dalı. `_counts_meta`:
`run_featurecounts(bams, catalog_annotation, ..., feature_type=...)` → `_write_count_outputs`
(counts/tpm/fpkm). assignment_rate DİAGNOSTİK (permissive profil, FAIL yok). `_counts_short`
gövdesini şablon al; GFF yerine `catalog_annotation`. (KO'ya toplama: katalog anotasyonu KO
taşıyorsa `feature_type`/attribute ile KO düzeyinde; aksi halde gen düzeyi + m10 KEGG köprüsü.)

- [ ] **Step 4: Run tests + regresyon** — `pytest tests/test_m05_meta.py -q && pytest -q` → PASS.

- [ ] **Step 5: Commit**

```bash
git commit -am "feat(m05): metatranscriptome katalog-sayım dalı (featureCounts, diagnostik)"
```

---

### Task 9: Bağımlılık zinciri + `rnaforge run` orkestrasyonu

**Files:**
- Modify: `rnaforge/modules/m04_quant.py` (metatranscriptome'da m04, m_rrna_deplete'i bekler),
  `rnaforge/cli.py` veya orkestratör (`rnaforge run` sırası)
- Test: `tests/test_run_orchestration_meta.py`

**Interfaces:**
- Produces: metatranscriptome için `rnaforge run` sırası: `validate → qc → trim → rrna-deplete →
  taxonomy → quant → counts → de → (report)`. prokaryot/ökaryot sırası DEĞİŞMEZ.

- [ ] **Step 1: Write the failing test** — organism_type=metatranscriptome config ile orkestratör
  aşama listesini üretir; rrna-deplete ve taxonomy quant'tan ÖNCE; prokaryot listesi değişmemiş.

- [ ] **Step 2: Run to verify it fails** — `pytest tests/test_run_orchestration_meta.py -v`.

- [ ] **Step 3: Implement** — orkestratörün aşama seçimini organism_type'a göre genişlet
  (mevcut `rnaforge run` --from/--to mantığını izle). `_quant_meta` girişinde m_rrna_deplete
  state guard'ı (yoksa dürüstçe dur, resume mesajı).

- [ ] **Step 4: Run tests + regresyon** — `pytest -q` → PASS.

- [ ] **Step 5: Commit**

```bash
git commit -am "feat(run): metatranscriptome orkestrasyon sırası (rrna-deplete+taxonomy)"
```

---

### Task 10: Rapor — taksonomi bölümü + read/organism rozeti + yöntem/atıf

**Files:**
- Modify: `rnaforge/report_html.py` (`section_dataset` rozet, yeni `section_taxonomy`,
  `_SOFTWARE`, `_METHODS_TEXT_META`, `_REFERENCES_META`, `render_report` akışı)
- Test: `tests/test_report_meta.py`

**Interfaces:**
- Consumes: `taxonomy/abundance_matrix.tsv`, `statistics/rrna_depletion.json`.
- Produces: rapora "Topluluk Kompozisyonu (Taksonomi)" bölümü (Bracken top-N + figür),
  `organism_type=metatranscriptome` rozeti, TR+EN metatranskriptom yöntem metni, SortMeRNA/
  Kraken2/Bracken atıfları. Kullanılmayan aracı atıflamaz.

- [ ] **Step 1: Write the failing test** — metatranscriptome verisiyle render; HTML'de
  "Topluluk Kompozisyonu" bölümü + Kraken2/Bracken atıfı var; prokaryot raporunda bu bölüm YOK;
  permissive damgası mevcut.

- [ ] **Step 2: Run to verify it fails** — `pytest tests/test_report_meta.py -v`.

- [ ] **Step 3: Implement** — `section_taxonomy` (abundance_matrix'ten top-N tür tablosu + best-effort
  bar figürü, `qcplots` deseni); `render_report` organism_type=metatranscriptome'da bu bölümü ekler;
  `_METHODS_TEXT_META`/`_REFERENCES_META` (SortMeRNA Kopylova 2012, Kraken2 Wood 2019, Bracken Lu 2017);
  software tablosuna kraken2/bracken/sortmerna cond'lu satırlar. read_type rozeti mevcut deseni izler.

- [ ] **Step 4: Run tests + regresyon** — `pytest -q` → PASS.

- [ ] **Step 5: Commit**

```bash
git commit -am "feat(report): metatranscriptome taksonomi bölümü + yöntem/atıf"
```

---

### Task 11: Ortam (env) + referans hazırlığı

**Files:**
- Create: `envs/rnaforge-meta.yml` (kraken2, bracken)
- Modify: `prepare_references.sh` (Kraken2 DB + Bracken + gen kataloğu + rRNA DB, parametreli),
  `install.sh` (yeni env), `rnaforge doctor` (yeni env kontrolü), `README` (EN+TR)
- Test: `tests/test_prepare_references_meta.py` (script'in metatranscriptome bayraklarını
  tanıdığını/yardım metnini doğrulayan hafif test; ağır indirmeler CI-dışı)

**Interfaces:**
- Produces: `rnaforge-meta` env; `prepare_references.sh --kraken2-db <name> --gene-catalog <url>
  --rrna-db <url>` seçenekleri; `.sha256` + gitignore.

- [ ] **Step 1:** `envs/rnaforge-meta.yml` yaz (bioconda kraken2 + bracken; sürümler pinli).
- [ ] **Step 2:** `prepare_references.sh`'e metatranscriptome bölümü (Kraken2 DB indir/derle,
  `bracken-build`, gen kataloğu + KO anotasyonu, rRNA DB); parametreli; `.sha256`.
- [ ] **Step 3:** `install.sh` + `rnaforge doctor`'a `rnaforge-meta` ekle.
- [ ] **Step 4:** README (EN+TR) metatranscriptome kolu + komutları belgele.
- [ ] **Step 5:** Testler + regresyon → `pytest -q` PASS.
- [ ] **Step 6: Commit**

```bash
git add envs/rnaforge-meta.yml prepare_references.sh install.sh README* tests/test_prepare_references_meta.py
git commit -m "chore(meta): rnaforge-meta env + referans hazırlığı + doktor/README"
```

---

### Task 12: Uçtan uca smoke (sentetik topluluk)

**Files:**
- Create: `tests/test_e2e_meta_smoke.py`
- Test: kendisi

**Interfaces:**
- Consumes: tüm metatranscriptome zinciri.
- Produces: sentetik küçük topluluk (2 koşul × replika, bilinen sinyal geni) → `rnaforge run`
  (rrna→tax→quant→counts→de→report); conda/araç yoksa skip.

- [ ] **Step 1:** Sentetik katalog (birkaç gen, biri koşullar arası farklı ekspresyonlu) + rRNA DB
  + küçük Kraken2 DB (veya kraken2/bracken monkeypatch); rRNA'lı+non-rRNA okumalar üret.
- [ ] **Step 2:** `rnaforge run --to report` çağır; assert: counts.tsv üretildi, m06 DE en az 1
  anlamlı gen, report.html üretildi + taksonomi bölümü + permissive damga; FAIL kapısı yok.
- [ ] **Step 3:** conda yoksa test skip (mevcut `test_e2e_smoke.py` deseni).
- [ ] **Step 4:** Run + regresyon → `pytest -q` PASS.
- [ ] **Step 5: Commit**

```bash
git add tests/test_e2e_meta_smoke.py
git commit -m "test(meta): uçtan uca sentetik topluluk smoke testi"
```

---

## Self-Review (spec kapsamı)

- **Spec §2 üçüncü organism_type** → Task 1 (config), Task 7/8 (router dalları). ✓
- **Spec §4.1 rRNA depletion işlem** → Task 3 (runner), Task 4 (modül+CLI). ✓
- **Spec §4.2 taksonomi** → Task 5 (kraken2/bracken), Task 6 (modül+CLI). ✓
- **Spec §4.3 katalog-kuant** → Task 7 (m04), Task 8 (m05). ✓
- **Spec §4.4 downstream değişmez** → hiçbir görev m06/m09-m12'yi değiştirmez; Task 10 yalnız
  rapora bölüm ekler. ✓
- **Spec §5 profil** → Task 2. ✓
- **Spec §6 referans/config** → Task 1 (config), Task 11 (prepare_references). ✓
- **Spec §7 CLI/orkestrasyon** → Task 4/6 (subcommand'lar), Task 9 (run sırası). ✓
- **Spec §8 test** → her görev TDD; Task 12 e2e. ✓
- **Spec §9 biyolojik doğrulama** → plan sonrası ayrı bilimsel adım (kod-dışı; DURUM'a not). ✓
- **Placeholder taraması:** kod adımları gerçek kod içerir; §"..." işaretli yerler mevcut desen
  referansıdır (m03/m05 gövdesi şablon), yeni tip/isim tanımlamaz. Tip tutarlılığı: `run_sortmerna_deplete`,
  `run_kraken2/run_bracken/parse_bracken`, `_quant_meta/_counts_meta`, `run_taxonomy/run_rrna_deplete`
  görevler arası tutarlı. ✓

## Not (biyolojik doğrulama — kod sonrası)
Kod tamamlanınca aile disiplini: 2 paralel ajan (ENA+literatür) yayınlanmış bilinen-DE'li kısa-okuma
metatranskriptom seti seçer (aday iHMP/HMP2 IBD bağırsak). Uçtan uca + konkordans. Veri seçimi o
aşamada birlikte; kör indirme yok. DURUM.md + bellek güncellenir.
