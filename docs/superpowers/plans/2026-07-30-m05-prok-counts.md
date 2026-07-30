# m05 Prokaryot Count Matrix (featureCounts) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** m04'ün BAM'lerini anotasyona göre sayıp `quantification/counts.tsv` gen×örnek matrisini (ortak sözleşme) üreten, `assignment_rate` FAIL kapısıyla denetleyen prokaryot `m05` yolu + `rnaforge counts` CLI komutu.

**Architecture:** m04 desenini izler. Saf parserlar (`parse_counts`/`parse_summary`) featureCounts çıktısını okur; `run_featurecounts` tüm BAM'lere tek çağrı (native matris) yapar; `run_counts` sütunları KONUMLA sample_id'ye eşler, temiz `counts.tsv` yazar, stats→gates→raise sırasıyla kapı üretir. featureCounts parametreleri (`feature_type`/`attribute`) config-driven.

**Tech Stack:** Python 3.11, pytest, conda (`rnaforge-quant-prok`: featureCounts 2.1.1), mevcut `rnaforge.{gates,quality,state,metadata,config}`.

## Global Constraints
- **Test komutu:** `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest -q`.
- **Sessiz devam YOK (Kural 7);** gerçek/müşteri verisi ASLA (Kural 8, sentetik fixture).
- **assignment_rate FAIL kapısı:** oran < `profile.threshold("assignment_rate")` (prok 0.50) → FAIL.
- **Girdi = m04 BAM'ler** (`quantification/<id>/aligned.sorted.bam`). **Ön koşul m04 done.**
- **featureCounts TEK çağrı, tüm BAM'ler;** sütun→sample_id KONUMLA (isimle değil).
- **n_genes==0 → net ValueError** (yanlış feature_type/attribute; sessiz boş matris değil).
- **Yazma sırası (m04 deseni):** counts.tsv + count_statistics.json → `write_gate_results` → `raise_if_failed`; `mark_done` yalnız FAIL yoksa.
- **Config sertleştirme:** `quantification` `KNOWN_TOP_LEVEL_KEYS`'e EKLENMELİ, aksi halde reddedilir.
- **Env:** `conda run -n rnaforge-quant-prok featureCounts ...`. Yoksa entegrasyon testi skip.

**Mevcut API imzaları (referans):**
- `rnaforge.gates`: `GateResult(...)`, `PASS/FAIL`, `write_gate_results`, `raise_if_failed`, `GateFailure`.
- `rnaforge.quality`: `load_profile`, `Profile.threshold`, `.overrides()`, `.name`.
- `rnaforge.state`: `RunState`, `.is_done`, `.mark_done`, `.heartbeat`; `resolve_run_dir`.
- `rnaforge.metadata`: `load_metadata -> list[Sample]`; `Sample(sample_id, condition, fastq_1, fastq_2=None, ...)`.
- `rnaforge.config`: `Config.reference.annotation_gff: Path`, `Config.resources.threads`, `Config.organism_type`, `Config.quality`. YENİ: `Config.quantification.feature_type/attribute`.
- featureCounts `counts.txt`: `#` yorumlar, sonra `Geneid Chr Start End Strand Length <bam...>` başlığı, sonra gen satırları (sayımlar 6. sütundan). `.summary`: `Status <bam...>`, `Assigned` + `Unassigned_*` satırları.

---

### Task 1: featureCounts saf parserlar (`parse_counts`, `parse_summary`)

**Files:**
- Create: `rnaforge/featurecounts.py`
- Test: `tests/test_featurecounts.py`

**Interfaces:**
- Produces:
  - `FeatureCountsResult` frozen dataclass: `gene_ids: list[str]`, `counts: dict[str, list[int]]` (sütun→sayımlar), `assignment_rates: dict[str, float]` (sütun→oran).
  - `parse_counts(counts_text: str) -> tuple[list[str], dict[str, list[int]]]`.
  - `parse_summary(summary_text: str) -> dict[str, float]`.
  - `FeatureCountsParseError(ValueError)`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_featurecounts.py
from __future__ import annotations

import pytest

from rnaforge.featurecounts import (
    FeatureCountsParseError, FeatureCountsResult, parse_counts, parse_summary,
)

_COUNTS = """# Program:featureCounts
Geneid\tChr\tStart\tEnd\tStrand\tLength\ts1.bam\ts2.bam
geneA\tchr1\t101\t1100\t+\t1000\t150\t80
geneB\tchr1\t2101\t3100\t+\t1000\t150\t60
"""

_SUMMARY = """Status\ts1.bam\ts2.bam
Assigned\t300\t140
Unassigned_Unmapped\t0\t0
Unassigned_NoFeatures\t0\t60
"""


def test_parse_counts_reads_genes_and_columns():
    genes, counts = parse_counts(_COUNTS)
    assert genes == ["geneA", "geneB"]
    assert counts["s1.bam"] == [150, 150]
    assert counts["s2.bam"] == [80, 60]
    assert list(counts.keys()) == ["s1.bam", "s2.bam"]   # insertion order = BAM sirasi


def test_parse_counts_rejects_missing_header():
    with pytest.raises(FeatureCountsParseError, match="Geneid"):
        parse_counts("# only a comment\n")


def test_parse_summary_computes_assignment_rate():
    rates = parse_summary(_SUMMARY)
    assert rates["s1.bam"] == pytest.approx(1.0)          # 300/300
    assert rates["s2.bam"] == pytest.approx(140 / 200)    # 140/(140+0+60)


def test_parse_summary_zero_total_is_zero():
    assert parse_summary("Status\tx.bam\nAssigned\t0\n")["x.bam"] == 0.0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_featurecounts.py -q`
Expected: FAIL (ModuleNotFoundError: rnaforge.featurecounts)

- [ ] **Step 3: Write minimal implementation**

```python
# rnaforge/featurecounts.py
"""featureCounts çıktısını parse eder ve çalıştırır. Parserlar saftır."""
from __future__ import annotations

from dataclasses import dataclass


class FeatureCountsParseError(ValueError):
    """featureCounts çıktısı beklenen biçimde değil."""


@dataclass(frozen=True)
class FeatureCountsResult:
    gene_ids: list[str]
    counts: dict[str, list[int]]           # sütun (BAM) -> sayımlar
    assignment_rates: dict[str, float]     # sütun (BAM) -> atama oranı


def parse_counts(counts_text: str) -> tuple[list[str], dict[str, list[int]]]:
    header = None
    gene_ids: list[str] = []
    columns: list[str] = []
    counts: dict[str, list[int]] = {}
    for line in counts_text.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        fields = line.split("\t")
        if header is None:
            if fields[0] != "Geneid":
                raise FeatureCountsParseError(
                    f"featureCounts counts file has no 'Geneid' header (got {fields[0]!r})"
                )
            header = fields
            columns = fields[6:]              # Geneid Chr Start End Strand Length <bam...>
            counts = {c: [] for c in columns}
            continue
        gene_ids.append(fields[0])
        for col, value in zip(columns, fields[6:]):
            counts[col].append(int(value))
    if header is None:
        raise FeatureCountsParseError("featureCounts counts file has no 'Geneid' header line")
    return gene_ids, counts


def parse_summary(summary_text: str) -> dict[str, float]:
    lines = [ln for ln in summary_text.splitlines() if ln.strip()]
    if not lines or not lines[0].startswith("Status"):
        raise FeatureCountsParseError("featureCounts summary has no 'Status' header")
    columns = lines[0].split("\t")[1:]
    assigned = {c: 0 for c in columns}
    totals = {c: 0 for c in columns}
    for line in lines[1:]:
        fields = line.split("\t")
        status = fields[0]
        for col, value in zip(columns, fields[1:]):
            v = int(value)
            totals[col] += v
            if status == "Assigned":
                assigned[col] += v
    return {c: (assigned[c] / totals[c] if totals[c] > 0 else 0.0) for c in columns}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_featurecounts.py -q`
Expected: PASS (4 passed)

- [ ] **Step 5: Commit**

```bash
git add rnaforge/featurecounts.py tests/test_featurecounts.py
git commit -m "feat(m05): featureCounts saf parserlar (counts + assignment_rate summary)"
```

---

### Task 2: gerçek featureCounts runner (`run_featurecounts`)

**Files:**
- Modify: `rnaforge/featurecounts.py`
- Test: `tests/test_featurecounts.py` (ekle)

**Interfaces:**
- Consumes: `parse_counts`, `parse_summary`, `FeatureCountsResult` (Task 1).
- Produces:
  - `run_featurecounts(bams: list[Path], gff: Path, out_dir: Path, feature_type: str, attribute: str, paired: bool = False, threads: int = 4, env: str = "rnaforge-quant-prok") -> FeatureCountsResult`.
  - `FeatureCountsRunError(RuntimeError)`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_featurecounts.py (ekle)
import shutil
from pathlib import Path

from rnaforge.bowtie2 import Bowtie2RunError
from rnaforge.featurecounts import run_featurecounts, FeatureCountsRunError


def _genome_gtf_bam(tmp_path):
    """Sentetik genom + 2-gen GTF + o bölgelerden okuma → BAM (gerçek bowtie2+samtools)."""
    import random
    from rnaforge.bowtie2 import build_index, run_bowtie2
    random.seed(8)
    genome = "".join(random.choice("ACGT") for _ in range(6000))
    (tmp_path / "genome.fa").write_text(">chr1\n" + genome + "\n")
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(
        'chr1\tsrc\texon\t101\t1100\t.\t+\t.\tgene_id "geneA";\n'
        'chr1\tsrc\texon\t2101\t3100\t.\t+\t.\tgene_id "geneB";\n'
    )
    reads = tmp_path / "reads.fastq"
    with reads.open("w") as f:
        for i in range(300):
            lo, hi = (101, 1100) if i % 2 == 0 else (2101, 3100)
            p = random.randint(lo - 1, hi - 100)
            f.write(f"@r{i}\n{genome[p:p+100]}\n+\n{'I'*100}\n")
    prefix = build_index(tmp_path / "genome.fa", tmp_path / "idx")
    result = run_bowtie2(prefix, tmp_path / "aln", reads)
    return gtf, result.bam


@pytest.mark.skipif(shutil.which("conda") is None, reason="conda yok")
def test_run_featurecounts_counts_genes(tmp_path):
    """Entegrasyon: gerçek featureCounts 2 geni sayar, yüksek atama. Env yoksa skip."""
    try:
        gtf, bam = _genome_gtf_bam(tmp_path)
        result = run_featurecounts([bam], gtf, tmp_path / "fc",
                                   feature_type="exon", attribute="gene_id")
    except (FeatureCountsRunError, Bowtie2RunError) as exc:
        pytest.skip(f"featureCounts/bowtie2 çalıştırılamadı: {exc}")
    assert set(result.gene_ids) == {"geneA", "geneB"}
    col = list(result.counts.keys())[0]
    assert sum(result.counts[col]) > 0
    assert result.assignment_rates[col] > 0.9
```

- [ ] **Step 2: Run test to verify it fails**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_featurecounts.py -q -k run_featurecounts`
Expected: FAIL (ImportError: run_featurecounts)

- [ ] **Step 3: Write minimal implementation**

```python
# rnaforge/featurecounts.py (ekle — importlara ekle)
import subprocess
from pathlib import Path


class FeatureCountsRunError(RuntimeError):
    """featureCounts çalıştırılamadı ya da beklenen çıktıyı üretmedi."""


def run_featurecounts(bams: list[Path], gff: Path, out_dir: Path, feature_type: str,
                      attribute: str, paired: bool = False, threads: int = 4,
                      env: str = "rnaforge-quant-prok") -> FeatureCountsResult:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    counts_path = out_dir / "counts.txt"
    cmd = ["conda", "run", "-n", env, "featureCounts",
           "-a", str(gff), "-o", str(counts_path),
           "-t", feature_type, "-g", attribute, "-T", str(threads)]
    if paired:
        cmd += ["-p", "--countReadPairs"]
    cmd += [str(b) for b in bams]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise FeatureCountsRunError(
            f"featureCounts failed (exit {r.returncode})\ncmd: {' '.join(cmd)}\n"
            f"stderr: {r.stderr.strip()}"
        )
    summary_path = counts_path.with_name(counts_path.name + ".summary")
    if not counts_path.exists() or not summary_path.exists():
        raise FeatureCountsRunError(
            f"featureCounts reported success but output missing at {counts_path}"
        )
    gene_ids, counts = parse_counts(counts_path.read_text())
    rates = parse_summary(summary_path.read_text())
    return FeatureCountsResult(gene_ids=gene_ids, counts=counts, assignment_rates=rates)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_featurecounts.py -q`
Expected: PASS (entegrasyon ya geçer ya skip)

- [ ] **Step 5: Commit**

```bash
git add rnaforge/featurecounts.py tests/test_featurecounts.py
git commit -m "feat(m05): gercek featureCounts runner (cok-ornekli, entegrasyon testli)"
```

---

### Task 3: config `quantification` bölümü + profil `assignment_rate`

**Files:**
- Modify: `rnaforge/config.py`
- Modify: `rnaforge/profiles/prokaryote.yml`, `rnaforge/profiles/eukaryote.yml`
- Test: `tests/test_config.py` (ekle)

**Interfaces:**
- Produces: `Config.quantification: Quantification` (`feature_type: str = "exon"`, `attribute: str = "gene_id"`).

- [ ] **Step 1: Write the failing test**

```python
# tests/test_config.py (ekle)
def test_quantification_defaults(tmp_path):
    cfg = load_config(_write(tmp_path, PROK_BODY))
    assert cfg.quantification.feature_type == "exon"
    assert cfg.quantification.attribute == "gene_id"


def test_quantification_overrides(tmp_path):
    cfg = load_config(_write(tmp_path, PROK_BODY +
                             "\nquantification:\n  feature_type: CDS\n  attribute: locus_tag\n"))
    assert cfg.quantification.feature_type == "CDS"
    assert cfg.quantification.attribute == "locus_tag"


def test_quantification_is_known_top_level_key(tmp_path):
    # 'quantification' reddedilmemeli (KNOWN_TOP_LEVEL_KEYS'te)
    load_config(_write(tmp_path, PROK_BODY + "\nquantification:\n  feature_type: gene\n"))
```

- [ ] **Step 2: Run test to verify it fails**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_config.py -q -k quantification`
Expected: FAIL (AttributeError: Config has no 'quantification' / ConfigError: unknown key)

- [ ] **Step 3: Write minimal implementation**

`rnaforge/config.py`:

1. `KNOWN_TOP_LEVEL_KEYS`'e `"quantification"` ekle (frozenset içine).
2. Yeni dataclass (diğer `@dataclass(frozen=True)` bloklarının yanına):
```python
@dataclass(frozen=True)
class Quantification:
    feature_type: str = "exon"
    attribute: str = "gene_id"
```
3. `Config` dataclass'ına alan ekle (defaults'lu alanların yanına, sona):
```python
    quantification: Quantification = field(default_factory=Quantification)
```
4. `load_config` içinde section'ı parse et (`resources_raw`'dan sonra):
```python
    quantification_raw = _section(raw, "quantification")
```
   ve `Config(...)` çağrısında `quality=...`'dan sonra ekle:
```python
        quantification=Quantification(
            feature_type=str(quantification_raw.get("feature_type", "exon")),
            attribute=str(quantification_raw.get("attribute", "gene_id")),
        ),
```

`rnaforge/profiles/prokaryote.yml` ve `eukaryote.yml`: `thresholds:` altına ekle:
```yaml
  assignment_rate: 0.50
```

- [ ] **Step 4: Run test to verify it passes**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_config.py -q`
Expected: PASS (quantification testleri + mevcut config testleri yeşil)

- [ ] **Step 5: Commit**

```bash
git add rnaforge/config.py rnaforge/profiles/prokaryote.yml rnaforge/profiles/eukaryote.yml tests/test_config.py
git commit -m "feat(m05): config quantification bolumu (feature_type/attribute) + profil assignment_rate"
```

---

### Task 4: assignment_rate kapısı (`build_count_gates`)

**Files:**
- Create: `rnaforge/modules/m05_counts.py`
- Test: `tests/test_m05_counts.py`

**Interfaces:**
- Produces: `MODULE_NAME = "m05_counts"`; `build_count_gates(assignment_rates: dict[str, float], profile: Profile) -> list[GateResult]` (anahtar = sample_id).

- [ ] **Step 1: Write the failing test**

```python
# tests/test_m05_counts.py
from __future__ import annotations

from rnaforge.gates import FAIL, PASS
from rnaforge.modules.m05_counts import build_count_gates
from rnaforge.quality import load_profile


def test_all_above_threshold_passes():
    profile = load_profile("prokaryote")  # assignment_rate = 0.50
    gates = build_count_gates({"s1": 0.95, "s2": 0.80}, profile)
    assert len(gates) == 1
    assert gates[0].name == "assignment_rate"
    assert gates[0].module == "m05_counts"
    assert gates[0].status == PASS


def test_below_threshold_fails():
    profile = load_profile("prokaryote")
    gates = build_count_gates({"s1": 0.95, "s2": 0.10}, profile)
    g = gates[0]
    assert g.status == FAIL
    assert g.samples == ("s2",)
    assert g.measured == 0.10
    assert g.threshold == 0.50


def test_override_marks_overridden():
    profile = load_profile("prokaryote", {"assignment_rate": 0.05})
    gates = build_count_gates({"s1": 0.10}, profile)
    assert gates[0].status == PASS
    assert gates[0].overridden is True
```

- [ ] **Step 2: Run test to verify it fails**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_m05_counts.py -q`
Expected: FAIL (ModuleNotFoundError: rnaforge.modules.m05_counts)

- [ ] **Step 3: Write minimal implementation**

```python
# rnaforge/modules/m05_counts.py
"""m05 — Count Matrix (prokaryot: featureCounts).

m04 BAM'lerini anotasyona göre sayıp gen×örnek count matrisi (ortak sözleşme,
PLAN §5) üretir. Veri kapısı `assignment_rate`: featureCounts'un gene atadığı
okuma oranı profil eşiğinin altındaysa FAIL — çok düşük atama yanlış anotasyon/
tür demektir, sayımlar güvenilmez."""
from __future__ import annotations

from rnaforge.gates import FAIL, PASS, GateResult
from rnaforge.quality import Profile

MODULE_NAME = "m05_counts"
_GATE = "assignment_rate"


def build_count_gates(assignment_rates: dict[str, float],
                      profile: Profile) -> list[GateResult]:
    threshold = profile.threshold(_GATE)
    offenders = sorted(sid for sid, r in assignment_rates.items() if r < threshold)
    lowest = min(assignment_rates.values(), default=1.0)
    overridden = _GATE in profile.overrides()
    if offenders:
        status = FAIL
        message = (
            f"gene atama oranı eşiğin altında ({len(offenders)} örnek: "
            f"{', '.join(offenders)}); en düşük {lowest:.2f} < {threshold:.2f}. "
            "Düşük atama yanlış anotasyon/tür → güvenilmez sayımlar."
        )
    else:
        status = PASS
        message = f"tüm örnekler assignment ≥ {threshold:.2f} (en düşük {lowest:.2f})."
    return [GateResult(
        name=_GATE, module=MODULE_NAME, status=status, message=message,
        remedy=("Anotasyon (GFF/GTF) ile referans genomun eşleştiğini ve feature_type/"
                "attribute config'inin anotasyon formatına uyduğunu doğrulayın."),
        measured=lowest, threshold=threshold, overridden=overridden,
        samples=tuple(offenders),
    )]
```

- [ ] **Step 4: Run test to verify it passes**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_m05_counts.py -q`
Expected: PASS (3 passed)

- [ ] **Step 5: Commit**

```bash
git add rnaforge/modules/m05_counts.py tests/test_m05_counts.py
git commit -m "feat(m05): assignment_rate FAIL kapisi (profil esigine karsi)"
```

---

### Task 5: `run_counts` orkestrasyonu (+ counts.tsv sözleşmesi)

**Files:**
- Modify: `rnaforge/modules/m05_counts.py`
- Test: `tests/test_m05_counts.py` (ekle)

**Interfaces:**
- Consumes: `run_featurecounts` (Task 2), `build_count_gates` (Task 4), `RunState`, `write_gate_results`, `raise_if_failed`, `load_metadata`, `load_profile`, `Config`.
- Produces: `run_counts(config: Config, metadata_path: Path, run_dir: Path, force: bool = False) -> dict`. Yan çıktı `quantification/counts.tsv`: satır1 `gene\t<sample_id...>`, sonra gen×sayım.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_m05_counts.py (ekle)
import json
import textwrap
from pathlib import Path

import pytest

from rnaforge.config import load_config
from rnaforge.featurecounts import FeatureCountsResult
from rnaforge.gates import GateFailure
from rnaforge.modules import m05_counts
from rnaforge.modules.m05_counts import run_counts
from tests.conftest import write_fastq


def _setup(tmp_path):
    (tmp_path / "ref").mkdir()
    (tmp_path / "ref" / "genome.fa").write_text(">c1\n" + "ACGT" * 25 + "\n")
    (tmp_path / "ref" / "genes.gtf").write_text('c1\ts\texon\t1\t80\t.\t+\t.\tgene_id "g1";\n')
    for n in ("c1.fastq", "c2.fastq", "t1.fastq", "t2.fastq"):
        write_fastq(tmp_path / n, 200, 150, "I")
    config_path = tmp_path / "config.yaml"
    config_path.write_text(textwrap.dedent(f"""
        organism: "E. coli"
        organism_type: "prokaryote"
        reference:
          genome_fasta: "{tmp_path / 'ref' / 'genome.fa'}"
          annotation_gff: "{tmp_path / 'ref' / 'genes.gtf'}"
    """))
    metadata_path = tmp_path / "samples.tsv"
    metadata_path.write_text(
        "sample_id\tcondition\tfastq_1\n"
        "s1\tcontrol\tc1.fastq\n" "s2\tcontrol\tc2.fastq\n"
        "s3\ttreated\tt1.fastq\n" "s4\ttreated\tt2.fastq\n"
    )
    return config_path, metadata_path


def _prep_through_m04(config_path, metadata_path, run_dir, monkeypatch):
    """m01(gerçek)+m03(fake fastp)+m04(fake bowtie2) done state + BAM üretir."""
    from rnaforge.modules.m01_validate import run_validation
    from rnaforge.modules import m03_trim, m04_quant
    from rnaforge.modules.m03_trim import run_trim, trimmed_name
    from rnaforge.modules.m04_quant import run_quant
    from rnaforge.fastp import FastpResult
    from rnaforge.bowtie2 import AlignmentResult
    run_validation(load_config(config_path), metadata_path, run_dir)

    def fake_fastp(fastq_1, out_dir, min_length, fastq_2=None, aggressive_quality=False, env="rnaforge-qc"):
        out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
        out1 = out_dir / trimmed_name(Path(fastq_1)); out1.write_text("@r\nACGT\n+\nIIII\n")
        (out_dir / "fastp.json").write_text("{}")
        return FastpResult(200, 196, 0.98, out1=out1)
    monkeypatch.setattr(m03_trim, "run_fastp", fake_fastp)
    run_trim(load_config(config_path), metadata_path, run_dir)

    monkeypatch.setattr(m04_quant, "build_index", lambda g, i, env="rnaforge-quant-prok": Path(i) / "genome")
    def fake_align(index_prefix, out_dir, fastq_1, fastq_2=None, threads=4, env="rnaforge-quant-prok"):
        out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
        bam = out_dir / "aligned.sorted.bam"; bam.write_bytes(b"BAM")
        return AlignmentResult(bam=bam, alignment_rate=0.95)
    monkeypatch.setattr(m04_quant, "run_bowtie2", fake_align)
    run_quant(load_config(config_path), metadata_path, run_dir)


def _fake_featurecounts(monkeypatch, rate=0.9, n_genes=3):
    def fake_run(bams, gff, out_dir, feature_type, attribute, paired=False, threads=4, env="rnaforge-quant-prok"):
        Path(out_dir).mkdir(parents=True, exist_ok=True)
        genes = [f"g{i}" for i in range(n_genes)]
        # sütunlar BAM sirasinda; her BAM icin n_genes sayim
        counts = {str(b): [10 + i for i in range(n_genes)] for b in bams}
        rates = {str(b): rate for b in bams}
        return FeatureCountsResult(gene_ids=genes, counts=counts, assignment_rates=rates)
    monkeypatch.setattr(m05_counts, "run_featurecounts", fake_run)


def test_run_counts_requires_m04_done(tmp_path, monkeypatch):
    _fake_featurecounts(monkeypatch)
    config_path, metadata_path = _setup(tmp_path)
    with pytest.raises(ValueError, match="m04"):
        run_counts(load_config(config_path), metadata_path, tmp_path / "run")


def test_run_counts_writes_matrix_and_passes(tmp_path, monkeypatch):
    config_path, metadata_path = _setup(tmp_path)
    run_dir = tmp_path / "run"
    _prep_through_m04(config_path, metadata_path, run_dir, monkeypatch)
    _fake_featurecounts(monkeypatch, rate=0.9, n_genes=3)
    summary = run_counts(load_config(config_path), metadata_path, run_dir)

    assert summary["n_samples"] == 4
    assert summary["n_genes"] == 3
    matrix = (run_dir / "quantification" / "counts.tsv").read_text().splitlines()
    assert matrix[0] == "gene\ts1\ts2\ts3\ts4"        # sample_id basliklari (BAM yollari degil)
    assert matrix[1].split("\t")[0] == "g0"
    gates = json.loads((run_dir / "quality" / "gates.json").read_text())["gates"]
    assert any(g["module"] == "m05_counts" and g["status"] == "PASS" for g in gates)
    assert any(g["module"] == "m04_quant" for g in gates)   # onceki kapilar korundu


def test_run_counts_empty_matrix_raises(tmp_path, monkeypatch):
    config_path, metadata_path = _setup(tmp_path)
    run_dir = tmp_path / "run"
    _prep_through_m04(config_path, metadata_path, run_dir, monkeypatch)
    _fake_featurecounts(monkeypatch, n_genes=0)   # yanlis feature_type senaryosu
    with pytest.raises(ValueError, match="no genes"):
        run_counts(load_config(config_path), metadata_path, run_dir)


def test_run_counts_low_assignment_fails(tmp_path, monkeypatch):
    config_path, metadata_path = _setup(tmp_path)
    run_dir = tmp_path / "run"
    _prep_through_m04(config_path, metadata_path, run_dir, monkeypatch)
    _fake_featurecounts(monkeypatch, rate=0.10)
    with pytest.raises(GateFailure):
        run_counts(load_config(config_path), metadata_path, run_dir)
    gates = json.loads((run_dir / "quality" / "gates.json").read_text())["gates"]
    assert any(g["module"] == "m05_counts" and g["status"] == "FAIL" for g in gates)


def test_run_counts_resumes(tmp_path, monkeypatch):
    config_path, metadata_path = _setup(tmp_path)
    run_dir = tmp_path / "run"
    _prep_through_m04(config_path, metadata_path, run_dir, monkeypatch)
    _fake_featurecounts(monkeypatch, rate=0.9)
    run_counts(load_config(config_path), metadata_path, run_dir)
    calls = []
    monkeypatch.setattr(m05_counts, "run_featurecounts", lambda *a, **k: calls.append(1))
    summary = run_counts(load_config(config_path), metadata_path, run_dir)
    assert summary.get("resumed") is True
    assert calls == []
```

- [ ] **Step 2: Run test to verify it fails**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_m05_counts.py -q -k run_counts`
Expected: FAIL (ImportError: run_counts)

- [ ] **Step 3: Write minimal implementation**

```python
# rnaforge/modules/m05_counts.py (ekle importlar + fonksiyon)
import json
from collections import Counter
from pathlib import Path

from rnaforge.config import Config
from rnaforge.featurecounts import run_featurecounts
from rnaforge.gates import raise_if_failed, write_gate_results
from rnaforge.metadata import load_metadata
from rnaforge.quality import load_profile
from rnaforge.state import RunState


def run_counts(config: Config, metadata_path: Path, run_dir: Path,
               force: bool = False) -> dict:
    run_dir = Path(run_dir)
    quant_dir = run_dir / "quantification"
    stats_dir = run_dir / "statistics"
    logs_dir = run_dir / "logs"
    for d in (quant_dir, stats_dir, logs_dir):
        d.mkdir(parents=True, exist_ok=True)
    state = RunState(run_dir)
    stats_path = stats_dir / "count_statistics.json"

    if not force and state.is_done(MODULE_NAME) and stats_path.exists():
        summary = json.loads(stats_path.read_text())
        summary["resumed"] = True
        return summary

    if not state.is_done("m04_quant"):
        raise ValueError(
            "m05 (counts) requires m04 (quant) to have completed in this run directory "
            f"first: {run_dir}. Run `rnaforge quant` with the same --run-id, then re-run counts."
        )

    profile = load_profile(config.organism_type, config.quality)
    log_path = logs_dir / "counts.log"
    with log_path.open("w") as log_file:
        def log(msg: str) -> None:
            log_file.write(msg + "\n")
            log_file.flush()

        samples = load_metadata(metadata_path)
        bams = [quant_dir / s.sample_id / "aligned.sorted.bam" for s in samples]
        paired = any(s.fastq_2 is not None for s in samples)
        log(f"m05 featureCounts: {len(samples)} sample(s), "
            f"feature_type={config.quantification.feature_type}, "
            f"attribute={config.quantification.attribute}, paired={paired}")
        result = run_featurecounts(
            bams, config.reference.annotation_gff, quant_dir / "_featurecounts",
            feature_type=config.quantification.feature_type,
            attribute=config.quantification.attribute,
            paired=paired, threads=config.resources.threads,
        )
        if not result.gene_ids:
            raise ValueError(
                "featureCounts assigned reads to no genes (empty matrix). Likely the "
                f"feature_type ({config.quantification.feature_type!r}) or attribute "
                f"({config.quantification.attribute!r}) does not match the annotation."
            )
        state.heartbeat()

        # Sütun→sample_id KONUMLA (bams örnek sırasında verildi).
        columns = list(result.counts.keys())
        sample_ids = [s.sample_id for s in samples]
        assignment_by_sample = {
            sid: result.assignment_rates[col] for sid, col in zip(sample_ids, columns)
        }

        # counts.tsv sözleşmesi: gene\t<sample_id...>
        matrix_path = quant_dir / "counts.tsv"
        with matrix_path.open("w") as fh:
            fh.write("gene\t" + "\t".join(sample_ids) + "\n")
            for i, gene in enumerate(result.gene_ids):
                row = [str(result.counts[col][i]) for col in columns]
                fh.write(gene + "\t" + "\t".join(row) + "\n")
        log(f"count matrix written: {matrix_path} ({len(result.gene_ids)} genes)")

        gates = build_count_gates(assignment_by_sample, profile)
        summary = {
            "n_samples": len(samples), "n_genes": len(result.gene_ids),
            "samples": {sid: {"assignment_rate": assignment_by_sample[sid]} for sid in sample_ids},
            "gate_counts": dict(Counter(g.status for g in gates)),
        }
        stats_path.write_text(json.dumps(summary, indent=2))
        write_gate_results(run_dir, gates)
        for g in gates:
            log(f"gate {g.name}: {g.status} — {g.message}")
        raise_if_failed(gates)

    state.mark_done(MODULE_NAME, [str(stats_path), str(log_path)])
    return summary
```

- [ ] **Step 4: Run test to verify it passes**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_m05_counts.py -q`
Expected: PASS (tümü)

- [ ] **Step 5: Commit**

```bash
git add rnaforge/modules/m05_counts.py tests/test_m05_counts.py
git commit -m "feat(m05): run_counts orkestrasyonu (counts.tsv sozlesmesi, sutun->sample konumla)"
```

---

### Task 6: CLI `counts` subcommand

**Files:**
- Modify: `rnaforge/cli.py`
- Test: `tests/test_m05_counts.py` (ekle)

**Interfaces:**
- Consumes: `run_counts` (Task 5).
- Produces: `rnaforge counts --config ... --metadata ... [--runs-dir ...] [--run-id ...] [--force]`; exit 0 başarıda (+verdict), exit 1 assignment FAIL'de.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_m05_counts.py (ekle)
def test_cli_counts_returns_zero_and_prints_verdict(tmp_path, monkeypatch, capsys):
    from rnaforge.cli import main
    config_path, metadata_path = _setup(tmp_path)
    common = ["--config", str(config_path), "--metadata", str(metadata_path),
              "--runs-dir", str(tmp_path / "runs"), "--run-id", "demo"]
    # validate→trim→quant hazirla (ayni run-id), sonra counts. m03/m04 araclarini fake'le.
    from rnaforge.modules import m03_trim, m04_quant
    from rnaforge.modules.m03_trim import trimmed_name
    from rnaforge.fastp import FastpResult
    from rnaforge.bowtie2 import AlignmentResult
    def fake_fastp(fastq_1, out_dir, min_length, fastq_2=None, aggressive_quality=False, env="rnaforge-qc"):
        out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
        out1 = out_dir / trimmed_name(Path(fastq_1)); out1.write_text("@r\nACGT\n+\nIIII\n")
        (out_dir / "fastp.json").write_text("{}")
        return FastpResult(200, 196, 0.98, out1=out1)
    monkeypatch.setattr(m03_trim, "run_fastp", fake_fastp)
    monkeypatch.setattr(m04_quant, "build_index", lambda g, i, env="rnaforge-quant-prok": Path(i) / "genome")
    def fake_align(index_prefix, out_dir, fastq_1, fastq_2=None, threads=4, env="rnaforge-quant-prok"):
        out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
        bam = out_dir / "aligned.sorted.bam"; bam.write_bytes(b"BAM")
        return AlignmentResult(bam=bam, alignment_rate=0.95)
    monkeypatch.setattr(m04_quant, "run_bowtie2", fake_align)
    _fake_featurecounts(monkeypatch, rate=0.9, n_genes=3)

    assert main(["validate", *common]) == 0
    assert main(["trim", *common]) == 0
    assert main(["quant", *common]) == 0
    capsys.readouterr()
    assert main(["counts", *common]) == 0
    assert "quality verdict" in capsys.readouterr().out
```

- [ ] **Step 2: Run test to verify it fails**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_m05_counts.py -q -k cli_counts`
Expected: FAIL (argparse: invalid choice 'counts')

- [ ] **Step 3: Write minimal implementation**

`build_parser()`'a (`quant` parser'ından sonra) ekle:

```python
    counts = sub.add_parser("counts", help="build gene x sample count matrix (m05)")
    counts.add_argument("--config", required=True, type=Path)
    counts.add_argument("--metadata", required=True, type=Path)
    counts.add_argument("--runs-dir", type=Path, default=Path("runs"))
    counts.add_argument("--run-id", default="run")
    counts.add_argument("--force", action="store_true",
                        help="re-run even if m05 already completed in this run directory")
```

importlara: `from rnaforge.modules.m05_counts import run_counts`

Handler (m04 `_cmd_quant` deseni):

```python
def _cmd_counts(args) -> int:
    config = load_config(args.config)
    run_dir = resolve_run_dir(args.runs_dir, args.run_id)
    profile = load_profile(config.organism_type, config.quality)
    try:
        summary = run_counts(config, args.metadata, run_dir, force=args.force)
    except GateFailure:
        write_confidence_card(run_dir, profile)
        raise
    if summary.get("resumed"):
        print("m05_counts already completed in this run directory — reusing its result "
              "(use --force to re-run).")
    print(f"count matrix OK: {summary['n_genes']} genes x {summary['n_samples']} sample(s)")
    print(f"run directory: {run_dir}")
    card_path = write_confidence_card(run_dir, profile)
    card = json.loads(card_path.read_text())
    print(f"quality verdict: {card['verdict']} "
          f"(PASS={card['counts']['PASS']} WARN={card['counts']['WARN']} "
          f"FAIL={card['counts']['FAIL']}, profile={profile.name})")
    return 0
```

`main()` yönlendirmesine ekle:

```python
        if args.command == "counts":
            return _cmd_counts(args)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_m05_counts.py -q`
Expected: PASS

- [ ] **Step 5: Tüm suite + commit**

```bash
conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest -q
git add rnaforge/cli.py tests/test_m05_counts.py
git commit -m "feat(m05): 'counts' CLI subcommand (+FAIL guvence karti)"
```

---

### Task 7: Canlı doğrulama + DURUM + merge

- [ ] **Step 1: Canlı smoke (gerçek featureCounts).** Sentetik genom + GTF (birkaç gen, exon/gene_id) + o bölgelerden okuma ile `validate→trim→quant→counts` (gerçek fastp+bowtie2+featureCounts). Doğrula: `quantification/counts.tsv` (gene + sample_id sütunları, gerçek sayımlar), `count_statistics.json` (assignment yüksek), verdict, önceki kapılar korunmuş. FAIL yolu: uyumsuz anotasyon/rastgele okuma → düşük assignment → exit 1 INVALID.
- [ ] **Step 2: DURUM.md** — m05 prok BİTTİ, prokaryot uçtan-uca count matrisine ulaştı; sırada m06 DESeq2 (counts.tsv tüketir) veya m04-euk/m05-euk. Test sayısı güncelle.
- [ ] **Step 3: merge + push**

```bash
git checkout main && git merge --no-ff feat/m05-prok-counts -m "merge: m05 prokaryot count matrisi (featureCounts, assignment_rate FAIL kapisi, counts subcommand)"
conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest -q
git push origin main && git push origin feat/m05-prok-counts
```

---

## Notlar (uygulayıcı için)
- m05 zinciri: `validate → trim → quant → counts` (aynı `--run-id`); ön koşul **m04**.
- featureCounts tüm BAM'lere tek çağrı → native matris; sütun→sample_id KONUMLA (BAM adıyla değil).
- `counts.tsv` = ORTAK SÖZLEŞME (m06 DESeq2 girdisi); başlık `gene\t<sample_id...>`.
- config `quantification` bölümü `KNOWN_TOP_LEVEL_KEYS`'e eklenmezse config sertleştirme reddeder.
- `python -m rnaforge.cli` ÇALIŞMAZ; canlı testte `rnaforge` entry point + referans göreli yolları için doğru CWD.
