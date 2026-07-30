# m04 Prokaryot Hizalama (bowtie2) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** m03'ün trimlenmiş okumalarını bowtie2 ile referans genoma hizalayan, `alignment_rate` FAIL kapısıyla denetleyen ve sıralı/indeksli BAM üreten prokaryot `m04` yolu + `rnaforge quant` CLI komutu (router; ökaryot yolu net NotImplementedError).

**Architecture:** m02/m03 desenini izler. Saf parser (`parse_bowtie2_summary`) bowtie2 stderr özetinden oranı çıkarır; `build_index`/`run_bowtie2` gerçek araçları `rnaforge-quant-prok` env'de shell'ler (bowtie2 → samtools sort → index, üç ayrı `conda run`); `run_quant` router `organism_type`'a dallanır, m03'ün trimlenmiş okumalarını tüketir, stats→gates→raise sırasıyla kapı üretir.

**Tech Stack:** Python 3.11, pytest, conda (`rnaforge-quant-prok`: bowtie2 2.5.5, samtools 1.24), mevcut `rnaforge.{gates,quality,state,metadata,config}` + m03 `trimmed_reads`.

## Global Constraints
- **Test komutu:** `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest -q`.
- **Sessiz devam YOK (Kural 7);** gerçek/müşteri verisi ASLA (Kural 8, sentetik fixture).
- **alignment_rate FAIL kapısı:** oran < `profile.threshold("alignment_rate")` (prok 0.70) → FAIL, koşu durur.
- **Girdi = m03 trimlenmiş okumalar** (`trimmed/<id>/`), ham FASTQ değil. **Ön koşul m03 done.**
- **eukaryote → `NotImplementedError`** (net; sessiz/stub değil).
- **Yazma sırası (m03 deseni):** tüm örnekler hizala → `alignment_statistics.json` → `write_gate_results` → EN SON `raise_if_failed`; `mark_done` yalnız FAIL yoksa.
- **Çıktı (PLAN §14):** `quantification/_index/` (bir kez), `quantification/<id>/aligned.sorted.bam`+`.bai`+`bowtie2.log`; `statistics/alignment_statistics.json`; `quality/gates.json`; `logs/quant.log`.
- **Env:** `conda run -n rnaforge-quant-prok {bowtie2-build,bowtie2,samtools}`. Yoksa entegrasyon testi skip.

**Mevcut API imzaları (referans):**
- `rnaforge.gates`: `GateResult(...)`, `PASS/FAIL`, `write_gate_results`, `raise_if_failed`, `GateFailure`.
- `rnaforge.quality`: `load_profile(organism_type, overrides=None) -> Profile`; `Profile.threshold(gate)`, `.overrides()`, `.name`.
- `rnaforge.state`: `RunState(run_dir)`, `.is_done`, `.mark_done`, `.heartbeat`; `resolve_run_dir`.
- `rnaforge.metadata`: `load_metadata(path) -> list[Sample]`; `Sample(sample_id, condition, fastq_1, fastq_2=None, ...)`.
- `rnaforge.config`: `Config.reference.genome_fasta: Path`, `Config.organism_type`, `Config.quality`, `Config.resources.threads`.
- bowtie2 stderr özeti: son satır `NN.NN% overall alignment rate`.

---

### Task 1: bowtie2 saf özet parser (`parse_bowtie2_summary` + `AlignmentResult`)

**Files:**
- Create: `rnaforge/bowtie2.py`
- Test: `tests/test_bowtie2.py`

**Interfaces:**
- Produces:
  - `AlignmentResult` frozen dataclass: `bam: Path`, `alignment_rate: float`.
  - `parse_bowtie2_summary(stderr_text: str) -> float` — oran [0,1].
  - `Bowtie2ParseError(ValueError)`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_bowtie2.py
from __future__ import annotations

import pytest

from rnaforge.bowtie2 import AlignmentResult, Bowtie2ParseError, parse_bowtie2_summary

_STDERR = """300 reads; of these:
  300 (100.00%) were unpaired; of these:
    12 (4.00%) aligned 0 times
    288 (96.00%) aligned exactly 1 time
96.00% overall alignment rate
"""


def test_parse_reads_overall_rate():
    assert parse_bowtie2_summary(_STDERR) == pytest.approx(0.96)


def test_parse_zero_rate():
    assert parse_bowtie2_summary("0.00% overall alignment rate\n") == 0.0


def test_parse_rejects_missing_summary():
    with pytest.raises(Bowtie2ParseError, match="overall alignment rate"):
        parse_bowtie2_summary("some unrelated bowtie2 chatter\n")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_bowtie2.py -q`
Expected: FAIL (ModuleNotFoundError: rnaforge.bowtie2)

- [ ] **Step 3: Write minimal implementation**

```python
# rnaforge/bowtie2.py
"""bowtie2 hizalama: çalıştırır ve özetini parse eder. Parser saftır."""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

_RATE_RE = re.compile(r"([0-9]+(?:\.[0-9]+)?)%\s+overall alignment rate")


class Bowtie2ParseError(ValueError):
    """bowtie2 özeti beklenen biçimde değil."""


@dataclass(frozen=True)
class AlignmentResult:
    bam: Path
    alignment_rate: float


def parse_bowtie2_summary(stderr_text: str) -> float:
    match = None
    for m in _RATE_RE.finditer(stderr_text):
        match = m
    if match is None:
        raise Bowtie2ParseError(
            "bowtie2 stderr has no 'overall alignment rate' line — run may have failed"
        )
    return float(match.group(1)) / 100.0
```

- [ ] **Step 4: Run test to verify it passes**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_bowtie2.py -q`
Expected: PASS (3 passed)

- [ ] **Step 5: Commit**

```bash
git add rnaforge/bowtie2.py tests/test_bowtie2.py
git commit -m "feat(m04): bowtie2 saf ozet parser (alignment rate)"
```

---

### Task 2: gerçek index + hizalama runnerları (`build_index`, `run_bowtie2`)

**Files:**
- Modify: `rnaforge/bowtie2.py`
- Test: `tests/test_bowtie2.py` (ekle)

**Interfaces:**
- Consumes: `parse_bowtie2_summary`, `AlignmentResult` (Task 1).
- Produces:
  - `build_index(genome_fasta: Path, index_dir: Path, env: str = "rnaforge-quant-prok") -> Path` — `bowtie2-build`; index prefix (`index_dir/genome`) döner.
  - `run_bowtie2(index_prefix: Path, out_dir: Path, fastq_1: Path, fastq_2: Path | None = None, threads: int = 4, env: str = "rnaforge-quant-prok") -> AlignmentResult` — bowtie2 → samtools sort → index; `aligned.sorted.bam`+`.bai`+`bowtie2.log`.
  - `Bowtie2RunError(RuntimeError)`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_bowtie2.py (ekle)
import shutil
from pathlib import Path

from rnaforge.bowtie2 import build_index, run_bowtie2, Bowtie2RunError


def _synthetic_genome_and_reads(tmp_path, aligned=True):
    import random
    random.seed(11)
    genome = "".join(random.choice("ACGT") for _ in range(5000))
    (tmp_path / "genome.fa").write_text(">chr1\n" + genome + "\n")
    reads = tmp_path / "reads.fastq"
    with reads.open("w") as f:
        for i in range(200):
            if aligned:
                p = random.randint(0, len(genome) - 100)
                seq = genome[p:p + 100]
            else:
                seq = "".join(random.choice("ACGT") for _ in range(100))
            f.write(f"@r{i}\n{seq}\n+\n{'I' * 100}\n")
    return tmp_path / "genome.fa", reads


@pytest.mark.skipif(shutil.which("conda") is None, reason="conda yok")
def test_build_index_and_align_reports_high_rate(tmp_path):
    """Entegrasyon: genomdan türetilmiş okumalar yüksek hizalanır. rnaforge-quant-prok
    yoksa skip."""
    genome, reads = _synthetic_genome_and_reads(tmp_path, aligned=True)
    try:
        prefix = build_index(genome, tmp_path / "idx")
        result = run_bowtie2(prefix, tmp_path / "aln", reads)
    except Bowtie2RunError as exc:
        pytest.skip(f"bowtie2 çalıştırılamadı (env yok?): {exc}")
    assert result.bam.exists()
    assert Path(str(result.bam) + ".bai").exists()   # samtools index -> aligned.sorted.bam.bai
    assert result.alignment_rate > 0.95


@pytest.mark.skipif(shutil.which("conda") is None, reason="conda yok")
def test_random_reads_align_poorly(tmp_path):
    genome, reads = _synthetic_genome_and_reads(tmp_path, aligned=False)
    try:
        prefix = build_index(genome, tmp_path / "idx")
        result = run_bowtie2(prefix, tmp_path / "aln", reads)
    except Bowtie2RunError as exc:
        pytest.skip(f"bowtie2 çalıştırılamadı: {exc}")
    assert result.alignment_rate < 0.70   # genom-dışı okumalar eşiğin altında
```

- [ ] **Step 2: Run test to verify it fails**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_bowtie2.py -q -k "align or random"`
Expected: FAIL (ImportError: build_index / run_bowtie2)

- [ ] **Step 3: Write minimal implementation**

```python
# rnaforge/bowtie2.py (ekle — importlara subprocess ekle)
import subprocess


class Bowtie2RunError(RuntimeError):
    """bowtie2/samtools çalıştırılamadı ya da beklenen çıktıyı üretmedi."""


def _run(cmd: list[str]) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, capture_output=True, text=True)


def build_index(genome_fasta: Path, index_dir: Path,
                env: str = "rnaforge-quant-prok") -> Path:
    genome_fasta = Path(genome_fasta)
    if not genome_fasta.exists():
        raise Bowtie2RunError(f"genome FASTA does not exist: {genome_fasta}")
    index_dir = Path(index_dir)
    index_dir.mkdir(parents=True, exist_ok=True)
    prefix = index_dir / "genome"
    r = _run(["conda", "run", "-n", env, "bowtie2-build", "-q",
              str(genome_fasta), str(prefix)])
    if r.returncode != 0:
        raise Bowtie2RunError(
            f"bowtie2-build failed (exit {r.returncode}) for {genome_fasta}\n"
            f"stderr: {r.stderr.strip()}"
        )
    return prefix


def run_bowtie2(index_prefix: Path, out_dir: Path, fastq_1: Path,
                fastq_2: Path | None = None, threads: int = 4,
                env: str = "rnaforge-quant-prok") -> AlignmentResult:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    sam = out_dir / "aligned.sam"
    bam = out_dir / "aligned.sorted.bam"
    log = out_dir / "bowtie2.log"

    bt = ["conda", "run", "-n", env, "bowtie2", "-x", str(index_prefix),
          "-p", str(threads), "-S", str(sam)]
    if fastq_2 is not None:
        bt += ["-1", str(fastq_1), "-2", str(fastq_2)]
    else:
        bt += ["-U", str(fastq_1)]
    r = _run(bt)
    log.write_text(r.stderr)
    if r.returncode != 0:
        raise Bowtie2RunError(
            f"bowtie2 failed (exit {r.returncode})\ncmd: {' '.join(bt)}\n"
            f"stderr: {r.stderr.strip()}"
        )
    rate = parse_bowtie2_summary(r.stderr)

    sort = _run(["conda", "run", "-n", env, "samtools", "sort", "-o", str(bam), str(sam)])
    if sort.returncode != 0:
        raise Bowtie2RunError(f"samtools sort failed: {sort.stderr.strip()}")
    idx = _run(["conda", "run", "-n", env, "samtools", "index", str(bam)])
    if idx.returncode != 0:
        raise Bowtie2RunError(f"samtools index failed: {idx.stderr.strip()}")
    sam.unlink(missing_ok=True)
    return AlignmentResult(bam=bam, alignment_rate=rate)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_bowtie2.py -q`
Expected: PASS (entegrasyon ya geçer ya skip)

- [ ] **Step 5: Commit**

```bash
git add rnaforge/bowtie2.py tests/test_bowtie2.py
git commit -m "feat(m04): gercek bowtie2 index + hizalama runner (samtools sort/index)"
```

---

### Task 3: trimlenmiş okuma yol sözleşmesi (`trimmed_reads`)

**Files:**
- Modify: `rnaforge/fastp.py` (`_trimmed_name` → public `trimmed_name`)
- Modify: `rnaforge/modules/m03_trim.py` (ekle `trimmed_reads`)
- Test: `tests/test_m03_trim.py` (ekle)

**Interfaces:**
- Produces: `rnaforge.modules.m03_trim.trimmed_reads(run_dir: Path, sample: Sample) -> tuple[Path, Path | None]` — m03'ün `trimmed/<id>/<stem>.trimmed.fastq` çıktı yolları (adlandırma kuralının TEK kaynağı).
- Değişiklik: `rnaforge.fastp.trimmed_name(fastq: Path) -> str` (eski `_trimmed_name`; `run_fastp` da bunu kullanır).

- [ ] **Step 1: Write the failing test**

```python
# tests/test_m03_trim.py (ekle — importlara Sample ekle)
from rnaforge.metadata import Sample
from rnaforge.modules.m03_trim import trimmed_reads


def test_trimmed_reads_single_end(tmp_path):
    sample = Sample("s1", "control", tmp_path / "c1.fastq")
    out1, out2 = trimmed_reads(tmp_path / "run", sample)
    assert out1 == tmp_path / "run" / "trimmed" / "s1" / "c1.trimmed.fastq"
    assert out2 is None


def test_trimmed_reads_paired_end(tmp_path):
    sample = Sample("s1", "control", tmp_path / "c1_R1.fastq", tmp_path / "c1_R2.fastq")
    out1, out2 = trimmed_reads(tmp_path / "run", sample)
    assert out1.name == "c1_R1.trimmed.fastq"
    assert out2.name == "c1_R2.trimmed.fastq"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_m03_trim.py -q -k trimmed_reads`
Expected: FAIL (ImportError: trimmed_reads)

- [ ] **Step 3: Write minimal implementation**

`rnaforge/fastp.py`: `_trimmed_name` → `trimmed_name` (public) yeniden adlandır; `run_fastp` içindeki iki çağrıyı da güncelle (`_trimmed_name(fastq_1)` → `trimmed_name(fastq_1)`, `_trimmed_name(fastq_2)` → `trimmed_name(fastq_2)`).

`rnaforge/modules/m03_trim.py`'ye ekle:

```python
from rnaforge.fastp import trimmed_name
from rnaforge.metadata import Sample


def trimmed_reads(run_dir, sample: Sample) -> tuple[Path, Path | None]:
    """m03'ün bir örnek için ürettiği trimlenmiş FASTQ yol(lar)ı. Adlandırma
    kuralının TEK kaynağı: m03 buraya yazar, m04 buradan okur (drift önlenir)."""
    d = Path(run_dir) / "trimmed" / sample.sample_id
    out1 = d / trimmed_name(sample.fastq_1)
    out2 = d / trimmed_name(sample.fastq_2) if sample.fastq_2 else None
    return out1, out2
```

- [ ] **Step 4: Run test to verify it passes**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_m03_trim.py -q`
Expected: PASS (trimmed_reads + mevcut m03 testleri hâlâ yeşil)

- [ ] **Step 5: Commit**

```bash
git add rnaforge/fastp.py rnaforge/modules/m03_trim.py tests/test_m03_trim.py
git commit -m "feat(m04): trimmed_reads yol sozlesmesi (adlandirma tek kaynak); trimmed_name public"
```

---

### Task 4: alignment_rate kapısı (`build_alignment_gates`)

**Files:**
- Create: `rnaforge/modules/m04_quant.py`
- Test: `tests/test_m04_quant.py`

**Interfaces:**
- Consumes: `AlignmentResult` (Task 1); `rnaforge.gates.{GateResult,PASS,FAIL}`; `rnaforge.quality.Profile`.
- Produces: `MODULE_NAME = "m04_quant"`; `build_alignment_gates(results: dict[str, AlignmentResult], profile: Profile) -> list[GateResult]`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_m04_quant.py
from __future__ import annotations

from pathlib import Path

from rnaforge.bowtie2 import AlignmentResult
from rnaforge.gates import FAIL, PASS
from rnaforge.modules.m04_quant import build_alignment_gates
from rnaforge.quality import load_profile


def _res(rate: float) -> AlignmentResult:
    return AlignmentResult(bam=Path("x.bam"), alignment_rate=rate)


def test_all_above_threshold_passes():
    profile = load_profile("prokaryote")  # alignment_rate = 0.70
    gates = build_alignment_gates({"s1": _res(0.98), "s2": _res(0.85)}, profile)
    assert len(gates) == 1
    assert gates[0].name == "alignment_rate"
    assert gates[0].module == "m04_quant"
    assert gates[0].status == PASS


def test_below_threshold_fails_and_lists_offenders():
    profile = load_profile("prokaryote")
    gates = build_alignment_gates({"s1": _res(0.98), "s2": _res(0.30)}, profile)
    g = gates[0]
    assert g.status == FAIL
    assert g.samples == ("s2",)
    assert g.measured == 0.30
    assert g.threshold == 0.70


def test_override_marks_gate_overridden():
    profile = load_profile("prokaryote", {"alignment_rate": 0.20})
    gates = build_alignment_gates({"s1": _res(0.30)}, profile)
    assert gates[0].status == PASS
    assert gates[0].overridden is True
    assert gates[0].threshold == 0.20
```

- [ ] **Step 2: Run test to verify it fails**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_m04_quant.py -q`
Expected: FAIL (ModuleNotFoundError: rnaforge.modules.m04_quant)

- [ ] **Step 3: Write minimal implementation**

```python
# rnaforge/modules/m04_quant.py
"""m04 — Quantification ROUTER (prokaryot: bowtie2 genom hizalama).

organism_type'a göre dallanır (PLAN §5). Şimdilik yalnız prokaryot yolu bağlı;
eukaryote (Salmon) net NotImplementedError verir. Veri kapısı `alignment_rate`:
bowtie2 overall alignment rate profil eşiğinin altındaysa FAIL — düşük hizalama
güvenilmez count matrisi üretir (PLAN §3)."""
from __future__ import annotations

from rnaforge.bowtie2 import AlignmentResult
from rnaforge.gates import FAIL, PASS, GateResult
from rnaforge.quality import Profile

MODULE_NAME = "m04_quant"
_GATE = "alignment_rate"


def build_alignment_gates(results: dict[str, AlignmentResult],
                          profile: Profile) -> list[GateResult]:
    threshold = profile.threshold(_GATE)
    offenders = sorted(sid for sid, r in results.items() if r.alignment_rate < threshold)
    lowest = min((r.alignment_rate for r in results.values()), default=1.0)
    overridden = _GATE in profile.overrides()
    if offenders:
        status = FAIL
        message = (
            f"hizalama oranı eşiğin altında ({len(offenders)} örnek: "
            f"{', '.join(offenders)}); en düşük {lowest:.2f} < {threshold:.2f}. "
            "Düşük hizalama güvenilmez count matrisi üretir."
        )
    else:
        status = PASS
        message = f"tüm örnekler alignment ≥ {threshold:.2f} (en düşük {lowest:.2f})."
    return [GateResult(
        name=_GATE, module=MODULE_NAME, status=status, message=message,
        remedy=("Referans genomu, kütüphane kimyasını ve tür kimliğini doğrulayın; "
                "yanlış referans veya kontaminasyon hizalamayı düşürür."),
        measured=lowest, threshold=threshold, overridden=overridden,
        samples=tuple(offenders),
    )]
```

- [ ] **Step 4: Run test to verify it passes**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_m04_quant.py -q`
Expected: PASS (3 passed)

- [ ] **Step 5: Commit**

```bash
git add rnaforge/modules/m04_quant.py tests/test_m04_quant.py
git commit -m "feat(m04): alignment_rate FAIL kapisi (profil esigine karsi)"
```

---

### Task 5: `run_quant` router orkestrasyonu

**Files:**
- Modify: `rnaforge/modules/m04_quant.py`
- Test: `tests/test_m04_quant.py` (ekle)

**Interfaces:**
- Consumes: `build_index`/`run_bowtie2` (Task 2), `trimmed_reads` (Task 3), `build_alignment_gates` (Task 4), `RunState`, `write_gate_results`, `raise_if_failed`, `load_metadata`, `load_profile`, `Config`.
- Produces: `run_quant(config: Config, metadata_path: Path, run_dir: Path, force: bool = False) -> dict`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_m04_quant.py (ekle)
import json
import textwrap

import pytest

from rnaforge.config import load_config
from rnaforge.gates import GateFailure
from rnaforge.modules import m04_quant
from rnaforge.modules.m04_quant import run_quant
from tests.conftest import write_fastq


def _setup(tmp_path, organism_type="prokaryote"):
    (tmp_path / "ref").mkdir()
    (tmp_path / "ref" / "genome.fa").write_text(">c1\n" + "ACGT" * 25 + "\n")
    (tmp_path / "ref" / "genes.gff").write_text("##gff-version 3\n")
    ref = ("genome_fasta" if organism_type == "prokaryote" else "transcriptome_fasta")
    extra = ("annotation_gff" if organism_type == "prokaryote" else "tx2gene")
    for n in ("c1.fastq", "c2.fastq", "t1.fastq", "t2.fastq"):
        write_fastq(tmp_path / n, 200, 150, "I")
    config_path = tmp_path / "config.yaml"
    config_path.write_text(textwrap.dedent(f"""
        organism: "E. coli"
        organism_type: "{organism_type}"
        reference:
          {ref}: "{tmp_path / 'ref' / 'genome.fa'}"
          {extra}: "{tmp_path / 'ref' / 'genes.gff'}"
    """))
    metadata_path = tmp_path / "samples.tsv"
    metadata_path.write_text(
        "sample_id\tcondition\tfastq_1\n"
        "s1\tcontrol\tc1.fastq\n" "s2\tcontrol\tc2.fastq\n"
        "s3\ttreated\tt1.fastq\n" "s4\ttreated\tt2.fastq\n"
    )
    return config_path, metadata_path


def _prep_m01_m03(config_path, metadata_path, run_dir, monkeypatch, survival=0.98):
    """m01 (gerçek) + m03 (fastp monkeypatch) hazırlar — m04 için trimlenmiş okuma
    ve done state gerekir."""
    from rnaforge.modules.m01_validate import run_validation
    from rnaforge.modules import m03_trim
    from rnaforge.modules.m03_trim import run_trim, trimmed_name
    from rnaforge.fastp import FastpResult
    run_validation(load_config(config_path), metadata_path, run_dir)

    def fake_fastp(fastq_1, out_dir, min_length, fastq_2=None,
                   aggressive_quality=False, env="rnaforge-qc"):
        from pathlib import Path as P
        out_dir = P(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
        out1 = out_dir / trimmed_name(P(fastq_1))
        out1.write_text("@r\nACGT\n+\nIIII\n")
        (out_dir / "fastp.json").write_text("{}")
        return FastpResult(reads_before=200, reads_after=int(200 * survival),
                           survival_rate=survival, out1=out1)
    monkeypatch.setattr(m03_trim, "run_fastp", fake_fastp)
    run_trim(load_config(config_path), metadata_path, run_dir)


def _fake_bowtie2(monkeypatch, rate=0.95):
    from rnaforge.bowtie2 import AlignmentResult
    from pathlib import Path as P
    monkeypatch.setattr(m04_quant, "build_index", lambda genome, index_dir, env=...: P(index_dir) / "genome")
    def fake_align(index_prefix, out_dir, fastq_1, fastq_2=None, threads=4, env=...):
        out_dir = P(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
        bam = out_dir / "aligned.sorted.bam"; bam.write_bytes(b"BAM")
        return AlignmentResult(bam=bam, alignment_rate=rate)
    monkeypatch.setattr(m04_quant, "run_bowtie2", fake_align)


def test_run_quant_requires_m03_done(tmp_path, monkeypatch):
    _fake_bowtie2(monkeypatch)
    config_path, metadata_path = _setup(tmp_path)
    from rnaforge.modules.m01_validate import run_validation
    run_dir = tmp_path / "run"
    run_validation(load_config(config_path), metadata_path, run_dir)  # yalniz m01
    with pytest.raises(ValueError, match="m03"):
        run_quant(load_config(config_path), metadata_path, run_dir)


def test_run_quant_eukaryote_not_implemented(tmp_path, monkeypatch):
    _fake_bowtie2(monkeypatch)
    config_path, metadata_path = _setup(tmp_path, organism_type="eukaryote")
    run_dir = tmp_path / "run"
    _prep_m01_m03(config_path, metadata_path, run_dir, monkeypatch)
    with pytest.raises(NotImplementedError, match="eukaryote"):
        run_quant(load_config(config_path), metadata_path, run_dir)


def test_run_quant_writes_bam_and_passes(tmp_path, monkeypatch):
    config_path, metadata_path = _setup(tmp_path)
    run_dir = tmp_path / "run"
    _prep_m01_m03(config_path, metadata_path, run_dir, monkeypatch)
    _fake_bowtie2(monkeypatch, rate=0.95)
    summary = run_quant(load_config(config_path), metadata_path, run_dir)

    assert summary["n_samples"] == 4
    assert (run_dir / "quantification" / "s1" / "aligned.sorted.bam").exists()
    stats = json.loads((run_dir / "statistics" / "alignment_statistics.json").read_text())
    assert set(stats["samples"]) == {"s1", "s2", "s3", "s4"}
    gates = json.loads((run_dir / "quality" / "gates.json").read_text())["gates"]
    assert any(g["module"] == "m04_quant" and g["status"] == "PASS" for g in gates)
    assert any(g["module"] == "m03_trim" for g in gates)   # onceki kapilar korundu


def test_run_quant_low_alignment_fails(tmp_path, monkeypatch):
    config_path, metadata_path = _setup(tmp_path)
    run_dir = tmp_path / "run"
    _prep_m01_m03(config_path, metadata_path, run_dir, monkeypatch)
    _fake_bowtie2(monkeypatch, rate=0.10)
    with pytest.raises(GateFailure):
        run_quant(load_config(config_path), metadata_path, run_dir)
    gates = json.loads((run_dir / "quality" / "gates.json").read_text())["gates"]
    assert any(g["module"] == "m04_quant" and g["status"] == "FAIL" for g in gates)


def test_run_quant_resumes(tmp_path, monkeypatch):
    config_path, metadata_path = _setup(tmp_path)
    run_dir = tmp_path / "run"
    _prep_m01_m03(config_path, metadata_path, run_dir, monkeypatch)
    _fake_bowtie2(monkeypatch, rate=0.95)
    run_quant(load_config(config_path), metadata_path, run_dir)
    calls = []
    monkeypatch.setattr(m04_quant, "run_bowtie2", lambda *a, **k: calls.append(1))
    summary = run_quant(load_config(config_path), metadata_path, run_dir)
    assert summary.get("resumed") is True
    assert calls == []
```

- [ ] **Step 2: Run test to verify it fails**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_m04_quant.py -q -k run_quant`
Expected: FAIL (ImportError: run_quant)

- [ ] **Step 3: Write minimal implementation**

```python
# rnaforge/modules/m04_quant.py (ekle importlar + fonksiyon)
import json
from collections import Counter
from pathlib import Path

from rnaforge.bowtie2 import build_index, run_bowtie2
from rnaforge.config import Config
from rnaforge.gates import raise_if_failed, write_gate_results
from rnaforge.metadata import load_metadata
from rnaforge.modules.m03_trim import trimmed_reads
from rnaforge.quality import load_profile
from rnaforge.state import RunState


def run_quant(config: Config, metadata_path: Path, run_dir: Path,
              force: bool = False) -> dict:
    run_dir = Path(run_dir)
    quant_dir = run_dir / "quantification"
    stats_dir = run_dir / "statistics"
    logs_dir = run_dir / "logs"
    for d in (quant_dir, stats_dir, logs_dir):
        d.mkdir(parents=True, exist_ok=True)
    state = RunState(run_dir)
    stats_path = stats_dir / "alignment_statistics.json"

    if not force and state.is_done(MODULE_NAME) and stats_path.exists():
        summary = json.loads(stats_path.read_text())
        summary["resumed"] = True
        return summary

    if config.organism_type == "eukaryote":
        raise NotImplementedError(
            "m04 eukaryote (Salmon) path not yet implemented; prokaryote only for now."
        )
    if not state.is_done("m03_trim"):
        raise ValueError(
            "m04 (quant) requires m03 (trim) to have completed in this run directory "
            f"first: {run_dir}. Run `rnaforge trim` with the same --run-id, then re-run quant."
        )

    profile = load_profile(config.organism_type, config.quality)
    log_path = logs_dir / "quant.log"
    with log_path.open("w") as log_file:
        def log(msg: str) -> None:
            log_file.write(msg + "\n")
            log_file.flush()

        samples = load_metadata(metadata_path)
        index_prefix = build_index(config.reference.genome_fasta, quant_dir / "_index")
        log(f"m04 bowtie2: index built, {len(samples)} sample(s)")
        results = {}
        per_sample = {}
        for sample in samples:
            state.heartbeat()
            t1, t2 = trimmed_reads(run_dir, sample)
            result = run_bowtie2(index_prefix, quant_dir / sample.sample_id, t1,
                                 fastq_2=t2, threads=config.resources.threads)
            results[sample.sample_id] = result
            per_sample[sample.sample_id] = {
                "alignment_rate": result.alignment_rate, "bam": str(result.bam),
            }
            log(f"{sample.sample_id}: alignment_rate={result.alignment_rate:.3f}")

        gates = build_alignment_gates(results, profile)
        summary = {
            "n_samples": len(samples), "samples": per_sample,
            "gate_counts": dict(Counter(g.status for g in gates)),
        }
        stats_path.write_text(json.dumps(summary, indent=2))
        write_gate_results(run_dir, gates)
        for g in gates:
            log(f"gate {g.name}: {g.status} — {g.message}")
        raise_if_failed(gates)
        log(f"alignment statistics written: {stats_path}")

    state.mark_done(MODULE_NAME, [str(stats_path), str(log_path)])
    return summary
```

- [ ] **Step 4: Run test to verify it passes**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_m04_quant.py -q`
Expected: PASS (tümü)

- [ ] **Step 5: Commit**

```bash
git add rnaforge/modules/m04_quant.py tests/test_m04_quant.py
git commit -m "feat(m04): run_quant router (bowtie2 prok yolu, m03 on kosul, euk NotImplemented)"
```

---

### Task 6: CLI `quant` subcommand

**Files:**
- Modify: `rnaforge/cli.py`
- Test: `tests/test_m04_quant.py` (ekle)

**Interfaces:**
- Consumes: `run_quant` (Task 5). `main()` GateFailure'ı zaten yakalıyor (exit 1).
- Produces: `rnaforge quant --config ... --metadata ... [--runs-dir ...] [--run-id ...] [--force]`; exit 0 başarıda (+verdict), exit 1 alignment FAIL'de.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_m04_quant.py (ekle)
from rnaforge.cli import main


def test_cli_quant_returns_zero_and_prints_verdict(tmp_path, monkeypatch, capsys):
    config_path, metadata_path = _setup(tmp_path)
    common = ["--config", str(config_path), "--metadata", str(metadata_path),
              "--runs-dir", str(tmp_path / "runs"), "--run-id", "demo"]
    # Aynı run-id ile m01 (validate) + m03 (trim) hazırla, sonra quant.
    # Not: resolve_run_dir ayni run-id'yi ayni dizine cozer.
    from rnaforge.modules import m03_trim
    from rnaforge.modules.m03_trim import trimmed_name
    from rnaforge.fastp import FastpResult
    def fake_fastp(fastq_1, out_dir, min_length, fastq_2=None, aggressive_quality=False, env="rnaforge-qc"):
        from pathlib import Path as P
        out_dir = P(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
        out1 = out_dir / trimmed_name(P(fastq_1)); out1.write_text("@r\nACGT\n+\nIIII\n")
        (out_dir / "fastp.json").write_text("{}")
        return FastpResult(200, int(200 * 0.98), 0.98, out1=out1)
    monkeypatch.setattr(m03_trim, "run_fastp", fake_fastp)
    _fake_bowtie2(monkeypatch, rate=0.95)

    assert main(["validate", *common]) == 0
    assert main(["trim", *common]) == 0
    capsys.readouterr()
    assert main(["quant", *common]) == 0
    assert "quality verdict" in capsys.readouterr().out
```

- [ ] **Step 2: Run test to verify it fails**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_m04_quant.py -q -k cli_quant`
Expected: FAIL (argparse: invalid choice 'quant')

- [ ] **Step 3: Write minimal implementation**

`build_parser()`'a (`trim` parser'ından sonra) ekle:

```python
    quant = sub.add_parser("quant", help="align reads to reference (m04)")
    quant.add_argument("--config", required=True, type=Path)
    quant.add_argument("--metadata", required=True, type=Path)
    quant.add_argument("--runs-dir", type=Path, default=Path("runs"))
    quant.add_argument("--run-id", default="run")
    quant.add_argument("--force", action="store_true",
                       help="re-run even if m04 already completed in this run directory")
```

importlara: `from rnaforge.modules.m04_quant import run_quant`

Handler (m03 `_cmd_trim` deseni: FAIL'de de güvence kartı):

```python
def _cmd_quant(args) -> int:
    config = load_config(args.config)
    run_dir = resolve_run_dir(args.runs_dir, args.run_id)
    profile = load_profile(config.organism_type, config.quality)
    try:
        summary = run_quant(config, args.metadata, run_dir, force=args.force)
    except GateFailure:
        write_confidence_card(run_dir, profile)
        raise
    if summary.get("resumed"):
        print("m04_quant already completed in this run directory — reusing its result "
              "(use --force to re-run).")
    print(f"alignment OK: {summary['n_samples']} sample(s)")
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
        if args.command == "quant":
            return _cmd_quant(args)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_m04_quant.py -q`
Expected: PASS

- [ ] **Step 5: Tüm suite + commit**

```bash
conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest -q
git add rnaforge/cli.py tests/test_m04_quant.py
git commit -m "feat(m04): 'quant' CLI subcommand (+FAIL guvence karti)"
```

---

### Task 7: Canlı doğrulama + DURUM + merge

- [ ] **Step 1: Canlı smoke (gerçek bowtie2).** Sentetik genom + O GENOMDAN türetilmiş okumalar (yüksek hizalama) ile `validate → trim → quant`. Doğrula: `quantification/<id>/aligned.sorted.bam`+`.bai`, `alignment_statistics.json` (rate >0.95), verdict, önceki kapılar korunmuş. Ayrıca genom-DIŞI okumalarla FAIL yolu (alignment < 0.70 → exit 1, INVALID).
- [ ] **Step 2: DURUM.md** — m04 prok BİTTİ, sırada m05 (featureCounts count matrisi) + m04 euk (salmon). Test sayısı güncelle.
- [ ] **Step 3: merge + push**

```bash
git checkout main && git merge --no-ff feat/m04-prok-alignment -m "merge: m04 prokaryot bowtie2 hizalama (alignment_rate FAIL kapisi, quant subcommand)"
conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest -q
git push origin main && git push origin feat/m04-prok-alignment
```

---

## Notlar (uygulayıcı için)
- m04 zinciri: `validate` → `trim` → `quant` (aynı `--run-id`); m04 ön koşulu **m03** (trimlenmiş okuma gerekir).
- Trimlenmiş yol TEK kaynaktan (`m03_trim.trimmed_reads`); m04 dosya adı uydurmaz.
- Router yapısı korunur: euk armı `NotImplementedError` (sonraki spec salmon yolunu doldurur).
- BAM ara üründür; gen×örnek count matrisi m05'in (featureCounts) işidir.
- `python -m rnaforge.cli` ÇALIŞMAZ; canlı testte `rnaforge` entry point + referans göreli yolları için doğru CWD.
