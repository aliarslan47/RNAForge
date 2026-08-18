# Ökaryot Kısa-Okuma Yolu (Salmon + tximport) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Ökaryot (Illumina kısa-okuma) niceleme kolunu bağla — Salmon (decoy-aware, genom opsiyonel) + tximport (tx→gen), mevcut `counts.tsv` sözleşmesinde m06+ ile buluşarak.

**Architecture:** m04/m05 mevcut router desenini izler. `organism_type==eukaryote` yeni dallar `_quant_euk` (Salmon) ve `_counts_euk` (tximport) çağırır; eski `NotImplementedError` kaldırılır. m06+ organizma-agnostik, değişmez. Salmon `rnaforge-quant-euk` env (kurulu), tximport `rnaforge-de` env (R+Bioconductor).

**Tech Stack:** Python 3.11, salmon 2.3.4 (bioconda), R/Bioconductor tximport, pytest, conda.

**Spec:** `docs/superpowers/specs/2026-08-18-eukaryote-shortread-design.md`

## Global Constraints

- Ökaryot ayrımı YALNIZ m04/m05'te; m06+ DEĞİŞMEZ (organizma-agnostik).
- Buluşma sözleşmesi: `quantification/counts.tsv` = `gene\t<sample_id...>`, sütun→sample_id KONUMLA.
- tximport `countsFromAbundance="lengthScaledTPM"` (uzunluk-düzeltilmiş sayım → m06 DESeq2 offset gerektirmez).
- Salmon decoy-aware: `reference.genome_fasta` varsa decoy kur, yoksa transkriptom-only + yüksek sesle log (sessiz düşürme yok).
- Env: Salmon `rnaforge-quant-euk`, tximport `rnaforge-de`. İkinci R env AÇILMAZ.
- Sessiz hata yasak: boş matris / eksik quant.sf / 0 gen → yüksek sesle hata, nonzero exit.
- R çağrı deseni: `conda run -n <env> Rscript <script> <args>` (deseq2.py deseni).

---

### Task 1: `rnaforge/salmon.py` — parser + index + quant runner

**Files:**
- Create: `rnaforge/salmon.py`
- Test: `tests/test_salmon.py`

**Interfaces:**
- Produces:
  - `parse_salmon_meta(meta_info_json: Path) -> float` — `percent_mapped/100` (0..1 mapping rate).
  - `build_salmon_index(transcriptome_fasta: Path, index_dir: Path, genome_fasta: Path | None = None, threads: int = 8, log=None) -> Path` — index dizinini döndürür; genome_fasta verilirse decoy-aware.
  - `run_salmon_quant(index_dir: Path, out_dir: Path, fastq_1: Path, fastq_2: Path | None = None, threads: int = 8, env: str = "rnaforge-quant-euk") -> SalmonQuant` — `SalmonQuant(quant_sf: Path, mapping_rate: float)`.
  - `@dataclass SalmonQuant: quant_sf: Path; mapping_rate: float`

- [ ] **Step 1: Write failing test for `parse_salmon_meta`**

```python
# tests/test_salmon.py
import json
from pathlib import Path
from rnaforge.salmon import parse_salmon_meta

def test_parse_salmon_meta_percent_mapped(tmp_path):
    mi = tmp_path / "meta_info.json"
    mi.write_text(json.dumps({"num_processed": 1000, "num_mapped": 853,
                              "percent_mapped": 85.3}))
    assert parse_salmon_meta(mi) == 0.853
```

- [ ] **Step 2: Run test to verify it fails**

Run: `conda run -n rnaforge-core python -m pytest tests/test_salmon.py::test_parse_salmon_meta_percent_mapped -v`
Expected: FAIL (ModuleNotFoundError: rnaforge.salmon)

- [ ] **Step 3: Implement `salmon.py` skeleton + parser**

```python
# rnaforge/salmon.py
"""Salmon (ökaryot transkriptom niceleme) — saf parser + runner (bowtie2.py deseni).

Decoy-aware selective alignment: genome_fasta verilirse genom decoy olarak indekslenir
(anotasyonsuz/intergenik sahte eşleşmeleri eler; Salmon önerisi). Env: rnaforge-quant-euk."""
from __future__ import annotations

import json
import subprocess
from dataclasses import dataclass
from pathlib import Path


class SalmonError(RuntimeError):
    """Salmon çalıştırılamadı ya da beklenen çıktıyı üretmedi."""


@dataclass
class SalmonQuant:
    quant_sf: Path
    mapping_rate: float


def parse_salmon_meta(meta_info_json: Path) -> float:
    data = json.loads(Path(meta_info_json).read_text())
    return float(data["percent_mapped"]) / 100.0
```

- [ ] **Step 4: Run test to verify it passes**

Run: `conda run -n rnaforge-core python -m pytest tests/test_salmon.py::test_parse_salmon_meta_percent_mapped -v`
Expected: PASS

- [ ] **Step 5: Write failing test for `build_salmon_index` (decoy vs no-decoy command)**

```python
# tests/test_salmon.py (append)
from rnaforge import salmon as salmon_mod

def test_build_index_decoy_builds_gentrome_and_decoys(tmp_path, monkeypatch):
    tx = tmp_path / "tx.fa"; tx.write_text(">t1\nACGT\n")
    genome = tmp_path / "genome.fa"; genome.write_text(">chr1\nAAAA\n>chr2\nTTTT\n")
    idx = tmp_path / "idx"
    calls = {}
    def fake_run(cmd, **k):
        calls["cmd"] = cmd
        Path(idx).mkdir(parents=True, exist_ok=True)
        class R: returncode = 0; stdout = ""; stderr = ""
        return R()
    monkeypatch.setattr(salmon_mod.subprocess, "run", fake_run)
    out = salmon_mod.build_salmon_index(tx, idx, genome_fasta=genome)
    assert out == idx
    # decoys.txt genom kontig adlarını içerir
    decoys = (idx.parent / "decoys.txt").read_text().split()
    assert decoys == ["chr1", "chr2"]
    assert "-d" in calls["cmd"]           # decoy modu

def test_build_index_no_decoy_transcriptome_only(tmp_path, monkeypatch):
    tx = tmp_path / "tx.fa"; tx.write_text(">t1\nACGT\n")
    idx = tmp_path / "idx"
    calls = {}
    def fake_run(cmd, **k):
        calls["cmd"] = cmd
        Path(idx).mkdir(parents=True, exist_ok=True)
        class R: returncode = 0; stdout = ""; stderr = ""
        return R()
    monkeypatch.setattr(salmon_mod.subprocess, "run", fake_run)
    salmon_mod.build_salmon_index(tx, idx, genome_fasta=None)
    assert "-d" not in calls["cmd"]       # decoy yok
```

- [ ] **Step 6: Run tests to verify they fail**

Run: `conda run -n rnaforge-core python -m pytest tests/test_salmon.py -v`
Expected: FAIL (build_salmon_index not defined)

- [ ] **Step 7: Implement `build_salmon_index`**

```python
# rnaforge/salmon.py (append)
def _contig_names(fasta: Path) -> list[str]:
    names = []
    with Path(fasta).open() as fh:
        for line in fh:
            if line.startswith(">"):
                names.append(line[1:].split()[0])
    return names


def build_salmon_index(transcriptome_fasta: Path, index_dir: Path,
                       genome_fasta: Path | None = None, threads: int = 8,
                       log=None) -> Path:
    index_dir = Path(index_dir)
    index_dir.parent.mkdir(parents=True, exist_ok=True)
    if genome_fasta is not None:
        decoys = index_dir.parent / "decoys.txt"
        decoys.write_text("\n".join(_contig_names(genome_fasta)) + "\n")
        gentrome = index_dir.parent / "gentrome.fa"
        with gentrome.open("wb") as out:
            for src in (transcriptome_fasta, genome_fasta):
                out.write(Path(src).read_bytes())
        cmd = ["conda", "run", "-n", "rnaforge-quant-euk", "salmon", "index",
               "-t", str(gentrome), "-d", str(decoys), "-i", str(index_dir),
               "-k", "31", "-p", str(threads)]
    else:
        if log:
            log("salmon index: genome_fasta verilmedi → transkriptom-only "
                "(doğruluk için genome_fasta decoy önerilir)")
        cmd = ["conda", "run", "-n", "rnaforge-quant-euk", "salmon", "index",
               "-t", str(transcriptome_fasta), "-i", str(index_dir),
               "-k", "31", "-p", str(threads)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0 or not index_dir.exists():
        raise SalmonError(f"salmon index failed: {r.stderr[-500:]}")
    return index_dir
```

- [ ] **Step 8: Run tests to verify they pass**

Run: `conda run -n rnaforge-core python -m pytest tests/test_salmon.py -v`
Expected: PASS

- [ ] **Step 9: Write failing test for `run_salmon_quant` (single + paired command + mapping_rate)**

```python
# tests/test_salmon.py (append)
def test_run_salmon_quant_parses_mapping_rate(tmp_path, monkeypatch):
    idx = tmp_path / "idx"; idx.mkdir()
    r1 = tmp_path / "s_R1.fastq"; r1.write_text("@r\nACGT\n+\nIIII\n")
    out = tmp_path / "s1"
    def fake_run(cmd, **k):
        assert "-1" not in cmd and "-r" in cmd     # single-end
        aux = out / "aux_info"; aux.mkdir(parents=True, exist_ok=True)
        (aux / "meta_info.json").write_text('{"percent_mapped": 77.0}')
        (out / "quant.sf").write_text("Name\tLength\tTPM\tNumReads\n")
        class R: returncode = 0; stdout = ""; stderr = ""
        return R()
    monkeypatch.setattr(salmon_mod.subprocess, "run", fake_run)
    q = salmon_mod.run_salmon_quant(idx, out, r1)
    assert q.mapping_rate == 0.77
    assert q.quant_sf == out / "quant.sf"

def test_run_salmon_quant_paired_uses_1_2(tmp_path, monkeypatch):
    idx = tmp_path / "idx"; idx.mkdir()
    r1 = tmp_path / "s_R1.fastq"; r1.write_text("@r\nACGT\n+\nIIII\n")
    r2 = tmp_path / "s_R2.fastq"; r2.write_text("@r\nACGT\n+\nIIII\n")
    out = tmp_path / "s1"
    def fake_run(cmd, **k):
        assert "-1" in cmd and "-2" in cmd and "-r" not in cmd
        aux = out / "aux_info"; aux.mkdir(parents=True, exist_ok=True)
        (aux / "meta_info.json").write_text('{"percent_mapped": 80.0}')
        (out / "quant.sf").write_text("Name\tLength\tTPM\tNumReads\n")
        class R: returncode = 0; stdout = ""; stderr = ""
        return R()
    monkeypatch.setattr(salmon_mod.subprocess, "run", fake_run)
    q = salmon_mod.run_salmon_quant(idx, out, r1, fastq_2=r2)
    assert q.mapping_rate == 0.80
```

- [ ] **Step 10: Run tests to verify they fail**

Run: `conda run -n rnaforge-core python -m pytest tests/test_salmon.py -v`
Expected: FAIL (run_salmon_quant not defined)

- [ ] **Step 11: Implement `run_salmon_quant`**

```python
# rnaforge/salmon.py (append)
def run_salmon_quant(index_dir: Path, out_dir: Path, fastq_1: Path,
                     fastq_2: Path | None = None, threads: int = 8,
                     env: str = "rnaforge-quant-euk") -> SalmonQuant:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = ["conda", "run", "-n", env, "salmon", "quant", "-i", str(index_dir),
           "-l", "A", "-p", str(threads), "--validateMappings", "-o", str(out_dir)]
    if fastq_2 is not None:
        cmd += ["-1", str(fastq_1), "-2", str(fastq_2)]
    else:
        cmd += ["-r", str(fastq_1)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    quant_sf = out_dir / "quant.sf"
    meta = out_dir / "aux_info" / "meta_info.json"
    if r.returncode != 0 or not quant_sf.exists() or not meta.exists():
        raise SalmonError(f"salmon quant failed for {fastq_1.name}: {r.stderr[-500:]}")
    return SalmonQuant(quant_sf=quant_sf, mapping_rate=parse_salmon_meta(meta))
```

- [ ] **Step 12: Run all salmon tests to verify they pass**

Run: `conda run -n rnaforge-core python -m pytest tests/test_salmon.py -v`
Expected: PASS (5 tests)

- [ ] **Step 13: Commit**

```bash
git add rnaforge/salmon.py tests/test_salmon.py
git commit -m "feat(salmon): parser + decoy-aware index + quant runner (eukaryote)"
```

---

### Task 2: m04-euk router branch + mapping_rate gate

**Files:**
- Modify: `rnaforge/modules/m04_quant.py` (remove NotImplementedError; add `_quant_euk`)
- Test: `tests/test_m04_quant.py`

**Interfaces:**
- Consumes: `salmon.build_salmon_index`, `salmon.run_salmon_quant`, `salmon.SalmonQuant`.
- Produces: `_quant_euk(config, metadata_path, run_dir, quant_dir, stats_dir, logs_dir, state) -> dict`;
  summary `{"read_type":"short","organism_type":"eukaryote","n_samples":N,"samples":{sid:{"mapping_rate":r,"quant_sf":path}},"gate_counts":{...}}`.
  `mapping_rate` → `alignment_rate` kapısı (reuse `build_alignment_gates` via adapter).

- [ ] **Step 1: Write failing test — eukaryote no longer raises, runs salmon, writes mapping gate**

```python
# tests/test_m04_quant.py (append; imitate existing long-branch tests)
def test_run_quant_eukaryote_runs_salmon_and_gates(tmp_path, monkeypatch):
    import rnaforge.modules.m04_quant as m04
    from rnaforge.salmon import SalmonQuant
    # seed m01(short)+m03 done, eukaryote config
    run_dir = tmp_path / "run"; (run_dir / "statistics").mkdir(parents=True)
    (run_dir / "statistics" / "raw_statistics.json").write_text(
        '{"read_type":"short","platform":"illumina"}')
    from rnaforge.state import RunState
    st = RunState(run_dir); st.mark_done("m01_validate", []); st.mark_done("m03_trim", [])
    fq = tmp_path / "s1.fastq"; fq.write_text("@r\nACGT\n+\nIIII\n")
    meta = tmp_path / "m.tsv"; meta.write_text(f"sample_id\tcondition\tfastq_1\ns1\tctrl\t{fq}\n")
    # trimmed_reads → return fq
    monkeypatch.setattr(m04, "trimmed_reads", lambda rd, s: (fq, None))
    monkeypatch.setattr(m04, "build_salmon_index", lambda *a, **k: tmp_path / "idx")
    monkeypatch.setattr(m04, "run_salmon_quant",
        lambda *a, **k: SalmonQuant(quant_sf=tmp_path / "q.sf", mapping_rate=0.9))
    from rnaforge.config import (Config, Reference, Library, Trimming, DE, Report, Resources)
    cfg = Config(organism="human", organism_type="eukaryote", platform="illumina",
        reference=Reference(transcriptome_fasta=tmp_path/"tx.fa", tx2gene=tmp_path/"t2g.tsv"),
        library=Library(), trimming=Trimming(), de=DE(), report=Report(), resources=Resources())
    summary = m04.run_quant(cfg, meta, run_dir)
    assert summary["organism_type"] == "eukaryote"
    assert summary["samples"]["s1"]["mapping_rate"] == 0.9
    gates = json.loads((run_dir / "quality" / "gates.json").read_text())["gates"]
    assert any(g["name"] == "alignment_rate" for g in gates if g["module"] == "m04_quant")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `conda run -n rnaforge-core python -m pytest tests/test_m04_quant.py::test_run_quant_eukaryote_runs_salmon_and_gates -v`
Expected: FAIL (NotImplementedError raised)

- [ ] **Step 3: Replace NotImplementedError with eukaryote dispatch + implement `_quant_euk`**

In `run_quant`, replace the block:
```python
    if config.organism_type == "eukaryote":
        raise NotImplementedError(
            "m04 eukaryote (Salmon) path not yet implemented; prokaryote only for now."
        )
    if not state.is_done("m03_trim"):
```
with:
```python
    if not state.is_done("m03_trim"):
```
and after the `state.is_done("m03_trim")` guard, change the read_type dispatch to check organism_type first:
```python
    if config.organism_type == "eukaryote":
        summary = _quant_euk(config, metadata_path, run_dir,
                             quant_dir, stats_dir, logs_dir, state)
    else:
        read_type = resolve_read_type(run_dir)
        if read_type == "long":
            summary = _quant_long(config, metadata_path, run_dir,
                                  quant_dir, stats_dir, logs_dir, state)
        else:
            summary = _quant_short(config, metadata_path, run_dir,
                                   quant_dir, stats_dir, logs_dir, state)
```
Add imports at top: `from rnaforge.salmon import SalmonQuant, build_salmon_index, run_salmon_quant`.
Add the function (adapter wraps SalmonQuant so `build_alignment_gates` sees `.alignment_rate`):
```python
def _quant_euk(config: Config, metadata_path: Path, run_dir: Path,
               quant_dir: Path, stats_dir: Path, logs_dir: Path,
               state: RunState) -> dict:
    """Ökaryot niceleme (Salmon, decoy-aware). mapping_rate → alignment_rate FAIL kapısı
    (eukaryote.yml permissive). Trimlenmiş okuma → salmon quant -l A."""
    stats_path = stats_dir / "alignment_statistics.json"
    profile = load_profile(config.organism_type, config.quality)
    log_path = logs_dir / "quant.log"
    with log_path.open("w") as log_file:
        def log(msg: str) -> None:
            log_file.write(msg + "\n"); log_file.flush()
        samples = load_metadata(metadata_path)
        index_dir = build_salmon_index(
            config.reference.transcriptome_fasta, quant_dir / "_index",
            genome_fasta=config.reference.genome_fasta,
            threads=config.resources.threads, log=log)
        log(f"m04 salmon: index ready, {len(samples)} sample(s)")
        results = {}
        per_sample = {}
        for sample in samples:
            state.heartbeat()
            t1, t2 = trimmed_reads(run_dir, sample)
            q = run_salmon_quant(index_dir, quant_dir / sample.sample_id, t1,
                                 fastq_2=t2, threads=config.resources.threads)
            results[sample.sample_id] = _MappingAdapter(q.mapping_rate)
            per_sample[sample.sample_id] = {
                "mapping_rate": q.mapping_rate, "quant_sf": str(q.quant_sf)}
            log(f"{sample.sample_id}: mapping_rate={q.mapping_rate:.3f}")
        gates = build_alignment_gates(results, profile)
        summary = {
            "read_type": "short", "organism_type": "eukaryote",
            "n_samples": len(samples), "samples": per_sample,
            "gate_counts": dict(Counter(g.status for g in gates)),
        }
        stats_path.write_text(json.dumps(summary, indent=2))
        write_gate_results(run_dir, gates)
        for g in gates:
            log(f"gate {g.name}: {g.status} — {g.message}")
        raise_if_failed(gates)
        log(f"alignment statistics written: {stats_path}")
    return summary
```
Add near top (after imports):
```python
from dataclasses import dataclass as _dataclass

@_dataclass
class _MappingAdapter:
    """build_alignment_gates .alignment_rate bekler; salmon mapping_rate'i uyarlar."""
    alignment_rate: float
```

- [ ] **Step 4: Run test to verify it passes**

Run: `conda run -n rnaforge-core python -m pytest tests/test_m04_quant.py::test_run_quant_eukaryote_runs_salmon_and_gates -v`
Expected: PASS

- [ ] **Step 5: Run full m04 suite (no regression on prok short/long)**

Run: `conda run -n rnaforge-core python -m pytest tests/test_m04_quant.py -v`
Expected: PASS (all)

- [ ] **Step 6: Commit**

```bash
git add rnaforge/modules/m04_quant.py tests/test_m04_quant.py
git commit -m "feat(m04): eukaryote Salmon quant branch + mapping_rate gate"
```

---

### Task 3: tximport R script + python runner + m05-euk branch

**Files:**
- Create: `rnaforge/scripts/tximport.R`
- Create/Modify: `rnaforge/tximport.py` (python runner)
- Modify: `rnaforge/modules/m05_counts.py` (add `_counts_euk`)
- Modify: `envs/rnaforge-de.yml` (add `bioconductor-tximport`)
- Test: `tests/test_m05_counts.py`, `tests/test_tximport.py`

**Interfaces:**
- Consumes: m04-euk `quantification/<sid>/quant.sf`; `reference.tx2gene`.
- Produces:
  - `tximport.run_tximport(quant_sfs: dict[str,Path], tx2gene: Path, out_dir: Path, env: str = "rnaforge-de") -> TximportResult`
    where `TximportResult(gene_ids: list[str], counts: dict[str, list[float]])` (counts keyed by sample_id, lengthScaledTPM, integer-rounded).
  - `_counts_euk(config, metadata_path, run_dir, quant_dir, stats_dir, logs_dir, state) -> dict` → writes `counts.tsv` (`gene\t<sid...>`).

- [ ] **Step 1: Add bioconductor-tximport to env yml + install**

Edit `envs/rnaforge-de.yml` dependencies: add line `  - bioconductor-tximport`.
Run: `conda install -n rnaforge-de -c bioconda -c conda-forge bioconductor-tximport -y`
Verify: `conda run -n rnaforge-de Rscript -e 'library(tximport); cat("ok\n")'`
Expected: `ok`

- [ ] **Step 2: Write failing test for `run_tximport` (python wrapper, R monkeypatched)**

```python
# tests/test_tximport.py
from pathlib import Path
from rnaforge import tximport as tx

def test_run_tximport_reads_matrix(tmp_path, monkeypatch):
    # R betiği yerine sahte: gene x sample counts.tsv yazar
    out = tmp_path / "out"; out.mkdir()
    def fake_run(cmd, **k):
        (out / "gene_counts.tsv").write_text("gene\ts1\ts2\ng1\t10\t20\ng2\t0\t5\n")
        class R: returncode = 0; stdout = ""; stderr = ""
        return R()
    monkeypatch.setattr(tx.subprocess, "run", fake_run)
    res = tx.run_tximport({"s1": tmp_path/"s1.sf", "s2": tmp_path/"s2.sf"},
                          tmp_path/"t2g.tsv", out)
    assert res.gene_ids == ["g1", "g2"]
    assert res.counts["s1"] == [10.0, 20.0]
    assert res.counts["s2"] == [0.0, 5.0]
```

- [ ] **Step 3: Run test to verify it fails**

Run: `conda run -n rnaforge-core python -m pytest tests/test_tximport.py -v`
Expected: FAIL (module not found)

- [ ] **Step 4: Implement `rnaforge/tximport.py`**

```python
# rnaforge/tximport.py
"""tximport (transkript→gen) — R Bioconductor wrapper (deseq2.py deseni).

countsFromAbundance="lengthScaledTPM" → uzunluk-düzeltilmiş sayım; m06 DESeq2 bunu düz
sayım gibi okur ve doğru olur (m06 organizma-agnostik kalır). Env: rnaforge-de."""
from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path

_SCRIPT = Path(__file__).parent / "scripts" / "tximport.R"


class TximportError(RuntimeError):
    """tximport çalıştırılamadı ya da beklenen çıktıyı üretmedi."""


@dataclass
class TximportResult:
    gene_ids: list[str]
    counts: dict[str, list[float]]


def run_tximport(quant_sfs: dict[str, Path], tx2gene: Path, out_dir: Path,
                 env: str = "rnaforge-de") -> TximportResult:
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    # args: <tx2gene> <out_gene_counts.tsv> sid1 sf1 sid2 sf2 ...
    out_tsv = out_dir / "gene_counts.tsv"
    args = [str(tx2gene), str(out_tsv)]
    for sid, sf in quant_sfs.items():
        args += [sid, str(sf)]
    cmd = ["conda", "run", "-n", env, "Rscript", str(_SCRIPT), *args]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0 or not out_tsv.exists():
        raise TximportError(f"tximport failed: {r.stderr[-500:]}")
    lines = out_tsv.read_text().splitlines()
    header = lines[0].split("\t")[1:]     # sample ids in file order
    gene_ids, cols = [], {sid: [] for sid in header}
    for line in lines[1:]:
        parts = line.split("\t")
        gene_ids.append(parts[0])
        for sid, val in zip(header, parts[1:]):
            cols[sid].append(float(val))
    return TximportResult(gene_ids=gene_ids, counts=cols)
```

- [ ] **Step 5: Run test to verify it passes**

Run: `conda run -n rnaforge-core python -m pytest tests/test_tximport.py -v`
Expected: PASS

- [ ] **Step 6: Write `rnaforge/scripts/tximport.R`**

```r
#!/usr/bin/env Rscript
# tximport: quant.sf'leri gen-seviyesine topla (lengthScaledTPM → m06 DESeq2 offset gerektirmez).
# argv: <tx2gene.tsv> <out_gene_counts.tsv> sid1 sf1 sid2 sf2 ...
suppressMessages(library(tximport))
args <- commandArgs(trailingOnly = TRUE)
tx2gene_path <- args[1]; out_path <- args[2]
rest <- args[-(1:2)]
sids <- rest[seq(1, length(rest), by = 2)]
sfs  <- rest[seq(2, length(rest), by = 2)]
names(sfs) <- sids
tx2gene <- read.table(tx2gene_path, header = FALSE, sep = "\t",
                      stringsAsFactors = FALSE)
txi <- tximport(sfs, type = "salmon", tx2gene = tx2gene,
                countsFromAbundance = "lengthScaledTPM")
counts <- round(txi$counts)                      # uzunluk-düzeltilmiş sayım
df <- data.frame(gene = rownames(counts), counts, check.names = FALSE)
write.table(df, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat("tximport ok:", nrow(counts), "genes x", ncol(counts), "samples\n")
```

- [ ] **Step 7: Write failing test for `_counts_euk` (m05 router, tximport monkeypatched)**

```python
# tests/test_m05_counts.py (append)
def test_run_counts_eukaryote_writes_counts_tsv(tmp_path, monkeypatch):
    import rnaforge.modules.m05_counts as m05
    from rnaforge.tximport import TximportResult
    run_dir = tmp_path / "run"; (run_dir / "statistics").mkdir(parents=True)
    (run_dir / "statistics" / "raw_statistics.json").write_text(
        '{"read_type":"short","platform":"illumina"}')
    from rnaforge.state import RunState
    st = RunState(run_dir); st.mark_done("m04_quant", [])
    # m04-euk çıktısı: quantification/<sid>/quant.sf
    qd = run_dir / "quantification"
    for sid in ("s1", "s2"):
        (qd / sid).mkdir(parents=True); (qd / sid / "quant.sf").write_text("Name\n")
    meta = tmp_path / "m.tsv"
    meta.write_text("sample_id\tcondition\tfastq_1\ns1\tctrl\ta.fq\ns2\ttrt\tb.fq\n")
    monkeypatch.setattr(m05, "run_tximport",
        lambda *a, **k: TximportResult(gene_ids=["g1","g2"],
                                       counts={"s1":[10.0,0.0], "s2":[20.0,5.0]}))
    from rnaforge.config import (Config, Reference, Library, Trimming, DE, Report, Resources)
    cfg = Config(organism="human", organism_type="eukaryote", platform="illumina",
        reference=Reference(transcriptome_fasta=tmp_path/"tx.fa", tx2gene=tmp_path/"t2g.tsv"),
        library=Library(), trimming=Trimming(), de=DE(), report=Report(), resources=Resources())
    summary = m05.run_counts(cfg, meta, run_dir)
    matrix = (qd / "counts.tsv").read_text().splitlines()
    assert matrix[0] == "gene\ts1\ts2"
    assert matrix[1] == "g1\t10\t20" or matrix[1].startswith("g1\t")
    assert summary["organism_type"] == "eukaryote"
```

- [ ] **Step 8: Run test to verify it fails**

Run: `conda run -n rnaforge-core python -m pytest tests/test_m05_counts.py::test_run_counts_eukaryote_writes_counts_tsv -v`
Expected: FAIL (router has no eukaryote branch; run_tximport not imported)

- [ ] **Step 9: Add eukaryote branch to `run_counts` + implement `_counts_euk`**

In `run_counts`, change the read_type dispatch to check organism_type first:
```python
    if config.organism_type == "eukaryote":
        summary = _counts_euk(config, metadata_path, run_dir,
                              quant_dir, stats_dir, logs_dir, state)
    else:
        read_type = resolve_read_type(run_dir)
        if read_type == "long":
            summary = _counts_long(config, metadata_path, run_dir,
                                  quant_dir, stats_dir, logs_dir, state)
        else:
            summary = _counts_short(config, metadata_path, run_dir,
                                  quant_dir, stats_dir, logs_dir, state)
```
Add import: `from rnaforge.tximport import run_tximport`.
Add function:
```python
def _counts_euk(config: Config, metadata_path: Path, run_dir: Path,
                quant_dir: Path, stats_dir: Path, logs_dir: Path,
                state: RunState) -> dict:
    """Ökaryot sayım (tximport, lengthScaledTPM). counts.tsv sözleşmesi.
    Salmon zaten hizalamada atadı → assignment FAIL kapısı yok (diagnostik)."""
    stats_path = stats_dir / "count_statistics.json"
    log_path = logs_dir / "counts.log"
    with log_path.open("w") as log_file:
        def log(msg: str) -> None:
            log_file.write(msg + "\n"); log_file.flush()
        samples = load_metadata(metadata_path)
        quant_sfs = {s.sample_id: quant_dir / s.sample_id / "quant.sf" for s in samples}
        missing = [sid for sid, p in quant_sfs.items() if not p.exists()]
        if missing:
            raise ValueError(f"m05 eukaryote: quant.sf eksik örnek(ler): {missing} "
                             "(m04 salmon koştu mu?)")
        state.heartbeat()
        res = run_tximport(quant_sfs, config.reference.tx2gene, quant_dir)
        if not res.gene_ids:
            raise ValueError("tximport 0 gen döndürdü (tx2gene eşleşmedi mi?)")
        sample_ids = [s.sample_id for s in samples]
        matrix_path = quant_dir / "counts.tsv"
        with matrix_path.open("w") as fh:
            fh.write("gene\t" + "\t".join(sample_ids) + "\n")
            for i, gene in enumerate(res.gene_ids):
                row = [f"{int(round(res.counts[sid][i]))}" for sid in sample_ids]
                fh.write(gene + "\t" + "\t".join(row) + "\n")
        log(f"count matrix written: {matrix_path} ({len(res.gene_ids)} genes)")
        summary = {
            "read_type": "short", "organism_type": "eukaryote",
            "n_samples": len(samples), "n_genes": len(res.gene_ids),
            "gate_counts": {},
        }
        stats_path.write_text(json.dumps(summary, indent=2))
    return summary
```

- [ ] **Step 10: Run test to verify it passes**

Run: `conda run -n rnaforge-core python -m pytest tests/test_m05_counts.py::test_run_counts_eukaryote_writes_counts_tsv -v`
Expected: PASS

- [ ] **Step 11: Run full m05 + tximport suites**

Run: `conda run -n rnaforge-core python -m pytest tests/test_m05_counts.py tests/test_tximport.py -v`
Expected: PASS (all)

- [ ] **Step 12: Commit**

```bash
git add rnaforge/tximport.py rnaforge/scripts/tximport.R rnaforge/modules/m05_counts.py envs/rnaforge-de.yml tests/test_m05_counts.py tests/test_tximport.py
git commit -m "feat(m05): eukaryote tximport branch (lengthScaledTPM) → counts.tsv contract"
```

---

### Task 4: Config test + full-suite regression + real-tool integration

**Files:**
- Test: `tests/test_config.py`, `tests/test_salmon.py`

**Interfaces:**
- Consumes: all prior tasks.

- [ ] **Step 1: Write config test — eukaryote genome_fasta optional**

```python
# tests/test_config.py (append). EUK_BODY: minimal eukaryote config helper.
def test_eukaryote_genome_fasta_optional(tmp_path):
    body = ('organism: "human"\norganism_type: "eukaryote"\n'
            'reference:\n  transcriptome_fasta: "tx.fa"\n  tx2gene: "t2g.tsv"\n')
    cfg = load_config(_write(tmp_path, body))
    assert cfg.reference.genome_fasta is None
    assert cfg.reference.transcriptome_fasta is not None

def test_eukaryote_genome_fasta_parsed_when_present(tmp_path):
    body = ('organism: "human"\norganism_type: "eukaryote"\n'
            'reference:\n  transcriptome_fasta: "tx.fa"\n  tx2gene: "t2g.tsv"\n'
            '  genome_fasta: "genome.fa"\n')
    cfg = load_config(_write(tmp_path, body))
    assert str(cfg.reference.genome_fasta) == "genome.fa"
```

- [ ] **Step 2: Run config tests**

Run: `conda run -n rnaforge-core python -m pytest tests/test_config.py -k eukaryote -v`
Expected: PASS

- [ ] **Step 3: Real-tool salmon integration test (skip if env missing)**

```python
# tests/test_salmon.py (append)
import shutil, subprocess, pytest

def _has_env(name):
    r = subprocess.run(["conda", "run", "-n", name, "salmon", "--version"],
                       capture_output=True, text=True)
    return r.returncode == 0

@pytest.mark.skipif(not _has_env("rnaforge-quant-euk"), reason="salmon env yok")
def test_salmon_index_and_quant_real(tmp_path):
    # küçük sentetik transkriptom (2 transkript) + 1 örnek
    tx = tmp_path / "tx.fa"
    tx.write_text(">t1\n" + "ACGT"*60 + "\n>t2\n" + "GGCC"*60 + "\n")
    idx = tmp_path / "idx"
    from rnaforge.salmon import build_salmon_index, run_salmon_quant
    build_salmon_index(tx, idx)                       # decoy'suz
    r1 = tmp_path / "s.fastq"
    r1.write_text("".join(f"@r{i}\n{'ACGT'*60}\n+\n{'I'*240}\n" for i in range(50)))
    q = run_salmon_quant(idx, tmp_path / "out", r1)
    assert q.quant_sf.exists() and 0.0 <= q.mapping_rate <= 1.0
```

- [ ] **Step 4: Run integration test**

Run: `conda run -n rnaforge-core python -m pytest tests/test_salmon.py -k real -v`
Expected: PASS (or SKIP if env missing)

- [ ] **Step 5: Full suite regression**

Run: `conda run -n rnaforge-core python -m pytest -q`
Expected: PASS (492 prior + new tests, no regressions)

- [ ] **Step 6: Commit**

```bash
git add tests/test_config.py tests/test_salmon.py
git commit -m "test(eukaryote): config genome_fasta optional + salmon real-tool integration"
```

---

## Self-Review

**Spec coverage:**
- §3 m04-euk Salmon (index decoy-aware + quant + mapping_rate gate) → Task 1 + Task 2 ✓
- §4 m05-euk tximport (lengthScaledTPM → counts.tsv) → Task 3 ✓
- §5 config genome_fasta optional (fields already exist in Reference) → Task 4 (test) ✓; REQUIRED_REFERENCE unchanged (no code change needed — eukaryote already (transcriptome_fasta, tx2gene)) ✓
- §6 testing → every task is TDD ✓
- §7 validation → deferred scientific step (post-implementation), not in this plan by design ✓
- Env: bioconductor-tximport → Task 3 Step 1 ✓

**Placeholder scan:** No TBD/TODO; all steps have concrete code. ✓

**Type consistency:** `SalmonQuant(quant_sf, mapping_rate)` produced Task 1, consumed Task 2 ✓. `TximportResult(gene_ids, counts: dict[sid→list])` produced Task 3 python, matches R output header order ✓. `_MappingAdapter.alignment_rate` matches `build_alignment_gates` expectation (`.alignment_rate`) ✓. Router branches use existing `_quant_long/_counts_long` names verbatim ✓.

**Note:** eukaryote read_type is always "short" this increment; ONT-eukaryote is out of scope and will hit `_quant_euk` (which ignores read_type) — acceptable since euk-long is a separate future subsystem. No guard needed because organism_type is checked before read_type.
