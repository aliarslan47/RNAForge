# Ökaryot Uzun-Okuma (gen düzeyi, transkriptom) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development or superpowers:executing-plans. Steps use `- [ ]`.

**Goal:** Ökaryot ONT/PacBio uzun-okumayı gen×örnek `counts.tsv`'ye getir (minimap2→transkriptom → primer-sayım → tx2gene topla), m06+ değişmeden.

**Architecture:** m04/m05 eukaryote dalı içinde read_type alt-dallanması. long → minimap2 transkriptoma (run_minimap2 reuse, preset platformdan, diagnostik) + primer-alignment sayımı → tx2gene aggregate. Referans = transcriptome_fasta+tx2gene (mevcut).

**Tech Stack:** Python 3.11, minimap2/samtools (rnaforge-longread), pytest.

**Spec:** `docs/superpowers/specs/2026-08-18-eukaryote-longread-design.md`

## Global Constraints

- Ayrım YALNIZ m04/m05; m06+ DEĞİŞMEZ. counts.tsv = `gene\t<sid...>`, sütun→sample_id KONUMLA.
- Long = DIAGNOSTİK (FAIL kapısı yok). Referans transcriptome_fasta+tx2gene; genom/GTF yok.
- Primer sayım = `samtools view -F 2308` (unmapped+secondary+supplementary hariç).
- Boş matris → yüksek sesle hata.

---

### Task 1: sayım + tx2gene yardımcıları

**Files:** Modify `rnaforge/minimap2.py`, `rnaforge/tximport.py`; Test `tests/test_minimap2.py`, `tests/test_tximport.py`

**Interfaces produced:**
- `minimap2.count_primary_alignments(bam_path, env="rnaforge-longread") -> dict[str,int]` (transcript_id → primer okuma sayısı)
- `tximport.parse_tx2gene(path) -> dict[str,str]` (tx_id → gene_id)

- [ ] **Step 1: parse_tx2gene testi (failing)**
```python
# tests/test_tximport.py (append)
def test_parse_tx2gene(tmp_path):
    from rnaforge.tximport import parse_tx2gene
    p = tmp_path / "t2g.tsv"; p.write_text("ENST1.1\tENSG1.1\nENST2.2\tENSG1.1\nENST3.1\tENSG2.4\n")
    assert parse_tx2gene(p) == {"ENST1.1":"ENSG1.1","ENST2.2":"ENSG1.1","ENST3.1":"ENSG2.4"}
```
- [ ] **Step 2: fails** — `conda run -n rnaforge-core python -m pytest tests/test_tximport.py::test_parse_tx2gene -v`
- [ ] **Step 3: implement parse_tx2gene**
```python
# rnaforge/tximport.py (append)
def parse_tx2gene(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    for line in Path(path).read_text().splitlines():
        if not line.strip(): continue
        parts = line.split("\t")
        if len(parts) >= 2:
            out[parts[0]] = parts[1]
    return out
```
- [ ] **Step 4: passes**
- [ ] **Step 5: count_primary_alignments — real-tool test (skip if env yok)**
```python
# tests/test_minimap2.py (append)
import subprocess, pytest
def _has_lr(): 
    return subprocess.run(["conda","run","-n","rnaforge-longread","samtools","--version"],
                          capture_output=True,text=True).returncode == 0
@pytest.mark.skipif(not _has_lr(), reason="rnaforge-longread yok")
def test_count_primary_alignments(tmp_path):
    from rnaforge.minimap2 import count_primary_alignments
    sam = tmp_path / "a.sam"
    sam.write_text(
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:tx1\tLN:100\n@SQ\tSN:tx2\tLN:100\n"
        "r1\t0\ttx1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n"       # primary tx1
        "r2\t0\ttx1\t5\t60\t4M\t*\t0\t0\tACGT\tIIII\n"       # primary tx1
        "r3\t0\ttx2\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n"       # primary tx2
        "r3\t256\ttx1\t1\t0\t4M\t*\t0\t0\tACGT\tIIII\n"      # secondary -> sayılmaz
        "r4\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n")          # unmapped -> sayılmaz
    env="rnaforge-longread"
    bam = tmp_path / "a.bam"
    subprocess.run(["conda","run","-n",env,"bash","-c",
                    f"samtools sort -o {bam} {sam} && samtools index {bam}"],check=True)
    assert count_primary_alignments(bam) == {"tx1": 2, "tx2": 1}
```
- [ ] **Step 6: fails**
- [ ] **Step 7: implement count_primary_alignments**
```python
# rnaforge/minimap2.py (append)
def count_primary_alignments(bam_path: Path, env: str = "rnaforge-longread") -> dict[str, int]:
    """Primer hizalama sayımı hedef (transkript) başına: -F 2308 = unmapped(4)+
    secondary(256)+supplementary(2048) hariç → okuma başına tek satır. Sütun 3 = hedef."""
    r = _run(["conda", "run", "-n", env, "samtools", "view", "-F", "2308", str(bam_path)])
    if r.returncode != 0:
        raise Minimap2RunError(f"samtools view failed: {r.stderr.strip()[-500:]}")
    counts: dict[str, int] = {}
    for line in r.stdout.splitlines():
        if not line: continue
        ref = line.split("\t")[2]
        if ref and ref != "*":
            counts[ref] = counts.get(ref, 0) + 1
    return counts
```
- [ ] **Step 8: passes** (or skip)
- [ ] **Step 9: Commit + push**
```bash
git add rnaforge/minimap2.py rnaforge/tximport.py tests/test_minimap2.py tests/test_tximport.py
git commit -m "feat(euk-long): count_primary_alignments + parse_tx2gene helpers"
git push origin main
```

---

### Task 2: m04 `_quant_euk_long` + read_type alt-dispatch

**Files:** Modify `rnaforge/modules/m04_quant.py`; Test `tests/test_m04_quant.py`

**Interfaces:** Consumes run_minimap2, minimap2_preset, resolve_platform. Produces `_quant_euk_long(...) -> dict` (read_type=long, organism_type=eukaryote, samples[sid].mapping_rate/bam, gate_counts boş — diagnostik).

- [ ] **Step 1: test (failing)** — eukaryote+long → minimap2, mapping_rate diagnostik (FAIL yok)
```python
# tests/test_m04_quant.py (append)
def test_run_quant_eukaryote_long_minimap2_diagnostic(tmp_path, monkeypatch):
    import rnaforge.modules.m04_quant as m04
    from rnaforge.bowtie2 import AlignmentResult
    from rnaforge.state import RunState
    run_dir = tmp_path / "run"; (run_dir / "statistics").mkdir(parents=True)
    (run_dir / "statistics" / "raw_statistics.json").write_text(
        '{"read_type":"long","platform":"ont","chemistry":"cdna"}')
    st = RunState(run_dir); st.mark_done("m01_validate", []); st.mark_done("m03_trim", [])
    fq = tmp_path / "s1.fastq"; fq.write_text("@r\nACGT\n+\nIIII\n")
    meta = tmp_path / "m.tsv"; meta.write_text(f"sample_id\tcondition\tfastq_1\ns1\tctrl\t{fq}\n")
    monkeypatch.setattr(m04, "trimmed_reads", lambda rd, s: (fq, None))
    bam = tmp_path / "b.bam"; bam.write_text("")
    monkeypatch.setattr(m04, "run_minimap2",
        lambda *a, **k: AlignmentResult(alignment_rate=0.55, bam=bam))
    from rnaforge.config import (Config, Reference, Library, Trimming, DE, Report, Resources)
    cfg = Config(organism="human", organism_type="eukaryote", platform="ont",
        reference=Reference(transcriptome_fasta=tmp_path/"tx.fa", tx2gene=tmp_path/"t2g.tsv"),
        library=Library(chemistry="cdna"), trimming=Trimming(), de=DE(), report=Report(),
        resources=Resources())
    summary = m04.run_quant(cfg, meta, run_dir)
    assert summary["organism_type"] == "eukaryote" and summary["read_type"] == "long"
    assert summary["samples"]["s1"]["mapping_rate"] == 0.55
    gates = __import__("json").loads((run_dir/"quality"/"gates.json").read_text())["gates"] if (run_dir/"quality"/"gates.json").exists() else []
    assert all(g.get("status") != "FAIL" for g in gates)   # diagnostik, FAIL yok
```
- [ ] **Step 2: fails** (eukaryote → _quant_euk salmon çağrılır, hata)
- [ ] **Step 3: implement** — eukaryote dalını read_type'a böl + `_quant_euk_long`
Replace eukaryote dispatch:
```python
    if config.organism_type == "eukaryote":
        read_type = resolve_read_type(run_dir)
        if read_type == "long":
            summary = _quant_euk_long(config, metadata_path, run_dir,
                                      quant_dir, stats_dir, logs_dir, state)
        else:
            summary = _quant_euk(config, metadata_path, run_dir,
                                 quant_dir, stats_dir, logs_dir, state)
    else:
        read_type = resolve_read_type(run_dir)
        if read_type == "long":
            summary = _quant_long(...)   # mevcut
        else:
            summary = _quant_short(...)  # mevcut
```
Add function:
```python
def _quant_euk_long(config, metadata_path, run_dir, quant_dir, stats_dir, logs_dir, state) -> dict:
    """Ökaryot uzun-okuma: minimap2 → transkriptom (transkript=hedef, splice yok).
    Preset platformdan. DIAGNOSTİK — FAIL kapısı yok (tüm long yolları gibi)."""
    platform = resolve_platform(run_dir)
    preset = minimap2_preset(platform)
    stats_path = stats_dir / "alignment_statistics.json"
    log_path = logs_dir / "quant.log"
    with log_path.open("w") as log_file:
        def log(m): log_file.write(m + "\n"); log_file.flush()
        samples = load_metadata(metadata_path)
        log(f"m04 eukaryote long minimap2 ({platform}→-ax {preset}) → transcriptome: {len(samples)} sample(s)")
        per_sample = {}
        for sample in samples:
            state.heartbeat()
            t1, _ = trimmed_reads(run_dir, sample)
            result = run_minimap2(config.reference.transcriptome_fasta,
                                  quant_dir / sample.sample_id, t1, preset,
                                  threads=config.resources.threads)
            per_sample[sample.sample_id] = {
                "mapping_rate": result.alignment_rate, "bam": str(result.bam)}
            log(f"{sample.sample_id}: mapping_rate={result.alignment_rate:.3f} (diagnostik)")
        summary = {"read_type": "long", "organism_type": "eukaryote",
                   "n_samples": len(samples), "samples": per_sample, "gate_counts": {}}
        stats_path.write_text(json.dumps(summary, indent=2))
        log(f"alignment statistics written: {stats_path}")
    return summary
```
- [ ] **Step 4: passes**
- [ ] **Step 5: full m04 suite** — no regression
- [ ] **Step 6: Commit + push**

---

### Task 3: m05 `_counts_euk_long` + read_type alt-dispatch

**Files:** Modify `rnaforge/modules/m05_counts.py`; Test `tests/test_m05_counts.py`

**Interfaces:** Consumes count_primary_alignments, parse_tx2gene. Produces `_counts_euk_long(...) -> dict` writing counts.tsv.

- [ ] **Step 1: test (failing)** — eukaryote+long → transcript counts aggregate to gene
```python
# tests/test_m05_counts.py (append)
def test_run_counts_eukaryote_long_aggregates_tx_to_gene(tmp_path, monkeypatch):
    import rnaforge.modules.m05_counts as m05
    from rnaforge.state import RunState
    run_dir = tmp_path/"run"; (run_dir/"statistics").mkdir(parents=True)
    (run_dir/"statistics"/"raw_statistics.json").write_text('{"read_type":"long","platform":"ont"}')
    st = RunState(run_dir); st.mark_done("m04_quant", [])
    qd = run_dir/"quantification"
    for sid in ("s1","s2"):
        (qd/sid).mkdir(parents=True); (qd/sid/"aligned.sorted.bam").write_text("")
    t2g = tmp_path/"t2g.tsv"; t2g.write_text("ENST1.1\tENSG1\nENST2.2\tENSG1\nENST3.1\tENSG2\n")
    a=tmp_path/"a.fq"; a.write_text("@r\nA\n+\nI\n"); b=tmp_path/"b.fq"; b.write_text("@r\nA\n+\nI\n")
    meta=tmp_path/"m.tsv"; meta.write_text(f"sample_id\tcondition\tfastq_1\ns1\tctrl\t{a}\ns2\ttrt\t{b}\n")
    # s1: ENST1=3,ENST3=1 -> ENSG1=3,ENSG2=1 ; s2: ENST2=5 -> ENSG1=5
    cnts = {"s1":{"ENST1.1":3,"ENST3.1":1}, "s2":{"ENST2.2":5}}
    monkeypatch.setattr(m05, "count_primary_alignments", lambda bam, **k: cnts["s1" if "s1" in str(bam) else "s2"])
    from rnaforge.config import (Config, Reference, Library, Trimming, DE, Report, Resources)
    cfg = Config(organism="human", organism_type="eukaryote", platform="ont",
        reference=Reference(transcriptome_fasta=tmp_path/"tx.fa", tx2gene=t2g),
        library=Library(), trimming=Trimming(), de=DE(), report=Report(), resources=Resources())
    summary = m05.run_counts(cfg, meta, run_dir)
    lines = (qd/"counts.tsv").read_text().splitlines()
    assert lines[0] == "gene\ts1\ts2"
    d = {l.split("\t")[0]: l.split("\t")[1:] for l in lines[1:]}
    assert d["ENSG1"] == ["3","5"] and d["ENSG2"] == ["1","0"]
    assert summary["read_type"]=="long" and summary["organism_type"]=="eukaryote"
```
- [ ] **Step 2: fails**
- [ ] **Step 3: implement** — eukaryote dalını read_type'a böl + `_counts_euk_long`
Import: `from rnaforge.minimap2 import count_primary_alignments`; `from rnaforge.tximport import run_tximport, parse_tx2gene`.
Replace eukaryote dispatch:
```python
    if config.organism_type == "eukaryote":
        read_type = resolve_read_type(run_dir)
        if read_type == "long":
            summary = _counts_euk_long(config, metadata_path, run_dir, quant_dir, stats_dir, logs_dir, state)
        else:
            summary = _counts_euk(config, metadata_path, run_dir, quant_dir, stats_dir, logs_dir, state)
    else:
        read_type = resolve_read_type(run_dir)
        if read_type == "long":
            summary = _counts_long(...)
        else:
            summary = _counts_short(...)
```
Add:
```python
def _counts_euk_long(config, metadata_path, run_dir, quant_dir, stats_dir, logs_dir, state) -> dict:
    """Ökaryot uzun-okuma sayımı: BAM primer-hizalama → transkript sayımı → tx2gene ile gen.
    Diagnostik (kapı yok). counts.tsv ortak sözleşme."""
    stats_path = stats_dir / "count_statistics.json"
    log_path = logs_dir / "counts.log"
    with log_path.open("w") as log_file:
        def log(m): log_file.write(m + "\n"); log_file.flush()
        samples = load_metadata(metadata_path)
        tx2gene = parse_tx2gene(config.reference.tx2gene)
        per_sample_gene: dict[str, dict[str, int]] = {}
        genes_seen: set[str] = set()
        for sample in samples:
            state.heartbeat()
            bam = quant_dir / sample.sample_id / "aligned.sorted.bam"
            if not bam.exists():
                raise ValueError(f"m05 eukaryote-long: BAM eksik: {bam} (m04 koştu mu?)")
            txc = count_primary_alignments(bam)
            gc: dict[str, int] = {}
            for tx, c in txc.items():
                g = tx2gene.get(tx) or tx2gene.get(tx.split(".")[0])
                if g is None: continue
                gc[g] = gc.get(g, 0) + c
            per_sample_gene[sample.sample_id] = gc
            genes_seen.update(gc)
            log(f"{sample.sample_id}: {sum(txc.values())} primer okuma → {len(gc)} gen")
        if not genes_seen:
            raise ValueError("m05 eukaryote-long: 0 gen (tx2gene eşleşmedi / hizalama boş?)")
        sample_ids = [s.sample_id for s in samples]
        genes = sorted(genes_seen)
        matrix_path = quant_dir / "counts.tsv"
        with matrix_path.open("w") as fh:
            fh.write("gene\t" + "\t".join(sample_ids) + "\n")
            for g in genes:
                fh.write(g + "\t" + "\t".join(str(per_sample_gene[s].get(g, 0)) for s in sample_ids) + "\n")
        log(f"count matrix written: {matrix_path} ({len(genes)} genes)")
        summary = {"read_type": "long", "organism_type": "eukaryote",
                   "n_samples": len(samples), "n_genes": len(genes), "gate_counts": {}}
        stats_path.write_text(json.dumps(summary, indent=2))
    return summary
```
- [ ] **Step 4: passes**
- [ ] **Step 5: full m05 suite**
- [ ] **Step 6: Commit + push**

---

### Task 4: full suite regression + push

- [ ] **Step 1:** `conda run -n rnaforge-core python -m pytest -q` → tümü yeşil (510 + yeni)
- [ ] **Step 2:** DURUM + memory güncelle, commit + push

## Self-Review
- Spec §3 m04-euk-long → Task 2 ✓ · §4 m05-euk-long (count+tx2gene aggregate) → Task 1+3 ✓ · §6 test → her task TDD ✓ · §7 validation → deferred ✓
- Types: `count_primary_alignments -> dict[str,int]` (Task1) consumed Task3 ✓; `parse_tx2gene -> dict[str,str]` ✓; `AlignmentResult(.alignment_rate,.bam)` reused ✓; router uses existing `_quant_long/_counts_long/_quant_euk/_counts_euk` verbatim ✓
- tx2gene versiyon eşleşmesi: transkriptomdan türetilmiş tx2gene versiyonlu; BAM ref adları = transkriptom FASTA adları (versiyonlu) → doğrudan eşleşir. Yedek: `tx.split(".")[0]` versiyonsuz deneme.
