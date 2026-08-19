"""Gerçek-araç uçtan-uca smoke testi (CI-default, conda varsa koşar).

512 birim/mock testinin YAKALAYAMADIĞI şeyi kanıtlar: pipeline BÜTÜN olarak, gerçek
fastp→bowtie2→featureCounts→DESeq2 zinciriyle, `rnaforge run` orkestratörü altında,
küçük sentetik ama sinyalli bir veri üzerinde çalışıyor mu. Araç sürümü çıktı formatını
değiştirirse bu test kırılır (parser'lar sabit fixture yerine gerçek çıktı görür).

Sentetik prokaryot: 12 gen, ilk 3'ü treated'da ~4× yukarı → DESeq2 bu 3'ü anlamlı
DEG bulmalı. conda ya da gerekli env'ler yoksa skip (Ali'nin makinesinde koşar)."""
from __future__ import annotations

import random
import shutil
import subprocess
import textwrap
from pathlib import Path

import pytest

from rnaforge.cli import main

_REQUIRED_ENVS = ("rnaforge-qc", "rnaforge-quant-prok", "rnaforge-de")
_GENE_LEN = 600
_SPACER = 100
_READ_LEN = 100
_N_GENES = 12
_SIGNAL_GENES = 3          # ilk 3 gen treated'da yukarı
_BASE_DEPTH = 80


def _have_required_envs() -> bool:
    if shutil.which("conda") is None:
        return False
    out = subprocess.run(["conda", "env", "list"], capture_output=True, text=True)
    names = {line.split()[0] for line in out.stdout.splitlines()
             if line.strip() and not line.startswith("#")}
    return all(env in names for env in _REQUIRED_ENVS)


def _build_fixture(tmp_path: Path):
    rng = random.Random(2024)
    genes, gtf_lines, pos = [], [], 1
    genome = []
    for i in range(_N_GENES):
        seq = "".join(rng.choice("ACGT") for _ in range(_GENE_LEN))
        start = pos
        end = pos + _GENE_LEN - 1
        genes.append((f"g{i}", seq, start))
        gtf_lines.append(f'contig1\tsynthetic\texon\t{start}\t{end}\t.\t+\t.\tgene_id "g{i}";')
        genome.append(seq)
        genome.append("".join(rng.choice("ACGT") for _ in range(_SPACER)))
        pos = end + _SPACER + 1

    ref = tmp_path / "ref"
    ref.mkdir()
    (ref / "genome.fa").write_text(">contig1\n" + "".join(genome) + "\n")
    (ref / "genes.gtf").write_text("\n".join(gtf_lines) + "\n")

    samples = [("c1", "control"), ("c2", "control"), ("t1", "treated"), ("t2", "treated")]
    for sid, cond in samples:
        reads = []
        for gi, (_gname, seq, _start) in enumerate(genes):
            depth = _BASE_DEPTH + rng.randint(-5, 5)
            if gi < _SIGNAL_GENES and cond == "treated":
                depth *= 4
            for r in range(depth):
                off = rng.randint(0, _GENE_LEN - _READ_LEN)
                frag = seq[off:off + _READ_LEN]
                reads.append(f"@{sid}_g{gi}_{r}\n{frag}\n+\n{'I' * _READ_LEN}")
        (tmp_path / f"{sid}.fastq").write_text("\n".join(reads) + "\n")

    config_path = tmp_path / "config.yaml"
    config_path.write_text(textwrap.dedent(f"""
        organism: "Synthetic prokaryote"
        organism_type: "prokaryote"
        reference:
          genome_fasta: "{ref / 'genome.fa'}"
          annotation_gff: "{ref / 'genes.gtf'}"
        de:
          design: "~condition"
          reference: control
    """))
    metadata_path = tmp_path / "samples.tsv"
    metadata_path.write_text(
        "sample_id\tcondition\tfastq_1\n"
        + "".join(f"{sid}\t{cond}\t{sid}.fastq\n" for sid, cond in samples)
    )
    return config_path, metadata_path


@pytest.mark.skipif(not _have_required_envs(),
                    reason=f"conda veya gerekli env'ler yok: {', '.join(_REQUIRED_ENVS)}")
def test_run_pipeline_end_to_end_real_tools(tmp_path):
    import json
    config_path, metadata_path = _build_fixture(tmp_path)
    runs_dir = tmp_path / "runs"
    rc = main(["run", "--config", str(config_path), "--metadata", str(metadata_path),
               "--runs-dir", str(runs_dir), "--run-id", "smoke", "--to", "de"])
    assert rc == 0, "gerçek-araç pipeline (validate→de) sıfır dönmeli"

    run_dir = sorted(runs_dir.glob("*_smoke"))[-1]
    counts = (run_dir / "quantification" / "counts.tsv").read_text().splitlines()
    assert len(counts) - 1 == _N_GENES, "count matrisi 12 gen içermeli"

    de = json.loads((run_dir / "statistics" / "de_statistics.json").read_text())
    assert de["n_genes"] == _N_GENES
    assert de["n_up"] >= _SIGNAL_GENES, "treated'da yukarı 3 sinyal geni anlamlı olmalı"

    results = {r.split("\t")[0]: r for r in
               (run_dir / "differential_expression" / "deseq2_results.tsv").read_text().splitlines()[1:]}
    for gi in range(_SIGNAL_GENES):
        fields = results[f"g{gi}"].split("\t")     # gene baseMean log2FC lfcSE stat pvalue padj
        assert float(fields[2]) > 1, f"g{gi} treated'da yukarı olmalı (log2FC>1)"


def _build_multifactor_fixture(tmp_path: Path):
    """Faz 3 canlı doğrulama: keyfi kovaryat (sex) + >2-faktör design + condition-dışı
    kontrast. 6 örnek (3 control/3 treated), sex dengeli (confounded değil)."""
    rng = random.Random(2025)
    genes, gtf_lines, pos, genome = [], [], 1, []
    for i in range(_N_GENES):
        seq = "".join(rng.choice("ACGT") for _ in range(_GENE_LEN))
        start, end = pos, pos + _GENE_LEN - 1
        genes.append((f"g{i}", seq, start))
        gtf_lines.append(f'contig1\tsynthetic\texon\t{start}\t{end}\t.\t+\t.\tgene_id "g{i}";')
        genome.append(seq)
        genome.append("".join(rng.choice("ACGT") for _ in range(_SPACER)))
        pos = end + _SPACER + 1
    ref = tmp_path / "ref"; ref.mkdir()
    (ref / "genome.fa").write_text(">contig1\n" + "".join(genome) + "\n")
    (ref / "genes.gtf").write_text("\n".join(gtf_lines) + "\n")

    # sex condition ile confounded DEĞİL (her condition'da hem M hem F var).
    samples = [("c1", "control", "M"), ("c2", "control", "F"), ("c3", "control", "M"),
               ("t1", "treated", "F"), ("t2", "treated", "M"), ("t3", "treated", "F")]
    for sid, cond, _sex in samples:
        reads = []
        for gi, (_g, seq, _s) in enumerate(genes):
            depth = _BASE_DEPTH + rng.randint(-5, 5)
            if gi < _SIGNAL_GENES and cond == "treated":
                depth *= 4
            for r in range(depth):
                off = rng.randint(0, _GENE_LEN - _READ_LEN)
                reads.append(f"@{sid}_g{gi}_{r}\n{seq[off:off + _READ_LEN]}\n+\n{'I' * _READ_LEN}")
        (tmp_path / f"{sid}.fastq").write_text("\n".join(reads) + "\n")

    config_path = tmp_path / "config.yaml"
    config_path.write_text(textwrap.dedent(f"""
        organism: "Synthetic prokaryote"
        organism_type: "prokaryote"
        reference:
          genome_fasta: "{ref / 'genome.fa'}"
          annotation_gff: "{ref / 'genes.gtf'}"
        de:
          design: "~sex + condition"
          reference: control
          contrasts:
            - [treated, control]
            - [sex, M, F]
    """))
    metadata_path = tmp_path / "samples.tsv"
    metadata_path.write_text(
        "sample_id\tcondition\tsex\tfastq_1\n"
        + "".join(f"{sid}\t{cond}\t{sex}\t{sid}.fastq\n" for sid, cond, sex in samples)
    )
    return config_path, metadata_path


@pytest.mark.skipif(not _have_required_envs(),
                    reason=f"conda veya gerekli env'ler yok: {', '.join(_REQUIRED_ENVS)}")
def test_run_pipeline_multifactor_covariate_to_report(tmp_path):
    """Faz 3 uçtan-uca (rapora kadar): keyfi kovaryat design + condition-dışı kontrast
    gerçek DESeq2'de koşar, condition sinyalini bulur, ayrı sex kontrast dosyası üretir,
    ve rapor çökmeden oluşur (A6 doğrulaması)."""
    import json
    config_path, metadata_path = _build_multifactor_fixture(tmp_path)
    runs_dir = tmp_path / "runs"
    rc = main(["run", "--config", str(config_path), "--metadata", str(metadata_path),
               "--runs-dir", str(runs_dir), "--run-id", "mf", "--to", "report"])
    assert rc == 0, "çok-faktör pipeline (validate→report) sıfır dönmeli"

    run_dir = sorted(runs_dir.glob("*_mf"))[-1]
    # coldata çok-sütunlu (sex yazıldı)
    coldata = (run_dir / "differential_expression" / "coldata.tsv").read_text().splitlines()
    assert coldata[0].split("\t") == ["sample", "condition", "sex"]
    # condition sinyali bulundu (birincil kontrast)
    de = json.loads((run_dir / "statistics" / "de_statistics.json").read_text())
    assert de["n_up"] >= _SIGNAL_GENES
    # condition-dışı faktör kontrastı ayrı dosya üretti
    assert (run_dir / "differential_expression" / "deseq2_results.sex.M_vs_F.tsv").exists()
    # rapor çökmeden oluştu
    reports = list(run_dir.glob("**/report*.html"))
    assert reports, "rapor HTML üretilmeli (çok-faktör + condition-dışı kontrast ile çökmeden)"
