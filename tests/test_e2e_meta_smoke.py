"""Task 12 — Metatranskriptom uçtan-uca smoke testi (gerçek araç zinciri, conda varsa koşar).

Metatranskriptom KOLUNU bütün olarak kanıtlar: `rnaforge run` organism_type=metatranscriptome'da
`trim`'den sonra rrna-deplete + taxonomy aşamalarını otomatik ekler, sonra rRNA'sı çıkarılmış
okumaları gen kataloğuna hizalar (Bowtie2) → sayar (featureCounts) → DESeq2 → rapor. Rapor
"Topluluk Kompozisyonu (Taksonomi)" bölümünü ve permissive damgayı taşımalı; hiçbir aşama FAIL
vermemeli (permissive metatranscriptome profili).

İki DB-bağımlı ARAÇ çağrısı monkeypatch'lenir (plan Task 12 §1 izin verir): SortMeRNA rRNA
depletion (gerçek rRNA DB gerektirir) ve Kraken2/Bracken (GB'larca DB gerektirir). Bunların
DIŞINDA her şey GERÇEK: fastp (m03) → Bowtie2 (m04 gen kataloğu) → featureCounts (m05) → DESeq2
(m06) → rapor (m08). Böylece orkestrasyon sırası, aşama bağımlılıkları, meta yönlendirmesi
(_quant_meta/_counts_meta) ve gerçek DE→rapor zinciri kanıtlanır; iki dış DB aracının kendisi
kendi birim testlerinde (test_m_rrna_deplete, test_m_taxonomy, test_kraken2) doğrulanır.

conda ya da gerekli env'ler (qc/quant-prok/de) yoksa skip."""
from __future__ import annotations

import json
import random
import shutil
import subprocess
import textwrap
from pathlib import Path

import pytest

from rnaforge.cli import main

# SortMeRNA/Kraken2 monkeypatch'li → seqqc/meta env'i GEREKMEZ; yalnız gerçek koşan araçların
# env'leri gerekir: fastp (qc), bowtie2+featureCounts (quant-prok), DESeq2 (de).
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


def _build_meta_fixture(tmp_path: Path):
    """Sentetik gen kataloğu (tek contig, 12 gen; ilk 3'ü treated'da ~4× yukarı) + katalog
    anotasyonu (GTF exon/gene_id) + dummy kraken2 DB dizini + dummy rRNA FASTA. Okumalar
    genlerden üretilir (gerçek Bowtie2 hizalar)."""
    rng = random.Random(2024)
    genes, gtf_lines, pos, catalog = [], [], 1, []
    for i in range(_N_GENES):
        seq = "".join(rng.choice("ACGT") for _ in range(_GENE_LEN))
        start, end = pos, pos + _GENE_LEN - 1
        genes.append((f"g{i}", seq))
        gtf_lines.append(f'contig1\tsynthetic\texon\t{start}\t{end}\t.\t+\t.\tgene_id "g{i}";')
        catalog.append(seq)
        catalog.append("".join(rng.choice("ACGT") for _ in range(_SPACER)))
        pos = end + _SPACER + 1

    ref = tmp_path / "ref"; ref.mkdir()
    (ref / "catalog.fa").write_text(">contig1\n" + "".join(catalog) + "\n")
    (ref / "catalog.gtf").write_text("\n".join(gtf_lines) + "\n")
    # Monkeypatch'lenen araçların DB'leri: var-olmaları ŞART değil (m01 yalnız gene_catalog_fasta
    # + catalog_annotation varlığını denetler) ama rapor yazılım tablosu için dummy koyalım.
    kdb = tmp_path / "kraken2_db"; kdb.mkdir()
    rrna = tmp_path / "rrna.fasta"; rrna.write_text(">rRNA_dummy\nACGT\n")

    samples = [("c1", "control"), ("c2", "control"), ("t1", "treated"), ("t2", "treated")]
    for sid, cond in samples:
        reads = []
        for gi, (_g, seq) in enumerate(genes):
            depth = _BASE_DEPTH + rng.randint(-5, 5)
            if gi < _SIGNAL_GENES and cond == "treated":
                depth *= 4
            for r in range(depth):
                off = rng.randint(0, _GENE_LEN - _READ_LEN)
                reads.append(f"@{sid}_g{gi}_{r}\n{seq[off:off + _READ_LEN]}\n+\n{'I' * _READ_LEN}")
        (tmp_path / f"{sid}.fastq").write_text("\n".join(reads) + "\n")

    config_path = tmp_path / "config.yaml"
    config_path.write_text(textwrap.dedent(f"""
        organism: "Synthetic community"
        organism_type: "metatranscriptome"
        reference:
          gene_catalog_fasta: "{ref / 'catalog.fa'}"
          catalog_annotation: "{ref / 'catalog.gtf'}"
        taxonomy:
          kraken2_db: "{kdb}"
        rrna:
          db_fasta: "{rrna}"
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


def _install_db_tool_stubs(monkeypatch):
    """İki DB-bağımlı ARAÇ sarmalayıcısını değiştir (plan Task 12 §1). Modüller bunları kendi
    ad alanlarına import ettiği için orada patch'lenir."""
    from rnaforge.modules import m_rrna_deplete, m_taxonomy

    def fake_sortmerna(reads, db_fasta, workdir, paired=False, threads=1, env=None):
        # rRNA depletion pass-through: "rRNA'sız" (other) okumalar = trimlenmiş girdi.
        # Gerçek Bowtie2 bunları gen kataloğuna hizalar. depletion_rate gate için (WARN-only).
        return {"other": [reads[0]], "depletion_rate": 0.42}

    def fake_kraken2(reads, db, prefix, paired=False, threads=1, env=None):
        return Path(f"{prefix}.kreport")            # run_bracken monkeypatch'li → içeriği önemsiz

    def fake_bracken(report, db, out_path, read_len=100, level="S", env=None):
        # Gerçek parse_bracken'in okuyacağı 7-sütun Bracken formatı (header + satırlar).
        Path(out_path).write_text(
            "name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\t"
            "new_est_reads\tfraction_total_reads\n"
            "Escherichia coli\t562\tS\t120\t30\t150\t0.60\n"
            "Bacteroides fragilis\t817\tS\t60\t15\t75\t0.30\n"
            "Faecalibacterium prausnitzii\t853\tS\t20\t5\t25\t0.10\n"
        )
        return Path(out_path)

    monkeypatch.setattr(m_rrna_deplete, "run_sortmerna_deplete", fake_sortmerna)
    monkeypatch.setattr(m_taxonomy, "run_kraken2", fake_kraken2)
    monkeypatch.setattr(m_taxonomy, "run_bracken", fake_bracken)


@pytest.mark.skipif(not _have_required_envs(),
                    reason=f"conda veya gerekli env'ler yok: {', '.join(_REQUIRED_ENVS)}")
def test_metatranscriptome_run_end_to_end_to_report(tmp_path, monkeypatch):
    config_path, metadata_path = _build_meta_fixture(tmp_path)
    _install_db_tool_stubs(monkeypatch)
    runs_dir = tmp_path / "runs"

    rc = main(["run", "--config", str(config_path), "--metadata", str(metadata_path),
               "--runs-dir", str(runs_dir), "--run-id", "meta", "--to", "report"])
    assert rc == 0, "metatranskriptom pipeline (validate→report) sıfır dönmeli (FAIL kapısı yok)"

    run_dir = sorted(runs_dir.glob("*_meta"))[-1]

    # Orkestrasyon meta aşamalarını gerçekten koştu (state'te işaretli).
    state = json.loads((run_dir / "state.json").read_text())
    done = state.get("modules", state)
    assert "m_rrna_deplete" in json.dumps(done)
    assert "m_taxonomy" in json.dumps(done)

    # Gerçek Bowtie2 + featureCounts → 12 genli sayım matrisi.
    counts = (run_dir / "quantification" / "counts.tsv").read_text().splitlines()
    assert len(counts) - 1 == _N_GENES, "count matrisi 12 gen içermeli"

    # Taksonomi çıktısı (stub Bracken) birleştirildi.
    matrix = (run_dir / "taxonomy" / "abundance_matrix.tsv").read_text()
    assert "Escherichia coli" in matrix

    # Gerçek DESeq2 → treated'da yukarı 3 sinyal geni.
    de = json.loads((run_dir / "statistics" / "de_statistics.json").read_text())
    assert de["n_up"] >= _SIGNAL_GENES, "treated'da yukarı 3 sinyal geni anlamlı olmalı"

    # Rapor: taksonomi bölümü + permissive damga (metatranscriptome profili).
    report = (run_dir / "report" / "report.html").read_text()
    assert 'id="taxonomy"' in report
    assert "Escherichia coli" in report
    assert "metatranscriptome" in report            # permissive-profil damgası (yöntem/güven)
