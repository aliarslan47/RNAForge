from __future__ import annotations

import json
import textwrap

from rnaforge.config import load_config


def _config(tmp_path):
    (tmp_path / "ref").mkdir()
    (tmp_path / "ref" / "genome.fa").write_text(">c1\nACGT\n")
    (tmp_path / "ref" / "genes.gff").write_text("##gff-version 3\n")
    cfg = tmp_path / "config.yaml"
    cfg.write_text(textwrap.dedent(f"""
        organism: "E. coli"
        organism_type: "prokaryote"
        reference:
          genome_fasta: "{tmp_path / 'ref' / 'genome.fa'}"
          annotation_gff: "{tmp_path / 'ref' / 'genes.gff'}"
    """))
    return load_config(cfg)


def _seed_read_type(run_dir, read_type):
    (run_dir / "statistics").mkdir(parents=True)
    (run_dir / "statistics" / "raw_statistics.json").write_text(
        json.dumps({"platform": "ont", "read_type": read_type})
    )


def test_load_run_profile_long_stamps_long_profile(tmp_path):
    from rnaforge.cli import _load_run_profile
    run_dir = tmp_path / "run"
    _seed_read_type(run_dir, "long")
    profile = _load_run_profile(_config(tmp_path), run_dir)
    assert profile.name == "prokaryote_long"
    assert profile.permissive is True


def test_load_run_profile_short_uses_organism_profile(tmp_path):
    from rnaforge.cli import _load_run_profile
    run_dir = tmp_path / "run"
    _seed_read_type(run_dir, "short")
    profile = _load_run_profile(_config(tmp_path), run_dir)
    assert profile.name == "prokaryote"


def test_effective_metadata_prefers_resolved(tmp_path):
    """Basecall (m00) çözülmüş metadata yazdıysa tüm aşamalar onu kullanır."""
    from rnaforge.cli import _effective_metadata
    from rnaforge.basecall import basecalled_metadata_path
    run_dir = tmp_path / "run"
    orig = tmp_path / "orig.tsv"; orig.write_text("x")
    assert _effective_metadata(orig, run_dir) == orig      # çözülmüş yok → orijinal
    rm = basecalled_metadata_path(run_dir); rm.parent.mkdir(parents=True); rm.write_text("y")
    assert _effective_metadata(orig, run_dir) == rm        # çözülmüş var → onu kullan


def test_load_run_profile_pre_m01_falls_back_to_short(tmp_path):
    """raw_statistics yoksa (m01 öncesi) kısa profile düşer, çökmez."""
    from rnaforge.cli import _load_run_profile
    run_dir = tmp_path / "run"; run_dir.mkdir()
    profile = _load_run_profile(_config(tmp_path), run_dir)
    assert profile.name == "prokaryote"
