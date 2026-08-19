from __future__ import annotations

import rnaforge.versions as versions
from rnaforge.versions import capture_tool_versions, parse_conda_list_json


def test_parse_conda_list_json():
    text = (
        '[{"name": "fastp", "version": "1.3.6", "channel": "bioconda"},'
        ' {"name": "python", "version": "3.11.15", "channel": "conda-forge"}]'
    )
    assert parse_conda_list_json(text) == {"fastp": "1.3.6", "python": "3.11.15"}


def test_parse_conda_list_json_invalid_returns_empty():
    assert parse_conda_list_json("not json") == {}
    assert parse_conda_list_json("") == {}


def test_capture_tool_versions_reads_installed(monkeypatch):
    """Kurulu sürümleri env başına tek conda sorgusuyla toplar; kurulu olmayan araç atlanır."""
    monkeypatch.setattr(versions.shutil, "which", lambda _x: "/usr/bin/conda")
    fake_envs = {
        "rnaforge-qc": {"fastp": "1.3.6", "fastqc": "0.12.1"},
        "rnaforge-de": {"bioconductor-deseq2": "1.50.2"},
        # rnaforge-core sorgusu boş → Python/networkx atlanır (kurulu değil gibi)
    }
    seen = []

    def fake_conda_list(env):
        seen.append(env)
        return fake_envs.get(env, {})

    monkeypatch.setattr(versions, "_conda_list", fake_conda_list)
    sources = {
        "fastp": ("rnaforge-qc", "fastp"),
        "FastQC": ("rnaforge-qc", "fastqc"),
        "DESeq2 (R)": ("rnaforge-de", "bioconductor-deseq2"),
        "networkx": ("rnaforge-core", "networkx"),
    }
    got = capture_tool_versions(sources)
    assert got == {"fastp": "1.3.6", "FastQC": "0.12.1", "DESeq2 (R)": "1.50.2"}
    # env başına TEK sorgu (rnaforge-qc iki araç için tekrar sorgulanmaz)
    assert seen.count("rnaforge-qc") == 1


def test_capture_tool_versions_no_conda_returns_empty(monkeypatch):
    """conda yoksa boş dict → çağıran curated fallback'e düşer (rapor çökmden)."""
    monkeypatch.setattr(versions.shutil, "which", lambda _x: None)
    assert capture_tool_versions({"fastp": ("rnaforge-qc", "fastp")}) == {}
