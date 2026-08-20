"""Task 11 — prepare_references.sh metatranskriptom bayrakları (hafif; AĞIR İNDİRME YOK).
Yalnız yardım metni + bayrak tanıma + söz dizimi doğrulanır; gerçek Kraken2/rRNA DB
indirmeleri CI-dışıdır (ağ + GB'lar). Betik indirme bloklarına GİRMEDEN test edilir."""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parent.parent / "prepare_references.sh"


def _run(*args):
    return subprocess.run(["bash", str(SCRIPT), *args], capture_output=True, text=True)


def test_script_exists_and_syntax_valid():
    assert SCRIPT.exists()
    # bash -n: söz dizimi denetimi, hiçbir komut ÇALIŞTIRMAZ (indirme tetiklenmez).
    r = subprocess.run(["bash", "-n", str(SCRIPT)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr


def test_help_lists_metatranscriptome_flags():
    r = _run("--help")
    assert r.returncode == 0
    for flag in ("--kraken2-db-url", "--kraken2-db-name", "--rrna-db-url"):
        assert flag in r.stdout, f"{flag} yardım metninde yok"
    # var olan bayraklar da korunmuş olmalı (regresyon)
    assert "--kegg-org" in r.stdout and "--string-taxid" in r.stdout


def test_unknown_flag_rejected():
    r = _run("--definitely-not-a-flag")
    assert r.returncode == 2
    assert "bilinmeyen argüman" in r.stderr


@pytest.mark.skipif(shutil.which("bash") is None, reason="bash yok")
def test_meta_flags_are_recognized_not_unknown(tmp_path):
    # --kraken2-db-name tek başına (URL yok) → meta bloğu ATLANIR (indirme yok),
    # betik hatasız biter. Bayrağın "bilinmeyen argüman" ile reddedilmediğini kanıtlar.
    r = subprocess.run(["bash", str(SCRIPT), "--skip-obo", "--kraken2-db-name", "testdb"],
                       capture_output=True, text=True, cwd=tmp_path)
    assert r.returncode == 0, r.stderr
    assert "bilinmeyen argüman" not in r.stderr
    assert "Kraken2 DB atlandı" in r.stdout      # URL verilmedi → blok atlandı
    assert "rRNA DB atlandı" in r.stdout
