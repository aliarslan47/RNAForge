from __future__ import annotations

import pytest

from rnaforge.quality import ProfileError, load_profile


def test_prokaryote_profile_is_strict_about_alignment():
    profile = load_profile("prokaryote")
    assert profile.name == "prokaryote"
    assert profile.threshold("alignment_rate") == 0.70
    assert profile.permissive is False


def test_eukaryote_profile_is_marked_permissive():
    """Elde okaryot dogrulamasi YOK; gevsek esikler rapora damgalanmali."""
    profile = load_profile("eukaryote")
    assert profile.permissive is True
    assert profile.threshold("alignment_rate") == 0.50


def test_override_changes_threshold_and_is_recorded():
    profile = load_profile("prokaryote", overrides={"alignment_rate": 0.30})
    assert profile.threshold("alignment_rate") == 0.30
    assert profile.is_overridden("alignment_rate") is True
    assert profile.is_overridden("survival_rate") is False


def test_unknown_override_key_is_rejected():
    """Yazim hatasi sessizce yutulmamali; kullanici esigi ezdigini saniyor olabilir."""
    with pytest.raises(ProfileError, match="alignment_rat"):
        load_profile("prokaryote", overrides={"alignment_rat": 0.3})


def test_unknown_organism_type_is_rejected():
    with pytest.raises(ProfileError, match="no quality profile"):
        load_profile("archaea")
