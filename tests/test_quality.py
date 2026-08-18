from __future__ import annotations

import types

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


def test_prokaryote_long_profile_is_permissive():
    """ONT uzun-okuma profili bilinçli permissive + damgalı; yalnız katastrofik
    hizalama (yanlış referans) FAIL, düşük survival/assignment WARN."""
    profile = load_profile("prokaryote_long")
    assert profile.name == "prokaryote_long"
    assert profile.permissive is True
    assert profile.threshold("alignment_rate") == 0.50    # katastrofik floor → FAIL
    assert profile.threshold("survival_rate") == 0.20     # Pychopper doğal düşük → WARN
    assert profile.threshold("assignment_rate") == 0.05   # ONT CDS-only düşük → WARN


def test_eukaryote_long_profile_is_permissive():
    """Ökaryot uzun-okuma profili de permissive + damgalı (ökaryot + ONT gerekçeleri birleşir)."""
    profile = load_profile("eukaryote_long")
    assert profile.name == "eukaryote_long"
    assert profile.permissive is True
    assert profile.threshold("alignment_rate") == 0.50
    assert profile.threshold("survival_rate") == 0.20


def test_profile_name_for_long_and_short():
    from rnaforge.quality import profile_name_for
    assert profile_name_for("prokaryote", "long") == "prokaryote_long"
    assert profile_name_for("prokaryote", "short") == "prokaryote"
    assert profile_name_for("eukaryote", "short") == "eukaryote"
    assert profile_name_for("eukaryote", "long") == "eukaryote_long"


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


def test_non_numeric_override_is_rejected():
    """Sayisal olmayan bir ezme degeri "expected a number" ile reddedilmeli."""
    with pytest.raises(ProfileError, match="expected a number"):
        load_profile("prokaryote", overrides={"alignment_rate": "yuksek"})


def test_overrides_returns_only_overridden_gates_sorted():
    profile = load_profile(
        "prokaryote",
        overrides={"rrna_fraction": 0.15, "alignment_rate": 0.30},
    )
    assert profile.overrides() == {"alignment_rate": 0.30, "rrna_fraction": 0.15}
    assert isinstance(profile.overrides(), dict)


def test_profile_description_is_loaded_and_stripped():
    profile = load_profile("prokaryote")
    assert profile.description.startswith("Bakteriyel RNA-seq")
    assert profile.description == profile.description.strip()


def test_threshold_return_type_is_always_float():
    """read_depth ve base_quality YAML'da int; donen deger yine de float olmali."""
    profile = load_profile("prokaryote")
    assert isinstance(profile.threshold("read_depth"), float)
    assert isinstance(profile.threshold("base_quality"), float)


def test_profile_internals_cannot_be_mutated_after_construction():
    """frozen=True sadece attribute atamasini engeller; ic yapilar da
    gercekten degistirilemez olmali - yoksa is_overridden() sessizce
    yanlis sonuc dondurebilir ve musteri raporundaki guven karti bozulur.
    """
    profile = load_profile("prokaryote")

    assert isinstance(profile._thresholds, types.MappingProxyType)
    with pytest.raises(TypeError):
        profile._thresholds["alignment_rate"] = 0.0

    assert isinstance(profile._overridden, frozenset)
    with pytest.raises(AttributeError):
        profile._overridden.add("read_depth")
