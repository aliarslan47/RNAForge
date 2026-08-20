import pytest

from rnaforge.quality import load_profile, profile_name_for


def test_metatranscriptome_profile_loads_permissive():
    p = load_profile("metatranscriptome")
    assert p.permissive is True
    assert p.threshold("replicate_correlation") == 0.75
    # alignment düşük eşik (katalog eksikliği doğal): FAIL yerine tolere
    assert p.threshold("alignment_rate") <= 0.10


def test_profile_name_for_metatranscriptome_short():
    assert profile_name_for("metatranscriptome", "short") == "metatranscriptome"


def test_metatranscriptome_profile_has_survival_rate_gate():
    """m03 kısa-okuma trim (_trim_short) survival_rate eşiğini profile.threshold'dan
    okur; yoksa ProfileError fırlatır. Topluluk okumaları rRNA-zengin → permissive
    (prokaryote_long felsefesiyle tutarlı, 0.20)."""
    p = load_profile("metatranscriptome")
    assert p.threshold("survival_rate") == pytest.approx(0.20)
    # sıkı prokaryote (0.50) değil, permissive-uzun-okuma (0.20) gibi toleranslı
    assert p.threshold("survival_rate") <= 0.20
