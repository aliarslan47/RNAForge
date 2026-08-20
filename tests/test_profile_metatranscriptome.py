from rnaforge.quality import load_profile, profile_name_for


def test_metatranscriptome_profile_loads_permissive():
    p = load_profile("metatranscriptome")
    assert p.permissive is True
    assert p.threshold("replicate_correlation") == 0.75
    # alignment düşük eşik (katalog eksikliği doğal): FAIL yerine tolere
    assert p.threshold("alignment_rate") <= 0.10


def test_profile_name_for_metatranscriptome_short():
    assert profile_name_for("metatranscriptome", "short") == "metatranscriptome"
