from pathlib import Path

from rnaforge import tximport as tx


def test_run_tximport_reads_matrix(tmp_path, monkeypatch):
    out = tmp_path / "out"; out.mkdir()

    def fake_run(cmd, **k):
        (out / "gene_counts.tsv").write_text("gene\ts1\ts2\ng1\t10\t20\ng2\t0\t5\n")
        class R: returncode = 0; stdout = ""; stderr = ""
        return R()

    monkeypatch.setattr(tx.subprocess, "run", fake_run)
    res = tx.run_tximport({"s1": tmp_path / "s1.sf", "s2": tmp_path / "s2.sf"},
                          tmp_path / "t2g.tsv", out)
    assert res.gene_ids == ["g1", "g2"]
    # gene_counts.tsv: g1→(s1=10,s2=20), g2→(s1=0,s2=5) ⇒ sütun s1=[10,0], s2=[20,5]
    assert res.counts["s1"] == [10.0, 0.0]
    assert res.counts["s2"] == [20.0, 5.0]


def test_run_tximport_raises_on_missing_output(tmp_path, monkeypatch):
    out = tmp_path / "out"; out.mkdir()

    def fake_run(cmd, **k):
        class R: returncode = 1; stdout = ""; stderr = "boom"
        return R()

    monkeypatch.setattr(tx.subprocess, "run", fake_run)
    import pytest
    with pytest.raises(tx.TximportError):
        tx.run_tximport({"s1": tmp_path / "s1.sf"}, tmp_path / "t2g.tsv", out)


def test_parse_tx2gene(tmp_path):
    from rnaforge.tximport import parse_tx2gene
    p = tmp_path / "t2g.tsv"
    p.write_text("ENST1.1\tENSG1.1\nENST2.2\tENSG1.1\nENST3.1\tENSG2.4\n")
    assert parse_tx2gene(p) == {"ENST1.1": "ENSG1.1", "ENST2.2": "ENSG1.1", "ENST3.1": "ENSG2.4"}
