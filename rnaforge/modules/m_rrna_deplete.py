"""m_rrna_deplete — Metatranskriptom: SortMeRNA ile rRNA'yı ÇIKAR (--other).

`rrna_deplete.run_sortmerna_deplete` yalnız SortMeRNA'yı çalıştırır ve sonucu döndürür;
bu modül aile deseninin geri kalanını sağlar: per-sample döngü, atomic state, heartbeat,
örnek-başı resume, ve rRNA'sız FASTQ'ları sözleşme yoluna (`rrna_depleted/<sid>/`) yazma.

Kapı: `rrna_depletion_rate` — WARN-ONLY, ASLA FAIL (permissive metatranscriptome profili;
gen kataloğu/rRNA DB'si eksik olabilir, düşük depletion verimi ŞÜPHELİ ama geçersiz değil).
"""
from __future__ import annotations

import gzip
import json
import shutil
from collections import Counter
from pathlib import Path

from rnaforge.config import Config
from rnaforge.gates import FAIL, PASS, WARN, GateResult, write_gate_results
from rnaforge.metadata import Sample, load_metadata
from rnaforge.quality import Profile, load_profile
from rnaforge.rrna_deplete import run_sortmerna_deplete
from rnaforge.state import RunState

MODULE_NAME = "m_rrna_deplete"
_GATE = "rrna_depletion_rate"


def rrna_depleted_reads(run_dir: Path, sample: Sample) -> list[Path]:
    """m_rrna_deplete'nin bir örnek için ürettiği rRNA'sız FASTQ yol(lar)ı. Adlandırma
    kuralının TEK kaynağı: m_rrna_deplete buraya yazar, downstream (m_taxonomy, m04) buradan
    okur (m03'ün trimmed_reads() deseni — drift önlenir)."""
    d = Path(run_dir) / "rrna_depleted" / sample.sample_id
    return sorted(d.glob("other_*.fastq.gz"))


def build_rrna_gates(depletion_rates: dict[str, float], profile: Profile,
                     warn_only: bool = True) -> list[GateResult]:
    """rrna_depletion_rate kapısı. warn_only=True (varsayılan ve TEK kullanılan biçim) —
    metatranscriptome profili bilinçli olarak permissive: düşük depletion verimi şüpheli
    (rRNA DB'si eksik/uyumsuz olabilir) ama sonucu GEÇERSİZ kılmaz, asla FAIL üretmez."""
    threshold = profile.threshold(_GATE)
    offenders = sorted(sid for sid, r in depletion_rates.items() if r < threshold)
    lowest = min(depletion_rates.values(), default=1.0)
    overridden = _GATE in profile.overrides()
    if offenders:
        status = WARN if warn_only else FAIL
        message = (
            f"rRNA depletion verimi eşiğin altında ({len(offenders)} örnek: "
            f"{', '.join(offenders)}); en düşük {lowest:.2f} < {threshold:.2f}. "
            "Kütüphane beklenenden fazla rRNA içeriyor olabilir."
        )
    else:
        status = PASS
        message = f"tüm örnekler rRNA depletion verimi >= {threshold:.2f} (en düşük {lowest:.2f})."
    return [GateResult(
        name=_GATE, module=MODULE_NAME, status=status, message=message,
        remedy=("Düşük depletion veriminde rrna.db_fasta (SortMeRNA rRNA referansı) ve "
                "kütüphane hazırlama protokolünü (rRNA-depletion kiti) gözden geçirin."),
        measured=lowest, threshold=threshold, overridden=overridden,
        samples=tuple(offenders),
    )]


def _gzip_into(src: Path, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if str(src).endswith(".gz"):
        shutil.copyfile(src, dest)
    else:
        with open(src, "rb") as f_in, gzip.open(dest, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    return dest


_FASTQ_EXTS = (".fastq.gz", ".fq.gz", ".fastq", ".fq")


def _dest_name(src: Path, idx: int) -> str:
    """Sözleşme adı: `other_<tag>.fastq.gz`. SortMeRNA'nın kendi adlandırması (`other.fastq`,
    `other_fwd.fastq`/`other_rev.fastq`...) korunur; yoksa idx ile tekil kılınır."""
    name = src.name
    for ext in _FASTQ_EXTS:
        if name.endswith(ext):
            base = name[: -len(ext)]
            break
    else:
        base = Path(name).stem
    tag = base[len("other"):].lstrip("_") if base.startswith("other") else base
    return f"other_{tag or idx}.fastq.gz"


def run_rrna_deplete(config: Config, metadata_path: Path, run_dir: Path,
                     force: bool = False) -> dict:
    run_dir = Path(run_dir)
    out_dir = run_dir / "rrna_depleted"
    stats_dir = run_dir / "statistics"
    logs_dir = run_dir / "logs"
    for d in (out_dir, stats_dir, logs_dir):
        d.mkdir(parents=True, exist_ok=True)
    state = RunState(run_dir)
    stats_path = stats_dir / "rrna_depletion.json"
    profile = load_profile(config.organism_type, config.quality)

    if not force and state.is_done(MODULE_NAME) and stats_path.exists():
        per_sample = json.loads(stats_path.read_text())
        gates = build_rrna_gates({sid: v["depletion_rate"] for sid, v in per_sample.items()},
                                 profile)
        return {
            "n_samples": len(per_sample), "samples": per_sample,
            "gate_counts": dict(Counter(g.status for g in gates)),
            "resumed": True,
        }

    if not state.is_done("m01_validate"):
        raise ValueError(
            "m_rrna_deplete requires m01 (validate) to have completed in this run "
            f"directory first: {run_dir}. Run `rnaforge validate` with the same "
            "--run-id, then re-run rrna-deplete."
        )

    log_path = logs_dir / "rrna_deplete.log"
    with log_path.open("w") as log_file:
        def log(msg: str) -> None:
            log_file.write(msg + "\n")
            log_file.flush()

        samples = load_metadata(metadata_path)
        log(f"m_rrna_deplete: {len(samples)} sample(s), rrna_db={config.rrna.db_fasta}, "
            f"env={config.rrna.env}")
        per_sample: dict[str, dict] = {}
        for sample in samples:
            state.heartbeat()
            sid = sample.sample_id
            if (not force and state.is_item_done(MODULE_NAME, sid)
                    and rrna_depleted_reads(run_dir, sample)):
                per_sample[sid] = state.item_payload(MODULE_NAME, sid)
                log(f"{sid}: resumed (cached depletion_rate="
                    f"{per_sample[sid].get('depletion_rate')})")
                continue

            paired = sample.fastq_2 is not None
            reads = [sample.fastq_1, sample.fastq_2] if paired else [sample.fastq_1]
            workdir = run_dir / "rrna_work" / sid
            result = run_sortmerna_deplete(
                reads, config.rrna.db_fasta, workdir, paired=paired,
                threads=config.resources.threads, env=config.rrna.env,
            )
            sample_out = out_dir / sid
            out_paths = [
                _gzip_into(src, sample_out / _dest_name(Path(src), idx))
                for idx, src in enumerate(sorted(result["other"]), start=1)
            ]
            shutil.rmtree(workdir, ignore_errors=True)

            depletion_rate = round(float(result["depletion_rate"]), 4)
            payload = {
                "depletion_rate": depletion_rate,
                "other_fastq": [str(p) for p in out_paths],
            }
            per_sample[sid] = payload
            state.mark_item_done(MODULE_NAME, sid, payload)
            log(f"{sid}: depletion_rate={depletion_rate:.3f} ({len(out_paths)} rRNA-free file(s))")

        gates = build_rrna_gates(
            {sid: ps["depletion_rate"] for sid, ps in per_sample.items()}, profile)
        # Sözleşme (Task 4 brief): dosya İÇERİĞİ tam olarak {sid: {"depletion_rate": float}} —
        # m08 rapor (Task 10) ve teşhis bunu doğrudan bu şekilde okur.
        stats_path.write_text(json.dumps(
            {sid: {"depletion_rate": ps["depletion_rate"]} for sid, ps in per_sample.items()},
            indent=2,
        ))
        write_gate_results(run_dir, gates)
        for g in gates:
            log(f"gate {g.name}: {g.status} — {g.message}")
        log(f"rrna depletion statistics written: {stats_path}")

    summary = {
        "n_samples": len(samples),
        "samples": {sid: {"depletion_rate": ps["depletion_rate"]} for sid, ps in per_sample.items()},
        "gate_counts": dict(Counter(g.status for g in gates)),
    }
    state.mark_done(MODULE_NAME, [str(stats_path), str(log_path)])
    return summary
