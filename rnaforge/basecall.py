"""ONT ham-sinyal (FAST5/POD5) basecalling: dorado (GPU) + pod5 dönüşümü.

Pipeline'ın girdisi FASTQ'dur; ham sinyal geldiğinde m00 bunu FASTQ'ya çevirir.
Basecaller = dorado (GPU zorunlu — CPU pratik değil). Parser saftır."""
from __future__ import annotations

import re
import subprocess
from pathlib import Path


class BasecallRunError(RuntimeError):
    """dorado/pod5 çalıştırılamadı ya da beklenen çıktıyı üretmedi."""


def is_signal_input(path: Path | str) -> str | None:
    """Girdi ham sinyal mi? 'pod5' | 'fast5' | None (FASTQ/başka).

    Tek dosya (.pod5/.fast5) ya da bu uzantılı dosya içeren bir dizin sinyaldir.
    FASTQ (.fastq/.fq[.gz]) ve diğerleri None döner (basecalling gerekmez)."""
    p = Path(path)
    if p.is_dir():
        if any(p.glob("*.pod5")):
            return "pod5"
        if any(p.glob("*.fast5")):
            return "fast5"
        return None
    name = p.name.lower()
    if name.endswith(".pod5"):
        return "pod5"
    if name.endswith(".fast5"):
        return "fast5"
    return None


def convert_fast5_to_pod5(fast5: Path, out_pod5: Path,
                          env: str = "rnaforge-basecall") -> Path:
    """FAST5 (dosya ya da dizin) → tek POD5. dorado POD5'i güvenilir okur."""
    fast5 = Path(fast5)
    out_pod5 = Path(out_pod5)
    out_pod5.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["conda", "run", "-n", env, "pod5", "convert", "fast5",
           str(fast5), "--output", str(out_pod5), "--force-overwrite"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0 or not out_pod5.exists():
        raise BasecallRunError(
            f"pod5 convert fast5 failed (exit {r.returncode}) for {fast5}\n"
            f"stderr: {r.stderr.strip()[-800:]}"
        )
    return out_pod5


_BASECALLED_RE = re.compile(r"Simplex reads basecalled:\s*(\d+)")


def parse_basecalled_count(stderr_text: str) -> int | None:
    """dorado stderr'inden basecall edilen read sayısı; yoksa None."""
    m = _BASECALLED_RE.search(stderr_text)
    return int(m.group(1)) if m else None


def run_dorado(pod5: Path, out_fastq: Path, dorado_bin: Path | str = "dorado",
               model: str = "hac", device: str = "cuda:all",
               models_dir: Path | None = None) -> int:
    """dorado basecaller (GPU) → FASTQ. model 'hac' (kompleks) → dorado, POD5 run
    metadata'sından tam modeli otomatik seçer/indirir (DNA-cDNA vs RNA/kimya).
    Basecall edilen read sayısını döner. Sıfır-olmayan çıkış ya da boş çıktı →
    yüksek sesle hata (feedback_gurultulu_hata)."""
    pod5 = Path(pod5)
    out_fastq = Path(out_fastq)
    out_fastq.parent.mkdir(parents=True, exist_ok=True)
    log_path = out_fastq.with_suffix(".dorado.log")

    cmd = [str(dorado_bin), "basecaller", "--emit-fastq", "--device", device]
    if models_dir is not None:
        Path(models_dir).mkdir(parents=True, exist_ok=True)
        cmd += ["--models-directory", str(models_dir)]
    cmd += [model, str(pod5)]

    with out_fastq.open("w") as out_fh:
        r = subprocess.run(cmd, stdout=out_fh, stderr=subprocess.PIPE, text=True)
    log_path.write_text(r.stderr)
    if r.returncode != 0:
        raise BasecallRunError(
            f"dorado basecaller failed (exit {r.returncode}) on {pod5}\n"
            f"cmd: {' '.join(cmd)}\nstderr: {r.stderr.strip()[-1000:]}"
        )
    n = parse_basecalled_count(r.stderr)
    if n is None:
        n = sum(1 for _ in out_fastq.open()) // 4
    if n == 0 or out_fastq.stat().st_size == 0:
        raise BasecallRunError(
            f"dorado produced no reads from {pod5} (empty FASTQ). "
            f"stderr tail: {r.stderr.strip()[-600:]}"
        )
    return n
