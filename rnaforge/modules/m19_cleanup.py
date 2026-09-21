"""m19 — Koşu-sonu ara-dosya temizliği (opsiyonel, son adım).

Trimlenmiş FASTQ'lar (`trimmed/`) ve hizalama BAM'leri (`quantification/**/
aligned.sorted.bam[.bai]`) tümüyle yeniden üretilebilir ara ürünlerdir ve diski
şişirir (tipik bir koşuda toplam boyutun ~%95'i). m19 bunları TÜM downstream tüketici
(m16 seqqc + m17 alignqc dahil) bittikten SONRA, `rnaforge run`'ın en son adımı olarak
siler. Sayım matrisi (counts.tsv/tpm/fpkm/quant.sf/nanocount.tsv), DE, figürler, rapor
ve loglar ASLA silinmez.

Kapı üretmez. `config.cleanup.remove_intermediates=False` iken no-op (yalnız log).
`keep_bam=True` BAM'leri korur, yalnız trimmed'i siler. m04/m05'in bittiğini şart koşar;
tekrar çalıştırmada (resume) zaten silinmiş dosyalar sessizce atlanır — yüksek sesle
loglanır, sessiz yutma yok."""
from __future__ import annotations

import json
from pathlib import Path

from rnaforge.config import Config
from rnaforge.state import RunState

MODULE_NAME = "m19_cleanup"


def _size_of(path: Path) -> int:
    """Bir dosya ya da dizinin toplam bayt boyutu (semboller izlenmez)."""
    if path.is_file():
        return path.stat().st_size
    total = 0
    for p in path.rglob("*"):
        if p.is_file():
            total += p.stat().st_size
    return total


def _human(nbytes: int) -> str:
    size = float(nbytes)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if size < 1024 or unit == "TB":
            return f"{size:.1f} {unit}"
        size /= 1024
    return f"{size:.1f} TB"


def run_cleanup(config: Config, metadata_path: Path, run_dir: Path,
                force: bool = False) -> dict:
    run_dir = Path(run_dir)
    stats_dir = run_dir / "statistics"
    logs_dir = run_dir / "logs"
    for d in (stats_dir, logs_dir):
        d.mkdir(parents=True, exist_ok=True)
    state = RunState(run_dir)
    stats_path = stats_dir / "cleanup_statistics.json"

    if not force and state.is_done(MODULE_NAME) and stats_path.exists():
        summary = json.loads(stats_path.read_text())
        summary["resumed"] = True
        return summary

    log_path = logs_dir / "cleanup.log"
    with log_path.open("w") as log_file:
        def log(msg: str) -> None:
            log_file.write(msg + "\n")
            log_file.flush()

        if not config.cleanup.remove_intermediates:
            log("cleanup.remove_intermediates=false → ara dosyalar KORUNDU (no-op).")
            summary = {"removed": False, "freed_bytes": 0, "removed_paths": []}
            stats_path.write_text(json.dumps(summary, indent=2))
            state.mark_done(MODULE_NAME, [str(stats_path), str(log_path)])
            return summary

        keep_bam = config.cleanup.keep_bam
        freed = 0
        removed_paths: list[str] = []

        # 1) trimlenmiş FASTQ'lar (m03 çıktısı) — her zaman silinir.
        trimmed_dir = run_dir / "trimmed"
        if trimmed_dir.is_dir():
            state.heartbeat()
            n = _size_of(trimmed_dir)
            _rmtree(trimmed_dir)
            freed += n
            removed_paths.append(str(trimmed_dir))
            log(f"silindi: {trimmed_dir} ({_human(n)})")
        else:
            log(f"trimmed/ yok, atlandı: {trimmed_dir}")

        # 2) hizalama BAM'leri (m04 çıktısı) — keep_bam=False iken silinir. Sayım
        #    tabloları (counts.tsv/quant.sf/nanocount.tsv) aynı dizinde KALIR.
        quant_dir = run_dir / "quantification"
        if keep_bam:
            log("cleanup.keep_bam=true → BAM'ler korundu (yalnız trimmed silindi).")
        elif quant_dir.is_dir():
            bams = sorted(quant_dir.rglob("aligned.sorted.bam")) + \
                sorted(quant_dir.rglob("aligned.sorted.bam.bai"))
            for bam in bams:
                state.heartbeat()
                n = bam.stat().st_size
                bam.unlink()
                freed += n
                removed_paths.append(str(bam))
                log(f"silindi: {bam} ({_human(n)})")
            if not bams:
                log(f"BAM bulunamadı, atlandı: {quant_dir}")

        log(f"toplam geri kazanılan: {_human(freed)} ({len(removed_paths)} yol)")
        summary = {
            "removed": True,
            "keep_bam": keep_bam,
            "freed_bytes": freed,
            "freed_human": _human(freed),
            "removed_paths": removed_paths,
        }
        stats_path.write_text(json.dumps(summary, indent=2))
        state.mark_done(MODULE_NAME, [str(stats_path), str(log_path)])
    return summary


def _rmtree(path: Path) -> None:
    """shutil.rmtree yerine stdlib-yalın, sembolleri izlemeyen özyineli silme."""
    for child in path.iterdir():
        if child.is_dir() and not child.is_symlink():
            _rmtree(child)
        else:
            child.unlink()
    path.rmdir()
