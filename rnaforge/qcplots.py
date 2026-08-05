"""QC figürlerini üretmek için ince sarmalayıcı: spec dict'ini geçici JSON'a
yazar ve `qc_plot.py`'yi matplotlib'in bulunduğu ortamda (rnaforge-seqqc)
çalıştırır. Saf veri → figür; çağıran modüller matplotlib'e doğrudan bağımlı değil."""
from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path

_SCRIPT = Path(__file__).parent / "qc_plot.py"


class QCPlotError(RuntimeError):
    """QC figürü üretilemedi."""


def render_qc_figure(spec: dict, out_png: Path, env: str = "rnaforge-seqqc") -> Path:
    out_png = Path(out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as fh:
        json.dump(spec, fh)
        spec_path = fh.name
    cmd = ["conda", "run", "-n", env, "python", str(_SCRIPT), spec_path, str(out_png)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    Path(spec_path).unlink(missing_ok=True)
    if r.returncode != 0 or not out_png.exists():
        raise QCPlotError(
            f"qc_plot failed (exit {r.returncode})\ncmd: {' '.join(cmd)}\n"
            f"stderr: {r.stderr.strip()}"
        )
    return out_png
