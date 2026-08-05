"""Bağımsız QC figür üreticisi — yalnız matplotlib+json+sys kullanır, rnaforge
paketini import ETMEZ. Böylece matplotlib'in bulunduğu ortamda (rnaforge-seqqc)
`python qc_plot.py <spec.json> <out.png>` ile çalıştırılabilir; ana paketin
kurulu olması gerekmez.

Spec şeması (JSON):
  {"type": "lines", "title": str, "xlabel": str, "ylabel": str,
   "x": [..], "series": {"ad": [y..], ...}}
  {"type": "bars",  "title": str, "xlabel": str, "ylabel": str,
   "x": [labels..], "y": [values..]}
"""
from __future__ import annotations

import json
import sys


def render(spec: dict, out_png: str) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    kind = spec.get("type")
    fig, ax = plt.subplots(figsize=(spec.get("width", 8), spec.get("height", 3.5)))
    if kind == "lines":
        x = spec["x"]
        for name, ys in spec["series"].items():
            ax.plot(range(len(ys)), ys, label=name, linewidth=1.2)
        # x etiketleri seyreltilir (pozisyon çok olabilir)
        n = len(x)
        step = max(1, n // 12)
        ax.set_xticks(range(0, n, step))
        ax.set_xticklabels([x[i] for i in range(0, n, step)], rotation=45, ha="right", fontsize=7)
        ax.legend(loc="upper right", fontsize=8, ncol=len(spec["series"]))
        ax.set_ylim(bottom=0)
    elif kind == "bars":
        x = spec["x"]
        y = spec["y"]
        ax.bar(range(len(y)), y, color=spec.get("color", "#4C72B0"))
        step = max(1, len(x) // 20)
        ax.set_xticks(range(0, len(x), step))
        ax.set_xticklabels([x[i] for i in range(0, len(x), step)], rotation=45, ha="right", fontsize=7)
    else:
        raise SystemExit(f"qc_plot: unknown spec type {kind!r}")

    ax.set_title(spec.get("title", ""), fontsize=11)
    ax.set_xlabel(spec.get("xlabel", ""), fontsize=9)
    ax.set_ylabel(spec.get("ylabel", ""), fontsize=9)
    ax.grid(True, alpha=0.25, linewidth=0.5)
    fig.tight_layout()
    fig.savefig(out_png, dpi=120)
    plt.close(fig)


def main(argv: list[str]) -> int:
    if len(argv) != 3:
        print("usage: python qc_plot.py <spec.json> <out.png>", file=sys.stderr)
        return 2
    spec = json.loads(open(argv[1]).read())
    render(spec, argv[2])
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
