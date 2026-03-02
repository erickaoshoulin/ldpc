#!/usr/bin/env python3
import csv
import math
from collections import defaultdict
from pathlib import Path

RESULTS_CSV = Path("results/performance.csv")
OUT_PNG = Path("results/performance.png")
OUT_SVG = Path("results/performance.svg")


def read_rows(path: Path):
    rows = defaultdict(lambda: {"snr": [], "ber": [], "fer": [], "avg_iterations": []})
    with path.open("r", newline="") as fp:
        reader = csv.DictReader(fp)
        for row in reader:
            name = row["algorithm"]
            rows[name]["snr"].append(float(row["snr_db"]))
            rows[name]["ber"].append(float(row["ber"]))
            rows[name]["fer"].append(float(row["fer"]))
            rows[name]["avg_iterations"].append(float(row["avg_iterations"]))
    return rows


def try_matplotlib(rows):
    try:
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        return False

    fig, axes = plt.subplots(2, 2, figsize=(13, 8))
    axes = axes.flatten()
    styles = {"conventional": "o", "rmas1": "s", "rmas2": "^"}
    colors = {"conventional": "#1f77b4", "rmas1": "#ff7f0e", "rmas2": "#2ca02c"}

    for alg, values in rows.items():
        marker = styles.get(alg, "o")
        color = colors.get(alg, "#555")
        axes[0].semilogy(values["snr"], values["ber"], marker=marker, color=color, label=alg.upper(), linewidth=2.0)
        axes[1].semilogy(values["snr"], values["fer"], marker=marker, color=color, label=alg.upper(), linewidth=2.0)
        axes[2].plot(values["snr"], values["avg_iterations"], marker=marker, color=color, label=alg.upper(), linewidth=2.0)

    best_snr = max(max(v["snr"]) for v in rows.values())
    ordered = []
    for alg, values in rows.items():
        idx = values["snr"].index(best_snr)
        ordered.append((alg, values["ber"][idx], values["fer"][idx], values["avg_iterations"][idx]))
    ordered.sort(key=lambda x: x[1])

    table = axes[3]
    table.axis("off")
    table.set_title(f"Summary @ {best_snr:.1f} dB")
    txt = ["Algorithm   BER        FER        AvgIter"]
    for alg, ber, fer, it in ordered:
        txt.append(f"{alg.upper():<10} {ber:>8.2e} {fer:>8.2e} {it:>8.2f}")
    table.text(0.02, 0.95, "\n".join(txt), va="top", family="monospace", fontsize=10)

    axes[0].set_title("BER vs SNR")
    axes[1].set_title("FER vs SNR")
    axes[2].set_title("Average Iterations vs SNR")
    for ax in axes[:3]:
        ax.set_xlabel("SNR (dB)")
        ax.grid(True, linestyle="--", alpha=0.4)
        ax.legend()
    axes[0].set_ylabel("BER")
    axes[1].set_ylabel("FER")
    axes[2].set_ylabel("Avg Iterations")
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=180)
    return True


def write_svg(rows):
    width, height = 1260, 640
    pad = 50
    panel_w = (width - 4 * pad) // 3
    panel_h = 260

    colors = {"conventional": "#1f77b4", "rmas1": "#ff7f0e", "rmas2": "#2ca02c"}

    snr_vals = [x for v in rows.values() for x in v["snr"]]
    x_min, x_max = min(snr_vals), max(snr_vals)

    def panel(metric, title, idx, y_log=False):
        x0 = pad + idx * (panel_w + pad)
        y0 = pad
        lines = [f'<rect x="{x0}" y="{y0}" width="{panel_w}" height="{panel_h}" fill="white" stroke="#bbb"/>',
                 f'<text x="{x0 + panel_w/2:.1f}" y="{y0 - 12}" text-anchor="middle" font-size="14">{title}</text>']
        vals = [y for v in rows.values() for y in v[metric]]
        y_min = min(vals)
        y_max = max(vals)
        if y_log:
            y_min = max(y_min, 1e-8)

        def tx(x):
            return x0 + 30 + (x - x_min) / (x_max - x_min + 1e-9) * (panel_w - 40)

        def ty(y):
            if y_log:
                lo, hi = math.log10(y_min), math.log10(y_max + 1e-12)
                yy = (math.log10(max(y, y_min)) - lo) / (hi - lo + 1e-9)
            else:
                yy = (y - y_min) / (y_max - y_min + 1e-9)
            return y0 + panel_h - 20 - yy * (panel_h - 40)

        for alg, v in rows.items():
            pts = " ".join(f"{tx(x):.1f},{ty(y):.1f}" for x, y in zip(v["snr"], v[metric]))
            lines.append(f'<polyline fill="none" stroke="{colors.get(alg,"#333")}" stroke-width="2" points="{pts}"/>')
            if v["snr"]:
                lines.append(f'<text x="{tx(v["snr"][-1])+6:.1f}" y="{ty(v[metric][-1]):.1f}" font-size="10" fill="{colors.get(alg,"#333")}">{alg.upper()}</text>')
        return "\n".join(lines)

    summary_y = pad + panel_h + 70
    right_snr = max(max(v["snr"]) for v in rows.values())
    summary = []
    for alg, values in rows.items():
        idx = values["snr"].index(right_snr)
        summary.append((alg, values["ber"][idx], values["fer"][idx], values["avg_iterations"][idx]))
    summary.sort(key=lambda x: x[1])

    summary_lines = [
        f'<text x="{pad}" y="{summary_y}" font-size="18" font-weight="bold">Performance Summary @ {right_snr:.1f} dB</text>',
        f'<text x="{pad}" y="{summary_y + 28}" font-size="13" font-family="monospace">Algorithm    BER         FER         AvgIter</text>'
    ]
    for i, (alg, ber, fer, it) in enumerate(summary):
        color = colors.get(alg, "#333")
        y = summary_y + 55 + i * 24
        summary_lines.append(
            f'<text x="{pad}" y="{y}" font-size="13" fill="{color}" font-family="monospace">{alg.upper():<10} {ber:>10.3e} {fer:>10.3e} {it:>8.3f}</text>'
        )

    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
           '<rect width="100%" height="100%" fill="#f8f9fb"/>',
           panel("ber", "BER vs SNR", 0, y_log=True),
           panel("fer", "FER vs SNR", 1, y_log=True),
           panel("avg_iterations", "Average Iterations vs SNR", 2, y_log=False),
           *summary_lines,
           '</svg>']
    OUT_SVG.write_text("\n".join(svg))


def main():
    if not RESULTS_CSV.exists():
        raise SystemExit(f"missing input CSV: {RESULTS_CSV}")
    rows = read_rows(RESULTS_CSV)
    made_png = try_matplotlib(rows)
    write_svg(rows)
    if made_png:
        print(f"saved {OUT_PNG} and {OUT_SVG}")
    else:
        print(f"saved {OUT_SVG} (install matplotlib to also export PNG)")


if __name__ == "__main__":
    main()
