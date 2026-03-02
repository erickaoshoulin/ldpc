#!/usr/bin/env python3
import csv
import math
from collections import defaultdict
from pathlib import Path

IN_CSV = Path("results/matrix_performance.csv")
OUT_SVG = Path("results/matrix_performance.svg")


def load_rows():
    data = defaultdict(lambda: defaultdict(lambda: {"snr": [], "ber": []}))
    with IN_CSV.open("r", newline="") as fp:
        for row in csv.DictReader(fp):
            p = row["profile"]
            a = row["algorithm"]
            data[p][a]["snr"].append(float(row["snr_db"]))
            data[p][a]["ber"].append(float(row["ber"]))
    return data


def main():
    if not IN_CSV.exists():
        raise SystemExit(f"missing {IN_CSV}")

    data = load_rows()
    profiles = [p for p in data.keys() if p != "custom"][:3]
    algs = ["conventional", "rmas1", "rmas2", "as"]
    colors = {"conventional": "#1f77b4", "rmas1": "#ff7f0e", "rmas2": "#2ca02c", "as": "#d62728"}

    width, height = 1320, 500
    pad = 40
    panel_w = 390
    panel_h = 320
    snr_vals = [x for p in profiles for a in algs for x in data[p][a]["snr"]]
    ber_vals = [max(y, 1e-8) for p in profiles for a in algs for y in data[p][a]["ber"]]
    x_min, x_max = min(snr_vals), max(snr_vals)
    y_min, y_max = min(ber_vals), max(ber_vals)
    log_lo, log_hi = math.log10(y_min), math.log10(y_max)

    def tx(x, x0):
        return x0 + 35 + (x - x_min) / (x_max - x_min + 1e-9) * (panel_w - 60)

    def ty(y, y0):
        yy = (math.log10(max(y, y_min)) - log_lo) / (log_hi - log_lo + 1e-9)
        return y0 + panel_h - 30 - yy * (panel_h - 55)

    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
             '<rect width="100%" height="100%" fill="#f8f9fb"/>',
             '<text x="40" y="24" font-size="18" font-weight="bold">BER comparison across IEEE-like matrix profiles</text>']

    for i, profile in enumerate(profiles):
        x0 = pad + i * (panel_w + 35)
        y0 = 60
        parts.append(f'<rect x="{x0}" y="{y0}" width="{panel_w}" height="{panel_h}" fill="white" stroke="#bbb"/>')
        parts.append(f'<text x="{x0 + panel_w/2:.1f}" y="{y0 - 10}" text-anchor="middle" font-size="14">{profile}</text>')
        for alg in algs:
            d = data[profile][alg]
            if not d["snr"]:
                continue
            pts = " ".join(f"{tx(x, x0):.1f},{ty(y, y0):.1f}" for x, y in zip(d["snr"], d["ber"]))
            parts.append(f'<polyline points="{pts}" fill="none" stroke="{colors[alg]}" stroke-width="2"/>')

    lx, ly = 50, 430
    for i, alg in enumerate(algs):
        xx = lx + i * 200
        parts.append(f'<line x1="{xx}" y1="{ly}" x2="{xx+25}" y2="{ly}" stroke="{colors[alg]}" stroke-width="3"/>')
        parts.append(f'<text x="{xx+32}" y="{ly+4}" font-size="12">{alg.upper()}</text>')

    parts.append('</svg>')
    OUT_SVG.write_text("\n".join(parts))
    print(f"saved {OUT_SVG}")


if __name__ == "__main__":
    main()
