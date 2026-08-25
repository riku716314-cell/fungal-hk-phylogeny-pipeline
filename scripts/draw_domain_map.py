#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SMART の結果からドメイン構成図 (SVG) を描く。

タンパク質を横棒で表し、ドメインを色付き矩形で重ねる形式。
論文の Figure でよく使われる表現で、複数タンパク質のドメイン構成を
並べて比較するのに適している。

依存ライブラリ無しで SVG を直接組み立てている。SVG はベクタ形式なので
拡大しても劣化せず、Illustrator や PowerPoint で個別要素として編集できる。
図の再生成がコマンド 1 つで済むため、データ修正のたびに手作業で描き直す
必要がない。

SMART の出力には STATUS=hidden|threshold の候補も大量に含まれる。
これらは閾値に達していない予測なので、visible のもののみを描画する。

使用例:
    python draw_domain_map.py \\
        --smart-dir smart_results/ --fasta seqs.fasta \\
        --groups groups.tsv --out domain_map.svg
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

# 描画するドメインと配色。ここに無いものは描画しない。
DOMAIN_STYLE = {
    "PAS":                   ("#e8232a", "PAS"),
    "PAC":                   ("#c9cde8", "PAC"),
    "GAF":                   ("#111111", "GAF"),
    "HAMP":                  ("#8ab5e0", "HAMP"),
    "HisKA":                 ("#f5e400", "HisKA"),
    "HATPase_c":             ("#f0a01e", "ATPase"),
    "REC":                   ("#2b8f9e", "REC"),
    "STYKc":                 ("#3d5c2a", "Protein kinase"),
    "transmembrane_domain":  ("#a34fc4", "TM"),
}
GROUP_ORDER = ["I", "II", "III", "IV", "V", "VI", "VII",
               "VIII", "IX", "X", "XI", "OUTGROUP", "UNASSIGNED"]

# レイアウト (px)
LABEL_W, PLOT_W, LEGEND_W = 210, 700, 210
TOP, ROW_H, BAR_H, AXIS_GAP = 26, 34, 19, 26


def parse_smart(path: Path) -> list[dict]:
    """SMART の結果から visible なドメインのみ抽出する。"""
    records, cur = [], {}
    with path.open(encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if line.startswith("DOMAIN="):
                if cur:
                    records.append(cur)
                cur = {"DOMAIN": line.split("=", 1)[1]}
            elif "=" in line and cur:
                k, v = line.split("=", 1)
                cur[k] = v
        if cur:
            records.append(cur)
    doms = [{"name": r["DOMAIN"], "start": int(r["START"]), "end": int(r["END"])}
            for r in records
            if r.get("DOMAIN") in DOMAIN_STYLE
            and r.get("STATUS", "").startswith("visible")]
    return sorted(doms, key=lambda d: d["start"])


def read_lengths(path: Path) -> dict[str, int]:
    lengths, header, buf = {}, None, []
    with path.open(encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if header:
                    lengths[header] = len("".join(buf))
                header, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
        if header:
            lengths[header] = len("".join(buf))
    return lengths


def esc(s: str) -> str:
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--smart-dir", type=Path, required=True)
    p.add_argument("--fasta", type=Path, required=True)
    p.add_argument("--groups", type=Path, default=None,
                   help="グループ TSV。指定するとグループ順に並べ替える")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--max-aa", type=int, default=None,
                   help="軸の上限。既定は最長配列から自動決定")
    args = p.parse_args()

    lengths = read_lengths(args.fasta)

    groups: dict[str, str] = {}
    if args.groups:
        with args.groups.open(encoding="utf-8") as f:
            for row in csv.reader(f, delimiter="\t"):
                if len(row) >= 2 and not row[0].startswith("#") \
                        and row[0].strip().lower() != "sequence_id":
                    groups[row[0].strip()] = row[1].strip()

    entries = []
    for path in sorted(args.smart_dir.glob("*_SMART_results.txt")):
        name = path.name.replace("_SMART_results.txt", "")
        if name not in lengths:
            print(f"WARNING: {name} が FASTA にありません。スキップします", file=sys.stderr)
            continue
        entries.append({"id": name, "length": lengths[name],
                        "group": groups.get(name, ""),
                        "domains": parse_smart(path)})

    if not entries:
        print("ERROR: 描画対象がありません", file=sys.stderr)
        return 1

    def sort_key(e):
        g = e["group"]
        gi = GROUP_ORDER.index(g) if g in GROUP_ORDER else len(GROUP_ORDER)
        return (gi, -e["length"])
    entries.sort(key=sort_key)

    max_aa = args.max_aa or (max(e["length"] for e in entries) // 500 + 1) * 500
    scale = PLOT_W / max_aa
    svg_w = LABEL_W + PLOT_W + LEGEND_W + 40
    svg_h = TOP + ROW_H * len(entries) + AXIS_GAP + 56
    x0 = LABEL_W

    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{svg_w}" height="{svg_h}" '
           f'viewBox="0 0 {svg_w} {svg_h}" font-family="Arial, Helvetica, sans-serif">',
           f'<rect width="{svg_w}" height="{svg_h}" fill="white"/>']

    for i, e in enumerate(entries):
        y = TOP + ROW_H * i
        yc = y + BAR_H / 2 + 4
        out.append(f'<text x="{x0 - 62}" y="{yc}" font-size="12.5" '
                   f'text-anchor="end" fill="#111">{esc(e["id"])}</text>')
        if e["group"]:
            out.append(f'<text x="{x0 - 8}" y="{yc}" font-size="11.5" '
                       f'text-anchor="end" fill="#777">Group {esc(e["group"])}</text>')
        bar_w = e["length"] * scale
        out.append(f'<rect x="{x0}" y="{y}" width="{bar_w:.1f}" height="{BAR_H}" '
                   f'fill="white" stroke="#111" stroke-width="1.2"/>')
        for d in e["domains"]:
            dx = x0 + d["start"] * scale
            dw = max((d["end"] - d["start"]) * scale, 3.0)
            color, _ = DOMAIN_STYLE[d["name"]]
            out.append(f'<rect x="{dx:.1f}" y="{y}" width="{dw:.1f}" height="{BAR_H}" '
                       f'fill="{color}" stroke="#111" stroke-width="1.0"/>')

    # 目盛り軸
    ay = TOP + ROW_H * len(entries) + AXIS_GAP - 12
    out.append(f'<line x1="{x0}" y1="{ay}" x2="{x0 + PLOT_W}" y2="{ay}" '
               f'stroke="#111" stroke-width="2"/>')
    step = 500 if max_aa > 1500 else 250
    for tick in range(0, max_aa + 1, step):
        tx = x0 + tick * scale
        out.append(f'<line x1="{tx:.1f}" y1="{ay - 6}" x2="{tx:.1f}" y2="{ay + 6}" '
                   f'stroke="#111" stroke-width="2"/>')
        out.append(f'<text x="{tx:.1f}" y="{ay + 25}" font-size="12.5" '
                   f'text-anchor="middle" fill="#111">{tick}</text>')
    out.append(f'<text x="{x0 + PLOT_W / 2:.1f}" y="{ay + 46}" font-size="12" '
               f'text-anchor="middle" fill="#555">amino acid position</text>')

    # 凡例 (実際に使われたドメインのみ)
    lx, ly = x0 + PLOT_W + 34, TOP + 2
    for name, (color, disp) in DOMAIN_STYLE.items():
        if not any(d["name"] == name for e in entries for d in e["domains"]):
            continue
        out.append(f'<rect x="{lx}" y="{ly}" width="17" height="14" '
                   f'fill="{color}" stroke="#111" stroke-width="1"/>')
        out.append(f'<text x="{lx + 24}" y="{ly + 12}" font-size="13" fill="#111">{disp}</text>')
        ly += 22

    out.append("</svg>")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(out), encoding="utf-8")

    print(f"{len(entries)} 配列を描画: {args.out}\n")
    for e in entries:
        ds = " ".join(f'{DOMAIN_STYLE[d["name"]][1]}({d["start"]}-{d["end"]})'
                      for d in e["domains"])
        print(f'  {e["id"]:24} {e["length"]:5d} aa  {ds}')
    return 0


if __name__ == "__main__":
    sys.exit(main())
