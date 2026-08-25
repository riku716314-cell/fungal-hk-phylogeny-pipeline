#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""HMMER の tblout を統合し、2/3 ドメインルールで HK 候補を絞り込む。

ヒスチジンキナーゼ (HK) は HisKA / HATPase_c / REC の 3 ドメインで特徴づけられる。
Mina et al. (2024) は「3 つのうち 2 つ以上を持つ」配列を HK 候補とみなす基準を
採用しており、本スクリプトはその基準を実装したもの。

3 ドメイン全てを必須にすると、REC を持たない非ハイブリッド型 HK を取りこぼす。
逆に 1 ドメインでは他のシグナル伝達因子が大量に混入する。2/3 はその折衷点。

使用例:
    python filter_hk_candidates.py \\
        --hiska   hiska.tblout \\
        --hatpase hatpase.tblout \\
        --rec     rec.tblout \\
        --evalue  1e-5 \\
        --out     hk_candidates.tsv
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

DOMAINS = ("HisKA", "HATPase_c", "REC")
MIN_DOMAINS = 2  # 2/3 ルール


def parse_tblout(path: Path, evalue_threshold: float) -> dict[str, float]:
    """hmmsearch --tblout を読み、閾値を満たす protein_id -> E-value を返す。

    tblout は空白区切りだが、末尾の description 列に空白が含まれるため
    列数を固定した split はできない。maxsplit で必要な列だけ切り出す。
    """
    hits: dict[str, float] = {}
    with path.open(encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = re.split(r"\s+", line.strip(), maxsplit=18)
            if len(fields) < 5:
                continue
            protein_id, full_evalue = fields[0], float(fields[4])
            if full_evalue > evalue_threshold:
                continue
            # 同一配列に複数ヒットがある場合は最良 (最小) の E-value を保持
            if protein_id not in hits or full_evalue < hits[protein_id]:
                hits[protein_id] = full_evalue
    return hits


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--hiska", type=Path, required=True, help="HisKA の tblout")
    p.add_argument("--hatpase", type=Path, required=True, help="HATPase_c の tblout")
    p.add_argument("--rec", type=Path, required=True, help="REC の tblout")
    p.add_argument("--evalue", type=float, default=1e-5, help="E-value 閾値 (既定: 1e-5)")
    p.add_argument("--out", type=Path, required=True, help="出力 TSV")
    args = p.parse_args()

    tblouts = {"HisKA": args.hiska, "HATPase_c": args.hatpase, "REC": args.rec}

    for name, path in tblouts.items():
        if not path.exists():
            print(f"ERROR: {name} の tblout が見つかりません: {path}", file=sys.stderr)
            return 1

    # protein_id -> {domain: evalue}
    profile: dict[str, dict[str, float]] = {}
    for domain, path in tblouts.items():
        for pid, ev in parse_tblout(path, args.evalue).items():
            profile.setdefault(pid, {})[domain] = ev
        print(f"{domain:12} : {len(parse_tblout(path, args.evalue)):5d} hits "
              f"(E-value <= {args.evalue})")

    candidates = {pid: d for pid, d in profile.items() if len(d) >= MIN_DOMAINS}

    print(f"\n全ヒット配列数        : {len(profile)}")
    print(f"{MIN_DOMAINS}/3 ルール通過 (HK候補) : {len(candidates)}")

    # ドメイン構成の内訳を出す (全 3 ドメイン型 = ハイブリッド型 HK の目安)
    by_combo: dict[str, int] = {}
    for d in candidates.values():
        combo = "+".join(sorted(d, key=DOMAINS.index))
        by_combo[combo] = by_combo.get(combo, 0) + 1
    print("\nドメイン構成の内訳:")
    for combo, n in sorted(by_combo.items(), key=lambda kv: -kv[1]):
        print(f"  {combo:28} {n:4d}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8", newline="") as f:
        f.write("protein_id\tdomains\tn_domains\t" +
                "\t".join(f"evalue_{d}" for d in DOMAINS) + "\n")
        for pid, d in sorted(candidates.items()):
            found = sorted(d, key=DOMAINS.index)
            evs = "\t".join(f"{d[x]:.2e}" if x in d else "NA" for x in DOMAINS)
            f.write(f"{pid}\t{'+'.join(found)}\t{len(found)}\t{evs}\n")

    print(f"\n出力: {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
