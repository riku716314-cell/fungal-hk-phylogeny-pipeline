#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""マージした FASTA を系統解析にかける前に検証する。

このスクリプトは実際の失敗から生まれた。データセットから一部の配列を除く際に
「`>` で始まるヘッダ行だけを削除する」という処理を書いたことがあり、
配列本体の行が残って直前の配列の末尾に連結された。結果として、本来 1220 aa の
配列が 7920 aa になっていた。アライメントも系統樹も正常に完走してしまい、
エラーは一切出なかった。

系統解析のツール群は「生物学的にありえない配列」を弾いてくれない。
だから入力段階で機械的に検査する必要がある。

検査項目:
  1. 重複ヘッダ  … 同名配列があるとツールによって挙動が変わる
  2. 空配列      … アライメント後に全ギャップ列となり脱落しうる
  3. 異常な配列長 … 中央値から大きく外れるものは連結事故の兆候
  4. 不正な文字   … 想定外の記号 (数字など) が混入していないか
  5. NCBI 照合   … --verify-ncbi で実際の配列長と一致するか問い合わせる

使用例:
    python validate_fasta.py --fasta merged.fasta
    python validate_fasta.py --fasta merged.fasta --verify-ncbi --acc-prefix KAJ,EHA
"""
from __future__ import annotations

import argparse
import re
import statistics
import sys
import urllib.parse
import urllib.request
from pathlib import Path

VALID_AA = set("ACDEFGHIKLMNPQRSTVWYBXZJUO*-")
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
OUTLIER_FACTOR = 3.0  # 中央値の何倍を異常とみなすか


def read_fasta(path: Path) -> tuple[dict[str, str], list[str]]:
    """ヘッダ単位で厳密にパースする。順序と重複情報を保つため list も返す。"""
    seqs: dict[str, str] = {}
    order: list[str] = []
    header, buf = None, []
    with path.open(encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if header is not None:
                    seqs[header] = "".join(buf)
                header = line[1:].split()[0] if line[1:].strip() else ""
                order.append(header)
                buf = []
            else:
                buf.append(line.strip())
        if header is not None:
            seqs[header] = "".join(buf)
    return seqs, order


def fetch_ncbi_lengths(accessions: list[str]) -> dict[str, int]:
    params = {"db": "protein", "id": ",".join(accessions),
              "rettype": "fasta", "retmode": "text"}
    url = f"{EUTILS}?{urllib.parse.urlencode(params)}"
    with urllib.request.urlopen(url, timeout=120) as r:
        text = r.read().decode("utf-8", errors="replace")
    lengths, header, buf = {}, None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if header:
                lengths[header] = len("".join(buf))
            header = line[1:].split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if header:
        lengths[header] = len("".join(buf))
    return lengths


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fasta", type=Path, required=True)
    p.add_argument("--verify-ncbi", action="store_true",
                   help="NCBI に問い合わせて配列長を照合する")
    p.add_argument("--acc-prefix", default="",
                   help="NCBI 照合の対象とするヘッダの接頭辞 (カンマ区切り, 例: KAJ,EHA)")
    args = p.parse_args()

    seqs, order = read_fasta(args.fasta)
    problems = 0

    print(f"配列数: {len(order)} (ユニークなヘッダ: {len(seqs)})\n")

    # --- 1. 重複ヘッダ ---
    dups = {h for h in order if order.count(h) > 1}
    if dups:
        problems += len(dups)
        print(f"[NG] 重複ヘッダ {len(dups)} 件: {sorted(dups)[:5]}")
    else:
        print("[OK] 重複ヘッダなし")

    # --- 2. 空配列 ---
    empty = [h for h, s in seqs.items() if not s]
    if empty:
        problems += len(empty)
        print(f"[NG] 空配列 {len(empty)} 件: {empty[:5]}")
    else:
        print("[OK] 空配列なし")

    # --- 3. 異常な配列長 ---
    lengths = {h: len(s) for h, s in seqs.items() if s}
    if lengths:
        med = statistics.median(lengths.values())
        outliers = {h: n for h, n in lengths.items()
                    if n > med * OUTLIER_FACTOR or n < med / OUTLIER_FACTOR}
        print(f"     配列長: 中央値 {med:.0f} / 最短 {min(lengths.values())} "
              f"/ 最長 {max(lengths.values())}")
        if outliers:
            print(f"[要確認] 中央値の {OUTLIER_FACTOR} 倍を超えて外れる配列 {len(outliers)} 件:")
            for h, n in sorted(outliers.items(), key=lambda kv: -kv[1])[:10]:
                print(f"           {h:24} {n:6d} aa")
            print("           連結事故の可能性。NCBI 等の一次情報と照合してください。")
        else:
            print("[OK] 極端な配列長の外れ値なし")

    # --- 4. 不正な文字 ---
    bad = {}
    for h, s in seqs.items():
        invalid = set(s.upper()) - VALID_AA
        if invalid:
            bad[h] = invalid
    if bad:
        problems += len(bad)
        print(f"[NG] 不正な文字を含む配列 {len(bad)} 件:")
        for h, chars in list(bad.items())[:5]:
            print(f"       {h}: {sorted(chars)}")
    else:
        print("[OK] 不正な文字なし")

    # --- 5. NCBI 照合 ---
    if args.verify_ncbi:
        prefixes = tuple(x.strip() for x in args.acc_prefix.split(",") if x.strip())
        targets = [h for h in seqs
                   if (not prefixes or h.startswith(prefixes))
                   and re.match(r"^[A-Z]{2,3}_?\d+\.\d+$", h)]
        if not targets:
            print("\n[skip] NCBI 照合対象のアクセッション形式ヘッダがありません")
        else:
            print(f"\nNCBI 照合中 ({len(targets)} 件) ...")
            remote = fetch_ncbi_lengths(targets)
            mismatch = [(h, len(seqs[h]), remote[h])
                        for h in targets if h in remote and len(seqs[h]) != remote[h]]
            missing = [h for h in targets if h not in remote]
            if mismatch:
                problems += len(mismatch)
                print(f"[NG] 配列長不一致 {len(mismatch)} 件:")
                for h, local, rem in mismatch:
                    print(f"       {h:24} local={local:6d}  NCBI={rem:6d}")
            if missing:
                print(f"[要確認] NCBI から取得できなかった {len(missing)} 件: {missing[:5]}")
            if not mismatch and not missing:
                print(f"[OK] {len(targets)} 件すべて NCBI と配列長が一致")

    print(f"\n{'=' * 56}")
    if problems:
        print(f"問題あり: {problems} 件。系統解析の前に修正してください。")
        return 1
    print("検証を通過しました。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
