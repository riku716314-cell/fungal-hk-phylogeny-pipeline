#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""NCBI E-utilities でタンパク質配列とメタデータ (locus_tag / 配列長) を取得する。

配列ファイルを配布する代わりにアクセッション番号のリストを配布し、
取得は利用者の手元で行う方式にしている。データの出所が NCBI であることが
明示され、リポジトリ側で配列が古くなる問題も避けられる。

locus_tag の取得は、系統樹のラベルと先行研究の遺伝子名を突き合わせるために使う。
配列そのものは種名や自動命名 ("hypothetical protein" 等) しか持たないことが多く、
論文の Table と対応づけるには locus_tag が事実上の共通キーになる。

使用例:
    # FASTA を取得
    python fetch_ncbi_sequences.py --acc-list accessions.txt --out seqs.fasta

    # locus_tag と配列長の対応表を取得
    python fetch_ncbi_sequences.py --acc-list accessions.txt --metadata --out meta.tsv
"""
from __future__ import annotations

import argparse
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
BATCH_SIZE = 100          # 1 リクエストあたりの ID 数
SLEEP_SEC = 0.4           # NCBI 推奨: API キー無しは 3 req/sec 以下


def efetch(ids: list[str], rettype: str, email: str | None) -> str:
    params = {
        "db": "protein",
        "id": ",".join(ids),
        "rettype": rettype,
        "retmode": "text",
    }
    if email:
        params["email"] = email
    url = f"{EUTILS}?{urllib.parse.urlencode(params)}"
    with urllib.request.urlopen(url, timeout=120) as r:
        return r.read().decode("utf-8", errors="replace")


def parse_genpept(text: str) -> list[tuple[str, str, str]]:
    """GenPept から (accession, locus_tag, 配列長) を抽出する。"""
    rows = []
    for rec in text.split("//\n"):
        acc = re.search(r"VERSION\s+(\S+)", rec)
        if not acc:
            continue
        tag = re.search(r'/locus_tag="([^"]+)"', rec)
        length = re.search(r"LOCUS\s+\S+\s+(\d+) aa", rec)
        rows.append((acc.group(1),
                     tag.group(1) if tag else "NA",
                     length.group(1) if length else "NA"))
    return rows


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--acc-list", type=Path, required=True,
                   help="アクセッション番号を 1 行 1 件で書いたファイル (# 以降はコメント)")
    p.add_argument("--out", type=Path, required=True, help="出力先")
    p.add_argument("--metadata", action="store_true",
                   help="FASTA ではなく locus_tag/配列長の TSV を出力する")
    p.add_argument("--email", default=None,
                   help="NCBI に通知する連絡先メールアドレス (推奨)")
    args = p.parse_args()

    accessions = []
    for line in args.acc_list.read_text(encoding="utf-8").splitlines():
        line = line.split("#")[0].strip()
        if line:
            accessions.append(line)

    if not accessions:
        print("ERROR: アクセッション番号が 1 件も読み取れませんでした", file=sys.stderr)
        return 1
    print(f"{len(accessions)} 件のアクセッション番号を読み込みました")

    rettype = "gp" if args.metadata else "fasta"
    chunks = [accessions[i:i + BATCH_SIZE] for i in range(0, len(accessions), BATCH_SIZE)]

    results: list[str] = []
    for i, chunk in enumerate(chunks, 1):
        print(f"  取得中 {i}/{len(chunks)} ({len(chunk)} 件) ...")
        results.append(efetch(chunk, rettype, args.email))
        if i < len(chunks):
            time.sleep(SLEEP_SEC)

    args.out.parent.mkdir(parents=True, exist_ok=True)

    if args.metadata:
        rows = [r for text in results for r in parse_genpept(text)]
        with args.out.open("w", encoding="utf-8", newline="") as f:
            f.write("accession\tlocus_tag\tlength_aa\n")
            for acc, tag, length in rows:
                f.write(f"{acc}\t{tag}\t{length}\n")
        print(f"\n{len(rows)} 件のメタデータを出力: {args.out}")
        missing = len(accessions) - len(rows)
        if missing:
            print(f"WARNING: {missing} 件が取得できませんでした", file=sys.stderr)
    else:
        text = "".join(results)
        args.out.write_text(text, encoding="utf-8")
        n = text.count(">")
        print(f"\n{n} 配列を出力: {args.out}")
        if n != len(accessions):
            print(f"WARNING: 要求 {len(accessions)} 件に対し取得 {n} 件", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
