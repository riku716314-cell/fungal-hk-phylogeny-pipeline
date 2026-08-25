#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SMART (https://smart.embl.de) にバッチで配列を投げ、ドメイン予測を取得する。

SMART は Web UI しか公開しておらず、公式には Perl 製のバッチスクリプトが
配布されている。本スクリプトはその通信仕様を Python で再実装したもの。

実装上の要点:
  1. セッション Cookie が必須。いきなり解析エンドポイントに POST しても
     モード選択ページが返るだけで解析は走らない。先に change_mode.cgi へ
     アクセスして normal モードのセッションを確立する必要がある。
  2. 解析は非同期。POST のレスポンスに含まれる jobId をポーリングして
     完了を待つ。
  3. 取得済みの配列はスキップする。SMART は公共サーバであり、
     再実行のたびに全件投げ直すのは避けるべき。

Pfam ベースの HMMER と SMART では検出されるドメインが一致しないことがある。
どちらが正しいという話ではなく、プロファイルと閾値の設計が異なるため。
先行研究と比較するときは、その研究が使ったツールに揃える必要がある。

使用例:
    python smart_batch_query.py --fasta seqs.fasta --outdir smart_results/
"""
from __future__ import annotations

import argparse
import http.cookiejar
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

MODE_URL = "https://smart.embl.de/smart/change_mode.cgi?mode=normal"
SUBMIT_URL = "https://smart.embl.de/smart/show_motifs.pl"
STATUS_URL = "https://smart.embl.de/results.cgi"

POLL_INTERVAL = 10    # 秒
POLL_MAX = 60         # 最大待機 = POLL_INTERVAL * POLL_MAX
REQUEST_GAP = 3       # 配列間の待機。公共サーバへの負荷を抑える


def build_opener() -> urllib.request.OpenerDirector:
    jar = http.cookiejar.CookieJar()
    opener = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(jar))
    opener.addheaders = [("User-Agent", "smart-batch-python/1.0")]
    # normal モードのセッションを確立しておく (これが無いと解析が走らない)
    opener.open(MODE_URL, timeout=30).read()
    return opener


def read_fasta(path: Path) -> dict[str, str]:
    seqs, header, buf = {}, None, []
    with path.open(encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if header:
                    seqs[header] = "".join(buf)
                header, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
        if header:
            seqs[header] = "".join(buf)
    return seqs


def query_one(opener, seq: str) -> str | None:
    data = urllib.parse.urlencode({"SEQUENCE": seq, "TEXTONLY": "1"}).encode()
    with opener.open(SUBMIT_URL, data=data, timeout=120) as r:
        content = r.read().decode("utf-8", errors="replace")

    # 即座に結果が返る場合
    if "DOMAIN=" in content:
        return content

    # 非同期ジョブの場合は jobId をポーリング
    m = re.search(r"jobId\s*=\s*'(\d+)'", content)
    if not m:
        return None
    job_id = m.group(1)
    for _ in range(POLL_MAX):
        time.sleep(POLL_INTERVAL)
        url = f"{STATUS_URL}?jobid={job_id}"
        with opener.open(url, timeout=120) as r:
            body = r.read().decode("utf-8", errors="replace")
        if "DOMAIN=" in body:
            return body
    return None


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fasta", type=Path, required=True)
    p.add_argument("--outdir", type=Path, required=True)
    args = p.parse_args()

    seqs = read_fasta(args.fasta)
    args.outdir.mkdir(parents=True, exist_ok=True)
    print(f"{len(seqs)} 配列を処理します\n")

    opener = build_opener()
    ok = skipped = failed = 0

    for i, (name, seq) in enumerate(sorted(seqs.items()), 1):
        out = args.outdir / f"{name}_SMART_results.txt"
        if out.exists() and out.stat().st_size > 0:
            print(f"[{i}/{len(seqs)}] {name}: 取得済みのためスキップ")
            skipped += 1
            continue

        print(f"[{i}/{len(seqs)}] {name} ({len(seq)} aa) 問い合わせ中 ...", flush=True)
        try:
            result = query_one(opener, seq)
        except Exception as e:                      # noqa: BLE001
            print(f"    ERROR: {e}", file=sys.stderr)
            failed += 1
            continue

        if result:
            out.write_text(result, encoding="utf-8")
            n = result.count("DOMAIN=")
            print(f"    OK ({n} ドメイン記載) -> {out.name}")
            ok += 1
        else:
            print("    FAILED: 結果を取得できませんでした", file=sys.stderr)
            failed += 1

        time.sleep(REQUEST_GAP)

    print(f"\n成功 {ok} / スキップ {skipped} / 失敗 {failed}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
