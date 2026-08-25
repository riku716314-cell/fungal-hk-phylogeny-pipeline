#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""既知の分類を持つアンカー配列から、系統樹上の未分類配列へグループを伝播する。

糸状菌の HK は Catlett et al. (2003) 以来 Group I〜XI に分類されており、
先行研究でグループが確定している配列 (アンカー) を系統樹に含めておけば、
未知の配列がどのアンカーと単系統群を成すかでグループを推定できる。

アルゴリズム:
  各アンカーについて、そのアンカーを含むクレードを根の方向へ広げていき、
  「別グループのアンカーが入り込む直前」で停止する。得られたクレードの
  葉すべてにそのアンカーのグループを与える。
  同一グループの複数アンカーは衝突とみなさない (Group XI のように
  1 種が複数のパラログを持つ場合、それらは同じクレードに入るのが正常なため)。

外群の扱い (重要):
  外群はグループ分類の枠組みそのものの外にある。クレード判定に含めると、
  外群と隣接するクレードが不自然に広がり誤ったグループが伝播する。
  実際に本パイプラインでは大腸菌 EnvZ が Group I と判定される誤りが生じ、
  外群を判定対象から明示的に除外することで解消した。

どのアンカーのクレードにも入らなかった配列は、最近縁アンカーへの
patristic distance によるフォールバック判定とし、確度が低いことを
`fallback` フラグで明示する。推定を断定と区別して記録するための措置。

使用例:
    python assign_groups.py \\
        --tree tree.contree --anchors anchors.tsv \\
        --outgroup EcoliEnvZ --out groups.tsv
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


def _import_phylo():
    """Biopython は --help を出すだけなら不要なので、実行時に import する。

    ImportError を一律「未インストール」と扱うと、依存 DLL の読み込み失敗
    (NumPy 等) まで同じ案内になってしまい原因を見誤る。両者を区別する。
    """
    try:
        from Bio import Phylo
    except ModuleNotFoundError:
        sys.exit("ERROR: Biopython が見つかりません:  pip install biopython")
    except ImportError as e:
        sys.exit(f"ERROR: Biopython の読み込みに失敗しました: {e}\n"
                 "       インストールはされているが依存ライブラリ (NumPy 等) を "
                 "読み込めない状態です。")
    return Phylo


def load_anchors(path: Path) -> dict[str, str]:
    """TSV (sequence_id, group) を読む。# 始まりはコメント。"""
    anchors: dict[str, str] = {}
    with path.open(encoding="utf-8") as f:
        for row in csv.reader(f, delimiter="\t"):
            if not row or row[0].startswith("#") or len(row) < 2:
                continue
            if row[0].strip().lower() in ("sequence_id", "id"):
                continue
            anchors[row[0].strip()] = row[1].strip()
    return anchors


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tree", type=Path, required=True, help="Newick 形式の系統樹")
    p.add_argument("--anchors", type=Path, required=True,
                   help="アンカー定義 TSV (sequence_id<TAB>group)")
    p.add_argument("--outgroup", default=None, help="外群の配列名")
    p.add_argument("--out", type=Path, required=True, help="出力 TSV")
    args = p.parse_args()

    Phylo = _import_phylo()
    tree = Phylo.read(str(args.tree), "newick")
    anchors = load_anchors(args.anchors)
    print(f"アンカー {len(anchors)} 件を読み込みました")

    all_leaves = [t.name for t in tree.get_terminals()]
    missing = [a for a in anchors if a not in all_leaves]
    if missing:
        print(f"WARNING: 系統樹に存在しないアンカー: {missing}", file=sys.stderr)

    if args.outgroup:
        if args.outgroup not in all_leaves:
            print(f"ERROR: 外群 {args.outgroup} が系統樹にありません", file=sys.stderr)
            return 1
        tree.root_with_outgroup(args.outgroup)
        print(f"外群 {args.outgroup} で再 root 化しました")

    def leaves_of(clade) -> set[str]:
        """外群は常に除外する。理由は docstring 参照。"""
        return {t.name for t in clade.get_terminals() if t.name != args.outgroup}

    assignments: dict[str, tuple[str, str]] = {}  # id -> (group, method)

    # --- アンカーごとにクレードを拡張 ---
    for anchor, group in anchors.items():
        if anchor not in all_leaves:
            continue
        path = tree.get_path(anchor)
        best: set[str] = {anchor}
        # 葉から根に向かって遡る
        for clade in reversed(path[:-1] if len(path) > 1 else path):
            members = leaves_of(clade)
            others = [a for a in anchors
                      if a in members and a != anchor and anchors[a] != group]
            if others:
                break  # 別グループのアンカーが混入する直前で停止
            best = members
        for leaf in best:
            # 既に別グループが割り当て済みなら、より小さいクレードの判定を優先
            if leaf not in assignments:
                assignments[leaf] = (group, "clade")

    # --- フォールバック: 最近縁アンカーへの距離で判定 ---
    unassigned = [l for l in all_leaves
                  if l not in assignments and l != args.outgroup]
    for leaf in unassigned:
        best_anchor, best_dist = None, float("inf")
        for anchor in anchors:
            if anchor not in all_leaves:
                continue
            d = tree.distance(leaf, anchor)
            if d < best_dist:
                best_anchor, best_dist = anchor, d
        if best_anchor:
            assignments[leaf] = (anchors[best_anchor], "fallback")

    if args.outgroup:
        assignments[args.outgroup] = ("OUTGROUP", "outgroup")

    # --- 出力 ---
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8", newline="") as f:
        f.write("sequence_id\tgroup\tmethod\n")
        for leaf in all_leaves:
            g, m = assignments.get(leaf, ("UNASSIGNED", "none"))
            f.write(f"{leaf}\t{g}\t{m}\n")

    n_clade = sum(1 for g, m in assignments.values() if m == "clade")
    n_fb = sum(1 for g, m in assignments.values() if m == "fallback")
    print(f"\nクレード判定 : {n_clade}")
    print(f"フォールバック: {n_fb}  (確度が低いため要確認)")

    # --- 自己検証: アンカー自身が正しく再現されるか ---
    wrong = [(a, g, assignments[a][0]) for a, g in anchors.items()
             if a in assignments and assignments[a][0] != g]
    if wrong:
        print(f"\n[NG] アンカー自身の再現に失敗 {len(wrong)} 件:", file=sys.stderr)
        for a, exp, got in wrong:
            print(f"       {a}: 期待 {exp} / 実際 {got}", file=sys.stderr)
        return 1
    print(f"[OK] アンカー {len(anchors)} 件すべて自身のグループを再現")
    print(f"\n出力: {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
