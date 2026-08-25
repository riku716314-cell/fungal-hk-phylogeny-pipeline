# fungal-hk-phylogeny-pipeline

糸状菌のヒスチジンキナーゼ (histidine kinase, HK) を同定し、系統樹上で
既知の分類群に対応づけ、ドメイン構成図まで自動生成する解析パイプラインです。

大学院の研究で必要になり、既製ツールでは繋がらない部分を自作しました。
**特定の研究データに依存しない形で、手法とスクリプトのみを公開しています。**

---

## 何を解決したかったのか

糸状菌の HK はゲノム中に十数個あり、環境応答や病原性に関わることが知られています。
新しい菌種でこれを調べるには、次の 4 つを順に片付ける必要があります。

| 課題 | 単体のツールでは足りない理由 |
|---|---|
| ゲノムから HK 候補を絞り込む | HMMER は「どのドメインが当たったか」しか出さない。複数の検索結果を統合して判定する処理が要る |
| 系統樹を作る | MAFFT・IQ-TREE で作れるが、入力データの破損を検出してくれない |
| 各配列を既知の分類群に対応づける | 系統樹は木を描くだけ。「この配列は Group V」という判定は人手か自作コードが必要 |
| ドメイン構成図を描く | 手作業だと、データを直すたびに図を描き直すことになる |

この 4 つを繋ぐ部分を書いたのがこのリポジトリです。

---

## パイプライン全体像

```
[1] ゲノムのタンパク質配列
        │  hmmsearch (HisKA / HATPase_c / REC)
        ▼
    filter_hk_candidates.py     ← 2/3 ドメインルールで HK 候補を確定
        │
        │  fetch_ncbi_sequences.py   ← 参照種の配列を NCBI から取得
        ▼
    validate_fasta.py           ← 系統解析にかける前に入力を検証【重要】
        │
        │  MAFFT → trimAl → IQ-TREE
        ▼
    assign_groups.py            ← 系統樹上でグループ (I〜XI) を伝播
        │
        │  smart_batch_query.py      ← SMART でドメイン予測
        ▼
    draw_domain_map.py          ← ドメイン構成図 (SVG) を生成
```

---

## 工夫した点

### 1. 「エラーが出ない失敗」を検出する仕組みを入れた

このパイプラインで最も時間を使った部分です。

データセットから一部の配列を除外する処理で、**`>` で始まるヘッダ行だけを削除する**
コードを書いたことがありました。配列本体の行が残り、直前の配列の末尾に連結されて、
本来 1220 aa の配列が **7920 aa** になっていました。

問題は、**アライメントも系統樹推定も正常に完走してしまった**ことです。
警告もエラーも一切出ず、それらしい系統樹が出力されました。
異常に気づいたのは、系統樹の枝の長さを眺めていて違和感を持ったときです。

系統解析ツールは「生物学的にありえない配列」を弾いてくれません。
そこで入力段階で機械的に検査する [`validate_fasta.py`](scripts/validate_fasta.py)
を書きました。重複ヘッダ・空配列・異常な配列長・不正な文字を調べ、
`--verify-ncbi` を付けると **NCBI に問い合わせて配列長を照合**します。

実際にこの事故を再現して検査にかけると、こう出ます。

```
配列数: 9 (ユニークなヘッダ: 9)
[OK] 重複ヘッダなし
     配列長: 中央値 1412 / 最短 1271 / 最長 3195
[OK] 極端な配列長の外れ値なし          ← 長さの統計だけでは見逃す
NCBI 照合中 (9 件) ...
[NG] 配列長不一致 1 件:
       EHA56736.1               local=  3195  NCBI=  1994
問題あり: 1 件。系統解析の前に修正してください。
```

**統計的な外れ値検出だけでは見逃していた**という点が重要でした。
3195 aa は中央値の 3 倍に届かないため閾値に掛かりません。
一次情報 (NCBI) と突き合わせる検査を足して初めて確実に捕まえられます。

### 2. 外群の扱いをアルゴリズムから分離した

系統樹のグループ判定で、外群である大腸菌 EnvZ が **Group I** と判定される
誤りが出ました。外群はグループ分類の枠組みそのものの外にある配列なので、
どのグループにも属さないのが正解です。

原因は、クレードを広げる処理が外群も「クレードの一員」として数えていたことでした。
外群を判定対象から明示的に除外することで解消しています。

このとき **系統樹自体は一切変更していません**。修正したのはラベル付けのロジックだけで、
その判断根拠も [`docs/METHODS.md`](docs/METHODS.md) に記録しています。
解析結果を都合よく書き換えたのではないことを、後から辿れる形で残すためです。

### 3. 推定と断定を出力上で区別した

[`assign_groups.py`](scripts/assign_groups.py) の出力には `method` 列があり、
クレードから確定した判定 (`clade`) と、距離による推定 (`fallback`) を区別します。

```
sequence_id	group	method
PoryHK5	V	clade
XXXXXXXX	V	fallback     ← 確度が低いことが出力に残る
```

どちらも「Group V」という同じ結論ですが、根拠の強さが違います。
これを混ぜると、後で見返したときに区別がつかなくなります。

さらに、**アンカー配列自身が正しく再現されるかを毎回自己検証**し、
1 件でも合わなければ異常終了します。実際の系統樹 (104 配列) では
先行研究で分類が確定している 10 配列すべてを再現しました。

### 4. 先行研究との対応づけを機械的に行った

系統樹のラベルと論文の遺伝子名を突き合わせる必要がありましたが、
配列のタイトルは自動命名 (`hypothetical protein` 等) で当てにならず、
一見すると別遺伝子に見えるものもありました。

そこで NCBI E-utilities から **locus_tag** を取得し、論文の Table と
機械的に照合しました。目視の印象ではなく共通キーで対応づけたことで、
10 件すべての 1 対 1 対応を確定できました。

> この過程で、配列のタイトルから「先行研究の遺伝子とは別物では」と
> 一度誤った判断をしています。locus_tag で照合し直して誤りと分かりました。
> 名前ではなく識別子で照合すべきという教訓です。

### 5. 配列データを同梱せず、取得スクリプトを置いた

`example/` にあるのは**アクセッション番号のリストだけ**です。
配列本体は [`fetch_ncbi_sequences.py`](scripts/fetch_ncbi_sequences.py) で
利用者の手元に取得します。データの出所が NCBI であることが明示され、
リポジトリ内の配列が古くなる問題も避けられます。

---

## 動作確認

`example/` のデータで、NCBI への問い合わせを含めて実際に動作します。

```console
$ python scripts/fetch_ncbi_sequences.py \
      --acc-list example/accessions_moryzae_hk.txt \
      --metadata --out example/moryzae_hk_metadata.tsv
10 件のアクセッション番号を読み込みました
  取得中 1/1 (10 件) ...
10 件のメタデータを出力: example\moryzae_hk_metadata.tsv

$ python scripts/validate_fasta.py --fasta example/moryzae_hk.fasta --verify-ncbi
配列数: 10 (ユニークなヘッダ: 10)
[OK] 重複ヘッダなし
[OK] 空配列なし
     配列長: 中央値 1408 / 最短 1201 / 最長 2580
[OK] 極端な配列長の外れ値なし
[OK] 不正な文字なし
NCBI 照合中 (10 件) ...
[OK] 10 件すべて NCBI と配列長が一致
検証を通過しました。
```

---

## 使い方

### 必要なもの

| ソフトウェア | 用途 |
|---|---|
| Python 3.10 以降 | スクリプト全般 |
| [Biopython](https://biopython.org/) | `assign_groups.py` の系統樹操作 |
| [HMMER](http://hmmer.org/) 3.x | ドメイン検索 (`hmmsearch`) |
| [MAFFT](https://mafft.cbrc.jp/alignment/software/) | アライメント |
| [trimAl](http://trimal.cgenomics.org/) | アライメントのトリミング |
| [IQ-TREE](http://www.iqtree.org/) 2.x | 最尤法による系統樹推定 |

Python 側の依存は Biopython のみです。

```console
$ pip install biopython
```

HMMER・MAFFT・trimAl・IQ-TREE は Conda での導入が簡単です。

```console
$ conda install -c bioconda hmmer mafft trimal iqtree
```

### 実行手順

**1. HK 候補の絞り込み**

Pfam の HMM プロファイルで 3 ドメインを検索し、2/3 ルールで判定します。

```console
$ for dom in HisKA HATPase_c REC; do
      hmmsearch --tblout ${dom}.tblout ${dom}.hmm proteome.fasta
  done

$ python scripts/filter_hk_candidates.py \
      --hiska HisKA.tblout --hatpase HATPase_c.tblout --rec REC.tblout \
      --evalue 1e-5 --out hk_candidates.tsv
```

**2. 参照配列の取得と検証**

```console
$ python scripts/fetch_ncbi_sequences.py \
      --acc-list example/accessions_moryzae_hk.txt --out reference.fasta

$ cat hk_candidates.fasta reference.fasta > merged.fasta

$ python scripts/validate_fasta.py --fasta merged.fasta --verify-ncbi
```

検証を通してから次へ進みます。ここを飛ばすと、前述の「エラーが出ない失敗」を
系統樹まで持ち込むことになります。

**3. アライメントと系統樹推定**

```console
$ mafft --auto merged.fasta > aligned.fasta
$ trimal -in aligned.fasta -out trimmed.fasta -gt 0.2
$ iqtree -s trimmed.fasta -m LG+G4 -bb 1000 -nt AUTO -o EcoliEnvZ
```

> `trimal` は `-automated1` ではなく `-gt 0.2` を使っています。
> `-automated1` では全長にわたってギャップが多い配列が
> `Removing sequence composed only by gaps` として**削除される**ことがあり、
> 実際に 1 配列を失いました。詳細は [`docs/METHODS.md`](docs/METHODS.md) に記載しています。

**4. グループの割り当て**

```console
$ python scripts/assign_groups.py \
      --tree trimmed.fasta.contree --anchors anchors.tsv \
      --outgroup EcoliEnvZ --out groups.tsv
```

`anchors.tsv` は先行研究で分類が確定している配列の一覧です。

```
sequence_id	group
PoryHK5	V
PoryHK6	VI
```

**5. ドメイン構成図の生成**

```console
$ python scripts/smart_batch_query.py --fasta candidates.fasta --outdir smart_results/

$ python scripts/draw_domain_map.py \
      --smart-dir smart_results/ --fasta candidates.fasta \
      --groups groups.tsv --out domain_map.svg
```

---

## ファイル構成

```
scripts/
  filter_hk_candidates.py    HMMER の結果を統合し 2/3 ルールで HK 候補を判定
  fetch_ncbi_sequences.py    NCBI E-utilities で配列・locus_tag・配列長を取得
  validate_fasta.py          系統解析前の入力検証 (NCBI 照合を含む)
  assign_groups.py           系統樹上でグループを伝播 + 自己検証
  smart_batch_query.py       SMART へのバッチ問い合わせ (セッション管理・再開対応)
  draw_domain_map.py         ドメイン構成図 (SVG) を生成

docs/
  METHODS.md                 使用ツール・パラメータ・判断根拠・既知の限界

example/
  accessions_moryzae_hk.txt  動作確認用のアクセッション番号リスト
```

---

## 技術的な補足

**なぜ SMART と HMMER を併用するのか**

同じ配列でも、Pfam ベースの HMMER と SMART では検出されるドメインが一致しない
ことがあります。どちらかが誤りというより、プロファイルと閾値の設計が違うためです。
先行研究の基準と比較する場合は、**その研究が使ったツールに揃える**必要があります。
このパイプラインでは HMMER を一次スクリーニング、SMART を最終的なドメイン判定に
使い分けています。

**SMART のバッチ問い合わせについて**

SMART は Web UI しか公開しておらず、公式には Perl 製のバッチスクリプトが
配布されています。これを Python で再実装しました。解析エンドポイントに直接
POST してもモード選択ページが返るだけで、事前に `change_mode.cgi` へアクセスして
セッション Cookie を確立する必要があります。解析自体も非同期なので、
`jobId` をポーリングして完了を待ちます。公共サーバなので、取得済みの配列は
スキップし、リクエスト間に待機を入れています。

---

## 既知の限界

- グループ判定は**アンカー配列の正しさに依存**します。アンカーの分類が誤っていれば
  伝播先もすべて誤ります。アンカーは必ず査読済みの文献から取ってください
- `fallback` 判定は距離ベースの推定であり、確度は `clade` 判定より劣ります
- ブートストラップ値の低いノードでは、クレードの境界自体が信頼できません。
  判定結果は系統樹の支持率と併せて解釈する必要があります
- `validate_fasta.py` の NCBI 照合は、ヘッダがアクセッション番号形式の配列のみが対象です

---

## ライセンス

MIT License ([LICENSE](LICENSE))
