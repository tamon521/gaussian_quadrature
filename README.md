# gaussian_quadrature
Rustでガウス求積を行うプログラム．

## ガウス求積とは
[ガウス求積](https://ja.wikipedia.org/wiki/ガウス求積)(Gaussian quadrature)とは，
区間[−1, 1] で定義された実数値関数の定積分値を，比較的少ない演算で精度良く求めることができるアルゴリズムである．
このプログラムでは**ガウス・ルジャンドル公式による求積**を用いた．

これらの導出については[ガウス求積法の導出](https://slpr.sakura.ne.jp/qp/guass-legendre-quadrature-derivation/)等を参照．

## 実行方法

```
$ git clone https://github.com/tamon521/gaussian_quadrature.git
$ cd gaussian_quadrature
$ cargo run 
```