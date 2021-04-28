#ifndef REPORT_H
#define REPORT_H

#if defined(__cplusplus)
extern "C" {
#endif

  // double のスカラ値で結果OK確認、確認結果出力を行う。
  // result:        計算した結果
  // reference:     リファレンス
  // percent_error: reference +/- |reference*(percent_error/100)| を誤差範囲として許容する
  // 備考: Fortran インタフェースは report.c を参照
  //       数回呼んでもよいが、配列の全ての要素に対して使うようなことはしない
  void report_validation(double result, double reference, double percent_error);

  // 二乗和を使う場合 (ret = arr[0]^2 + arr[1]^2 + ... + arr[n-1]^2)
  double get_ss_r8(double *arr, int n);
  float get_ss_r4(float *arr, int n);

  // 乱数を使う場合
  double get_rand(double min, double max); // [min,max] の乱数
  void set_seed(unsigned int seed);

#ifdef __cplusplus
}
#endif

#endif
