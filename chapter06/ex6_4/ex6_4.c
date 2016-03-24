#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "wave.h"
#include "iir_filter.h"

int main(void)
{
  MONO_PCM pcm0, pcm1;
  int n, m, I, J, L, offset, frame, number_of_frame;
  double fc, Q, a[3], b[3], *x, *y;
  
  mono_wave_read(&pcm0, "sample04.wav"); /* WAVEファイルからモノラルの音データを入力する */
  
  pcm1.fs = pcm0.fs; /* 標本化周波数 */
  pcm1.bits = pcm0.bits; /* 量子化精度 */
  pcm1.length = pcm0.length; /* 音データの長さ */
  pcm1.s = calloc(pcm1.length, sizeof(double)); /* メモリの確保 */
  
  fc = 1000.0 / pcm0.fs; /* 遮断周波数 */
  Q = 1.0 / sqrt(2.0); /* クオリティファクタ */
  I = 2; /* 遅延器の数 */
  J = 2; /* 遅延器の数 */
  
  IIR_LPF(fc, Q, a, b); /* IIRフィルタの設計 */
  
  L = 256; /* フレームの長さ */
  
  number_of_frame = pcm0.length / L; /* フレームの数 */
  
  x = calloc((L + J), sizeof(double)); /* メモリの確保 */
  y = calloc((L + I), sizeof(double)); /* メモリの確保 */
  
  for (frame = 0; frame < number_of_frame; frame++)
  {
    offset = L * frame;
    
    /* 直前のフレームの後半のJサンプルをつけ加える */
    for (n = 0; n < L + J; n++)
    {
      if (offset - J + n < 0)
      {
        x[n] = 0.0;
      }
      else
      {
        x[n] = pcm0.s[offset - J + n];
      }
    }
    
    /* 直前のフィルタリング結果の後半のIサンプルをつけ加える */
    for (n = 0; n < L + I; n++)
    {
      if (offset - I + n < 0)
      {
        y[n] = 0.0;
      }
      else
      {
        y[n] = pcm1.s[offset - I + n];
      }
    }
    
    /* フィルタリング */
    for (n = 0; n < L; n++)
    {
      for (m = 0; m <= J; m++)
      {
        y[I + n] += b[m] * x[J + n - m];
      }
      for (m = 1; m <= I; m++)
      {
        y[I + n] += -a[m] * y[I + n - m];
      }
    }
    
    /* フィルタリング結果の連結 */
    for (n = 0; n < L; n++)
    {
      pcm1.s[offset + n] = y[I + n];
    }
  }
  
  mono_wave_write(&pcm1, "ex6_4.wav"); /* WAVEファイルにモノラルの音データを出力する */
  
  free(pcm0.s); /* メモリの解放 */
  free(pcm1.s); /* メモリの解放 */
  
  return 0;
}
