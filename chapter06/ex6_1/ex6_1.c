#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "wave.h"
#include "window_function.h"
#include "sinc.h"
#include "fir_filter.h"

int main(void)
{
  MONO_PCM pcm0, pcm1;
  int n, m, J;
  double fe, delta, *b, *w;
  
  mono_wave_read(&pcm0, "sample04.wav"); /* WAVEファイルからモノラルの音データを入力する */
  
  pcm1.fs = pcm0.fs; /* 標本化周波数 */
  pcm1.bits = pcm0.bits; /* 量子化精度 */
  pcm1.length = pcm0.length; /* 音データの長さ */
  pcm1.s = calloc(pcm1.length, sizeof(double)); /* メモリの確保 */
  
  fe = 1000.0 / pcm0.fs; /* エッジ周波数 */
  delta = 1000.0 / pcm0.fs; /* 遷移帯域幅 */
  
  J = (int)(3.1 / delta + 0.5) - 1; /* 遅延器の数 */
  if (J % 2 == 1)
  {
    J++; /* J+1が奇数になるように調整する */
  }
  
  b = calloc((J + 1), sizeof(double)); /* メモリの確保 */
  w = calloc((J + 1), sizeof(double)); /* メモリの確保 */
  
  Hanning_window(w, (J + 1)); /* ハニング窓 */
  
  FIR_LPF(fe, J, b, w); /* FIRフィルタの設計 */
  
  /* フィルタリング */
  for (n = 0; n < pcm1.length; n++)
  {
    for (m = 0; m <= J; m++)
    {
      if (n - m >= 0)
      {
        pcm1.s[n] += b[m] * pcm0.s[n - m];
      }
    }
  }
  
  mono_wave_write(&pcm1, "ex6_1.wav"); /* WAVEファイルにモノラルの音データを出力する */
  
  free(pcm0.s); /* メモリの解放 */
  free(pcm1.s); /* メモリの解放 */
  free(b); /* メモリの解放 */
  free(w); /* メモリの解放 */
  
  return 0;
}
