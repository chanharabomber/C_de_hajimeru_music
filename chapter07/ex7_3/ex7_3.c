#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "wave.h"
#include "iir_filter.h"

int main(void)
{
  MONO_PCM pcm0, pcm1;
  int n, m, d, I, J;
  double fc, Q, depth, rate, a[3], b[3];
  
  mono_wave_read(&pcm0, "sample04.wav"); /* WAVEファイルからモノラルの音データを入力する */
  
  pcm1.fs = pcm0.fs; /* 標本化周波数 */
  pcm1.bits = pcm0.bits; /* 量子化精度 */
  pcm1.length = pcm0.length; /* 音データの長さ */
  pcm1.s = calloc(pcm1.length, sizeof(double)); /* メモリの確保 */
  
  I = 2; /* 遅延器の数 */
  J = 2; /* 遅延器の数 */
  
  depth = 800.0; /* 800Hz */
  rate = 1.0; /* 1Hz */
  
  /* ワウ */
  for (n = 0; n < pcm1.length; n++)
  {
    fc = (1000.0 + depth * sin(2.0 * M_PI * rate * n / pcm1.fs)) / pcm1.fs; /* 中心周波数 */
    Q = 2.0; /* クオリティファクタ */
    
    IIR_resonator(fc, Q, a, b); /* IIRフィルタの設計 */
    
    for (m = 0; m <= J; m++)
    {
      if (n - m >= 0)
      {
        pcm1.s[n] += b[m] * pcm0.s[n - m];
      }
    }
    for (m = 1; m <= I; m++)
    {
      if (n - m >= 0)
      {
        pcm1.s[n] += -a[m] * pcm1.s[n - m];
      }
    }
  }
  
  mono_wave_write(&pcm1, "ex7_3.wav"); /* WAVEファイルにモノラルの音データを出力する */
  
  free(pcm0.s); /* メモリの解放 */
  free(pcm1.s); /* メモリの解放 */
  
  return 0;
}
