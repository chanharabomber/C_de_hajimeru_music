#include <stdio.h>
#include <stdlib.h>
#include "wave.h"

int main(void)
{
  MONO_PCM pcm0, pcm1;
  int n, m, template_size, pmin, pmax, p, q, offset0, offset1;
  double rate, max_of_r, *x, *y, *r;
  
  mono_wave_read(&pcm0, "sample12.wav"); /* WAVEファイルからモノラルの音データを入力する */
  
  rate = 0.66; /* rateの範囲は0.5<=rate<1.0 */
  
  pcm1.fs = pcm0.fs; /* 標本化周波数 */
  pcm1.bits = pcm0.bits; /* 量子化精度 */
  pcm1.length = (int)(pcm0.length / rate) + 1; /* 音データの長さ */
  pcm1.s = calloc(pcm1.length, sizeof(double)); /* メモリの確保 */
  
  template_size = (int)(pcm1.fs * 0.01); /* 10ms */
  pmin = (int)(pcm1.fs * 0.005); /* 5ms */
  pmax = (int)(pcm1.fs * 0.02); /* 20ms */
  
  x = calloc(template_size, sizeof(double)); /* メモリの確保 */
  y = calloc(template_size, sizeof(double)); /* メモリの確保 */
  r = calloc((pmax + 1), sizeof(double)); /* メモリの確保 */
  
  offset0 = 0;
  offset1 = 0;
  
  while (offset0 + pmax * 2 < pcm0.length)
  {
    for (n = 0; n < template_size; n++)
    {
      x[n] = pcm0.s[offset0 + n]; /* 本来の音データ */
    }
    
    max_of_r = 0.0;
    p = pmin;
    for (m = pmin; m <= pmax; m++)
    {
      for (n = 0; n < template_size; n++)
      {
        y[n] = pcm0.s[offset0 + m + n]; /* mサンプルずらした音データ */
      }
      r[m] = 0.0;
      for (n = 0; n < template_size; n++)
      {
        r[m] += x[n] * y[n]; /* 相関関数 */
      }
      if (r[m] > max_of_r)
      {
        max_of_r = r[m]; /* 相関関数のピーク */
        p = m; /* 音データの基本周期 */
      }
    }
    
    for (n = 0; n < p; n++)
    {
      pcm1.s[offset1 + n] = pcm0.s[offset0 + n];
    }
    for (n = 0; n < p; n++)
    {
      pcm1.s[offset1 + p + n] = pcm0.s[offset0 + p + n] * (p - n) / p; /* 単調減少の重みづけ */
      pcm1.s[offset1 + p + n] += pcm0.s[offset0 + n] * n / p; /* 単調増加の重みづけ */
    }
    
    q = (int)(p * rate / (1.0 - rate) + 0.5);
    for (n = p; n < q; n++)
    {
      if (offset0 + n >= pcm0.length)
      {
        break;
      }
      pcm1.s[offset1 + p + n] = pcm0.s[offset0 + n];
    }
    
    offset0 += q; /* offset0の更新 */
    offset1 += p + q; /* offset1の更新 */
  }
  
  mono_wave_write(&pcm1, "ex11_2.wav"); /* WAVEファイルにモノラルの音データを出力する */
  
  free(pcm0.s); /* メモリの解放 */
  free(pcm1.s); /* メモリの解放 */
  free(x); /* メモリの解放 */
  free(y); /* メモリの解放 */
  free(r); /* メモリの解放 */
  
  return 0;
}
