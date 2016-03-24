#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "wave.h"

int main(void)
{
  MONO_PCM pcm0, pcm1;
  int n, m;
  double d, depth, rate, t, tau, delta;
  
  mono_wave_read(&pcm0, "sample07.wav"); /* WAVEファイルからモノラルの音データを入力する */
  
  pcm1.fs = pcm0.fs; /* 標本化周波数 */
  pcm1.bits = pcm0.bits; /* 量子化精度 */
  pcm1.length = pcm0.length; /* 音データの長さ */
  pcm1.s = calloc(pcm1.length, sizeof(double)); /* メモリの確保 */
  
  d = pcm1.fs * 0.002; /* 2ms */
  depth = pcm1.fs * 0.002; /* 2ms */
  rate = 0.5; /* 0.5Hz */
  
  /* フランジャ */
  for (n = 0; n < pcm1.length; n++)
  {
    pcm1.s[n] = pcm0.s[n];
    
    tau = d + depth * sin(2.0 * M_PI * rate * n / pcm1.fs);
    t = (double)n - tau;
    m = (int)t;
    delta = t - (double)m;
    if (m >= 0 && m + 1 < pcm1.length)
    {
      pcm1.s[n] += delta * pcm0.s[m + 1] + (1.0 - delta) * pcm0.s[m]; 
    }
  }
  
  mono_wave_write(&pcm1, "ex9_2.wav"); /* WAVEファイルにモノラルの音データを出力する */
  
  free(pcm0.s); /* メモリの解放 */
  free(pcm1.s); /* メモリの解放 */
  
  return 0;
}
