#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "wave.h"
#include "sinc.h"

int main(void)
{
  MONO_PCM pcm0, pcm1;
  int n, m, J, offset;
  double t, pitch;
  
  mono_wave_read(&pcm0, "sample12.wav"); /* WAVEファイルからモノラルの音データを入力する */
  
  pitch = 1.5; /* 音の高さを1.5倍にする */
  
  pcm1.fs = pcm0.fs; /* 標本化周波数の変更 */
  pcm1.bits = pcm0.bits; /* 量子化精度 */
  pcm1.length = (int)(pcm0.length / pitch); /* 音データの長さ */
  pcm1.s = calloc(pcm1.length, sizeof(double)); /* メモリの確保 */
  
  J = 24;
  
  for (n = 0; n < pcm1.length; n++)
  {
    t = pitch * n;
    offset = (int)t;
    for (m = offset - J / 2; m <= offset + J / 2; m++)
    {
      if (m >= 0 && m < pcm0.length)
      {
        pcm1.s[n] += pcm0.s[m] * sinc(M_PI * (t - m));
      }
    }
  }
  
  mono_wave_write(&pcm1, "ex12_2.wav"); /* WAVEファイルにモノラルの音データを出力する */
  
  free(pcm0.s); /* メモリの解放 */
  free(pcm1.s); /* メモリの解放 */
  
  return 0;
}
