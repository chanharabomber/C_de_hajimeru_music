#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "wave.h"

int main(void)
{
  MONO_PCM pcm0, pcm1;
  int n, m, i, repeat;
  double a, d;
  
  mono_wave_read(&pcm0, "sample01.wav"); /* WAVEファイルからモノラルの音データを入力する */
  
  pcm1.fs = pcm0.fs; /* 標本化周波数 */
  pcm1.bits = pcm0.bits; /* 量子化精度 */
  pcm1.length = pcm0.length; /* 音データの長さ */
  pcm1.s = calloc(pcm1.length, sizeof(double)); /* メモリの確保 */
  
  a = 0.5; /* 減衰率 */
  d = pcm1.fs * 0.05; /* 遅延時間 */
  repeat = 10; /* 繰り返し回数 */
  
  for (n = 0; n < pcm1.length; n++)
  {
    pcm1.s[n] = pcm0.s[n]; /* 現在の時刻の音データ */
    
    for (i = 1; i <= repeat; i++)
    {
      m = (int)((double)n - (double)i * d);
      
      if (m >= 0)
      {
        pcm1.s[n] += pow(a, (double)i) * pcm0.s[m]; /* 過去の音データをミックスする */
      }
    }
  }
  
  mono_wave_write(&pcm1, "ex3_2.wav"); /* WAVEファイルにモノラルの音データを出力する */
  
  free(pcm0.s); /* メモリの解放 */
  free(pcm1.s); /* メモリの解放 */
  
  return 0;
}
