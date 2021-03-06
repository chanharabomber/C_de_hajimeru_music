#include <stdio.h>
#include <stdlib.h>
#include "wave.h"

int main(void)
{
  MONO_PCM pcm0, pcm1;
  int n;
  
  //mono_wave_read(&pcm0, "a.wav"); /* WAVEファイルからモノラルの音データを入力する */
  mono_wave_read(&pcm0, "/home/akihabara/C_de_hajimeru/chapter02/guitar_A4.wav");
  pcm1.fs = pcm0.fs; /* 標本化周波数 */
  pcm1.bits = pcm0.bits; /* 量子化精度 */
  pcm1.length = pcm0.length; /* 音データの長さ */
  pcm1.s = calloc(pcm1.length, sizeof(double)); /* メモリの確保 */
  for (n = 0; n < pcm1.length; n++)
  {
    pcm1.s[n] = pcm0.s[n]; /* 音データのコピー */
    printf("%d %f\n",n,pcm1.s[n]);
  }
  
  //mono_wave_write(&pcm1, "b.wav"); /* WAVEファイルにモノラルの音データを出力する */
  mono_wave_write(&pcm1,"b_another.wav");
  free(pcm0.s); /* メモリの解放 */
  free(pcm1.s); /* メモリの解放 */
  
  return 0;
}
