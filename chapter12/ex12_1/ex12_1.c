#include <stdio.h>
#include <stdlib.h>
#include "wave.h"

int main(void)
{
  MONO_PCM pcm0, pcm1;
  int n;
  double pitch;
  
  //mono_wave_read(&pcm0, "sample12.wav"); /* WAVEファイルからモノラルの音データを入力する */
  mono_wave_read(&pcm0, "/home/akihabara/C_de_hajimeru/chapter02/guitar_A4.wav");
  pitch = 1.5; /* 音の高さを1.5倍にする */
  
  /* 音データのコピー */
  pcm1.fs = (int)(pcm0.fs * pitch); /* 標本化周波数の変更 */
  pcm1.bits = pcm0.bits; /* 量子化精度 */
  pcm1.length = pcm0.length; /* 音データの長さ */
  pcm1.s = calloc(pcm1.length, sizeof(double)); /* メモリの確保 */
  for (n = 0; n < pcm1.length; n++)
  {
    pcm1.s[n] = pcm0.s[n]; /* 音データ */
  }
  
  mono_wave_write(&pcm1, "ex12_1_guitar.wav"); /* WAVEファイルにモノラルの音データを出力する */
  
  free(pcm0.s); /* メモリの解放 */
  free(pcm1.s); /* メモリの解放 */
  
  return 0;
}
