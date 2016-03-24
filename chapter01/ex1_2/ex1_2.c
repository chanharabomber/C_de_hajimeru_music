#include <stdio.h>
#include <stdlib.h>
#include "wave.h"

int main(void)
{
  STEREO_PCM pcm0, pcm1;
  int n;
  
  stereo_wave_read(&pcm0, "c.wav"); /* WAVEファイルからステレオの音データを入力する */
  
  pcm1.fs = pcm0.fs; /* 標本化周波数 */
  pcm1.bits = pcm0.bits; /* 量子化精度 */
  pcm1.length = pcm0.length; /* 音データの長さ */
  pcm1.sL = calloc(pcm1.length, sizeof(double)); /* メモリの確保 */
  pcm1.sR = calloc(pcm1.length, sizeof(double)); /* メモリの確保 */
  for (n = 0; n < pcm1.length; n++)
  {
    pcm1.sL[n] = pcm0.sL[n]; /* 音データのコピー */
    pcm1.sR[n] = pcm0.sR[n]; /* 音データのコピー */
  }
  
  stereo_wave_write(&pcm1, "d.wav"); /* WAVEファイルにステレオの音データを出力する */
  
  free(pcm0.sL); /* メモリの解放 */
  free(pcm0.sR); /* メモリの解放 */
  free(pcm1.sL); /* メモリの解放 */
  free(pcm1.sR); /* メモリの解放 */
  
  return 0;
}
