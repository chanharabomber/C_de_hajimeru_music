#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "wave.h"
#include "window_function.h"
#include "fft.h"

#include "sinc.h"
#include "morph_function.h"

int main(void)
{
  MONO_PCM  pcm0, pcm1;
  int n, k, N, offset, frame, number_of_frame;
  double *w,*x2_real,*x2_imag,*A2,*T2,*y2_real,*y2_imag;
  mono_wave_read(&pcm0,"sample1.wav");
  N = 4096;
  N = 1024;
  N = 2048;
  N = 512;
 // N= 8192;
  w = calloc(N, sizeof(double)); /* メモリの確保 */
  x2_real = calloc(N, sizeof(double)); /* メモリの確保 */
  x2_imag = calloc(N, sizeof(double)); /* メモリの確保 */
  A2 = calloc(N, sizeof(double)); /* メモリの確保 */
  T2 = calloc(N, sizeof(double)); /* メモリの確保 */
  y2_real = calloc(N, sizeof(double)); /* メモリの確保 */
  y2_imag = calloc(N, sizeof(double)); /* メモリの確保 */
  
  Hanning_window(w,N);
  
  number_of_frame = (pcm0.length - N / 2) / (N / 2);
  

  pcm1.fs = pcm0.fs;
  pcm1.bits = pcm0.bits; /* 量子化精度 */
  pcm1.length = pcm0.length; /* 音データの長さ */
  pcm1.s = calloc(pcm1.length, sizeof(double)); /* メモリの確保 */
 
for (frame = 0; frame < number_of_frame;frame++) //
{ //
    offset = N / 2 * frame; //
    for (n = 0; n < N; n++)
    {
      x2_real[n] = pcm0.s[n+offset] * w[n]; //
      //x2_real[n] = pcm0.s[n+offset];
      x2_imag[n] = 0.0;
    }
    FFT(x2_real,x2_imag,N);
    /* 振幅スペクトルと位相スペクトル */
    for (k = 0; k < N; k++)
    {
      A2[k] = sqrt(x2_real[k] * x2_real[k] + x2_imag[k] * x2_imag[k]);
      if (x2_imag[k] != 0.0 && x2_real[k] != 0.0)
      {
        T2[k] = atan2(x2_imag[k], x2_real[k]);
      }
    }
    for (k = 0; k < N; k++)
    {
      y2_real[k] = A2[k] * cos(T2[k]);
      y2_imag[k] = A2[k] * sin(T2[k]);
    }
    IFFT(y2_real,y2_imag,N); 
    
    for(k=0;k<N;k++) {
        pcm1.s[k+offset] += y2_real[k];//
    }
}//framesyorinotoki


  mono_wave_write(&pcm1, "ex10_2try.wav"); /* WAVEファイルにモノラルの音データを出力する */
  
  free(pcm0.s); /* メモリの解放 */
  free(pcm1.s); /* メモリの解放 */
  free(w); /* メモリの解放 */
  free(x2_real);
  free(x2_imag);
  free(A2);
  free(T2);
  free(y2_real);
  free(y2_imag);
  
  return 0;
}
