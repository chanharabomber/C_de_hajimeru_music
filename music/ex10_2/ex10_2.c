#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "wave.h"
#include "window_function.h"
#include "fft.h"

int main(void)
{
  MONO_PCM pcm0, pcm1;
  int n, k, N, offset, frame, number_of_frame;
  double threshold, *w, *x_real, *x_imag, *A, *T, *y_real, *y_imag;
  
  FILE *fp;
  char *fname = "data.txt";
  fp = fopen( fname, "w" );
  if( fp == NULL ){
    printf( "%sファイルが開けません\n", fname );
    return -1;
  }


  mono_wave_read(&pcm0, "sample10.wav"); /* WAVEファイルからモノラルの音データを入力する */
  //mono_wave_read(&pcm0, "../guitar_A5.wav");
  //mono_wave_read(&pcm0, "../resampling/ex12_2.wav");
  pcm1.fs = pcm0.fs; /* 標本化周波数 */
  pcm1.bits = pcm0.bits; /* 量子化精度 */
  pcm1.length = pcm0.length; /* 音データの長さ */
  pcm1.s = calloc(pcm1.length, sizeof(double)); /* メモリの確保 */
  
  threshold = 0.5; /* しきい値 */
  
  N = 1024; /* DFTのサイズ */
  
  number_of_frame = (pcm0.length - N / 2) / (N / 2); /* フレームの数 */
  
  w = calloc(N, sizeof(double)); /* メモリの確保 */
  x_real = calloc(N, sizeof(double)); /* メモリの確保 */
  x_imag = calloc(N, sizeof(double)); /* メモリの確保 */
  A = calloc(N, sizeof(double)); /* メモリの確保 */
  T = calloc(N, sizeof(double)); /* メモリの確保 */
  y_real = calloc(N, sizeof(double)); /* メモリの確保 */
  y_imag = calloc(N, sizeof(double)); /* メモリの確保 */
  
  Hanning_window(w, N); /* ハニング窓 */
  
  for (frame = 0; frame < number_of_frame; frame++)
  {
    offset = N / 2 * frame;
    
    /* x(n)のFFT */
    for (n = 0; n < N; n++)
    {
      x_real[n] = pcm0.s[offset + n] * w[n];
      x_imag[n] = 0.0;
    }
    FFT(x_real, x_imag, N);
    
    /* 振幅スペクトルと位相スペクトル */
    for (k = 0; k < N; k++)
    {
      A[k] = sqrt(x_real[k] * x_real[k] + x_imag[k] * x_imag[k]);
      if (x_imag[k] != 0.0 && x_real[k] != 0.0)
      {
        T[k] = atan2(x_imag[k], x_real[k]);
      }
      
      printf("%d %f %f\n",k,A[k],T[k]);
      
    }
    /* スペクトルサブトラクション */
    for (k = 0; k < N; k++)
    {
      A[k] -= threshold;
      if (A[k] < 0.0)
      {
        A[k] = 0.0;
     }
    }
    
    for (k = 0; k < N; k++)
    {
      y_real[k] = A[k] * cos(T[k]);
      y_imag[k] = A[k] * sin(T[k]);
    }
    IFFT(y_real, y_imag, N);
    
    /* 加工結果の連結 */
    for (n = 0; n < N; n++)
    {
      pcm1.s[offset + n] += y_real[n];
    }
  }
  
  mono_wave_write(&pcm1, "ex10_2.wav"); /* WAVEファイルにモノラルの音データを出力する */
  
  free(pcm0.s); /* メモリの解放 */
  free(pcm1.s); /* メモリの解放 */
  free(w); /* メモリの解放 */
  free(x_real); /* メモリの解放 */
  free(x_imag); /* メモリの解放 */
  free(A); /* メモリの解放 */
  free(T); /* メモリの解放 */
  free(y_real); /* メモリの解放 */
  free(y_imag); /* メモリの解放 */
  fclose(fp);
  return 0;
}
