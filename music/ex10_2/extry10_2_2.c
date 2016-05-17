#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "wave.h"
#include "window_function.h"
#include "fft.h"

#include "sinc.h"

int main(void)
{
  MONO_PCM pcm00, pcm0, pcm1;
  int n, k, N, offset, frame, number_of_frame;
  double threshold, *w, *x_real, *x_imag, *A, *T, *y_real, *y_imag;
 
  double *x2_real,*x2_imag,*A2,*T2,*y2_real,*y2_imag;

  FILE *fp;
  char *fname = "data.txt";
  fp = fopen( fname, "w" );
  if( fp == NULL ){
    printf( "%sファイルが開けません\n", fname );
    return -1;
  }

  //mono_wave_read(&pcm0, "grouwlvoice.wav"); /* WAVEファイルからモノラルの音データを入力する */
  mono_wave_read(&pcm0, "sample12.wav");
  mono_wave_read(&pcm00,"sample10.wav");
  //mono_wave_read(&pcm00, "sample12.wav");

  MONO_PCM pcmtemp;
  double t, pitch;
  int J,m;
  pcmtemp.fs = pcm0.fs; /* 標本化周波数の変更 */
  /* サンプリング周波数のピッチではなく、声の基本周波数のピッチ */
  pitch = 2.0;
  //pitch = 250.0 / 70.0; // 決め打ち
  //pitch = 70.0 / pcm00.fs;
  printf("pcm00fs = %d,pitch= %f\n",pcm00.fs,pitch);
  printf("pcm0fs = %d,pitch = %f\n",pcm0.fs,pitch);
  pcmtemp.bits = pcm0.bits; /* 量子化精度 */
  pcmtemp.length = (int)(pcm0.length / pitch); /* 音データの長さ */
  pcmtemp.s = calloc(pcmtemp.length, sizeof(double)); /* メモリの確保 */

  pcm1.fs = pcm0.fs;
  pcm1.bits = pcm0.bits; /* 量子化精度 */
  pcm1.length = (int)(pcm0.length / pitch); /* 音データの長さ */
  pcm1.s = calloc(pcmtemp.length, sizeof(double)); /* メモリの確保 */
  

  J = 24;
  
  for (n = 0; n < pcmtemp.length; n++)
  {
    t = pitch * n;
    offset = (int)t;
    for (m = offset - J / 2; m <= offset + J / 2; m++)
    {
      if (m >= 0 && m < pcm0.length)
      {
        pcmtemp.s[n] += pcm0.s[m] * sinc(M_PI * (t - m));
      }
    }
  }
  mono_wave_write(&pcmtemp, "test.wav");
 
  threshold = 0.5; /* しきい値 */
  
  N = 1024; /* DFTのサイズ */
  
  number_of_frame = (fmin(pcmtemp.length,pcm00.length) - N / 2) / (N / 2);
  printf("framenum = %d\n",number_of_frame);

  w = calloc(N, sizeof(double)); /* メモリの確保 */
  x_real = calloc(N, sizeof(double)); /* メモリの確保 */
  x_imag = calloc(N, sizeof(double)); /* メモリの確保 */
  A = calloc(N, sizeof(double)); /* メモリの確保 */
  T = calloc(N, sizeof(double)); /* メモリの確保 */
  y_real = calloc(N, sizeof(double)); /* メモリの確保 */
  y_imag = calloc(N, sizeof(double)); /* メモリの確保 */
  
  x2_real = calloc(N, sizeof(double)); /* メモリの確保 */
  x2_imag = calloc(N, sizeof(double)); /* メモリの確保 */
  A2 = calloc(N, sizeof(double)); /* メモリの確保 */
  T2 = calloc(N, sizeof(double)); /* メモリの確保 */
  y2_real = calloc(N, sizeof(double)); /* メモリの確保 */
  y2_imag = calloc(N, sizeof(double)); /* メモリの確保 */
  
  Hanning_window(w, N); /* ハニング窓 */

  for (frame = 0; frame < number_of_frame; frame++)
  {
    offset = N / 2 * frame;

    /* x(n)のFFT */
    for (n = 0; n < N; n++)
    {
      x_real[n] = pcmtemp.s[offset + n] * w[n];
      x_imag[n] = 0.0;
    }
    
    for (n = 0; n < N; n++)
    {
      x2_real[n] = pcm00.s[offset + n] * w[n];
      x2_imag[n] = 0.0;
    }
    FFT(x2_real,x2_imag,N);

    FFT(x_real, x_imag, N);
    
    /* 振幅スペクトルと位相スペクトル */
    for (k = 0; k < N; k++)
    {
      A[k] = sqrt(x_real[k] * x_real[k] + x_imag[k] * x_imag[k]);
      if (x_imag[k] != 0.0 && x_real[k] != 0.0)
      {
        T[k] = atan2(x_imag[k], x_real[k]);
      }
      
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

    /* 振幅スペクトルと位相スペクトル */
    for (k = 0; k < N; k++)
    {
      A2[k] = sqrt(x2_real[k] * x2_real[k] + x2_imag[k] * x2_imag[k]);
      if (x2_imag[k] != 0.0 && x2_real[k] != 0.0)
      {
        T2[k] = atan2(x2_imag[k], x2_real[k]);
      }
    }
    /* スペクトルサブトラクション */
    for (k = 0; k < N; k++)
    {
      A2[k] -= threshold;
      if (A2[k] < 0.0)
      {
        A2[k] = 0.0;
     }
    }
    
    for (k = 0; k < N; k++)
    {
      y2_real[k] = A2[k] * cos(T[k]);
      y2_imag[k] = A2[k] * sin(T[k]);
    }
   
    IFFT(y_real, y_imag, N);
    IFFT(y2_real, y2_imag, N);
    /* 加工結果の連結 */
    for (n = 0; n < N; n++)
    {
      pcm1.s[offset + n ] += y_real[n];
      pcm1.s[offset + n] += y2_real[n];
    }
  }
  
  mono_wave_write(&pcm1, "ex10_2try.wav"); /* WAVEファイルにモノラルの音データを出力する */
  
  free(pcm0.s); /* メモリの解放 */
  free(pcm00.s);
  free(pcm1.s); /* メモリの解放 */
  free(w); /* メモリの解放 */
  free(x_real); /* メモリの解放 */
  free(x_imag); /* メモリの解放 */
  free(A); /* メモリの解放 */
  free(T); /* メモリの解放 */
  free(y_real); /* メモリの解放 */
  free(y_imag); /* メモリの解放 */

  free(pcmtemp.s);
  
  fclose(fp);
  free(x2_real);
  free(x2_imag);
  free(A2);
  free(T2);
  free(y2_real);
  free(y2_imag);
  return 0;
}
