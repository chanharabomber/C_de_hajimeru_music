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
  MONO_PCM pcm00, pcm0, pcm1;
  int n, k, N, offset, frame, number_of_frame;
  double threshold, *w, *x_real, *x_imag, *A, *T, *y_real, *y_imag;
 
  double *x2_real,*x2_imag,*A2,*T2,*y2_real,*y2_imag;
  double *y3_real,*y3_imag;
  
  FILE *fp;
  char tring = k;
  char *fname = "result.dat";
  fp = fopen( fname, "w" );

  mono_wave_read(&pcm00,"grouwl.wav"); 
  mono_wave_read(&pcm0,"sample.wav");
  //mono_wave_read(&pcm0,"hajimemasite.wav");


  N = 4096;
  N = 16384;
  N = 32768;
  N = 44000;
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
 
  y3_real = calloc(N, sizeof(double)); /* メモリの確保 */
  y3_imag = calloc(N, sizeof(double)); /* メモリの確保 */

  Hanning_window(w,N);
  number_of_frame = (pcm0.length - N / 2) / (N / 2);
  
  pcm1.fs = pcm0.fs;
  pcm1.bits = pcm0.bits; /* 量子化精度 */
  pcm1.length = pcm0.length; /* 音データの長さ */
  pcm1.s = calloc(pcm1.length, sizeof(double)); /* メモリの確保 */

    for (n = 0; n < N; n++)
    {
      //x_real[n] = pcm00.s[n+33000];
      x_real[n] = pcm00.s[n];
      x_imag[n] = 0.0;
    }
    FFT(x_real,x_imag,N);
    
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
      IFFT(y_real,y_imag,N);
      
      for(k=0;k<N;k++) {
        pcm1.s[k] = y_real[k];
      }

  mono_wave_write(&pcm1, "ffttry.wav"); /* WAVEファイルにモノラルの音データを出力する */
  
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

  fclose(fp);
  free(x2_real);
  free(x2_imag);
  free(A2);
  free(T2);
  free(y2_real);
  free(y2_imag);

  free(y3_real);
  free(y3_imag);
  return 0;
}
