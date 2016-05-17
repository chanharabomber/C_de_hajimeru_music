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
       

  mono_wave_read(&pcm0, "grouwlvoice.wav"); /* WAVEファイルからモノラルの音データを入力する */
  //mono_wave_read(&pcm0, "sample12.wav");
  mono_wave_read(&pcm00,"hajimemasite.wav");
  //mono_wave_read(&pcm00, "sample12.wav");

  MONO_PCM pcmtemp;
  double t, pitch;
  int J,m;
  pcmtemp.fs = pcm0.fs; /* 標本化周波数の変更 */
  /* サンプリング周波数のピッチではなく、声の基本周波数のピッチ */
  pitch = 1.0;
  //pitch = 250.0 / 70.0; // 決め打ち
  //pitch = 70.0 / pcm00.fs;
  pcmtemp.bits = pcm0.bits; /* 量子化精度 */
  pcmtemp.length = (int)(pcm0.length / pitch); /* 音データの長さ */
  pcmtemp.s = calloc(pcmtemp.length, sizeof(double)); /* メモリの確保 */

  pcm1.fs = pcm0.fs;
  pcm1.bits = pcm0.bits; /* 量子化精度 */
  pcm1.length = (int)(pcm0.length / pitch); /* 音データの長さ */
  pcm1.s = calloc(pcmtemp.length, sizeof(double)); /* メモリの確保 */

  printf("#pcm00fs = %d,pitch= %f length=%d\n",pcm00.fs,pitch,pcm00.length);
  printf("#pcm0fs = %d,pitch = %f length=%d\n",pcm0.fs,pitch,pcm0.length);
 

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
  N = 2048;
  N = 4096;
  number_of_frame = (fmin(pcmtemp.length,pcm00.length) - N / 2) / (N / 2);
  printf("length = %d , %d framenum = %d\n",pcmtemp.length,pcm00.length,number_of_frame);

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
  


  Hanning_window(w, N); /* ハニング窓 */

  for (frame = 0; frame < number_of_frame;frame++)
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
      //printf("%d %f %f \n",k,A[k],A2[k]);
      printf("%f %f \n",A[k],A2[k]);
      //if( fp == NULL ){
      //  printf( "%sファイルが開けません¥n", fname );
      //  return -1;
      //}
      //fprintf(fp,"%d %f %f¥n",k,A[k],A2[k]);
  
    }
    printf("\n");//gnuplotyo
    
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
      //printf("%f %f %f %f ¥n",A[k],T[k],A2[k],T2[k]);
      y2_real[k] = A2[k] * cos(T2[k]);
      y2_imag[k] = A2[k] * sin(T2[k]);
    }
   
    //2の方がリサンプルじゃない方
    int mi,di;
    for(k=0;k<N;k++)
    {
       mi = mapping(k,1.0/*r=ratio of pitch*/);
       di = d_i_filter(50.0/*F0 of input*/,(double)k,(double)mi,(double)(N/2),pcm1.fs); 
       //printf("k=%d mi=%d di=%d¥n",k,mi,di);
       //printf("A[k+di]= %f A[k]= %f A2[mi]= %f T[k]= %f T2[mi]= %f¥n",A[k+di],A[k],A2[mi],T[k],T2[mi]);
       //y3_real[k] = A[k+di]*(A[k]/A2[mi])*cos(T[k]/T2[mi]);
       //y3_imag[k] = A[k+di]*(A[k]/A2[mi])*sin(T[k]/T2[mi]);
       y3_real[k] = A2[k]*cos(T2[k]);
       y3_imag[k] = A2[k]*sin(T2[k]);
    }
    
    
    IFFT(y_real, y_imag, N);
    //IFFT(y2_real, y2_imag, N);
    IFFT(y3_real, y3_imag, N);
    /* 加工結果の連結 */
    for (n = 0; n < N; n++)
    {
        pcm1.s[offset + n ] += y3_real[n];
        pcm1.s[offset + n ] += y_real[n];
      //pcm1.s[offset + n] += y2_real[n];
    }
  }
  printf("#len= %d\n",pcm1.length);
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

  free(y3_real);
  free(y3_imag);
  return 0;
}
