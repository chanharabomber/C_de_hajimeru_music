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
  mono_wave_read(&pcm0,"sample1.wav");
  //mono_wave_read(&pcm0,"hajimemasite.wav");


  N = 4096;
  N = 8192;

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

/****** growl ******/
    for (n = 0; n < N; n++)
    {
      x_real[n] = pcm00.s[n+33000];
      x_imag[n] = 0.0;
    }
    FFT(x_real,x_imag,N);
   /* 
    for (k = 0; k < N; k++)
    {
      A[k] = sqrt(x_real[k] * x_real[k] + x_imag[k] * x_imag[k]);
      printf("%d %f\n",k,A[k]);
    }*/
  
  double pitch = 430.0 / 1250.0;
  pcm1.fs = pcm0.fs;
  pcm1.bits = pcm0.bits; /* 量子化精度 */
  pcm1.length = pcm0.length; /* 音データの長さ */
  pcm1.s = calloc(pcm1.length, sizeof(double)); /* メモリの確保 */
 
  printf("#%d %d %d\n",pcm1.fs,pcm1.bits,pcm1.length);
  printf("#%d %d %d\n",pcm0.fs,pcm0.bits,pcm0.length);
  printf("#numberofframe= %d\n",number_of_frame); 

/****** hajimemasite *******/
//for(frame=0;frame<pcm1.length-N;frame=frame+N) {
for (frame = 0; frame < number_of_frame;frame++) //
{ //
    offset = N / 2 * frame; //

    for (n = 0; n < N; n++)
    {
      x2_real[n] = pcm0.s[n+offset] * w[n]; //
//      x2_real[n] = pcm0.s[n+28000];
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
    /* スペクトルサブトラクション */
    /*for (k = 0; k < N; k++)
    {
      A2[k] -= threshold;
      if (A2[k] < 0.0)
      {
        A2[k] = 0.0;
     }
    }*/
    
    for (k = 0; k < N; k++)
    {
      y2_real[k] = A2[k] * cos(T2[k]);
      y2_imag[k] = A2[k] * sin(T2[k]);
    }

    IFFT(y2_real,y2_imag,N); 

/** growl sample**/
    for (k = 0; k < N; k++)
    {
      A[k] = sqrt(x_real[k] * x_real[k] + x_imag[k] * x_imag[k]);
      if (x_imag[k] != 0.0 && x_real[k] != 0.0)
      {
        T[k] = atan2(x_imag[k], x_real[k]);
      }
    }
    // スペクトルサブトラクション
    /*for (k = 0; k < N; k++)
    {
      A[k] -= threshold;
      if (A[k] < 0.0)
      {
        A[k] = 0.0;
     }
    }*/

/*
  //mapping
   for(k=0;k<N;k++) {
        //printf("%d %d\n",k,mapping((double)k,0.33));
        if(mapping((double)k,0.33) < N ) {
            A[k] = A[mapping((double)k,0.33)];
            T[k] = T[mapping((double)k,0.33)];
        } else {
            A[k] = 0.0;
            T[k] = 0.0;
        }
   }

    int i,d;
    for(k=0;k<N;k++) {
       if(k>=0 && k < 47)
        i = 31;
       else if(k>=47 && k < 78)
        i = 62;
       else if(k>=78 && k < 109)
        i = 93;
       else if(k>=109 && k < 140)
        i = 124;
       else if(k>=140 && k < 171)
        i = 155;
       else if(k>=171 && k < 204)
        i = 186;
        else if(k>=204 && k < 235)
        i = 217;
        else if(k>=235 && k < 266)
        i = 248;
         else if(k>=266 && k < 297)
        i = 271;
        else
         i = 302;
        
*/       
//       d = d_i_filter((double)300.0/*p_v*/,(double)i /*i*/,(double)mapping((double)i,0.33) /*m_i*/,(double)N /*L*/,(double)44100.0 /*f_s*/);
//       y_real[k] = A[k+d]*(A2[i]/A[mapping((double)i,0.33)]) * cos(T[k+d]*(T2[i]/T[mapping((double)i,0.33)]));
//       y_imag[k] = A[k+d]*(A2[i]/A[mapping((double)i,0.33)]) * sin(T[k+d]*(T2[i]/T[mapping((double)i,0.33)]));

      // printf("%d %f %f\n",k, A[k], sqrt(y_real[k] * y_real[k] + y_imag[k] * y_imag[k]));

//    }
    
    for (k = 0; k < N; k++)
    {
      y_real[k] = A[k] * cos(T[k]);
      y_imag[k] = A[k] * sin(T[k]);
    }

   for(k=0;k<N;k++)
    printf("%d %f\n",k,A2[k]);

    IFFT(y_real, y_imag, N);
 
  for(k=0;k<N;k++) {
  //for(k=0;k<15*N;k++) {
      //pcm1.s[k] += 2.0*y_real[k%N];
      //pcm1.s[k] += 0.1*y2_real[k%N];
//      pcm1.s[k+frame] += 2.0*y_real[k%N];
//      pcm1.s[k+frame] += 0.1*y2_real[k%N];
//      printf("%d %f\n",k,pcm1.s[k]);
      //pcm1.s[k] = y_real[k];
    //pcm1.s[k+offset] += y_real[k];//
    pcm1.s[k+offset] += y2_real[k];//
  }


}//framesyorinotoki
  

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
