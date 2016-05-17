#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "wave.h"
#include "fft.h"

int main(void)
{
  MONO_PCM pcm0;
  int n, k, N;
  double *x_real, *x_imag;
  
  double *power,*temppower;
  int *index_sort;
  //mono_wave_read(&pcm0, "ex2_1.wav"); /* WAVEファイルからモノラルの音データを入力する */
  mono_wave_read(&pcm0, "../guitar_A4.wav");
  //mono_wave_read(&pcm0, "../guitar_E5.wav");

  printf("length=%d\n",pcm0.length);
  
  //N = 64;
  N = 4096;
  x_real = calloc(N, sizeof(double)); /* メモリの確保 */
  x_imag = calloc(N, sizeof(double)); /* メモリの確保 */
  
  power=calloc(N,sizeof(double));
  temppower=calloc(N,sizeof(double));
  index_sort=calloc(N,sizeof(int));

  for (n = 0; n < N; n++)
  {
    x_real[n] = pcm0.s[n]; /* x(n)の実数部 */
    x_imag[n] = 0.0; /* x(n)の虚数部 */
    power[n] = 0.0;
  }
  
  FFT(x_real, x_imag, N); /* FFTの計算結果はx_realとx_imagに上書きされる */
  
  /* 周波数特性 */
  for (k = 0; k < N; k++)
  {
    // double power[k] = sqrt(x_real[k]*x_real[k] + x_imag[k]*x_imag[k]);
    power[k] = sqrt(x_real[k]*x_real[k] + x_imag[k]*x_imag[k]);
    printf("%d %f %f %f\n", k, x_real[k], x_imag[k],power[k]);

  }
 
  /*
  int maxindex=0;
  double max=0.0;
  for(k=0;k<N;k++) {
    if(max<power[k]) {
        max = power[k];
        maxindex = k;
    }
  }
  printf("%d %f \n",maxindex,max);
  */
  int i,j;
  double temp=0.0;
  for(i=0;i<N;i++) {
    index_sort[i]=i;
    temppower[i]=power[i];
  }
  int tempindex=0;
  for(i=0;i<N-1;i++) {
    for(j=N-1;j>i;j--) {
        if(power[j]<power[j-1]) {
            temp=temppower[j];temppower[j]=temppower[j-1];temppower[j-1]=temp;
            tempindex=j;index_sort[j]=j-1;index_sort[j-1]=tempindex;
        }
    }
  }
  for(i=0;i<N;i++)
    printf("%d %f\n",i,temppower[i]);


  free(pcm0.s); /* メモリの解放 */
  free(x_real); /* メモリの解放 */
  free(x_imag); /* メモリの解放 */
  
  free(power);
  free(index_sort);
  free(temppower);
  return 0;
}
