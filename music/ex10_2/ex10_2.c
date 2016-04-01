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
    printf( "%s�t�@�C�����J���܂���\n", fname );
    return -1;
  }


  mono_wave_read(&pcm0, "sample10.wav"); /* WAVE�t�@�C�����烂�m�����̉��f�[�^����͂��� */
  //mono_wave_read(&pcm0, "../guitar_A5.wav");
  //mono_wave_read(&pcm0, "../resampling/ex12_2.wav");
  pcm1.fs = pcm0.fs; /* �W�{�����g�� */
  pcm1.bits = pcm0.bits; /* �ʎq�����x */
  pcm1.length = pcm0.length; /* ���f�[�^�̒��� */
  pcm1.s = calloc(pcm1.length, sizeof(double)); /* �������̊m�� */
  
  threshold = 0.5; /* �������l */
  
  N = 1024; /* DFT�̃T�C�Y */
  
  number_of_frame = (pcm0.length - N / 2) / (N / 2); /* �t���[���̐� */
  
  w = calloc(N, sizeof(double)); /* �������̊m�� */
  x_real = calloc(N, sizeof(double)); /* �������̊m�� */
  x_imag = calloc(N, sizeof(double)); /* �������̊m�� */
  A = calloc(N, sizeof(double)); /* �������̊m�� */
  T = calloc(N, sizeof(double)); /* �������̊m�� */
  y_real = calloc(N, sizeof(double)); /* �������̊m�� */
  y_imag = calloc(N, sizeof(double)); /* �������̊m�� */
  
  Hanning_window(w, N); /* �n�j���O�� */
  
  for (frame = 0; frame < number_of_frame; frame++)
  {
    offset = N / 2 * frame;
    
    /* x(n)��FFT */
    for (n = 0; n < N; n++)
    {
      x_real[n] = pcm0.s[offset + n] * w[n];
      x_imag[n] = 0.0;
    }
    FFT(x_real, x_imag, N);
    
    /* �U���X�y�N�g���ƈʑ��X�y�N�g�� */
    for (k = 0; k < N; k++)
    {
      A[k] = sqrt(x_real[k] * x_real[k] + x_imag[k] * x_imag[k]);
      if (x_imag[k] != 0.0 && x_real[k] != 0.0)
      {
        T[k] = atan2(x_imag[k], x_real[k]);
      }
      
      printf("%d %f %f\n",k,A[k],T[k]);
      
    }
    /* �X�y�N�g���T�u�g���N�V���� */
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
    
    /* ���H���ʂ̘A�� */
    for (n = 0; n < N; n++)
    {
      pcm1.s[offset + n] += y_real[n];
    }
  }
  
  mono_wave_write(&pcm1, "ex10_2.wav"); /* WAVE�t�@�C���Ƀ��m�����̉��f�[�^���o�͂��� */
  
  free(pcm0.s); /* �������̉�� */
  free(pcm1.s); /* �������̉�� */
  free(w); /* �������̉�� */
  free(x_real); /* �������̉�� */
  free(x_imag); /* �������̉�� */
  free(A); /* �������̉�� */
  free(T); /* �������̉�� */
  free(y_real); /* �������̉�� */
  free(y_imag); /* �������̉�� */
  fclose(fp);
  return 0;
}
