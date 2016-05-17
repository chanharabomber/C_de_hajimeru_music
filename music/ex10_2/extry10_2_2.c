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
    printf( "%s�t�@�C�����J���܂���\n", fname );
    return -1;
  }

  //mono_wave_read(&pcm0, "grouwlvoice.wav"); /* WAVE�t�@�C�����烂�m�����̉��f�[�^����͂��� */
  mono_wave_read(&pcm0, "sample12.wav");
  mono_wave_read(&pcm00,"sample10.wav");
  //mono_wave_read(&pcm00, "sample12.wav");

  MONO_PCM pcmtemp;
  double t, pitch;
  int J,m;
  pcmtemp.fs = pcm0.fs; /* �W�{�����g���̕ύX */
  /* �T���v�����O���g���̃s�b�`�ł͂Ȃ��A���̊�{���g���̃s�b�` */
  pitch = 2.0;
  //pitch = 250.0 / 70.0; // ���ߑł�
  //pitch = 70.0 / pcm00.fs;
  printf("pcm00fs = %d,pitch= %f\n",pcm00.fs,pitch);
  printf("pcm0fs = %d,pitch = %f\n",pcm0.fs,pitch);
  pcmtemp.bits = pcm0.bits; /* �ʎq�����x */
  pcmtemp.length = (int)(pcm0.length / pitch); /* ���f�[�^�̒��� */
  pcmtemp.s = calloc(pcmtemp.length, sizeof(double)); /* �������̊m�� */

  pcm1.fs = pcm0.fs;
  pcm1.bits = pcm0.bits; /* �ʎq�����x */
  pcm1.length = (int)(pcm0.length / pitch); /* ���f�[�^�̒��� */
  pcm1.s = calloc(pcmtemp.length, sizeof(double)); /* �������̊m�� */
  

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
 
  threshold = 0.5; /* �������l */
  
  N = 1024; /* DFT�̃T�C�Y */
  
  number_of_frame = (fmin(pcmtemp.length,pcm00.length) - N / 2) / (N / 2);
  printf("framenum = %d\n",number_of_frame);

  w = calloc(N, sizeof(double)); /* �������̊m�� */
  x_real = calloc(N, sizeof(double)); /* �������̊m�� */
  x_imag = calloc(N, sizeof(double)); /* �������̊m�� */
  A = calloc(N, sizeof(double)); /* �������̊m�� */
  T = calloc(N, sizeof(double)); /* �������̊m�� */
  y_real = calloc(N, sizeof(double)); /* �������̊m�� */
  y_imag = calloc(N, sizeof(double)); /* �������̊m�� */
  
  x2_real = calloc(N, sizeof(double)); /* �������̊m�� */
  x2_imag = calloc(N, sizeof(double)); /* �������̊m�� */
  A2 = calloc(N, sizeof(double)); /* �������̊m�� */
  T2 = calloc(N, sizeof(double)); /* �������̊m�� */
  y2_real = calloc(N, sizeof(double)); /* �������̊m�� */
  y2_imag = calloc(N, sizeof(double)); /* �������̊m�� */
  
  Hanning_window(w, N); /* �n�j���O�� */

  for (frame = 0; frame < number_of_frame; frame++)
  {
    offset = N / 2 * frame;

    /* x(n)��FFT */
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
    
    /* �U���X�y�N�g���ƈʑ��X�y�N�g�� */
    for (k = 0; k < N; k++)
    {
      A[k] = sqrt(x_real[k] * x_real[k] + x_imag[k] * x_imag[k]);
      if (x_imag[k] != 0.0 && x_real[k] != 0.0)
      {
        T[k] = atan2(x_imag[k], x_real[k]);
      }
      
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

    /* �U���X�y�N�g���ƈʑ��X�y�N�g�� */
    for (k = 0; k < N; k++)
    {
      A2[k] = sqrt(x2_real[k] * x2_real[k] + x2_imag[k] * x2_imag[k]);
      if (x2_imag[k] != 0.0 && x2_real[k] != 0.0)
      {
        T2[k] = atan2(x2_imag[k], x2_real[k]);
      }
    }
    /* �X�y�N�g���T�u�g���N�V���� */
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
    /* ���H���ʂ̘A�� */
    for (n = 0; n < N; n++)
    {
      pcm1.s[offset + n ] += y_real[n];
      pcm1.s[offset + n] += y2_real[n];
    }
  }
  
  mono_wave_write(&pcm1, "ex10_2try.wav"); /* WAVE�t�@�C���Ƀ��m�����̉��f�[�^���o�͂��� */
  
  free(pcm0.s); /* �������̉�� */
  free(pcm00.s);
  free(pcm1.s); /* �������̉�� */
  free(w); /* �������̉�� */
  free(x_real); /* �������̉�� */
  free(x_imag); /* �������̉�� */
  free(A); /* �������̉�� */
  free(T); /* �������̉�� */
  free(y_real); /* �������̉�� */
  free(y_imag); /* �������̉�� */

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
