#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "wave.h"

int main(void)
{
  MONO_PCM pcm0; /* ���m�����̉��f�[�^ */
  STEREO_PCM pcm1; /* �X�e���I�̉��f�[�^ */
  int n;
  double a, depth, rate;
  
  mono_wave_read(&pcm0, "sample08.wav"); /* WAVE�t�@�C�����烂�m�����̉��f�[�^����͂��� */
  
  pcm1.fs = pcm0.fs; /* �W�{�����g�� */
  pcm1.bits = pcm0.bits; /* �ʎq�����x */
  pcm1.length = pcm0.length; /* ���f�[�^�̒��� */
  pcm1.sL = calloc(pcm1.length, sizeof(double)); /* �������̊m�� */
  pcm1.sR = calloc(pcm1.length, sizeof(double)); /* �������̊m�� */
  
  depth = 1.0;
  rate = 0.2; /* 0.2Hz */
  
  /* �I�[�g�p�� */
  for (n = 0; n < pcm1.length; n++)
  {
    a = 1.0 + depth * sin(2.0 * M_PI * rate * n / pcm1.fs);
    pcm1.sL[n] = a * pcm0.s[n];
    
    a = 1.0 + depth * sin(2.0 * M_PI * rate * n / pcm1.fs + M_PI);
    pcm1.sR[n] = a * pcm0.s[n];
  }
  
  stereo_wave_write(&pcm1, "ex9_3.wav"); /* WAVE�t�@�C���ɃX�e���I�̉��f�[�^���o�͂��� */
  
  free(pcm0.s); /* �������̉�� */
  free(pcm1.sL); /* �������̉�� */
  free(pcm1.sR); /* �������̉�� */
  
  return 0;
}
