#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "wave.h"

int main(void)
{
  MONO_PCM pcm0, pcm1;
  int n;
  double a, depth, rate;
  
  mono_wave_read(&pcm0, "sample12.wav"); /* WAVE�t�@�C�����烂�m�����̉��f�[�^����͂��� */
  
  pcm1.fs = pcm0.fs; /* �W�{�����g�� */
  pcm1.bits = pcm0.bits; /* �ʎq�����x */
  pcm1.length = pcm0.length; /* ���f�[�^�̒��� */
  pcm1.s = calloc(pcm1.length, sizeof(double)); /* �������̊m�� */
  
  depth = 1.0;
  rate = 150.0; /* 150Hz */
  
  /* �{�C�X�`�F���W�� */
  for (n = 0; n < pcm1.length; n++)
  {
    a = depth * sin(2.0 * M_PI * rate * n / pcm1.fs);
    pcm1.s[n] = a * pcm0.s[n];
  }
  
  mono_wave_write(&pcm1, "ex12_5.wav"); /* WAVE�t�@�C���Ƀ��m�����̉��f�[�^���o�͂��� */
  
  free(pcm0.s); /* �������̉�� */
  free(pcm1.s); /* �������̉�� */
  
  return 0;
}
