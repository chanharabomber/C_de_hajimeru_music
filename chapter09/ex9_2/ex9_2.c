#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "wave.h"

int main(void)
{
  MONO_PCM pcm0, pcm1;
  int n, m;
  double d, depth, rate, t, tau, delta;
  
  mono_wave_read(&pcm0, "sample07.wav"); /* WAVE�t�@�C�����烂�m�����̉��f�[�^����͂��� */
  
  pcm1.fs = pcm0.fs; /* �W�{�����g�� */
  pcm1.bits = pcm0.bits; /* �ʎq�����x */
  pcm1.length = pcm0.length; /* ���f�[�^�̒��� */
  pcm1.s = calloc(pcm1.length, sizeof(double)); /* �������̊m�� */
  
  d = pcm1.fs * 0.002; /* 2ms */
  depth = pcm1.fs * 0.002; /* 2ms */
  rate = 0.5; /* 0.5Hz */
  
  /* �t�����W�� */
  for (n = 0; n < pcm1.length; n++)
  {
    pcm1.s[n] = pcm0.s[n];
    
    tau = d + depth * sin(2.0 * M_PI * rate * n / pcm1.fs);
    t = (double)n - tau;
    m = (int)t;
    delta = t - (double)m;
    if (m >= 0 && m + 1 < pcm1.length)
    {
      pcm1.s[n] += delta * pcm0.s[m + 1] + (1.0 - delta) * pcm0.s[m]; 
    }
  }
  
  mono_wave_write(&pcm1, "ex9_2.wav"); /* WAVE�t�@�C���Ƀ��m�����̉��f�[�^���o�͂��� */
  
  free(pcm0.s); /* �������̉�� */
  free(pcm1.s); /* �������̉�� */
  
  return 0;
}
