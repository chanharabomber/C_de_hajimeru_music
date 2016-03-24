#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "wave.h"
#include "sinc.h"

int main(void)
{
  MONO_PCM pcm0, pcm1;
  int n, m, J, offset;
  double t, pitch;
  
  mono_wave_read(&pcm0, "sample12.wav"); /* WAVE�t�@�C�����烂�m�����̉��f�[�^����͂��� */
  
  pitch = 1.5; /* ���̍�����1.5�{�ɂ��� */
  
  pcm1.fs = pcm0.fs; /* �W�{�����g���̕ύX */
  pcm1.bits = pcm0.bits; /* �ʎq�����x */
  pcm1.length = (int)(pcm0.length / pitch); /* ���f�[�^�̒��� */
  pcm1.s = calloc(pcm1.length, sizeof(double)); /* �������̊m�� */
  
  J = 24;
  
  for (n = 0; n < pcm1.length; n++)
  {
    t = pitch * n;
    offset = (int)t;
    for (m = offset - J / 2; m <= offset + J / 2; m++)
    {
      if (m >= 0 && m < pcm0.length)
      {
        pcm1.s[n] += pcm0.s[m] * sinc(M_PI * (t - m));
      }
    }
  }
  
  mono_wave_write(&pcm1, "ex12_2.wav"); /* WAVE�t�@�C���Ƀ��m�����̉��f�[�^���o�͂��� */
  
  free(pcm0.s); /* �������̉�� */
  free(pcm1.s); /* �������̉�� */
  
  return 0;
}
