#include <stdio.h>
#include <stdlib.h>
#include "wave.h"

int main(void)
{
  MONO_PCM pcm0; /* ���m�����̉��f�[�^ */
  STEREO_PCM pcm1; /* �X�e���I�̉��f�[�^ */
  int n, m;
  double d;
  
  mono_wave_read(&pcm0, "sample09.wav"); /* WAVE�t�@�C�����烂�m�����̉��f�[�^����͂��� */
  
  pcm1.fs = pcm0.fs; /* �W�{�����g�� */
  pcm1.bits = pcm0.bits; /* �ʎq�����x */
  pcm1.length = pcm0.length; /* ���f�[�^�̒��� */
  pcm1.sL = calloc(pcm1.length, sizeof(double)); /* �������̊m�� */
  pcm1.sR = calloc(pcm1.length, sizeof(double)); /* �������̊m�� */
  
  d = pcm1.fs * 0.005; /* 5ms */
  
  for (n = 0; n < pcm1.length; n++)
  {
    m = (int)((double)n - d);
    
    if (m >= 0)
    {
      pcm1.sL[n] = pcm0.s[n] + pcm0.s[m];
      pcm1.sR[n] = pcm0.s[n] - pcm0.s[m];
    }
  }
  
  stereo_wave_write(&pcm1, "ex9_4.wav"); /* WAVE�t�@�C���ɃX�e���I�̉��f�[�^���o�͂��� */
  
  free(pcm0.s); /* �������̉�� */
  free(pcm1.sL); /* �������̉�� */
  free(pcm1.sR); /* �������̉�� */
  
  return 0;
}
