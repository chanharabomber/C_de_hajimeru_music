#include <stdio.h>
#include <stdlib.h>
#include "wave.h"

int main(void)
{
  MONO_PCM pcm0, pcm1;
  int n;
  
  //mono_wave_read(&pcm0, "a.wav"); /* WAVE�t�@�C�����烂�m�����̉��f�[�^����͂��� */
  mono_wave_read(&pcm0, "/home/akihabara/C_de_hajimeru/chapter02/guitar_A4.wav");
  pcm1.fs = pcm0.fs; /* �W�{�����g�� */
  pcm1.bits = pcm0.bits; /* �ʎq�����x */
  pcm1.length = pcm0.length; /* ���f�[�^�̒��� */
  pcm1.s = calloc(pcm1.length, sizeof(double)); /* �������̊m�� */
  for (n = 0; n < pcm1.length; n++)
  {
    pcm1.s[n] = pcm0.s[n]; /* ���f�[�^�̃R�s�[ */
    printf("%d %f\n",n,pcm1.s[n]);
  }
  
  //mono_wave_write(&pcm1, "b.wav"); /* WAVE�t�@�C���Ƀ��m�����̉��f�[�^���o�͂��� */
  mono_wave_write(&pcm1,"b_another.wav");
  free(pcm0.s); /* �������̉�� */
  free(pcm1.s); /* �������̉�� */
  
  return 0;
}
