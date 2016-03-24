#include <stdio.h>
#include <stdlib.h>
#include "wave.h"

int main(void)
{
  MONO_PCM pcm0, pcm1;
  int n;
  double pitch;
  
  //mono_wave_read(&pcm0, "sample12.wav"); /* WAVE�t�@�C�����烂�m�����̉��f�[�^����͂��� */
  mono_wave_read(&pcm0, "/home/akihabara/C_de_hajimeru/chapter02/guitar_A4.wav");
  pitch = 1.5; /* ���̍�����1.5�{�ɂ��� */
  
  /* ���f�[�^�̃R�s�[ */
  pcm1.fs = (int)(pcm0.fs * pitch); /* �W�{�����g���̕ύX */
  pcm1.bits = pcm0.bits; /* �ʎq�����x */
  pcm1.length = pcm0.length; /* ���f�[�^�̒��� */
  pcm1.s = calloc(pcm1.length, sizeof(double)); /* �������̊m�� */
  for (n = 0; n < pcm1.length; n++)
  {
    pcm1.s[n] = pcm0.s[n]; /* ���f�[�^ */
  }
  
  mono_wave_write(&pcm1, "ex12_1_guitar.wav"); /* WAVE�t�@�C���Ƀ��m�����̉��f�[�^���o�͂��� */
  
  free(pcm0.s); /* �������̉�� */
  free(pcm1.s); /* �������̉�� */
  
  return 0;
}
