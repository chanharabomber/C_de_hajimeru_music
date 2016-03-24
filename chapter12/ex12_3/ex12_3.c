#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "wave.h"
#include "sinc.h"

int main(void)
{
  MONO_PCM pcm0, pcm1, pcm2;
  int n, m, template_size, pmin, pmax, p, q, offset0, offset1, J, offset;
  double rate, max_of_r, t, pitch, *x, *y, *r;
  
  mono_wave_read(&pcm0, "sample12.wav"); /* WAVE�t�@�C�����烂�m�����̉��f�[�^����͂��� */
  
  rate = 1.5; /* rate�͈̔͂�1.0<rate */
  
  pcm1.fs = pcm0.fs; /* �W�{�����g�� */
  pcm1.bits = pcm0.bits; /* �ʎq�����x */
  pcm1.length = (int)(pcm0.length / rate) + 1; /* ���f�[�^�̒��� */
  pcm1.s = calloc(pcm1.length, sizeof(double)); /* �������̊m�� */
  
  template_size = (int)(pcm1.fs * 0.01); /* 10ms */
  pmin = (int)(pcm1.fs * 0.005); /* 5ms */
  pmax = (int)(pcm1.fs * 0.02); /* 20ms */
  
  x = calloc(template_size, sizeof(double)); /* �������̊m�� */
  y = calloc(template_size, sizeof(double)); /* �������̊m�� */
  r = calloc((pmax + 1), sizeof(double)); /* �������̊m�� */
  
  offset0 = 0;
  offset1 = 0;
  
  while (offset0 + pmax * 2 < pcm0.length)
  {
    for (n = 0; n < template_size; n++)
    {
      x[n] = pcm0.s[offset0 + n]; /* �{���̉��f�[�^ */
    }
    
    max_of_r = 0.0;
    p = pmin;
    for (m = pmin; m <= pmax; m++)
    {
      for (n = 0; n < template_size; n++)
      {
        y[n] = pcm0.s[offset0 + m + n]; /* m�T���v�����炵�����f�[�^ */
      }
      r[m] = 0.0;
      for (n = 0; n < template_size; n++)
      {
        r[m] += x[n] * y[n]; /* ���֊֐� */
      }
      if (r[m] > max_of_r)
      {
        max_of_r = r[m]; /* ���֊֐��̃s�[�N */
        p = m; /* ���f�[�^�̊�{���� */
      }
    }
    
    for (n = 0; n < p; n++)
    {
      pcm1.s[offset1 + n] = pcm0.s[offset0 + n] * (p - n) / p; /* �P�������̏d�ݕt�� */
      pcm1.s[offset1 + n] += pcm0.s[offset0 + p + n] * n / p; /* �P�������̏d�ݕt�� */
    }
    
    q = (int)(p / (rate - 1.0) + 0.5);
    for (n = p; n < q; n++)
    {
      if (offset0 + p + n >= pcm0.length)
      {
        break;
      }
      pcm1.s[offset1 + n] = pcm0.s[offset0 + p + n];
    }
    
    offset0 += p + q; /* offset0�̍X�V */
    offset1 += q; /* offset1�̍X�V */
  }
  
  pitch = 0.66; /* ���̍�����0.66�{�ɂ��� */
  
  pcm2.fs = pcm0.fs; /* �W�{�����g���̕ύX */
  pcm2.bits = pcm0.bits; /* �ʎq�����x */
  pcm2.length = pcm0.length; /* ���f�[�^�̒��� */
  pcm2.s = calloc(pcm2.length, sizeof(double)); /* �������̊m�� */
  
  J = 24;
  
  for (n = 0; n < pcm2.length; n++)
  {
    t = pitch * n;
    offset = (int)t;
    for (m = offset - J / 2; m <= offset + J / 2; m++)
    {
      if (m >= 0 && m < pcm1.length)
      {
        pcm2.s[n] += pcm1.s[m] * sinc(M_PI * (t - m));
      }
    }
  }
  
  mono_wave_write(&pcm2, "ex12_3.wav"); /* WAVE�t�@�C���Ƀ��m�����̉��f�[�^���o�͂��� */
  
  free(pcm0.s); /* �������̉�� */
  free(pcm1.s); /* �������̉�� */
  free(pcm2.s); /* �������̉�� */
  free(x); /* �������̉�� */
  free(y); /* �������̉�� */
  free(r); /* �������̉�� */
  
  return 0;
}
