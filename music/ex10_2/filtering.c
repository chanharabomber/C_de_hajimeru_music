#include <stdio.h>
int mapping(double x,double y);

int main() {
  int i,r;
  double fs = 44100.0;
  double L = 44100.0;
  for(r=0;r<100;r++) {
    for(i=0;i<100;i++)
      printf("i=%d r=%d m_i=%d d_i=%d\n",i,r,mapping((double)i,(double)r),d_i_filter(250.0,(double)i,(double)mapping((double)i,(double)r),L,fs) );
  }
  return 0;
}

/*    m[i] = floor{ (i/r) + 0.5 }     */
int mapping(double i,double r) {
  double m_i;
  //avoid 0wari 
  if(r!=0) m_i = (i/r) + 0.5;
  else m_i = i+0.5;
  return (int)( m_i < 0.0 ? m_i-0.9 : m_i );
}

//d[i]
int d_i_filter(double p_v,double i,double m_i,double L,double f_s) {
  int d_i;
  d_i = p_v * (m_i-i)*(L/f_s) + 0.5;
  return  (int)( d_i < 0.0 ? d_i-0.9 : d_i );
}

double g_i(double Av_i,double As_mi) {
  return Av_i / As_mi;
}

double theta_i(double tv_i,double ts_mi) {
  return tv_i / ts_mi;
}

