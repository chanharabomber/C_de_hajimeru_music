#include <stdio.h>

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

