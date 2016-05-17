#include <stdio.h>
int mapping(double x,double y);

int main() {
  int i,r;
  for(r=0;r<100;r++) {
    for(i=0;i<100;i++)
      printf("i=%d j=%d %d\n",i,r,mapping((double)i,(double)r) );
  }
  return 0;
}

/*    m[i] = floor{ (i/r) + 0.5 }     */
int mapping(double i,double r) {
  double temp;
  //avoid 0wari 
  if(r!=0) temp = (i/r) + 0.5;
  else temp = i+0.5;
  return (int)( temp < 0.0 ? temp-0.9 : temp );
}
