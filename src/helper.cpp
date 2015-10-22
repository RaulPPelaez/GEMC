#include "helper.h"

double pow(double x, int n){ 
  double res = x;
  for(int i=0;i<n-1;i++)res*=x; 
  return res;
}
double square_mod(const double *a){
  double res = 0.0;
  double temp;
  fori(0,3){
    temp = a[i];
    res+= SQUARE(temp);
  }
  return abs(res);
}

//From -0.5*bound to 0.5*bound
void get_inbounds(double *a){
  fori(0, 3) a[i]-=int( ( (a[i]<0)?-0.5:0.5 ) + a[i]);
}
