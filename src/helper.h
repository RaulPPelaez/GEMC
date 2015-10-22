#ifndef HELPER_H
#define HELPER_H

#include<sstream>   
#include<vector>
#include<cmath>
using namespace std;
#define fori(x,y) for(int i=x; i<y;i++)
#define forj(x,y) for(int j=x; j<y;j++)
#define fork(x,y) for(int k=x; k<y;k++)
#define forl(x,y) for(int l=x; l<y;l++)

#define RANDESP ((double)rand()/(double)RAND_MAX) 
#define PI2 3.1415*2.0

#define SQUARE(x) x*x

double pow(double x, int n);
double square_mod(const double *a);
void get_inbounds(double *a);

class RRand{
 public:
  RRand(){}

  void init(){
    s[0] = rand();
    s[1] = rand();
  }
  void init(uint64_t a, uint64_t b){
    s[0] = a;
    s[1] = b;
  }
  uint64_t next(){
    uint64_t x = s[0];
    uint64_t const y = s[1];
    s[0] = y;
    x ^= x << 23; // a
    x ^= x >> 17; // b
    x ^= y ^ (y >> 26); // c
    s[1] = x;
    return x + y;
  }
  double nextd(){
    return next()/(double)max;
  }
  uint64_t max = 18446744073709551615ULL;

  template <typename T> T nextd(T a, T b){
    return nextd()*(a-b)+b;
  }

  bool flip_coin(){
    double temp = nextd();
    while(temp==0.5) temp = nextd();

    return temp >0.5f;
  }
 private:
  uint64_t s[2];
  
};

#endif
