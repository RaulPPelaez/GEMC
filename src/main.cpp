#include"Experiment.h"

#include<fstream>
#include<iostream>

void init();

void mcmove();
void mcvol();
void mcswap();
void sample();

void write_exp();

Experiment e1, e2;
RRand rng;

double V, vstep, beta;
double mu[2];
uint vtry = 0, vaccept = 0; 

ofstream out("data.out");

uint nsteps = 0;
ofstream out2;

int main(){
  init();
  
  uint npart = 3000;
  uint nvol = 15;
  uint nswap = 6000;
  
  float R;
  uint do_steps = 3e6;
  uint relax_steps = 0;
  
  cout<<"Relaxing..."<<endl;
  fori(0,relax_steps) mcmove();
  fori(0,relax_steps){
    R = rng.next()%(npart + nvol + nswap);
    if(R <= npart) mcmove();
    else if(R<=(npart+nvol)) mcvol();
    else mcswap();  
  }
  cout<<"Averaging..."<<endl;
  while(nsteps<do_steps){
    nsteps++;
    R = rng.next()%(npart + nvol + nswap);
    if(R <= npart) mcmove();
    else if(R<=(npart+nvol)) mcvol();
    else mcswap();  
    if(nsteps%500==0){sample();}
  }

  write_exp();
  return 0;
}
void write_exp(){
  //  cout<<"Write step: "<<nsteps<<endl;
  uint N2 = e2.ps.size()/3;
  out2<<"#L="<<0.5*e2.L<<";"<<endl;
  fori(0, N2){
    double ips[3];
    forj(0,3) ips[j] = e2.ps[3*i+j];
    get_inbounds(&ips[0]);
    forj(0,3)out2<<ips[j]*e2.L<<" ";
    out2<<"0.5 1"<<endl;
  }
  uint N1 = e1.ps.size()/3;
  fori(0, N1){
    double ips[3];
    forj(0,3) ips[j] = e1.ps[3*i+j];
    get_inbounds(&ips[0]);
    forj(0,3)out2<<ips[j]*e1.L+e2.L<<" ";

    out2<<"0.5 2"<<endl;
  }
}

void init(){
  srand(time(NULL));
  rng.init();
  out2.open("p.dat");

  ifstream in("read.in");

  double T, L1, L2, rc;
  uint N1, N2;

  in>>T;
  in>>N1>>N2;
  in>>L1>>L2;
  in>>rc;

  e1.init(T, N1, L1, rc);
  e2.init(T, N2, L2, rc);


  V = L1*L1*L1 + L2*L2*L2;
  vstep = 0.1;
  beta = 1.0/T;
  cout<<"beta "<<beta<<endl;
}

void mcmove(){
  //if(rng.flip_coin()) e1.update();
  e2.update();
}

void mcvol(){


  vtry++;
  double en1o = e1.total_ener();
  double en2o = e2.total_ener();
  
  double v1o = e1.get_V();
  double v2o = V-v1o;

  double N1 = e1.ps.size()/3;
  double N2 = e2.ps.size()/3;

  float lnvn = log(v1o/v2o) + (rng.nextd()-0.5f)*vstep;
  float elnvn = exp(lnvn);

  double v1n = V*elnvn/(1+elnvn);
  double v2n = V-v1n;
  
  double L1o = e1.L;
  double L2o = e2.L;
  e1.L = cbrt(v1n);
  e2.L = cbrt(v2n);

  double en1n = e1.total_ener();
  double en2n = e2.total_ener();

  float arg1 = -beta*( (en1n-en1o) + (N1+1)*log(v1n/v1o)/beta);
  float arg2 = -beta*( (en2n-en2o) + (N2+1)*log(v2n/v2o)/beta);

  if(rng.nextd()>exp(arg1+arg2)){ //rejected
    e1.L = L1o;
    e2.L = L2o;
  }
  else{
    vaccept++;
    e1.make_linked_list();
    e2.make_linked_list();
  }
  
  if(vtry%1000==0){
    float TRUST = 0.5;
    float ratio = (float)vaccept/(float)vtry;
    if(ratio<TRUST){if(vstep>0.0000000001f) vstep *=0.9;}
    else if(ratio>TRUST){vstep *= 1.1;}
    vtry = vaccept = 0;
  }

}

void mcswap(){
  Experiment *ex[2];

  uint in=0, out=1;
  in = rng.flip_coin()?0:1;
  out = !in;
  ex[0] = &e1;
  ex[1] = &e2;
  
  

 
  uint N1 = ex[in]->ps.size()/3;
  uint N2 = ex[out]->ps.size()/3;


  float v1 = ex[in]->get_V();
  float v2 = ex[out]->get_V();

  double nps[3];
  fori(0,3) nps[i] = rng.nextd()-0.5;
  
  double enn = ex[in]->enerNL(&nps[0]);
  
  mu[in] += v1* exp(-beta*enn)/(N1+1);
  if(N2==0) return;  
  uint premove = rng.next()%N2;
  
  uint icell;
  double eno = ex[out]->enerNL(premove, icell);
  
  float arg = exp(-beta*(enn-eno + log( v2*(N1+1)/(v1*N2) )/beta) );

  if(rng.nextd()<arg){ //accepted
    fori(0, 3) ex[in]->ps.push_back(nps[i]); //add particle
    ex[out]->ps.erase(ex[out]->ps.begin()+3*premove, ex[out]->ps.begin() + 3*premove + 3); //remove particle
    ex[in]->make_linked_list();
    ex[out]->make_linked_list();
  }





}


void sample(){
  
  out<<e1.get_rho()<<" "<<e2.get_rho()<<" "<<e1.get_N()<<" "<<e2.get_N()<<" "<<e1.get_V()<<" "<<e2.get_V()<<" "<<-log(mu[0]/nsteps)/beta<<" "<<-log(mu[1]/nsteps)/beta<<endl;
}
