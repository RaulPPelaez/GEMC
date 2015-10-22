#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include"helper.h"
#include<algorithm>
#include<fstream>
typedef unsigned int uint;

class Experiment{
public:
  Experiment(){}

  void init(double T, uint N, double L, double rcut);
  
  void update();

  double enerNL(const double *ips) const;
  double enerNL(uint index, uint &icell) const;
  double total_ener() const;
  
  void make_linked_list();

  double get_V() const{ return L*L*L;}
  double get_rho() const{
    float N = ps.size()/3;
    return N/get_V();
  }
  uint get_N(){ return ps.size()/3;}

  vector<double> ps;
  double L;

private:
  double LJ(const double *p1, const double *p2) const;
  void adjust_step();
  void heal_list(uint cella, uint cellb);
  
  RRand rng;
  double beta;
  double rc, erc;

  double js; //jump step
  uint ntry, naccept;

  vector<uint> head, list;
  float lcell;
  uint ncells;
  vector<uint> redo;
  ofstream outjs;
};




#endif
