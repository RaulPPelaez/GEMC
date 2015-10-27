#include"Experiment.h"
#define PBC 
void Experiment::init(double T, uint N, double L, double rcut){
  this-> beta = 1.0/T;
  this->L = L;
  this->rc = rcut;
  this->rc2 = rc*rc;
  double rc6 = pow(1.0/rc,6);
  this->erc = 4.0*rc6*(rc6-1);
  ps.resize(3*N*N*N);
  int pindex = 0;
  float fac = 1.0 - 2.0/L;
  fori(0,N)forj(0,N)fork(0,N){
    ps[3*pindex+0] = (fac*( i/float(N-1) ) -0.5)+1.0/L;
    ps[3*pindex+1] = (fac*( j/float(N-1) ) -0.5)+1.0/L;
    ps[3*pindex+2] = (fac*( k/float(N-1) ) -0.5)+1.0/L;
    ++pindex;
  }

  this->rng.init();

  this->ncells = uint(L/rc + 0.5);
  this->lcell = L/(double)ncells;
  this->js = 0.005;

  make_linked_list();
  redo.resize(N*N*N+1); 

  outjs.open("js.dat");

  this->ener = total_ener();
}

void Experiment::update(){
  uint N = ps.size()/3;
  if(N==0) return;
  ntry++;
  uint p = rng.next()%N;

  double eno, enn;
  uint cella, cellb;
  eno = enerNL(p, cella);

  double jump[3];
  float js_scaled = js*lcell;
  fori(0,3){
    jump[i] = rng.nextd(-js_scaled, js_scaled);
    //    jump[2] = -js_scaled;
    ps[3*p+i] += jump[i];
  }

  float P=0.0;
  #ifdef PBC
    enn = enerNL(p, cellb);
    P = exp(-beta*(enn-eno));
  #else
  if(is_inbounds(&ps[3*p])){
    enn = enerNL(p, cellb);
    P = exp(-beta*(enn-eno));
  }
  #endif
  if(rng.nextd()>P){//if reject
    fori(0,3) ps[3*p+i] -= jump[i];
  }
  else{
    if(cella!=cellb) heal_list(cella, cellb); //if accepted fix linked list
    naccept++;
    ener = ener-eno+enn;
  }
  
  if(ntry%1000 == 0) adjust_step();
}

void Experiment::adjust_step(){
  float TRUST = 0.5;
  float ratio = (float)naccept/(float)ntry;
  if(ratio<TRUST){if(js>1e-6)this->js *=0.99000000f;}
  else if(ratio>TRUST){if(js<1) this->js *= 1.01000000f;}
  ntry = naccept = 0;
  //  outjs<<js*lcell<<" "<<lcell<<" "<<ratio<<"\n";
}

double Experiment::LJ(const double *p1, const double *p2) const{
  double r12[3];
  fori(0,3) r12[i] = p1[i]-p2[i];
  #ifdef PBC
  get_inbounds(r12); //PBC are in here
  #endif
  double r = square_mod(r12)*L*L;
  if(r>rc2) return 0.0;
  double r6=1/(r*r*r); 
  return 4.0*r6*(r6 - 1)-erc;  
}

void Experiment::make_linked_list(){
  uint N = ps.size()/3;
  this->ncells = uint(L/rc + 0.5);
  head.resize(ncells*ncells*ncells+1,0);
  list.resize(N+1,0);
  this->lcell = L/(double)ncells;
  std::fill(head.begin(), head.end(), 0);
  if(N==0) return;
  uint icell;
  uint celli, cellj, cellk;
  double ipos[3];
  fori(1,N+1){
    fork(0,3) ipos[k] = ps[3*(i-1)+k];
    #ifdef PBC
    get_inbounds(&ipos[0]);
    #endif
    celli = (uint)( (0.5 + ipos[0])*ncells);
    cellj = (uint)( (0.5 + ipos[1])*ncells);
    cellk = (uint)( (0.5 + ipos[2])*ncells);
    icell=1+celli+(cellj+cellk*ncells)*ncells;
    list[i] = head[icell];
    head[icell] = i;
  }
}

void Experiment::heal_list(uint cella, uint cellb){
  uint celli, cellj, cellk, icell, to_redo=0;
  uint i = head[cella];
  double ipos[3];

  while(i!=0){
    redo[to_redo]=i;
    i = list[i];
    to_redo++;
  }
  i = head[cellb];
  while(i!=0){
    redo[to_redo]=i;
    i = list[i];
    to_redo++;
  }
  std::sort(redo.begin(), redo.begin()+to_redo);
  head[cella] = head[cellb] = 0;

  forj(0,to_redo){
    i = redo[j];
    fork(0,3) ipos[k] = ps[3*(i-1)+k];
    #ifdef PBC
    get_inbounds(&ipos[0]);
    #endif
    celli = (uint)( (0.5 + ipos[0])*ncells);
    cellj = (uint)( (0.5 + ipos[1])*ncells);
    cellk = (uint)( (0.5 + ipos[2])*ncells);
    icell=1+celli+(cellj+cellk*ncells)*ncells;
    list[i] = head[icell];
    head[icell] = i;
  }
}

double Experiment::enerNL(const int &index, uint &icell) const{

  double ips[3];
  fori(0,3) ips[i] = ps[3*index + i];
  #ifdef PBC
  get_inbounds(&ips[0]);  
  #endif
  return enerNL(&ips[0], icell, index);
}
double Experiment::enerNL(const double *ips) const{
  uint icell;
  return enerNL(ips, icell,-1);
}

double Experiment::enerNL(const double *ips, uint &icell, const int &index = -1) const{
  double H = 0.0;
  int cell[3];
  forj(0,3) cell[j] = 1+int( ( 0.5 + ips[j])*ncells );
  icell=cell[0]+(cell[1]-1+(cell[2]-1)*ncells)*ncells;

  int j, jcel, jcelx,jcely,jcelz;
  for(int jx=cell[0]-1; jx<=cell[0]+1;jx++)
    for(int jy=cell[1]-1; jy<=cell[1]+1;jy++)
      for(int jz=cell[2]-1; jz<=cell[2]+1;jz++){
        if(jx==0) jcelx=ncells;else if(jx==ncells+1) jcelx=1;else jcelx = jx;
        if(jy==0) jcely=ncells;else if(jy==ncells+1) jcely=1;else jcely = jy;
        if(jz==0) jcelz=ncells;else if(jz==ncells+1) jcelz=1;else jcelz = jz;
        jcel =jcelx + (jcely-1)*ncells+(jcelz-1)*ncells*ncells;
        j = head[jcel];
        if(j==0) continue;
        do{
	  if((j-1)!=index)  H += LJ(ips, &ps[3*(j-1)] );
	  j=list[j];
        }while(j!=0);
      }
  return H;
}


double Experiment::total_ener() const{
  uint N = ps.size()/3;
  if(N==0) return 0.0;
  double E = 0.0;
  uint icell;
  fori(0, N) E += enerNL(i, icell); 
  return 0.5*E;
}

