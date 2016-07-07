#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H

#include "../include/ckcomplex.h"
#include "states.h"
#include "sysinfo.h"
#include "matrix.h"

class INTERPOLATOR{

 public:
  int thisqIndex;    // q index

  complex** colG;   // column of G: G_z,k,r' = sum_c psi_k,c.conj * psi_k,c,r'/(E_c - z)
  double* z;         // energy values for interpolation

  int nzpt;          // number of energy points to be calculated not including EoccMax and EoccMin
  double EoccMax;    // maximum value of occupied state energy at this k+q
  double EoccMin;    // minimum value of occupied state energy at this k+q
  double EunoccMin;  // minimum value of unoccupied state energy at this k
  double Egap;       // energy gap between occupied
  
  STATES* psiR;      // psi R at this k point
  STATES* psiRq;     // psi R at this k+q point


  SYSINFO* sys;
  int nkpt;          // number of k points
  int nocc;          // number of occupied states
  int nunocc;        // number of unoccupied states
  int nr;            // number of real-space grid

  // constructor
  INTERPOLATOR(){};
  INTERPOLATOR(int q, int n, STATES* _psiR, STATES* _psiRq, SYSINFO& _sys){
    thisqIndex = q;
    nzpt = n;
    psiR = _psiR;
    psiRq = _psiRq;
    sys = &_sys;
    nkpt = sys->nkpt;
    nocc = sys->nocc;
    nunocc = sys->nunocc;
    EoccMax = psiRq->eig[nocc-1]; // psi_v,k+q
    EoccMin = psiRq->eig[0];      //psi_v,k+q
    EunoccMin = psiR->eig[nocc];  //psi_c,k
    Egap = EunoccMin - EoccMax;
    nr = psiR->ndata;

    // for now, we comment it out
    /*
    // allocate colG (dimension: nzpt x nr)
    colG = new complex* [nzpt+1];  // in make_zlist routine, nzpt will increment
    for(int iz=0; iz<nzpt+1; iz++){
      colG[iz] = new complex [nr];
    }

    // initialize
    for(int iz=0; iz<nzpt+1; iz++){
      for(int ir=0; ir<nr; ir++){
	colG[iz][ir] = 0;
      }
    }
    */

  }// end constructor

  // member functions
  void make_zlist();
  void initialize(int);
  void G_calculator(int, int);
  void interpolate(double, complex*);
  void LinearInterpolation(double, complex*);
  void QuadraticInterpolation(double, complex*);

#ifdef DEBUG
  void print_G(int);
  void print_values();
  void G_eachband(int,int);
#endif
  

};

#endif
