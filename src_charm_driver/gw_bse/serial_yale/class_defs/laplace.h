// class LAPLACE
// LAPLACE contains Gauss-Laugerre nodes and weights
// and also calculate gap and bandwidth of the system
// Sets number of windows

#ifndef LAPLACE_H
#define LAPLACE_H

#include "../include/ckcomplex.h"
#include "sysinfo.h"
#include "usrinput.h"
#include "matrix.h"
#include "states.h"

class LAPLACE{
  
 public:

  // system information
  SYSINFO* sys;

  // user input
  USRINPUT* usrin;
  
  // error tolerance (percentage tolerance) 
  double ptol;
  // number of nodes written in a file (GLagNodeWeight.dat)
  int *npt;
  // Gauss-Laguerre quadrature nodes
  double **n;
  // Gauss-Laguerre quadrature weights
  double **w;
  
  // number of windows in unoccupied states
  int nwunocc;
  // number of windows in occupied states
  int nwocc;
  // total number of windows = nwunocc x nwocc
  int totnw;

  // when nwocc and nwunocc are determined, we also need to know how to slice the energy
  // these values are determined in set_windows() member function
  int *slice_occ;
  int *slice_unocc;

  // to find optimal number of windows, we need to know gap and bandwidth of the system
  double maxEocc;
  double minEocc;
  double maxEunocc;
  double minEunocc;

  // gap of the system
  double gap;
  // bandwidth of the system
  double bandwidth;

  // energy bin
  double *bin_occ;
  double *bin_unocc;

  // number of nodes to be used for each window
  int **Nnodes;
  // optimal a
  double **opta;
  

  // set two parameter (epsc,epsv) : see the note
  // this may need to be adjusted for windowing. Ask Glenn.
  double epsc;
  double epsv;
 
  // counting number of calculations
  int ncounter;

 
  // constructor
  LAPLACE(){};
  LAPLACE(USRINPUT& _usrin, SYSINFO& _sys, STATES*** psi){

    usrin = &_usrin;
    sys = &_sys;
    ptol = usrin->ptol;
    // read Gauss-Laguerre nodes and weights from file
    get_nodes_weights();
    // initialize minE and maxE for occupied and unoccupied states
    initialize(psi[0][0], usrin->nocc, usrin->nunocc);
    // calculate band gap and band width
    for (int is=0; is < usrin->nspin; is++){
      for (int ik=0; ik < usrin->nkpt; ik++){
	set_GapBandwidth(psi[is][ik], usrin->nocc, usrin->nunocc);
      }
    }
    // set number of windows to be used
    set_windows();
    // if you do not want to use windows...
    set_noWindow();
    // set bins according to nwocc, nwunocc, slice_occ, and slice_unocc
    set_bin();
    // now we need to set up the number of nodes to be used at each window to satisfy the error tolerance
    // inside of this function, func_opta is called
    set_numnodes();

    // set epsc and epsv. this may need to be modified
    epsc = minEunocc - (minEunocc-maxEocc)/3;
    epsv = maxEocc + (minEunocc-maxEocc)/3;

    // initialize counter
    ncounter = 0;

    // add some padding at the edge of the bin
    bin_occ[0] -= 0.00001;
    bin_occ[nwocc] += 0.00001;
    bin_unocc[0] -= 0.00001;
    bin_unocc[nwunocc] += 0.00001; 

#ifdef DEBUG
    LAPLACE_class_out();
#endif
   
  }// end constructor


  // member functions
  void get_nodes_weights(); // get Gauss-Laguerre quadrature nodes and weights from the file
  void initialize(STATES*,int,int); // initialize maxEocc,maxEunocc,minEocc,minEunocc
  void set_GapBandwidth(STATES*,int,int); // find the gap and bandwidth of the system
  void set_windows();  // decide how many windows needed to minimize the number of computation
  void set_bin(); 
  void set_numnodes(); // set Nnodes for each windows pair
  void func_opta(double,double,int,double&); // function that finds optimal a for the quadrature
  void set_noWindow(); // do not use windowing


  // printing out all variables
#ifdef DEBUG
  void LAPLACE_class_out(){
    FILE *fp; fp = fopen("LAPLACE_class_out","w");
    fprintf(fp,"error tolerance (percentage): %lg\n",ptol);
    fprintf(fp,"occupied state energy range [ %lg , %lg ]\n",minEocc,maxEocc);
    fprintf(fp,"unoccupied state energy range [ %lg , %lg ]\n",minEunocc,maxEunocc);
    fprintf(fp,"gap and bandwidth: %lg, %lg\n",gap,bandwidth);
    fprintf(fp,"number of windows for occupied: %d and unoccupied: %d\n",nwocc,nwunocc);
    fprintf(fp,"total number of windows: %d\n",totnw);
    fprintf(fp,"slicing (max=10) - occupied\n");
    for(int i=0; i<nwocc; i++){
      fprintf(fp,"%d\n",slice_occ[i]);
    }
    fprintf(fp,"slicing (max=10) - unoccupied\n");
    for(int i=0; i<nwunocc; i++){
      fprintf(fp,"%d\n",slice_unocc[i]);
    }
    fprintf(fp,"number of nodes and optimal a values\n");
    for(int i=0; i<nwocc; i++){
      for(int j=0; j<nwunocc; j++){
	fprintf(fp,"window (%d %d): Nnodes %d, opta %lg\n",i,j,Nnodes[i][j],opta[i][j]);
      }
    }
    fprintf(fp,"energy bins - occupied\n");
    for(int i=0; i<nwocc; i++){
      fprintf(fp,"[ %lg  %lg ]\n",bin_occ[i],bin_occ[i+1]);
    }
    fprintf(fp,"energy bins - unoccupied\n");
    for(int i=0; i<nwunocc; i++){
      fprintf(fp,"[ %lg %lg ]\n",bin_unocc[i],bin_unocc[i+1]);
    }
    fclose(fp);
  }
#endif    
  




  
};

#endif
