#ifndef VCOULB_H
#define VCOULB_H

#include "../include/ckcomplex.h"
#include "gspace.h"
#include "sysinfo.h"

// Coulomb potential calculation

class VCOULB{

 public:

  // q index
  int qIndex; 

  // g vector information
  int ng;   // number of g vectors
  int *ga;   
  int *gb;  
  int *gc;

  // GSPACE is used
  GSPACE *gs;// = NULL;

  // SYSINFO is used
  SYSINFO sys;

  double *data;     // coulomb potential, 1D data

  // constructor 
  // set q index
  VCOULB(int q, GSPACE* _gs, SYSINFO& _sys){
    qIndex = q;
    gs = _gs;
    sys = _sys;
    construct();
    // TODO We should be able to call compute_spherical or compute_sheet
    // depends on the coulomb truncation option
    compute();
  }



  //====================
  // member functions
  //====================
  // set ng, ga, gb, gc and check pointers
  void construct();
  // computing coulomb potential data[iG] = sqrt(q+G)^-1 
  void compute(); 
  // spherical truncation for molecules, clusters, etc...
  void compute_spherical();
  // truncation for 2D structure
  void compute_sheet();
  // print
  void print();

};

#endif
