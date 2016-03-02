// TODO spherical truncation and sheet truncation

#include <cstdlib>
#include <iostream>
#include "include/ckcomplex.h"
#include "class_defs/vcoulb.h"
#include "util.h"
#include "coulomb.h"

// set ng, ga, gb, gc and check if GSPACE is set 
void VCOULB::construct(){

  if (gs == NULL){
    Die("class VCOULB requires GSPACE");
  }
  else{
    ng = gs->ng;
    ga = gs->ig;
    gb = gs->jg;
    gc = gs->kg;
    data = new double [ng];
  }

}


// TODO combine calc_vcoulb routine
void VCOULB::compute(){

  calc_vcoulb(qIndex, data, gs, sys);

}

void VCOULB::compute_spherical(){

}

void VCOULB::compute_sheet(){

}


void VCOULB::print(){

  for(int i=0; i<ng; i++){
    printf(" G: (%d, %d, %d), Vc: %lg \n", ga[i], gb[i], gc[i], data[i]);
  }

}
