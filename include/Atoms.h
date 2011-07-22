/*****************************************************************************
 * $Source$
 * $Author$
 * $Date$
 * $Revision$
 *****************************************************************************/

/** \file Atoms.h
 *
 */

#include "charm++.h"

#ifndef _atoms_h_
#define _atoms_h_

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class Atom {
//============================================================================
 public:
  Atom() {}
  Atom(double x_, double y_, double z_, double q_, double m_) {
     x = x_;   y = y_;   z = z_;  q = q_; m = m_; 
     xu =  yu = zu=0.0;
    vx = 0.0; vy = 0.0; vz = 0.0;
    fx = 0.0; fy = 0.0; fz = 0.0;
    fxu = 0.0; fyu = 0.0; fzu = 0.0;
    mvx2 = 0.0; mvy2 = 0.0; mvz2 = 0.0;
  }//end constructor
  Atom(double x_, double y_, double z_,double q_, double m_,
       double vx_, double vy_, double vz_) {
     x =  x_;  y =  y_;  z =  z_;  q = q_; m = m_; 
     xu =  yu = zu=0.0;
    vx = vx_; vy = vy_; vz = vz_; 
    fx = 0.0; fy = 0.0; fz = 0.0;
    fxu = 0.0; fyu = 0.0; fzu = 0.0;
    mvx2 = 0.0; mvy2 = 0.0; mvz2 = 0.0;
  }//end constructor
  void Init(double x_, double y_, double z_, double q_, double m_) {
    x = x_; y = y_; z = z_; q = q_; m = m_; 
    xu =  yu = zu=0.0;
    xcm =  ycm = zcm=0.0;
    vx = 0.0; vy = 0.0; vz = 0.0;
    fx = 0.0; fy = 0.0; fz = 0.0;
    fxu = 0.0; fyu = 0.0; fzu = 0.0;
    mvx2 = 0.0; mvy2 = 0.0; mvz2 = 0.0;
  }//end Init
  void Init(double x_, double y_, double z_, double q_, double m_,
            double vx_, double vy_, double vz_) {
     x =  x_;  y =  y_;  z =  z_;  q = q_; m = m_; 
     xu =  yu = zu=0.0;
    vx = vx_; vy = vy_; vz = vz_; 
    fx = 0.0; fy = 0.0; fz = 0.0;
    fxu = 0.0; fyu = 0.0; fzu = 0.0;
    mvx2 = 0.0; mvy2 = 0.0; mvz2 = 0.0;
  }//end Init
  Atom(Atom *a) { 
     x = a->x;   y = a->y;   z = a->z;  q = a->q; m = a->m; 
     xu = a->xu;   yu = a->yu;   zu = a->zu;  mu = a->mu; 
     vx = a->vx; vy = a->vy; vz = a->vz;
     fx = a->fx; fy = a->fy; fz = a->fz;
     fxu = a->fxu; fyu = a->fyu; fzu = a->fzu;
     mvx2 = 0.0; mvy2 = 0.0; mvz2 = 0.0;
  }//end constructor


  double x, y, z, q, m;
  double fx,fy,fz;
  double vxold,vyold,vzold;
  double xold,yold,zold;
  double vx,vy,vz,mvx2,mvy2,mvz2;

  double xu, yu, zu, mu;
  double xcm, ycm, zcm;
  double fxu,fyu,fzu;

  void pup(PUP::er &p) {
     p|x;  p|y;  p|z;  p|q; p|m; 
     p|fx; p|fy; p|fz; 
     p|vx; p|vy; p|vz;
     p|mvx2; p|mvy2; p|mvz2;
     p|xu;  p|yu;  p|zu; p|mu; 
     p|fxu; p|fyu; p|fzu; 
  }//pup routine

};//end Atom

PUPmarshall(Atom);
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class FastAtoms {
//============================================================================
 public:

  int natm;
  double *xu,*yu,*zu,*q;  // no masses or velocities cause these are working vectors
  double *fxu,*fyu,*fzu;
  double *x,*y,*z;  // no masses or velocities cause these are working vectors
  double *fx,*fy,*fz;
  FastAtoms(){}
 ~FastAtoms(){}
  
};//end FastAtom
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class PIMD_CM {
//============================================================================
 public:
  int natm;
  double *x,*y,*z;
  PIMD_CM(){x=NULL; y=NULL;z=NULL;}
  ~PIMD_CM(){}
  
};//end PIMD_CM
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class AtomNHC {
//============================================================================
 public:
   int len_nhc;
   double kT;
   double posKT;
   double fx[5],fy[5],fz[5];
   double vx[5],vy[5],vz[5];
   double m[5];
  AtomNHC() {}
  AtomNHC(int len_nhc_in,double kT_in, double *vx_in, double *vy_in, double *vz_in, 
          double *m_in) {
     len_nhc = len_nhc_in;
     checkLenNHC(len_nhc);
     kT      = kT_in;
     posKT = 0.0;
     for(int i=0;i<len_nhc;i++){
        m[i] =  m_in[i]; 
       vx[i] = vx_in[i]; vy[i] = vy_in[i]; vz[i] = vz_in[i];
       fx[i] = 0.0;      fy[i] = 0.0;      fz[i] = 0.0;
     }//endfor
  }//end constructor
  AtomNHC(AtomNHC *a) { 
     len_nhc = a->len_nhc;
     kT      = a->kT;
     posKT = 0.0;
     checkLenNHC(len_nhc);
     for(int i=0;i<len_nhc;i++){
        m[i] = a->m[i]; 
       vx[i] = a->vx[i]; vy[i] = a->vy[i]; vz[i] = a->vz[i];
       fx[i] = a->fx[i]; fy[i] = a->fy[i]; fz[i] = a->fz[i];
     }//endfor
  }//end constructor

  void Init(int len_nhc_in,double kT_in, double *vx_in, double *vy_in, double *vz_in, 
            double *m_in) {
     len_nhc = len_nhc_in;
     kT      = kT_in;
     posKT = 0.0;
     checkLenNHC(len_nhc);
     for(int i=0;i<len_nhc;i++){
        m[i] =  m_in[i]; 
       vx[i] = vx_in[i]; vy[i] = vy_in[i]; vz[i] = vz_in[i];
       fx[i] = 0.0;      fy[i] = 0.0;      fz[i] = 0.0;
     }//endfor
  }//end Init
  void Init(AtomNHC *a) { 
     len_nhc = a->len_nhc;
     kT      = a->kT;
     posKT = 0.0;
     checkLenNHC(len_nhc);
     for(int i=0;i<len_nhc;i++){
        m[i] = a->m[i]; 
       vx[i] = a->vx[i]; vy[i] = a->vy[i]; vz[i] = a->vz[i];
       fx[i] = a->fx[i]; fy[i] = a->fy[i]; fz[i] = a->fz[i];
     }//endfor
  }//end Init

  ~AtomNHC() { }// destructor

   void pup(PUP::er &p) {
     p|len_nhc; p|kT; p|posKT;
     p(m,len_nhc);
     p(vx,len_nhc); p(vy,len_nhc); p(vz,len_nhc);
     p(fx,len_nhc); p(fy,len_nhc); p(fz,len_nhc);
  }//routine

   void checkLenNHC(int mylen_nhc){
     if(mylen_nhc>5){
       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
       CkPrintf("The maximum allowed Atm NHC len is 5.\n");
       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
       CkExit();
     }//endif
   }//endif
};//end AtomNHC

PUPmarshall(AtomNHC);
//============================================================================

#endif
