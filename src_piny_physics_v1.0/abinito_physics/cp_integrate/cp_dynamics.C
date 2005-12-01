#include "standard_include.h"
#include "ckcomplex.h"

#include "../class_defs/allclass_gen.h"
#include "../class_defs/allclass_cp.h"
#include "../class_defs/allclass_mdintegrate.h"

#include "../class_defs/CP_OPERATIONS/class_cpintegrate.h"
#include "../../mathlib/proto_math.h"

#define ISOKIN_OPT 0

//============================================================================
// Dynamics controller
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CPINTEGRATE::CP_integrate_dyn(int ncoef, int istate,int iteration,
              complex *forces,complex *vpsi,complex *psi,double *cmass, 
              int *k_x, int *k_y,int *k_z,int len_nhc, int num_nhc,
              double **fNHC,double **vNHC,double *xNHC,double mNHC,double kTCP,
              double *fictEke,int nkx0_red,int nkx0_uni,int nkx0_zero)
//============================================================================
    { // Begin Function 
//----------------------------------------------------------------------------
// Local Variables

  GENERAL_DATA *general_data = GENERAL_DATA::get();
  CP *cp = CP::get();
#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_cp.h"

   double dt   = gentimeinfo->dt;
   double dt2  = dt*0.5;
   int istrt   = nkx0_red;
   int iend    = nkx0_uni+nkx0_red;
   int applyorder;

//============================================================================
// Update the velocities from time, t-dt/2, to time, t (if not the intial step).
// Compute the Fictious Kinetic Energy at time, t (e.g. t=0 on initial step).

   if(iteration>1){
     applyorder = 2;
     cp_evolve_vel(ncoef,forces,vpsi,cmass,len_nhc,num_nhc,fNHC,vNHC,
                   xNHC,mNHC,kTCP,nkx0_red,nkx0_uni,nkx0_zero,applyorder,iteration);
   }//endif

   get_fictKE(ncoef,vpsi,cmass,len_nhc,num_nhc,vNHC,mNHC,
              nkx0_red,nkx0_uni,nkx0_zero,fictEke);

//============================================================================
// Update the velocities from time, t, to time, t+dt/2. 

   applyorder = 1;
   cp_evolve_vel(ncoef,forces,vpsi,cmass,len_nhc,num_nhc,fNHC,vNHC,
                 xNHC,mNHC,kTCP,nkx0_red,nkx0_uni,nkx0_zero,applyorder,iteration);

//============================================================================
// Update the positions to the next step (t+dt)

   for(int i=istrt;i<ncoef;i++){
     psi[i] += vpsi[i]*dt;
   }//endfor

//============================================================================
  } // End function 
//============================================================================


//============================================================================
// Fancy Dan velocity integration under CP
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CPINTEGRATE::cp_evolve_vel(int ncoef, complex *forces, complex *vpsi,
                    double *cmass,int len_nhc, int num_nhc, double **fNHC,
                    double **vNHC,double *xNHC,double mNHC,double kTCP,
                    int nkx0_red,int nkx0_uni,int nkx0_zero,
                    int applyorder,int iteration)
//============================================================================
   { // Begin Function
//----------------------------------------------------------------------------
// Local Variables and error checking

  GENERAL_DATA *general_data = GENERAL_DATA::get();
  CP *cp = CP::get();
#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_cp.h"

   double dt      = gentimeinfo->dt;
   double dt2     = dt*0.5;
   int istrt      = nkx0_red;
   int iend       = nkx0_uni+nkx0_red;
   int igo        = 1;
   if(applyorder==2 && iteration==1){igo=0;}

//============================================================================
// Isokin-NHC apply first

#ifdef NEXT_STEP
   if(applyorder==1 && ISOKIN_OPT==1){
     cp_NHC_apply();
   }//endif
#endif

//============================================================================
// Evolve velocity with forces

   if(ISOKIN_OPT==0 && igo==1){
     for(int i=istrt;i<ncoef;i++){
       vpsi[i] += forces[i]*(dt2/cmass[i]);     
     }//endfor
   }//endif

#ifdef NEXT_STEP
   if(ISOKIN_OPT==1 && igo==1){
     cp_isokin_update();
   }//endif
#endif

//============================================================================
// Isokin-NHC apply second

#ifdef NEXT_STEP
   if(applyorder==2 && igo==1){
     cp_NHC_apply();
   }//endif
#endif

//============================================================================
  } // End function
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CPINTEGRATE::get_fictKE(int n,complex *v, double *m,int num, int len,
          double **vNHC,double mNHC,int nkx0_red,int nkx0_uni,int nkx0_zero,
          double *EKin_ret)
//============================================================================
  {// Begin Function 
//============================================================================

  int istrt0     = nkx0_red;
  int istrt      = nkx0_red+nkx0_zero;
  int iend       = nkx0_red+nkx0_uni;

  double EKin = 0.0;
  if(ISOKIN_OPT==1){
    for(int i=0;i<num;i++){
    for(int j=0;j<len;j++){
      EKin += mNHC*vNHC[i][j]*vNHC[i][j];
    }}//endfor
  }//endif
  for(int i=istrt0;i<istrt;i++){EKin += v[i].getMagSqr()*m[i];}
  for(int i=istrt;i<iend;i++)  {EKin += v[i].getMagSqr()*(2.0*m[i]);}
  for(int i=iend;i<n;i++)      {EKin += v[i].getMagSqr()*m[i];}
  EKin       *=0.5;

  (*EKin_ret) = EKin;

//============================================================================
  }//end routine
//============================================================================



