#include "standard_include.h"
#include "ckcomplex.h"

#include "../class_defs/allclass_gen.h"
#include "../class_defs/allclass_cp.h"
#include "../class_defs/allclass_mdintegrate.h"

#include "../class_defs/CP_OPERATIONS/class_cpintegrate.h"
#include "../../mathlib/proto_math.h"

#define ISOKIN_OPT 1

//============================================================================
// Dynamics controller
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CPINTEGRATE::CP_integrate_dyn(int ncoef, int istate,int iteration,
              complex *forces,complex *vpsi,complex *psi,double *cmass, 
              int *k_x, int *k_y,int *k_z,int len_nhc, int num_nhc,
              double **fNHC,double **vNHC,double *xNHC,double mNHC,double kTCP,
              double *fictEke,int nkx0_red,int nkx0_uni,int nkx0_zero,
              double *ekeNhc, double degfree,double degfreeNHC,
              double gammaNHC)
//============================================================================
    { // Begin Function 
//----------------------------------------------------------------------------
// Local Variables

  GENERAL_DATA *general_data = GENERAL_DATA::get();
  CP *cp = CP::get();
#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_cp.h"

   double dt      = gentimeinfo->dt;
   double dt2     = dt*0.5;
   int istrt      = nkx0_red;
   int applyorder;

//============================================================================
// Update the velocities from time, t-dt/2, to time, t (if not the intial step).
// Compute the Fictious Kinetic Energy at time, t (e.g. t=0 on initial step).

   if(iteration>1){
     applyorder = 2;
     cp_evolve_vel(ncoef,forces,vpsi,cmass,len_nhc,num_nhc,fNHC,vNHC,
                   xNHC,mNHC,kTCP,nkx0_red,nkx0_uni,nkx0_zero,applyorder,iteration,
                   degfree,degfreeNHC,gammaNHC);
   }//endif

   get_fictKE(ncoef,vpsi,cmass,len_nhc,num_nhc,vNHC,mNHC,nkx0_red,nkx0_uni,nkx0_zero,
              fictEke,ekeNhc,gammaNHC);

//============================================================================
// Update the velocities from time, t, to time, t+dt/2. 

   applyorder = 1;
   cp_evolve_vel(ncoef,forces,vpsi,cmass,len_nhc,num_nhc,fNHC,vNHC,
                 xNHC,mNHC,kTCP,nkx0_red,nkx0_uni,nkx0_zero,applyorder,iteration,
                 degfree,degfreeNHC,gammaNHC);

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
void CPINTEGRATE::cp_evolve_vel(int ncoef_full, complex *forces, complex *vpsi,
                    double *cmass,int len_nhc, int num_nhc, double **fNHC,
                    double **vNHC,double *xNHC,double mNHC,double kTCP,
                    int nkx0_red,int nkx0_uni,int nkx0_zero,
                    int applyorder,int iteration,double degfree,
                    double degfreeNHC,double gammaNHC)
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
   int istrt0     = nkx0_red;
   int istrt1     = nkx0_red+nkx0_zero;
   int iend       = nkx0_uni+nkx0_red;
   int igo        = 1;
   int ncoef      = ncoef_full-istrt0;
   if(applyorder==2 && iteration==1){igo=0;}

//============================================================================
// Make life simple

   for(int i=istrt0;i<istrt1;i++){forces[i].im=0.0;vpsi[i].im=0.0;}
   for(int i=istrt1;i<iend;i++){
     forces[i].re *= 2.0;
     forces[i].im *= 2.0;
     cmass[i]     *= 2.0;
   }//endif

//============================================================================
// Isokin-NHC apply first

   if(applyorder==1 && ISOKIN_OPT==1){
     cp_isoNHC_update(ncoef,&vpsi[istrt0],&cmass[istrt0],
                      len_nhc,num_nhc,xNHC,vNHC,fNHC,mNHC,kTCP,
                      degfree,degfreeNHC,gammaNHC);
   }//endif

//============================================================================
// Evolve velocity with forces : isokinetically or newtonianalyy

   if(ISOKIN_OPT==0 && igo==1){
     for(int i=istrt0;i<ncoef_full;i++){
       vpsi[i] += forces[i]*(dt2/cmass[i]);     
     }//endfor
   }//endif

   if(ISOKIN_OPT==1 && igo==1){
     cp_isoVel_update(ncoef,&vpsi[istrt0],&forces[istrt0],&cmass[istrt0],
                      num_nhc,vNHC,mNHC,kTCP,degfree,degfreeNHC,gammaNHC);
   }//endif

//============================================================================
// Isokin-NHC apply second

   if(ISOKIN_OPT==1 && applyorder==2 && igo==1){
     cp_isoNHC_update(ncoef,&vpsi[istrt0],&cmass[istrt0],
                      len_nhc,num_nhc,xNHC,vNHC,fNHC,mNHC,kTCP,
                      degfree,degfreeNHC,gammaNHC);
   }//endif

//============================================================================

   for(int i=istrt0;i<istrt1;i++){forces[i].im=0.0;vpsi[i].im=0.0;}
   for(int i=istrt1;i<iend;i++){
     forces[i].re /= 2.0;
     forces[i].im /= 2.0;
     cmass[i]     /= 2.0;
   }//endif

//============================================================================
  } // End function
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CPINTEGRATE::get_fictKE(int n,complex *v, double *m,int len, int num,
          double **vNHC,double mNHC,int nkx0_red,int nkx0_uni,int nkx0_zero,
          double *ekin_ret,double *ekinNHC_ret,double gamma)
//============================================================================
  {// Begin Function 
//============================================================================

  int istrt0     = nkx0_red;
  int istrt      = nkx0_red+nkx0_zero;
  int iend       = nkx0_red+nkx0_uni;

  double ekin = 0.0;
  double ekinNHC = 0.0;
  if(ISOKIN_OPT==1){
    for(int i=0;i<num;i++){
    for(int j=1;j<len;j++){
      ekinNHC += mNHC*vNHC[i][j]*vNHC[i][j];
    }}//endfor
    for(int i=0;i<num;i++){
      ekin += gamma*mNHC*vNHC[i][0]*vNHC[i][0];
    }//endfor
  }//endif
  for(int i=istrt0;i<istrt;i++){ekin += v[i].getMagSqr()*m[i];}
  for(int i=istrt;i<iend;i++)  {ekin += v[i].getMagSqr()*(2.0*m[i]);}
  for(int i=iend;i<n;i++)      {ekin += v[i].getMagSqr()*m[i];}
  ekin     *=0.5;
  ekinNHC  *=0.5;

  (*ekin_ret)    = ekin;
  (*ekinNHC_ret) = ekinNHC;

//============================================================================
  }//end routine
//============================================================================



