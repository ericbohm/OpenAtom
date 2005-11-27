#include "standard_include.h"
#include "ckcomplex.h"

#include "../class_defs/allclass_gen.h"
#include "../class_defs/allclass_cp.h"
#include "../class_defs/allclass_mdintegrate.h"

#include "../class_defs/CP_OPERATIONS/class_cpintegrate.h"
#include "../../mathlib/proto_math.h"


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
// Dynamics
//============================================================================
void CPINTEGRATE::CP_integrate_dyn(const int ncoef, const int istate,int iteration,
                                   complex *forces,complex *vpsi,complex *psi,
                                   const int *k_x, const int *k_y, const int *k_z,
                                   double *cmass, double *fictEke_ret)
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

//============================================================================
// Update the velcoties to full step (if not the intial step).
// Compute the Fictious Kinetic Energy at time t (e.g. t=0 on initial step)

   if(iteration>1){
      for(int i = 0; i < ncoef; i++){
        vpsi[i] += forces[i]*(dt2/cmass[i]);     
      }//endfor
      // apply nhc
   }//endif

   double fictEke = 0.0;
   for(int i = 0; i < ncoef; i++){
     fictEke += vpsi[i].getMagSqr()*(0.5*cmass[i]);
   }//endfor
   (*fictEke_ret) = fictEke;

//============================================================================
// Update the velocities to half step (t+dt/2): 
// Update the positions to full step

   // apply nhc
   for(int i = 0; i < ncoef; i++){
     vpsi[i] += forces[i]*(dt2/cmass[i]);     
   }//endfor

   for(int i = 0; i < ncoef; i++){
     psi[i] += vpsi[i]*dt;
   }//endfor

//============================================================================
} // End function 
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
// Dynamics
//============================================================================
void CPINTEGRATE::CP_integrate_half_vel(int ncoef, int iteration, 
                        complex *forces,complex *vpsi,complex *psi,
                        int *k_x,int *k_y,int *k_z,double *cmass)
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

//============================================================================

   if(iteration>1){
      for(int i = 0; i < ncoef; i++){
        vpsi[i] += forces[i]*(dt2/cmass[i]);     
      }//endfor
      // apply nhc
   }//endif

//============================================================================
} // End function
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CPINTEGRATE::initCPNHC(int npts,int maxLen,int *len_nhc_ret,double *kTCP_ret,
                            double *tau_ret,double *mNHC)
//============================================================================
  {// Begin Function 
//============================================================================

  GENERAL_DATA *general_data = GENERAL_DATA::get();
  CP *cp = CP::get();
#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_cp.h"

  int len_nhc   = cptherm_info->len_c_nhc;
  double Text   = cpopts->te_ext;
  double tau_in = cpopts->cp_tau_nhc;

  double kT     = Text/BOLTZ;
  double tau    = tau_in/TIME_CONV;

  if(len_nhc>maxLen){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Maximum NHC length is %d < %d\n",maxLen,len_nhc);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
  }//endif

//============================================================================

  (*len_nhc_ret) = len_nhc;
  (*kTCP_ret)    = kT;
  (*tau_ret)     = tau;
  
  double pre     = tau*tau*kT;
  mNHC[0]        = pre*((double)npts);
  for(int i=1;i<len_nhc;i++){mNHC[i] = pre;}

//============================================================================
   }// End function
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CPINTEGRATE::CPSmplVel(int n,double *m,complex *v,int nNhc,
                            double *mNhc,double *vNhc,double kT,int istart_cp)
//============================================================================
    {//Begin Function
//============================================================================
// Sample the velocities 

  if(istart_cp<3){
     double *vtmp = new double [n];
       sampl1DVelOneT(n,vtmp,m,kT);
       for(int i=0;i<n;i++){v[i].re = vtmp[i];}
       sampl1DVelOneT(n,vtmp,m,kT);
       for(int i=0;i<n;i++){v[i].im = vtmp[i];}
     delete [] vtmp;
  }//endif
  sampl1DVelOneT(nNhc,vNhc,mNhc,kT);

//============================================================================
// Scale the velocities

  double temp = mNhc[0]*vNhc[0]*vNhc[0];
  for(int i = 0;i<n;i++){temp += v[i].getMagSqr()*m[i];}
  double scale = sqrt( (((double)(2*n))*kT/temp) );

  vNhc[0] *= scale;
  for(int i = 0;i<n;i++){
    v[i].re *= scale;  v[i].im *= scale;
  }//endfor

//============================================================================
    }// End function
//============================================================================


//==================================================================== 
// Velocities  
//==================================================================== 
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==================================================================== 
void CPINTEGRATE::sampl1DVelOneT(int n, double* v,double* mass,
                                   double  kT)
//=================================================================== 
   {//begin routine 
//===================================================================

  MDINTEGRATE  *mdintegrate  = MDINTEGRATE::get();                  
#include "../class_defs/allclass_strip_mdintegrate.h"
   double *qseed  =  &(mdvel_samp->qseed);
   int    *iseed  =  &(mdvel_samp->iseed);
   int    *iseed2 =  &(mdvel_samp->iseed2);

//===================================================================
// Sample unit gaussian

   double *temp  = new double[n];
   double *ptemp = temp-1;

   gaussran(n,iseed,iseed2,qseed,ptemp);
   for(int i=0;i<n;i++){v[i] = temp[i]; }

   delete [] temp;

//===================================================================
// Apply the width

   for(int i=0;i<n;i++){
     double width = sqrt(kT/mass[i]); //kT has boltz
     v[i] *= width;
   }//endfor 

//------------------------------------------------------------------
 } //end routine 
//=================================================================== 

