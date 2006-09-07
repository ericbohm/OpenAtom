//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
#include "standard_include.h"
#include "ckcomplex.h"

#include "../class_defs/allclass_gen.h"
#include "../class_defs/allclass_cp.h"
#include "../class_defs/allclass_mdintegrate.h"

#include "../class_defs/CP_OPERATIONS/class_cpintegrate.h"
#include "../../mathlib/proto_math.h"

//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CPINTEGRATE::initCPNHC(int npts,int maxLen,int maxNum, int *len_nhc_ret,
                        int *num_nhc_ret, double *kTCP_ret,
                        double *tau_ret,double *mNHC,double *degfree_ret,
                        double *degfreeNHC_ret,double *gammaNHC_ret,int ncoef_true,
                        int ncoef_zero,double *v0,double *a2, double *a4)
//============================================================================
  {// Begin Function 
//============================================================================

  GENERAL_DATA *general_data = GENERAL_DATA::get();
  CP *cp = CP::get();
#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_cp.h"


  int len_nhc   = cptherm_info->len_c_nhc;
  int num_nhc   = cptherm_info->num_c_nhc_iso;
  double Text   = cpopts->te_ext;
  double tau_in = cpopts->cp_tau_nhc;

  double kT     = Text/BOLTZ;
  double tau    = tau_in/TIME_CONV;

//============================================================================

  if(len_nhc>maxLen){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Maximum NHC length is %d < %d\n",maxLen,len_nhc);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
  }//endif
  if(num_nhc>maxNum){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Maximum NHC number is %d < %d\n",maxNum,num_nhc);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
  }//endif
  if(len_nhc!=2){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Minimum NHC len_nhc is 2 > %d\n",maxNum,num_nhc);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
  }//endif

  double degfree    = (double)(2*ncoef_true-ncoef_zero);
  double degfreeNHC = (double)(num_nhc);     // 1 constraint for each NHC of length 2
  double gammaNHC   = degfree/(degfree+1.0); // not used anymore

//============================================================================

  (*len_nhc_ret)    = len_nhc;
  (*num_nhc_ret)    = num_nhc;
  (*kTCP_ret)       = kT;
  (*tau_ret)        = tau;
  
  double pre        = tau*tau*kT;
  mNHC[0]           = pre*degfree;  // NHC mass
  mNHC[1]           = 2.5*pre;      // Double well mass : pos is dimensionless

  (*degfree_ret)    = degfree;
  (*degfreeNHC_ret) = degfreeNHC;
  (*gammaNHC_ret)   = gammaNHC;

  for(int i=0;i<num_nhc;i++){
    double r1   = (double)(num_nhc+i+5);
    double r2   = (double)(num_nhc+5);
    double frac = r1/r2;
    v0[i]       = 0.8*kT*frac;     // dimensionfull barrier
    a2[i]       = 1.6*(0.8*frac);  // dimensionless `length^2'
    a4[i]       = a2[i]*a2[i];
  }//endfor

//============================================================================
  }// End function
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CPINTEGRATE::cp_isoNHC_update(int n,complex *v,
                  double *m,int len_nhc,int num_nhc,double **xNHC,double **xNHCP,
                  double **vNHC,double **fNHC,double *mNHC,
                  double *v0,double *a2,double *a4, double kT,
                  double degfree, double degfreeNHC, double gammaNHC)
//============================================================================
  {//begin routine
//============================================================================
// Get the yoshida stuff out

  GENERAL_DATA *general_data = GENERAL_DATA::get();
  CP *cp = CP::get();
#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_cp.h"

  double dt  = gentimeinfo->dt;
  int nresp  = cptherm_info->nres_c_nhc;
  int nyosh  = cptherm_info->nyosh_c_nhc;

  double dti = dt/((double)nresp);
  double gkt = degfree*kT;

  double wdti[20],wdti2[20],wdti4[20],wdti8[20];
  set_yosh(nyosh,dti,wdti,wdti2,wdti4,wdti8);

//============================================================================
// Precompute the kinetic energy and get initial forces

   double ekin = 0.0;
   for(int i=0;i<n;i++){
     ekin += m[i]*(v[i].re*v[i].re+v[i].im*v[i].im);
   }//endfor

   for(int i=0;i<num_nhc;i++){
     double ekinNHC = mNHC[0]*vNHC[i][0]*vNHC[i][0]
                    + mNHC[1]*vNHC[i][1]*vNHC[i][1];
   }

   for(int i=0;i<num_nhc;i++){xNHCP[i][0] = 0.0;}

   for(int i=0;i<num_nhc;i++){
     fNHC[i][0] = (ekin - gkt)/mNHC[0];
     fNHC[i][1] = -( 4.0*v0[i]/(a4[i]*mNHC[1]) )
                  *(xNHCP[i][1]*(xNHCP[i][1]*xNHCP[i][1]-a2[i]));
   }//endfor

//============================================================================
// (II) Yoshida-Suzuki step yourself to heaven

  double scale = 1.0;
  for(int iresn=1;iresn<=nresp;iresn++){
  for(int iyosh=1;iyosh<=nyosh;iyosh++){
   //-----------------------------------------------------------------------
   // Evolve extended system velocities
    double vsum = 0.0;
    for(int i=0;i<num_nhc;i++){
      double ffdot = (fNHC[i][0]*fNHC[i][0]*mNHC[0]
	             +fNHC[i][1]*fNHC[i][1]*mNHC[1])/kT;
      double pfdot = (fNHC[i][0]*vNHC[i][0]*mNHC[0]
		     +fNHC[i][1]*vNHC[i][1]*mNHC[1])/kT;
      double arg2  = ffdot*wdti4[iyosh]*wdti4[iyosh];
      double h, hdot;
      if(arg2 > 1.0e-10){
        double ffdot_r = sqrt(ffdot);
        double arg     = ffdot_r*wdti4[iyosh];
        double cosha   = cosh(arg);
        double sinha   = sinh(arg);
        h    = (pfdot/ffdot)*(cosha-1.0) + (1.0/ffdot_r)*sinha;
        hdot = (pfdot/ffdot_r)*sinha + cosha;
      }else{
        h    = (((ffdot*pfdot/24.0*wdti4[iyosh]+ffdot/6.0)*wdti4[iyosh]
                      +0.50*pfdot)*wdti4[iyosh]+1.0)*wdti4[iyosh];
        hdot = ((ffdot*pfdot/6.0*wdti4[iyosh]  +ffdot/2.0)*wdti4[iyosh]
  		            +pfdot)*wdti4[iyosh]+1.0;
      }//endif
      vNHC[i][0] = (vNHC[i][0]+fNHC[i][0]*h)/hdot;
      vNHC[i][1] = (vNHC[i][1]+fNHC[i][1]*h)/hdot;
      vsum += vNHC[i][0];
    }//endfor
   //-----------------------------------------------------------------------
   // Evolve coef velocities
    double arg  = -wdti4[iyosh]*vsum;
    double aa   = exp(arg);
    scale       *= aa;
    ekin        *= (aa*aa);
   //-----------------------------------------------------------------------
   // Evolve extended system positions
    for(int i=0;i<num_nhc;i++){
      xNHC[i][0]  += kT*vNHC[i][0]*wdti2[iyosh];
      xNHC[i][1]  += vNHC[i][0]*(ekin-gkt)*wdti2[iyosh];
      xNHCP[i][1] += vNHC[i][1]*wdti2[iyosh];
    }//endfor
   //-----------------------------------------------------------------------
   // Evolve coef velocities
    scale       *= aa;
    ekin        *= (aa*aa);
   //-----------------------------------------------------------------------
   // Compute new forces
    for(int i=0;i<num_nhc;i++){
      fNHC[i][0] = (ekin - gkt)/mNHC[0];
      fNHC[i][1] = -( 4.0*v0[i]/(a4[i]*mNHC[1]) )
                   *( xNHCP[i][1]*(xNHCP[i][1]*xNHCP[i][1]-a2[i] ));
    }//endfor
   //-----------------------------------------------------------------------
   // Evolve extended system velocities
    for(int i=0;i<num_nhc;i++){
      double ffdot = (fNHC[i][0]*fNHC[i][0]*mNHC[0]
	             +fNHC[i][1]*fNHC[i][1]*mNHC[1])/kT;
      double pfdot = (fNHC[i][0]*vNHC[i][0]*mNHC[0]
		     +fNHC[i][1]*vNHC[i][1]*mNHC[1])/kT;
      double arg2  = ffdot*wdti4[iyosh]*wdti4[iyosh];
      double h, hdot;
      if(arg2 > 1.0e-10){
        double ffdot_r = sqrt(ffdot);
        double arg     = ffdot_r*wdti4[iyosh];
        double cosha   = cosh(arg);
        double sinha   = sinh(arg);
        h    = (pfdot/ffdot)*(cosha-1.0) + (1.0/ffdot_r)*sinha;
        hdot = (pfdot/ffdot_r)*sinha + cosha;
      }else{
        h    = (((ffdot*pfdot/24.0*wdti4[iyosh]+ffdot/6.0)*wdti4[iyosh]
                      +0.50*pfdot)*wdti4[iyosh]+1.0)*wdti4[iyosh];
        hdot = ((ffdot*pfdot/6.0*wdti4[iyosh]  +ffdot/2.0)*wdti4[iyosh]
  		            +pfdot)*wdti4[iyosh]+1.0;
      }//endif
      vNHC[i][0] = (vNHC[i][0]+fNHC[i][0]*h)/hdot;
      vNHC[i][1] = (vNHC[i][1]+fNHC[i][1]*h)/hdot;
    }//endfor
  }}//endfor :: yoshida, respa loops

//==========================================================================
// Apply the scaling factor to the coef velocities 

   for(int i=0;i<n;i++){
     v[i].re *=scale; v[i].im *=scale;
   }//endfor

//==========================================================================
// Avoid roundoff error : Apply constraint exactly

   for(int i=0;i<num_nhc;i++){
     double ekinNHC = mNHC[0]*vNHC[i][0]*vNHC[i][0]
                    + mNHC[1]*vNHC[i][1]*vNHC[i][1];
#ifdef _CHECKME_CPDYN_
     PRINTF("bot: ekinNHC %g %g %g\n",ekinNHC,kT,ekin);
#endif
     double scNHC   = sqrt(kT/ekinNHC);
     vNHC[i][0]    *= scNHC;
     vNHC[i][1]    *= scNHC;
   }//endfor

//============================================================================
  }//end routine
//============================================================================

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void CPINTEGRATE::set_yosh(int nyosh,double dt,double *wdt,double *wdt2,
                             double *wdt4,double *wdt8)
//========================================================================
  {//begin routine
//========================================================================
// Find the yoshida steps you want

  double temp,p2,onethird;
  switch(nyosh){
    case 1: 
          wdt[1] = 1.0;
         break;
    case 3:
          onethird = 1.0/3.0;
          temp    = pow(2.0,onethird);
          wdt[1] =  1.0/(2.0-temp);
          wdt[2] = -temp/(2.0-temp);
          wdt[3] =  1.0/(2.0-temp);
         break;
    case 5:
          onethird = 1.0/3.0;
          temp = pow(4.0,onethird);
          p2    = 1.0/(4.0-temp);
          wdt[1]  = p2;
          wdt[2]  = p2;
          wdt[3]  = 1.0-4.0*p2;
          wdt[4]  = p2;
          wdt[5]  = p2;
         break;
    case 7:
          wdt[1] =  0.784513610477560;
          wdt[2] =  0.235573213359357;
          wdt[3] = -1.17767998417887;
          wdt[4] =  1.0 - 2.0*(wdt[1]+wdt[2]+wdt[3]);
          wdt[5] = -1.17767998417887;
          wdt[6] =  0.235573213359357;
          wdt[7] =  0.784513610477560;
         break;
    case 9:
          wdt[1] =  0.192;
          wdt[2] =  0.554910818409783619692725006662999;
          wdt[3] =  0.124659619941888644216504240951585;
          wdt[4] = -0.843182063596933505315033808282941;
          wdt[5] =  1.0 - 2.0*(wdt[1]+wdt[2]+wdt[3]+wdt[4]);
          wdt[6] = -0.843182063596933505315033808282941;
          wdt[7] =  0.124659619941888644216504240951585;
          wdt[8] =  0.554910818409783619692725006662999;
          wdt[9] =  0.192;
         break;
  }//switch

//========================================================================
// Scale the yosida steps by the time step

  for(int i=1;i<=nyosh;i++){
    wdt[i]  = wdt[i]*dt;
    wdt2[i] = wdt[i]*0.5;
    wdt4[i] = wdt[i]*0.25;
    wdt8[i] = wdt[i]*0.125;
  }//endfor

//-------------------------------------------------------------------------
   }// end routine
//==========================================================================





