#include "standard_include.h"
#include "ckcomplex.h"

#include "../class_defs/allclass_gen.h"
#include "../class_defs/CP_OPERATIONS/class_cpxcfnctls.h"
#include "../class_defs/allclass_cp.h"


//===========================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===========================================================================

void CPXCFNCTS::CP_exc_pw_calc(
              const int numFFT, const int nf1, const int nf2, const int nf3,
              double *density,double *result,double *exc_ret,double *muxc_ret,int nfreq_cmi_update)

//============================================================================
// Function:  Exchange-correlation functional
//
// NOTE FOR RAMKUMAR:  BOX VOLUME (vol) MUST BE PASSED IN AND 
//                     exc_ret AND muxc_ret MUST BE SENT OUT.
//                     I ALSO ASSUME result IS ZEROED SOMEWHERE SO THAT I
//                     CAN ACCUMULATE IT.
//============================================================================
{ /* Begin function */
//--------------------------------------------------------------------------
// Static variables        


  GENERAL_DATA *general_data = GENERAL_DATA::get();
  CP           *cp           = CP::get();

#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_cp.h"

//--------------------------------------------------------------------------
// Local variables 

  double vol            = gencell->vol;
  double vscale;
  double pi,rho_r,power,xfact,cfact,rs,fpin,fpi;
  double ex,ec,mufact,sqtrs,vxc;
  double cf1,cf2,logterm,opars,tc1,tc2;
  double rat1,rat3,rat4;
  double mtwoa;

  int cp_lyp   = cpopts->cp_lyp;

  /*pw_lda correlation parameters */
  static double a      = 0.03109070;
  static double alpha1 = 0.213700;
  static double beta1  = 7.59570;
  static double beta2  = 3.58760;
  static double beta3  = 1.63820;
  static double beta4  = 0.492940;

  if(cp_lyp==1){
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     PRINTF("Don't use LYP with PW_lda!\n");
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     EXIT(1);
  }//endif

//=================================================================
// EVALUATE THE FUNCTIONAL AND ITS DERIVATIVE
//-----------------------------------------------------------------
// I. Assign use constants                                         

   pi     = M_PI;
   fpi    = 4.0*pi;
   power    = 1.0/3.0;

   rat1  =  0.50*beta1;
   rat3  =  1.50*beta3;
   rat4  =  2.00*beta4;
   mtwoa = -2.00*a;

   ex     = 0.0;
   ec     = 0.0;
   vxc    = 0.0;

//-------------------------------------------------------------------------
// II.  Start loop over FFT grid and evaluate terms in the functional

   for(int y=0 ; y< nf2; y++){
     int i = y*(nf1+2) + nf1; 
     result[i]    =0.0;
     result[(i+1)]=0.0;
   }//endfor

   for(int y=0 ; y< nf2; y++){
   for(int x=0 ; x< nf1; x++){
      int i = y*(nf1+2) + x;
#ifdef CMK_BLUEGENEL
      if((i+1)%nfreq_cmi_update==0){CmiNetworkProgress();}
#endif   
//-------------------------------------------------------------------------
// III. Exchange part

      rho_r = density[i];
      fpin = fpi*rho_r;
      rs   = pow((3.0/fpin),power);
      sqtrs = sqrt(rs);

      xfact = pow((3.0*rho_r/pi),power);
      xfact = -xfact;

//-------------------------------------------------------------------------
// IV. Correlation part

      cf1 = 2.0*a*(beta1*sqtrs + beta2*rs +
		   beta3*pow(rs,1.50) + beta4*rs*rs);
      cf2 = 2.0*a*(rat1/sqtrs + beta2 +
		   rat3*sqtrs + rat4*rs);
      logterm = log(1.0 + 1.0/cf1);
      opars = 1.0 + alpha1*rs;
      cfact = mtwoa*opars*logterm;

      tc1 = mtwoa*(1.0 + 2.0*alpha1*rs/3.0);
      tc2 = mtwoa*rs*opars/(3.0*cf1*(1.0 + cf1));
      mufact = tc1*logterm + tc2*cf2;

//-------------------------------------------------------------------------
// V. Energies

    ex += xfact*rho_r;
    ec += cfact*rho_r;

//-------------------------------------------------------------------------
// VI. KS potential contributions

    result[i] = (xfact + mufact);
    vxc      += (xfact + mufact)*rho_r;

  }}/* endfor */

//-------------------------------------------------------------------------
// VII. Finish the energies

   vscale = vol/((double)numFFT);

   ex *= 0.750;

//-------------------------------------------------------------------------
// VIII.  Return Values

   (*exc_ret)  = (ex+ec)*vscale;
   (*muxc_ret) = vxc*vscale;

//============================================================================
} /* End */
//============================================================================
