#include "standard_include.h"
#include "ckcomplex.h"

#include "../class_defs/allclass_gen.h"
#include "../class_defs/CP_OPERATIONS/class_cpxcfnctls.h"
#include "../class_defs/allclass_cp.h"


//===========================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===========================================================================

void CPXCFNCTS::CP_exc_calc(
              const int numFFT, const int nf1, const int nf2, const int nf3,
              double *density,double *result,double *exc_ret,double *muxc_ret)

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

  int cp_lyp   = cpopts->cp_lyp;

  static const double gamma = -0.14230;
  static const double beta1 = 1.05290;
  static const double beta2 = 0.33340;
  static const double au    = 0.03110;
  static const double bu    = -0.0480;
  static const double cu    = 0.00200;
  static const double du    = -0.01160;

//--------------------------------------------------------------------------
// Local variables 

  double vol            = gencell->vol;
  double pvx,pvc,vscale;
  double pi,rho_r_val,power,xfact,cfact,rs,fpin,fpi;
  double ex,ec,mufact,rat1,rat2,cf1,cf2,sqtrs,vxc;
  double lnrs;

  double lyp_fact=1.0;
  if(cp_lyp==1)lyp_fact=0.0;


//=================================================================
// EVALUATE THE FUNCTIONAL AND ITS DERIVATIVE
//-----------------------------------------------------------------
// I. Assign use constants                                         

   pi     = M_PI;
   fpi    = 4.0*pi;
   rat1   = (7.0/6.0)*beta1;
   rat2   = (4.0/3.0)*beta2;
   power    = 1.0/3.0;
   ex     = 0.0;
   ec     = 0.0;
   pvx    = 0.0;
   pvc    = 0.0;
   vxc    = 0.0;

//-------------------------------------------------------------------------
// Expressions for the exchange and correlation energies
//
//     ex = -sum_i rho[i]*(3*rho[i]/pi)**{1/3}
//
//  for rs < 1 :
//     ec = sum_i rho[i]*(a*ln(rs) + b + c*rs*ln(rs) + d*rs)
//
//          where rs = [3/(4*pi*rho[i])]**{1/3}
//
//  for rs >= 1 :
//     ec = sum_i rho[i]*(gamma/(1 + beta1*(rs)**{1/2} + beta2*rs))
//
//
//-------------------------------------------------------------------------
// II.  Start loop over FFT grid and evaluate terms in the functional

   for(int y=0 ; y< nf2; y++){
     int i = y*(nf1+2) + nf1; 
     result[i]    =0.0;
     result[(i+1)]=0.0;
   }//endfor

   int nfreq_cmi_update = 8;
   for(int y=0 ; y< nf2; y++){
   for(int x=0 ; x< nf1; x++){
      int i = y*(nf1+2) + x;
#ifdef CMK_VERSION_BLUEGENE
      if((i+1)%nfreq_cmi_update==0){CmiNetworkProgress();}
#endif   
//-------------------------------------------------------------------------
// III. Exchange part
      rho_r_val = density[i];
      fpin = fpi*rho_r_val;
      rs   = pow((3.0/fpin),power);
      sqtrs = sqrt(rs);
      xfact = -pow(((3.0*rho_r_val/pi)),power);

//-------------------------------------------------------------------------
// IV. Correlation part

    if(rs >= 1.0){  /* rs > 1  */ 

       cf1 = 1.0 + beta1*sqtrs + beta2*rs;
       cf2 = 1.0 + rat1*sqtrs + rat2*rs;
       cfact = gamma/cf1;
       mufact = cfact*(cf2/cf1);

    }else{ /* rs < 1  */

      lnrs = log(rs);
      cfact = au*lnrs + bu + cu*rs*lnrs + du*rs;
      mufact = au*lnrs + (bu - power*au) + 2.0*power*cu*rs*lnrs
            + power*(2.0*du - cu)*rs;

    }/*endif*/
    
    cfact  *= lyp_fact;
    mufact *= lyp_fact;

//-------------------------------------------------------------------------
// V. Energies
    ex += xfact*rho_r_val;
    ec += cfact*rho_r_val;

//-------------------------------------------------------------------------
// VI. Potential contributions


    result[i] = (xfact + mufact);


    vxc += (xfact + mufact)*rho_r_val;


  }}/* endfor */

//-------------------------------------------------------------------------
// VII. Finish the energies

   ex *= 0.75;

   vscale = vol/((double)numFFT);

//-------------------------------------------------------------------------
// VIII.  Return Values

   (*exc_ret)  = (ex+ec)*vscale;
   (*muxc_ret) = vxc*vscale;
  
//============================================================================
  } /* End */
//============================================================================
