//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
/** \file name cp_pz_lda.C
 ** \brief Compute Purdue Zunger exchange correlation
 */
//==========================================================================
#include "standard_include.h"
#include "ckcomplex.h"

#include "../class_defs/allclass_gen.h"
#include "../class_defs/CP_OPERATIONS/class_cpxcfnctls.h"
#include "../class_defs/allclass_cp.h"

//==========================================================================
/* 
 ** \brief compute Purdue Zunger exchange correlation
 */
//===========================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===========================================================================

void CPXCFNCTS::CP_exc_calc(const int numFFT, const int nf1, const int nf2, 
    const int nf3, double *density, double *result, double *exc_ret, 
    double *muxc_ret, int nfreq_cmi_update)

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
  double pvx, pvc, vscale;
  double pi, rho_r_val, power, xfact, cfact, rs, fpin, fpi;
  double ex, ec, mufact, rat1, rat2, cf1, cf2, sqtrs, vxc;
  double lnrs;

  double lyp_fact = 1.0;
  if(cp_lyp == 1) lyp_fact = 0.0;

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

  int x_len = 2 * (nf1/2 + 1);
  for(int z = 0; z < nf3; z++) {
    for(int y = 0 ; y < nf2; y++) {
      int i = z * nf2 * x_len +  y * x_len + nf1; 
      result[i]     = 0.0;
      result[(i+1)] = 0.0;
    }
  }//endfor

  int jump = 2 - (nf1 & 1);
  int i = 0;
  for(int z = 0; z < nf3; z++) { 
    for(int y = 0 ; y < nf2; y++) {
      for(int x = 0 ; x < nf1; x++) {
        //-------------------------------------------------------------------------
        // III. Exchange part
        rho_r_val = density[i];
        fpin = fpi * rho_r_val;
        rs   = pow((3.0/fpin), power);
        sqtrs = sqrt(rs);
        xfact = -pow(((3.0*rho_r_val/pi)),power);

        //-------------------------------------------------------------------------
        // IV. Correlation part

        if(rs >= 1.0) {  /* rs > 1  */ 
          cf1 = 1.0 + beta1 * sqtrs + beta2 * rs;
          cf2 = 1.0 + rat1 * sqtrs + rat2 * rs;
          cfact = gamma/cf1;
          mufact = cfact * (cf2/cf1);
        } else { /* rs < 1  */
          lnrs = log(rs);
          cfact = au * lnrs + bu + cu * rs * lnrs + du * rs;
          mufact = au * lnrs + (bu - power * au) + 2.0 * power * cu * rs * lnrs
                   + power * (2.0 * du - cu) * rs;
        }/*endif*/

        cfact  *= lyp_fact;
        mufact *= lyp_fact;

        //-------------------------------------------------------------------------
        // V. Energies
        ex += xfact * rho_r_val;
        ec += cfact * rho_r_val;

        //-------------------------------------------------------------------------
        // VI. Potential contributions

        result[i] = (xfact + mufact);
        vxc += (xfact + mufact) * rho_r_val;
        i++;
      }
      i += jump;
    }
  }/* endfor */

  //-------------------------------------------------------------------------
  // VII. Finish the energies

  ex *= 0.75;
  vscale = vol / ((double)numFFT);

  //-------------------------------------------------------------------------
  // VIII.  Return Values

  (*exc_ret)  = (ex + ec) * vscale;
  (*muxc_ret) = vxc * vscale;

  //============================================================================
} /* End */
//============================================================================
//===========================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===========================================================================

void CPXCFNCTS::CP_exc_lsda_calc(const int nfft, const int nf1, 
                                 const int nf2,  const int nf3,
                                 double* density, double* density_dn, 
                                 double* Vks,     double* Vks_dn,
                                 double* exc_ret, double* muxc_ret, 
                                 int nfreq_cmi_update)

//============================================================================
// RAZ added routine:  this subroutine calculates the perdew-zunger lsda   
// exchange correlation function form the old prb      
// Grabbed from piny code on gmartyna@csmlogin:
// /farm/gmartyna/PINY_KPT/para3/energy/cp/old/xc_functionals.c
//============================================================================
// Function:  Exchange-correlation functional for spin based PZ
//
// NOTE FOR RAMKUMAR:  BOX VOLUME (vol) MUST BE PASSED IN AND 
//                     exc_ret AND muxc_ret MUST BE SENT OUT.
//                     I ALSO ASSUME result IS ZEROED SOMEWHERE SO THAT I
//                     CAN ACCUMULATE IT.
//============================================================================
{ /* Begin function */

  GENERAL_DATA *general_data = GENERAL_DATA::get();
  CP           *cp           = CP::get();

#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_cp.h"

  int cp_lyp   = cpopts->cp_lyp;
  int cp_lypm1 = cpopts->cp_lypm1;

//
//============================================================================
// Local Variables:

  double vol            = gencell->vol;

  double pi,rho_r,rho_ru,rho_rd,power,exc;
  double xfact,rs,fpi,fpin,zeta,xff,xfp,muxf,muxp,rat1,rat2;
  double ex,vxc,muxfact_up,muxfact_dn,ec,mucf,mucp,mucfact_up;
  double mucfact_dn,cfact,cff,cfp;
  double cff1,cfp1,cff2,cfp2,ratbf1,ratbf2;
  double ratbp1,ratbp2;
  double fzeta,fpzeta,opzeta,omzeta,rden,sqtrs;
  double pvx,pvc,vscale;
  double lyp_fact;
  int kk;
  int nfft2 = nfft/2;

//============================================================================
// lsda correlation parameters:

  static double cf = 0.57730;
  static double cpconst= 0.45820;
  static double gammaf = -0.08430;
  static double gammap = -0.14230;
  static double beta1f = 1.05290;
  static double beta1p = 1.39810;
  static double beta2f = 0.33340;
  static double beta2p = 0.26110;

//============================================================================
// I) useful constants:

  pi = M_PI;
  fpi = 4.0*pi;
  rat1 = 2.0/3.0;
  rat2 = 4.0/3.0;
  ratbf1 = (7.0/6.0)*beta1f;
  ratbf2 = (4.0/3.0)*beta2f;
  ratbp1 = (7.0/6.0)*beta1p;
  ratbp2 = (4.0/3.0)*beta2p;
  power = 1.0/3.0;
  rden = 1.0/( pow(2.0,rat2) - 2.0);
  lyp_fact = (double)(1-(cp_lyp || cp_lypm1));

  ex = 0.0;
  ec = 0.0;
  vxc = 0.0;
  pvx = 0.0;
  pvc = 0.0;

//============================================================================
// II) Loop over real space grid

  int jump = 2 - (nf1 & 1);
  int i = 0;
  for(int z = 0; z < nf3; z++) { 
    for(int y = 0 ; y < nf2; y++) {
      for(int x = 0 ; x < nf1; x++) {

#ifdef CMK_BLUEGENEL
      if((i+1)%nfreq_cmi_update==0){CmiNetworkProgress();}
#endif

      rho_ru = density[i];
      rho_rd = density_dn[i];
      rho_r = rho_ru + rho_rd;
      zeta = (rho_ru - rho_rd)/rho_r;
      opzeta = 1.0 + zeta;
      omzeta = 1.0 - zeta;
      fpin = fpi*rho_r;
      rs = pow((3.0/fpin),power);
      sqtrs = sqrt(rs);
      xff = -cf/rs;
      xfp = -cpconst/rs;
      muxf = rat2*xff;
      muxp = rat2*xfp;
      cff1 = 1.0 + beta1f*sqtrs + beta2f*rs; 
      cff2 = 1.0 + ratbf1*sqtrs + ratbf2*rs;
      cfp1 = 1.0 + beta1p*sqtrs + beta2p*rs;
      cfp2 = 1.0 + ratbp1*sqtrs + ratbp2*rs;
      cff = gammaf/cff1;
      cfp = gammap/cfp1;
      mucf = cff*(cff2/cff1);
      mucp = cfp*(cfp2/cfp1);
      fzeta = rden*(pow(opzeta,rat2) + pow(omzeta,rat2) - 2.0);
      fpzeta = rat2*rden*(pow(opzeta,power) - pow(omzeta,power));
      xfact = xfp + fzeta*(xff - xfp);
      cfact = cfp + fzeta*(cff - cfp);

      muxfact_up = muxp + fzeta*(muxf-muxp) +
        omzeta*fpzeta*(xff - xfp);
      muxfact_dn = muxp + fzeta*(muxf-muxp) -
        opzeta*fpzeta*(xff - xfp);

      mucfact_up = mucp + fzeta*(mucf-mucp) +
        omzeta*fpzeta*(cff - cfp);
      mucfact_dn = mucp + fzeta*(mucf-mucp) -
        opzeta*fpzeta*(cff - cfp);

      cfact = cfact*lyp_fact;
      ex += xfact*rho_r;
      ec += cfact*rho_r;
      mucfact_up *= lyp_fact;
      mucfact_dn *= lyp_fact;

      pvx +=  (muxfact_up*rho_ru + muxfact_dn*rho_rd - xfact*rho_r);
      pvc +=  (mucfact_up*rho_ru + mucfact_dn*rho_rd - cfact*rho_r);

      Vks[i]     = (muxfact_up + mucfact_up);
      Vks_dn[i]  = (muxfact_dn + mucfact_dn);
      vxc += (muxfact_up + mucfact_up)*rho_ru + (muxfact_dn + mucfact_dn)*rho_rd;
      i++;
      }//endfor x
      i+=jump;
    }//endfor y
  }// endfor z
//============================================================================
// III) Finish return energies: 
//

  vscale = vol/((double)nfft);
  
  *exc_ret  = (ex+ec)*vscale;
  *muxc_ret = vxc*vscale;

/*** Nuke this pressure tensor stuff for now:
  pvx *= vscale;
  pvc *= vscale;

  if(cp_ptens_calc==1){
    pvten_cp[1] += (pvx + pvc);
    pvten_cp[5] += (pvx + pvc);
    pvten_cp[9] += (pvx + pvc);
  }//endif
***/

//============================================================================
  } /* End */
//============================================================================

