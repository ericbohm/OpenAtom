#include "standard_include.h"
#include "ckcomplex.h"

#include "../class_defs/CP_OPERATIONS/class_cpxcfnctls.h"
#include "../class_defs/allclass_gen.h"
#include "../class_defs/allclass_cp.h"

//============================================================================
/* The following function takes the gradient of the density as well as the 
*  density, itself, and evaluates the GGA functional (BLYP in this case).
*  
*     E_GGA = sum_r f(rho(r),|grad rho(r)|)
*
*     dE_GGA/drho(r) = df/drho(r)
*
*     dE_GGA/d|grad rho(r)| = df/d|grad rho(r)|
*  
*  It, then, takes df/d|grad rho(r)| and multiplies by grad rho(r)/|grad rho(r)|
*  the result of which is stored in rhoIRX, rhoIRY and rhoIRZ.  The df/drho(r) term
*  is stored in gradientCorrection
*/
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CPXCFNCTS::CP_getGGAFunctional(
          const int npts, const int nf1,const int nf2,const int nf3,
          double *density,double *rhoIRX, double *rhoIRY, double *rhoIRZ, 
          double *Vks,int iyPlane_ind,double *exc_gga_ret,int nfreq_cmi_update)
//============================================================================
   {/* Begin function */
//============================================================================
// Local variables

  GENERAL_DATA *general_data = GENERAL_DATA::get();
  CP           *cp           = CP::get();

#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_cp.h"

    double srho,dsrho;            /* Switching function and its derivative */
    double rho,rho_cut,g_rho2,g_rhoi,rsw;
    int m;
    int cp_becke = cpopts->cp_becke;
    int cp_lyp   = cpopts->cp_lyp;

    double ex,ec,dfx_drho,fx,fc,dfx_dgrho,dfc_drho,dfc_dgrho,dfxc_dgrho;
    double unit_gx,unit_gy,unit_gz;

    static const double beta = 0.0042;

    double gc_cut = cppseudo->gga_cut;
    double vol    = gencell->vol;
    double vscale = vol/((double)npts);
    
//============================================================================
// Initialization and useful costants

    ex = ec = 0.0;
/* Inverse of healing length on switching fn */
    static const double rho_heali = 1.0/(3.9*gc_cut);
    rho_cut   = 0.1*gc_cut;

#ifdef _CPDEBUG_XC_FUNCTIONALS_
    char fname[1000];
    sprintf(fname,"gradcorr.%d.out",iyPlane_ind);
    FILE *fp = fopen(fname,"w");
    fprintf(fp,"vscale %g rho_cut %g\n",vscale,rho_cut);
#endif

//============================================================================
// Start Loop over part of grid I have 

   for(int y = 0; y < nf2; y++){
   for(int x = 0; x < nf1; x++){

     int i = y*(nf1+2) + x;
//----------------------------------------------------------
// Get density and decide if the loop should be done

     rho   =  density[i];
#ifdef _CPDEBUG_XC_FUNCTIONALS_
     fprintf(fp,"rho %g",rho);
#endif

     if(rho > rho_cut){

//----------------------------------------------------------
// Calculate magnitude square of gradient of density

        g_rho2 = rhoIRX[i]*rhoIRX[i] 
               + rhoIRY[i]*rhoIRY[i] 
               + rhoIRZ[i]*rhoIRZ[i];

#ifdef _CPDEBUG_XC_FUNCTIONALS_
       fprintf(fp," rhoi %g %g %g grho2 %g",rhoIRX[i],rhoIRY[i],rhoIRZ[i],g_rho2);
#endif

//----------------------------------------------------------
// Calculate the switching function and its derivative

        rsw = (rho - rho_cut)*rho_heali;
        rsw = (rsw < 1.0 ? rsw : 1.0);
        rsw = (rsw > 0.0 ? rsw : 0.0);
        srho = rsw*rsw*(3.0-2.0*rsw);
        dsrho = 6.0*rsw*(1.0-rsw)*rho_heali;

//----------------------------------------------------------
// Calculate the exchange functional

       fx  = 0.0;  dfx_drho  = 0.0;  dfx_dgrho = 0.0;
       if(cp_becke==1){
         becke_gcx_lda(rho,g_rho2,&fx,&dfx_drho,&dfx_dgrho,beta);
       }//endif

       fc = 0.0;  dfc_drho  = 0.0;  dfc_dgrho = 0.0;
       if(cp_lyp==1){ 
         lyp_gcc(rho,g_rho2,&fc,&dfc_drho,&dfc_dgrho);
       }//endif

//----------------------------------------------------------
// Process the output via the switching function

       dfx_drho = (dfx_drho*srho + fx*dsrho);
       dfx_dgrho *= srho;
       fx *= srho;

       dfc_drho = (dfc_drho*srho + fc*dsrho);
       dfc_dgrho *= srho;
       fc *= srho;

//----------------------------------------------------------
// Finish the energy

#ifdef _CPDEBUG_XC_FUNCTIONALS_
       fprintf(fp," fx %g fc %g",fx,fc);
#endif

       ex += fx;
       ec += fc;

//----------------------------------------------------------
// Construct output for the FFT

      Vks[i] += (dfx_drho + dfc_drho);

      g_rhoi = 1.0/sqrt(g_rho2);
      unit_gx = rhoIRX[i]*g_rhoi;
      unit_gy = rhoIRY[i]*g_rhoi;
      unit_gz = rhoIRZ[i]*g_rhoi;

      dfxc_dgrho = dfx_dgrho + dfc_dgrho;
      rhoIRX[i] = unit_gx*dfxc_dgrho;
      rhoIRY[i] = unit_gy*dfxc_dgrho;
      rhoIRZ[i] = unit_gz*dfxc_dgrho;

   } else {  /* Density is less than cutoff */

      rhoIRX[i] = 0.0;
      rhoIRY[i] = 0.0;
      rhoIRZ[i] = 0.0;

   } /* Endif density is greater than cutoff */
#ifdef _CPDEBUG_XC_FUNCTIONALS_
    fprintf(fp,"\n");
#endif
 }}//endfor

//============================================================================
// Store the energy

   (*exc_gga_ret) = (ex+ec)*vscale;

#ifdef _CPDEBUG_XC_FUNCTIONALS_
   fprintf(fp,"done %g\n",exc_gga_ret[0]);
   fclose(fp);
#endif

//============================================================================
   }/* end function */
//============================================================================
