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

void CPXCFNCTS::CP_getGGAFunctional(
          const int npts, const int nf1,const int nf2,const int nf3,
          complex *density,complex *rhoIRX, complex *rhoIRY, complex *rhoIRZ, 
          complex *gradientCorrection,int iyPlane_ind,double *exc_gga_ret)

//============================================================================
{/* Begin function */
//============================================================================
// Local variables

  GENERAL_DATA *general_data = GENERAL_DATA::get();
  CP           *cp           = CP::get();

#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_cp.h"
    static double gc_cut = 1.0e-6;/* Cutoff for evaluating (grad rho/rho) */
    static double rho_heali;      /* Inverse of healing length on switching fn */
    double srho,dsrho;            /* Switching function and its derivative */
    double rho,rho_cut,g_rho2,g_rhoi,rsw;
    int m;

    double ex,ec,dfx_drho,fx,fc,dfx_dgrho,dfc_drho,dfc_dgrho,dfxc_dgrho;
    double dfx_drho_tmp;
    static double beta = 0.0042;
    double unit_gx,unit_gy,unit_gz;

    double vol    = gencell->vol;
    double vscale = vol/((double)(nf1*nf2*nf3));
    
//============================================================================
// Initialization

    ex = ec = 0.0;

//============================================================================
// Some useful constants

    rho_heali = 1.0/(3.9*gc_cut);
    rho_cut = 0.1*gc_cut;

//============================================================================
// Start Loop over part of grid I have 

   for(int i = 0; i < npts; i++){

//----------------------------------------------------------
// Get density and decide if the loop should be done

        rho   =  density[i].re;
        if(rho > rho_cut){

//----------------------------------------------------------
// Calculate magnitude square of gradient of density
         //FIX later in fft!!!

         rhoIRX[i].re /=  vol;
         rhoIRY[i].re /=  vol;
         rhoIRZ[i].re /=  vol;

          g_rho2 = rhoIRX[i].re*rhoIRX[i].re 
                 + rhoIRY[i].re*rhoIRY[i].re 
                 + rhoIRZ[i].re*rhoIRZ[i].re;

//----------------------------------------------------------
// Calculate the switching function and its derivative

         rsw = (rho - rho_cut)*rho_heali;
         rsw = (rsw < 1.0 ? rsw : 1.0);
         rsw = (rsw > 0.0 ? rsw : 0.0);
         srho = rsw*rsw*(3.0-2.0*rsw);
         dsrho = 6.0*rsw*(1.0-rsw)*rho_heali;

//----------------------------------------------------------
// Calculate the exchange functional

         becke_gcx_lda(rho,g_rho2,&fx,&dfx_drho,&dfx_dgrho,beta);

//----------------------------------------------------------
// Process the output via the switching function

         dfx_drho = (dfx_drho*srho + fx*dsrho);
         dfx_dgrho *= srho;
         fx *= srho;

//----------------------------------------------------------
// Calculate the correlation functional

         /*  Nothing for now */

         dfc_drho = 0.0;
         dfc_dgrho = 0.0;
         fc = 0.0;

//----------------------------------------------------------
// Finish the energy

         ex += fx;
         ec += fc;

//----------------------------------------------------------
// Construct output for the FFT

         gradientCorrection[i] = complex((dfx_drho + dfc_drho));

         g_rhoi = 1.0/sqrt(g_rho2);
         unit_gx = rhoIRX[i].re*g_rhoi;
         unit_gy = rhoIRY[i].re*g_rhoi;
         unit_gz = rhoIRZ[i].re*g_rhoi;

         dfxc_dgrho = dfx_dgrho + dfc_dgrho;
         rhoIRX[i] = complex(unit_gx*dfxc_dgrho);
         rhoIRY[i] = complex(unit_gy*dfxc_dgrho);
         rhoIRZ[i] = complex(unit_gz*dfxc_dgrho);

        } else {  /* Density is less than cutoff */

          gradientCorrection[i] = complex(0.0,0.0);
          rhoIRX[i] = complex(0.0,0.0);
          rhoIRY[i] = complex(0.0,0.0);
          rhoIRZ[i] = complex(0.0,0.0);

        } /* Endif density is greater than cutoff */
         

//============================================================================
// End Loop over grid

}
//============================================================================
// Store the energy

   (*exc_gga_ret) = (ex+ec)*vscale;

//============================================================================
}/* end function */
//============================================================================



