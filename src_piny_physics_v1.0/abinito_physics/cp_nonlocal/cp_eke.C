#include "standard_include.h"
#include "ckcomplex.h"

#include "../../../include/debug_flags.h"
#include "../../../include/Atoms.h"

#include "../class_defs/allclass_cp.h"
#include "../class_defs/allclass_gen.h"

#include "../class_defs/CP_OPERATIONS/class_cpnonlocal.h"
#include "../class_defs/allclass_gen.h"
#include "../class_defs/allclass_cp.h"

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================

void
CPNONLOCAL::CP_eke_calc
         (const int ncoef, const int istate,complex *forces,complex *psi_g,
          const int *k_x, const int *k_y, const int *k_z,double *eke_ret,
          int mydoublePack)

//============================================================================
  { /* Begin Function */
//----------------------------------------------------------------------------
// Local Variables

  GENERAL_DATA *general_data = GENERAL_DATA::get();
  CP           *cp           = CP::get();

#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_cp.h"

   double eke;
   double gx, gy, gz, g2;
   double *hmati    = gencell->hmati;
   double ecut      = cpcoeffs_info->ecut_psi; // KS-state cutoff in Ryd
   double tpi       = 2.0*M_PI; 
   double wght,wght_now;

   FILE *fp;

//----------------------------------------------------------------------------
// I. Finish the forces and energy (kinetic energy contribution)

#ifdef _CP_DEBUG_OLDFORCE_
   if(istate==0 && ncoef>0){
     fp = fopen("force_old_s0.out", "a+");
   }
#endif

#ifdef GJM_DEBUG_SIZE
   PRINTF(" %d : coefs in CP_eke_calc\n",ncoef);
#endif

   wght = 1.0;  if(mydoublePack==1){wght = 2.0;}

   eke = 0.0;
   for(int i = 0; i < ncoef; i++){
     
     gx = tpi*(k_x[i]*hmati[1] + k_y[i]*hmati[2] + k_z[i]*hmati[3]);
     gy = tpi*(k_x[i]*hmati[4] + k_y[i]*hmati[5] + k_z[i]*hmati[6]);
     gz = tpi*(k_x[i]*hmati[7] + k_y[i]*hmati[8] + k_z[i]*hmati[9]);
     g2 = gx*gx + gy*gy + gz*gz;

     if(g2<=ecut){
       wght_now = (k_x[i]==0 ? 1.0 : wght);
       forces[i] -= (psi_g[i]*(g2*wght_now)); 
       eke       += (wght_now*g2)*psi_g[i].getMagSqr();
#ifdef _CP_DEBUG_OLDFORCE_
       if(istate==0){
         fprintf(fp,"old force H+Ext+Exc+Eke+Enl : is=%d %d %d %d : %g %g\n",
		istate,k_x[i],k_y[i],k_z[i],forces[i].re,forces[i].im);
       }
#endif
     }else{
       forces[i].re = 0.0;
       forces[i].im = 0.0;
     }//endif

   }/* endfor */

#ifdef _CP_DEBUG_OLDFORCE_
   if(istate==0 && ncoef>0){
     fclose(fp);
   }
#endif

   (*eke_ret) = eke/2.0;

//============================================================================
} /* End function */
//============================================================================
