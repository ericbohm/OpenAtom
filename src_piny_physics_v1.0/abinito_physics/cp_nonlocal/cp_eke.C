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

void CPNONLOCAL::CP_eke_calc(int ncoef, int istate,complex *forces,complex *psi_g,
                  int *k_x, int *k_y, int *k_z, double *g2, double *eke_ret,
                  int mydoublePack,int nkx0)

//============================================================================
  { /* Begin Function */
//----------------------------------------------------------------------------
// Local Variables

  GENERAL_DATA *general_data = GENERAL_DATA::get();
  CP           *cp           = CP::get();

#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_cp.h"

   double ecut      = cpcoeffs_info->ecut_psi; // KS-state cutoff in Ryd
   int nfreq = 400;

//============================================================================
// Debugging

#ifdef GJM_DEBUG_SIZE
   PRINTF(" %d : coefs in CP_eke_calc\n",ncoef);
#endif

#ifdef _CP_DEBUG_OLDFORCE_
   FILE *fp;
   if(istate==0 && ncoef>0){fp = fopen("force_old_s0.out", "a+");}
#endif

//============================================================================
// I. Forces and energy (kinetic energy contribution) : gx=0

   double eke = 0.0;
   for(int i = 0; i < nkx0; i++){
     if(g2[i]<=ecut){
       forces[i] -= (psi_g[i]*(g2[i])); 
       eke       += (g2[i])*psi_g[i].getMagSqr();
#ifdef _CP_DEBUG_OLDFORCE_
       if(istate==0){
         fprintf(fp,"old force H+Ext+Exc+Eke+Enl : is=%d %d %d %d : %g %g\n",
		istate,k_x[i],k_y[i],k_z[i],forces[i].re,forces[i].im);
       }//endif
#endif
     }else{
       forces[i].re = 0.0;
       forces[i].im = 0.0;
     }//endif
#ifdef CMK_VERSION_BLUEGENE
     if(i%nfreq==0){CmiNetworkProgress();}
#endif
   }/* endfor */
#ifdef CMK_VERSION_BLUEGENE
       CmiNetworkProgress();
#endif

//============================================================================
// I. Forces and energy (kinetic energy contribution) : gx>=0

   double wght = 1.0;  if(mydoublePack==1){wght = 2.0;}
   for(int i = nkx0; i < ncoef; i++){
     if(g2[i]<=ecut){
       forces[i] -= (psi_g[i]*(wght*g2[i])); 
       eke       += (wght*g2[i])*psi_g[i].getMagSqr();
#ifdef _CP_DEBUG_OLDFORCE_
       if(istate==0){
         fprintf(fp,"old force H+Ext+Exc+Eke+Enl : is=%d %d %d %d : %g %g\n",
		istate,k_x[i],k_y[i],k_z[i],forces[i].re,forces[i].im);
       }//endif
#endif
     }else{
       forces[i].re = 0.0;
       forces[i].im = 0.0;
     }//endif
#ifdef CMK_VERSION_BLUEGENE
     if(i%nfreq==0){CmiNetworkProgress();}
#endif
   }/* endfor */
#ifdef CMK_VERSION_BLUEGENE
       CmiNetworkProgress();
#endif

//============================================================================
// Energy return value

   eke_ret[0] = eke/2.0;

//============================================================================
// Complete Debugging

#ifdef _CP_DEBUG_OLDFORCE_
   if(istate==0 && ncoef>0){fclose(fp);}
#endif

//============================================================================
  } /* End function */
//============================================================================
