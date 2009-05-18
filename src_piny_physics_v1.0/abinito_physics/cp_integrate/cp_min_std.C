#include "standard_include.h"
#include "ckcomplex.h"

#include "../../../include/debug_flags.h"
#include "../class_defs/allclass_gen.h"
#include "../class_defs/allclass_cp.h"

#include "../class_defs/CP_OPERATIONS/class_cpintegrate.h"


//============================================================================
// Steepest decent minimization
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================

void CPINTEGRATE::CP_integrate_min_STD
              (int ncoef,int istate,complex *forces,complex *psi_g,
               int *k_x, int *k_y, int *k_z,
               double *cmass, int nfreq_cmi_update)

//============================================================================
   { /* Begin Function */
//----------------------------------------------------------------------------
// Local Variables

  GENERAL_DATA *general_data = GENERAL_DATA::get();
  CP *cp = CP::get();
#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_cp.h"

   double gx, gy, gz, g2;

   double *hmati    = gencell->hmati;
   double ecut      = cpcoeffs_info->ecut_psi;
   double tpi       = 2.0*M_PI; 
   double dt        = gentimeinfo->dt;
//----------------------------------------------------------------------------
// I. Perform a step of steepest descent minimization
#ifdef _CP_DEBUG_NEWFORCE_
  FILE *fp  = fopen("force_new.out","a+");
  FILE *fp2 = fopen("coef_before_std.out","a+");
  FILE *fp3 = fopen("coef_after_std.out","a+");
  FILE *fp4;
  FILE *fp5;
  if(istate==0 && ncoef > 0){
    fp4 = fopen("coef_after_std_s0.out","a+");
    fp5 = fopen("force_new_s0.out","a+");
  }
#endif

#ifdef GJM_DEBUG_SIZE
   PRINTF(" %d : coefs in CP_integrate_min_STD\n",ncoef);
#endif

   for(int i = 0; i < ncoef; i++){

     gx = tpi*(k_x[i]*hmati[1] + k_y[i]*hmati[2] + k_z[i]*hmati[3]);
     gy = tpi*(k_x[i]*hmati[4] + k_y[i]*hmati[5] + k_z[i]*hmati[6]);
     gz = tpi*(k_x[i]*hmati[7] + k_y[i]*hmati[8] + k_z[i]*hmati[9]);
     g2 = gx*gx + gy*gy + gz*gz;

#ifdef _CP_DEBUG_NEWFORCE_
     if(k_x[i]==0 && k_y[i]==1 && k_z[i]==4){
        fprintf(fp,"new force H+Ext+Exc+Eke+Enl : is=%d %d %d %d : %g %g\n",
               istate,k_x[i],k_y[i],k_z[i],forces[i].re,forces[i].im);
        fprintf(fp2,"is=%d : %d %d %d : coef_before_std %.10g %.10g\n",
               istate,k_x[i],k_y[i],k_z[i],psi_g[i].re,psi_g[i].im);
     }//endif
#endif

       psi_g[i] += (forces[i]*(dt/cmass[i]));

#ifdef _CP_DEBUG_NEWFORCE_
     if(k_x[i]==0 && k_y[i]==1 && k_z[i]==4){
        fprintf(fp3,"is=%d : %d %d %d : coef_after_std %.10g %.10g\n",
 		istate,k_x[i],k_y[i],k_z[i],psi_g[i].re,psi_g[i].im);
     }//endif
     if(istate==0 && g2<=ecut){
        fprintf(fp4,"is=%d : %d %d %d : coef_after_std %.10g %.10g %.10g : dt %g\n",
 		istate,k_x[i],k_y[i],k_z[i],psi_g[i].re,psi_g[i].im,cmass[i],dt);
        fprintf(fp5,"%d %d %d : forc_after_lmatrix %.10g %.10g \n",
 		  k_x[i],k_y[i],k_z[i],forces[i].re,forces[i].im);
     }
#endif
#ifdef CMK_BLUEGENEL
     if(i%nfreq_cmi_update==0){CmiNetworkProgress();}
#endif
   } /* endfor */

#ifdef _CP_DEBUG_NEWFORCE_
    fclose(fp);
    fclose(fp2);
    fclose(fp3);
    if(istate==0 && ncoef > 0){
      fclose(fp4);
      fclose(fp5);
    }
#endif

//============================================================================
  } /* End function */
//============================================================================
