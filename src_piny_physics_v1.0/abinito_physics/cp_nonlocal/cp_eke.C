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
                  int *k_x, int *k_y, int *k_z, double **g2, double *eke_ret,
   	          int mydoublePack,int nkx0, int kpt, int nfreq_cmi_update)

//============================================================================
  { /* Begin Function */
//----------------------------------------------------------------------------
// Local Variables

  GENERAL_DATA *general_data = GENERAL_DATA::get();
  CP           *cp           = CP::get();

#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_cp.h"

   double ecut      = cpcoeffs_info->ecut_psi; // KS-state cutoff in Ryd
   double *occ      = cpcoeffs_info->occ_up;   // 
   double *wght_kpt = cpcoeffs_info->wght_kpt;
					       // Occupation numbers from file
   double *hmati   = gencell->hmati;
   double tpi = 2.0*M_PI;

//============================================================================
// Debugging

#ifdef GJM_DEBUG_SIZE
   PRINTF(" %d : coefs in CP_eke_calc\n",ncoef);
#endif

#ifdef _CP_DEBUG_OLDFORCE_
   FILE *fp;
   if(istate==0 && ncoef>0){fp = fopen("force_old_s0.out", "a+");}
#endif

   if(mydoublePack==0&&nkx0!=0){
     PRINTF("non-doublePack eke is broken!!\n");
     EXIT(1);
   }//endif

//============================================================================
// I. Forces and energy (kinetic energy contribution) : gx=0

   double occ_now      = occ[istate+1];
   double wght_kpt_now = wght_kpt[kpt+1];
   double wght_tot     = occ_now*wght_kpt_now;

   double eke = 0.0;
   for(int i = 0; i < nkx0; i++){
     if(g2[kpt][i]<=ecut){
       forces[i] -= (psi_g[i]*(wght_tot*g2[kpt][i])); 
       eke       += (wght_tot*g2[kpt][i])*psi_g[i].getMagSqr();
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
#ifdef CMK_BLUEGENEL
     if(i%nfreq_cmi_update==0){CmiNetworkProgress();}
#endif
   }/* endfor */
#ifdef CMK_BLUEGENEL
       CmiNetworkProgress();
#endif

//============================================================================
// I. Forces and energy (kinetic energy contribution) : gx>=0

   double wght = 1.0;  if(mydoublePack==1){wght = 2.0;}
   wght_tot *= wght;
   for(int i = nkx0; i < ncoef; i++){
#ifdef _DEBUG_CACHE_
     double aka = (double)(k_x[i]);
     double akb = (double)(k_y[i]);
     double akc = (double)(k_z[i]);
     double xk  = (aka*hmati[1]+akb*hmati[2]+akc*hmati[3])*tpi;
     double yk  = (aka*hmati[4]+akb*hmati[5]+akc*hmati[6])*tpi;
     double zk  = (aka*hmati[7]+akb*hmati[8]+akc*hmati[9])*tpi;
     double g2p = xk*xk+yk*yk+zk*zk;
     if(fabs(g2p-g2[kpt][i])>1.e-10){
       PRINTF("oops %g %g %d %d\n",g2p,g2[kpt][i],i,kpt);
       EXIT(1);
     }//endif
#endif
     if(g2[kpt][i]<=ecut){
       forces[i] -= (psi_g[i]*(wght_tot*g2[kpt][i])); 
       eke       += (wght_tot*g2[kpt][i])*psi_g[i].getMagSqr();
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
#ifdef CMK_BLUEGENEL
     if(i%nfreq_cmi_update==0){CmiNetworkProgress();}
#endif
   }/* endfor */
#ifdef CMK_BLUEGENEL
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
