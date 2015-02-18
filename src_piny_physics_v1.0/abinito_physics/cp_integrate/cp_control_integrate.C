//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
/** \file name cp_control_integrate
 ** \brief The physics routines that evolve the wavefunctions (psi)
 for both minimization and CPAIMD
 */
//==========================================================================

#include "standard_include.h"
#include "ckcomplex.h"

#include "../class_defs/allclass_gen.h"
#include "../class_defs/allclass_cp.h"

#include "../class_defs/CP_OPERATIONS/class_cpintegrate.h"

//============================================================================
/** \brief Generic call for the charm_driver to invoke psi integration/evolution
 */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CPINTEGRATE::CP_integrate(int ncoef, int istate,int iteration,
    complex *forces,complex *forcesold,complex *psi_g,double *cmass,
    int *k_x,int *k_y, int *k_z, int len_nhc, int num_nhc, int nck_nhc,
    double ***fNHC,double***vNHC,double ***xNHC,double ***xNHCP,
    double *mNHC,double *v0NHC, double *a2NHC, double *a4NHC, double kTCP,
    double gamma_conj_grad,double *fictEke,int nkx0_red,int nkx0_uni,
    int nkx0_zero,double *ekeNHC,double *potNHC, double degfree,double degfreeNHC,
    double *degFreeSplt, int *istrNHC,int *iendNHC,int halfStepEvolve, int nfreq_cmi_update)
  //============================================================================
{ /* Begin Function */
  //---------------------------------------------------------------------------

  GENERAL_DATA *general_data = GENERAL_DATA::get();
#include "../class_defs/allclass_strip_gen.h"
  complex *vpsi_g = forcesold; //they share memory

  int ifound=0;

  int cp_min_on = gensimopts->cp_min
    + gensimopts->cp_wave_min_pimd 
    + gensimopts->cp_wave_min;

  int cp_on     = gensimopts->cp     
    + gensimopts->cp_pimd
    + gensimopts->cp_wave
    + gensimopts->cp_wave_pimd;

#define DEBUG_PARAMS_OFF
#ifdef DEBUG_PARAMS
  PRINTF("Integrate Stuff : %d \n",gs->size);
#endif

  //============================================================================
  // I) CP Minimization :

  if(cp_min_on==1){

    if(genminopts->cp_min_std==1){
      ifound++;
      CP_integrate_min_STD(ncoef,istate,forces,psi_g,k_x,k_y,k_z,cmass,nfreq_cmi_update);
    }//endif

    if(genminopts->cp_min_cg==1){
      ifound++;
      CP_integrate_min_CG(ncoef,istate,forces,forcesold,psi_g,k_x,k_y,k_z,cmass,
          gamma_conj_grad,nfreq_cmi_update);
      //      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      //      PRINTF("CP-CG wave function minimization not implemented");
      //      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      //      EXIT(1);
    }//endif

    if((gensimopts->cp_wave_min_pimd)==1){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("CP-PIMD wave function minimization not implemented");
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT(1);
    }//endif

    if((gensimopts->cp_min)==1){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("Full cp minimization not implemented");
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT(1);
    }//endif

    if((gensimopts->cp_bomd)==1){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("CP bomd not yet implemented");
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT(1);
    }//endif

  }//endif

  //============================================================================
  // II) Car Parrinello MD:

  if(cp_on==1){
    ifound++;
    CP_integrate_dyn(ncoef,istate,iteration,forces,vpsi_g,psi_g,cmass,
        k_x,k_y,k_z,len_nhc,num_nhc,nck_nhc,fNHC,vNHC,xNHC,xNHCP,
        mNHC,v0NHC,a2NHC,a4NHC,kTCP,fictEke,nkx0_red,nkx0_uni,nkx0_zero,
        ekeNHC,potNHC,degfree,degfreeNHC,degFreeSplt,
        istrNHC,iendNHC,halfStepEvolve);
  }//endif

  //============================================================================
  // Error Handling

  if(ifound==0){

    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Unknown Integration option\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    EXIT();

  }//endif

  //---------------------------------------------------------------------------
}// end routine
//============================================================================



//============================================================================
/** \brief Create fictitious psi masses for CPAIMD and minimization
 */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================

void CPINTEGRATE::CP_create_mass(int nktot,int *k_x,int *k_y,int *k_z,
    double *cmass, int mydoublePack)

  //============================================================================
{ /* Begin Function */
  //---------------------------------------------------------------------------

  GENERAL_DATA *general_data = GENERAL_DATA::get();
  CP           *cp           = CP::get();
#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_cp.h"

  double *hmati        = gencell->hmati;
  double cmass_tau_def = cpcoeffs_info->tau_mass;
  double cmass_cut_def = cpcoeffs_info->ecut_mass;
  double tpi           = 2.0*M_PI;

  double cmass0,pre,pre1,rat,swit;
  double gx,gy,gz,g2;
  double wght,wght_now;

#define DEBUG_PARAMS_OFF
#ifdef DEBUG_PARAMS
  PRINTF("mass stuff : %g %g %g %d\n",ecut,cmass_tau_def,cmass_cut_def,nktot);
  PRINTF("box stuff  : %g %g %g\n",(1.0/hmati[1]),(1.0/hmati[5]),(1.0/hmati[9]));
#endif

  //============================================================================
  // Compute Constants

  cmass_tau_def /= TIME_CONV;
  cmass0         = (2.0*cmass_tau_def*cmass_tau_def)*CP_EMAGIC;
  pre            = 10.0*cmass_cut_def;
  pre1           = 1.0-pre;

  wght           = 1.0;
  if(mydoublePack==1){wght=2.0;}

  //   PRINTF("mydoublePack %d %g %g %g\n",mydoublePack,cmass0,cmass_cut_def,wght);

  //============================================================================
  // Create masses 

  for(int i = 0; i < nktot; i++){

    gx = tpi*(k_x[i]*hmati[1] + k_y[i]*hmati[2] + k_z[i]*hmati[3]);
    gy = tpi*(k_x[i]*hmati[4] + k_y[i]*hmati[5] + k_z[i]*hmati[6]);
    gz = tpi*(k_x[i]*hmati[7] + k_y[i]*hmati[8] + k_z[i]*hmati[9]);
    g2 = gx*gx + gy*gy + gz*gz;

    wght_now = (k_x[i]==0 ? 1.0 : wght);
    if(g2 > cmass_cut_def) {
      rat      = 1.0/g2;
      swit     = pre1*rat + pre;  //=1 at g2=cmass_cut : = 10*cmass_cut_def at g2=infty
      cmass[i] = (wght_now*g2)*(cmass0/cmass_cut_def);
#ifdef DEBUG_MASS
      if(k_x[i]==3 && k_y[i]==3 && k_z[i]==3){
        PRINTF("3 3 3 %g %g %g %g %g\n",wght_now,cmass0,g2,cmass_cut_def,cmass[i]);
      }
#endif
    } else {
      cmass[i] = wght_now*cmass0;
    }//endif

  }//end for

  //---------------------------------------------------------------------------
}// end routine
//============================================================================




//============================================================================
/** \brief Print a warning if the sum of squares of the F=H\psi is too large.
  under CPAIMD
 */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================

void CPINTEGRATE::CheckCoefGradMag(double fovlap)

  //============================================================================
{ /* Begin Function */
  //---------------------------------------------------------------------------

  GENERAL_DATA *general_data = GENERAL_DATA::get();
  CP *cp = CP::get();
#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_cp.h"

  int ifound=0;

  int cp_min_on = gensimopts->cp_min
    + gensimopts->cp_wave_min_pimd 
    + gensimopts->cp_wave_min;

  int cp_on     = gensimopts->cp     
    + gensimopts->cp_pimd
    + gensimopts->cp_wave
    + gensimopts->cp_wave_pimd;


  if(cp_on == 1){
    if(sqrt(fovlap) > cpopts->tol_coef){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("Coef force magnitude greater than tolerance set\n");
      PRINTF("in the input file.  Please reconsider your choices\n");
      PRINTF("of adiabaticity parameters.\n");
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT();
    }// endif 
  }// endif


  //---------------------------------------------------------------------------
}// end routine
//============================================================================

