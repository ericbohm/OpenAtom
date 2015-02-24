//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
/** \file name cp_integrate_dyn
 ** \brief The physics routines that evolve the wavefunctions (psi)
 for CPAIMD
 */
//==========================================================================
#include "standard_include.h"
#include "ckcomplex.h"

#include "../class_defs/allclass_gen.h"
#include "../class_defs/allclass_cp.h"
#include "../class_defs/allclass_mdintegrate.h"

#include "../class_defs/CP_OPERATIONS/class_cpintegrate.h"
#include "../../mathlib/proto_math.h"

#define ISOKIN_OPT 1

//============================================================================
/* 
 ** \brief The physics routine that controls wavefunctions (psi)
 evoluation under CPAIMD
 */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CPINTEGRATE::CP_integrate_dyn(int ncoef, int istate,int iteration,
    complex *forces,complex *vpsi,complex *psi,double *cmass, 
    int *k_x, int *k_y,int *k_z,int len_nhc, int num_nhc,int nck_nhc,
    double ***fNHC,double ***vNHC,double ***xNHC,double ***xNHCP,
    double *mNHC,double *v0NHC, double *a2NHC, double *a4NHC, double kTCP,
    double *fictEke,int nkx0_red,int nkx0_uni,int nkx0_zero,
    double *ekeNhc, double *potNHC,double degfree,double degfreeNHC,
    double *degFreeSplt,int *istrNHC, int *iendNHC,int halfStepEvolve)
  //============================================================================
{ // Begin Function 
  //----------------------------------------------------------------------------
  // Local Variables

  GENERAL_DATA *general_data = GENERAL_DATA::get();
  CP *cp = CP::get();
#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_cp.h"

  double dt      = gentimeinfo->dt * gentimeinfo->bomd_scale;
  double dt2     = dt*0.5;
  int istrt      = nkx0_red;
  int applyorder;
  int fwdFlag    = 1;

  //============================================================================
  // Update the velocities from time, t-dt/2, to time, t (if not the intial step).
  // Compute the Fictious Kinetic Energy at time, t (e.g. t=0 on initial step).

  if(iteration>1 && halfStepEvolve==1){
    applyorder = 2;
    cp_evolve_vel(ncoef,forces,vpsi,cmass,len_nhc,num_nhc,nck_nhc,fNHC,vNHC,xNHC,xNHCP,
        mNHC,v0NHC,a2NHC,a4NHC,kTCP,nkx0_red,nkx0_uni,nkx0_zero,
        applyorder,iteration,degfree,degfreeNHC,degFreeSplt,
        istrNHC,iendNHC,fwdFlag);
  }//endif

  get_fictKE(ncoef,vpsi,cmass,len_nhc,num_nhc,nck_nhc,vNHC,xNHC,xNHCP,mNHC,
      nkx0_red,nkx0_uni,nkx0_zero,fictEke,ekeNhc,potNHC);

  //============================================================================
  // Update the velocities from time, t, to time, t+dt/2. 

  applyorder = 1;
  cp_evolve_vel(ncoef,forces,vpsi,cmass,len_nhc,num_nhc,nck_nhc,fNHC,vNHC,xNHC,xNHCP,
      mNHC,v0NHC,a2NHC,a4NHC,kTCP,nkx0_red,nkx0_uni,nkx0_zero,
      applyorder,iteration,degfree,degfreeNHC,degFreeSplt,
      istrNHC,iendNHC,fwdFlag);

  //============================================================================
  // Update the positions to the next step (t+dt)

  for(int i=istrt;i<ncoef;i++){
    psi[i] += vpsi[i]*dt;
  }//endfor

  //============================================================================
} // End function 
//============================================================================


//============================================================================
/* 
 ** \brief Evolve the CP faux velocities to dt/2
 */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CPINTEGRATE::cp_evolve_vel(int ncoef_full, complex *forces, complex *vpsi,
    double *cmass,int len_nhc, int num_nhc, int nck_nhc, 
    double ***fNHC,double ***vNHC,double ***xNHC,double ***xNHCP,
    double *mNHC,double *v0, double *a2, double *a4, double kTCP,
    int nkx0_red,int nkx0_uni,int nkx0_zero,
    int applyorder,int iteration,double degfree,
    double degfreeNHC,double *degFreeSplt,
    int *istrNHC,int *iendNHC,int fwdFlag)
  //============================================================================
{ // Begin Function
  //----------------------------------------------------------------------------
  // Local Variables and error checking

  GENERAL_DATA *general_data = GENERAL_DATA::get();
  CP *cp = CP::get();
#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_cp.h"

  double dt      = gentimeinfo->dt * gentimeinfo->bomd_scale;
  if(fwdFlag==-1){dt=-dt;}

  double dt2     = dt*0.5;
  int istrt0     = nkx0_red;
  int istrt1     = nkx0_red+nkx0_zero;
  int iend       = nkx0_uni+nkx0_red;
  int ncoef      = ncoef_full-istrt0;

  if(applyorder==2 && iteration==1){
    PRINTF("@@@@@@@@@_ERROR_@@@@@@@@@@@\n");
    PRINTF("applyorder==2 when iteration==1 !\n");
    PRINTF("Internal error in cp_dynamics.C.\n");
    PRINTF("@@@@@@@@@_ERROR_@@@@@@@@@@@\n");
    EXIT(1);
  }//endif


  //============================================================================
  // Make life simple for gx=0 and g=0 at Gamma: This does nothing for kpoints

  for(int i=istrt0;i<istrt1;i++){forces[i].im=0.0;vpsi[i].im=0.0;}
  for(int i=istrt1;i<iend;i++){
    forces[i].re *= 2.0;
    forces[i].im *= 2.0;
    cmass[i]     *= 2.0;
  }//endif

  //============================================================================
  // Isokin-NHC apply first : start a new step

  if(applyorder==1 && ISOKIN_OPT==1){
    cp_isoNHC_update(ncoef,&vpsi[istrt0],&cmass[istrt0],
        len_nhc,num_nhc,nck_nhc,xNHC,xNHCP,vNHC,fNHC,mNHC,
        v0,a2,a4,kTCP,degfree,degfreeNHC,degFreeSplt,
        istrNHC,iendNHC,fwdFlag);
  }//endif

  //============================================================================
  // Evolve velocity with forces : apply==1 its a new step : apply==2 finish old step

  for(int i=istrt0;i<ncoef_full;i++){
    vpsi[i] += forces[i]*(dt2/cmass[i]);     
  }//endfor

  //============================================================================
  // Isokin-NHC apply second : finish old step

  if(ISOKIN_OPT==1 && applyorder==2){
    cp_isoNHC_update(ncoef,&vpsi[istrt0],&cmass[istrt0],
        len_nhc,num_nhc,nck_nhc,xNHC,xNHCP,vNHC,fNHC,mNHC,
        v0,a2,a4,kTCP,degfree,degfreeNHC,degFreeSplt,
        istrNHC,iendNHC,fwdFlag);
  }//endif

  //============================================================================
  // fix up gx=0 ; g=0 at gamma : This does nothing for k-point sampling

  for(int i=istrt0;i<istrt1;i++){forces[i].im=0.0;vpsi[i].im=0.0;}
  for(int i=istrt1;i<iend;i++){
    forces[i].re /= 2.0;
    forces[i].im /= 2.0;
    cmass[i]     /= 2.0;
  }//endif

  //============================================================================
} // End function
//============================================================================


//============================================================================
/* 
 ** \brief Compute the fictitious electronic kinetic energy
 */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CPINTEGRATE::get_fictKE(int n,complex *v, double *m,int len, int num,int nck,
    double ***vNHC,double ***xNHC,double ***xNHCP,double *mNHC,
    int nkx0_red,int nkx0_uni,int nkx0_zero,
    double *ekin_ret,double *ekinNHC_ret,double *potNHC_ret)
  //============================================================================
{// Begin Function 
  //============================================================================

  int istrt0     = nkx0_red;
  int istrt      = nkx0_red+nkx0_zero;
  int iend       = nkx0_red+nkx0_uni;

  double ekin    = 0.0;
  double ekinNHC = 0.0;
  double potNHC  = 0.0;

  if(ISOKIN_OPT==1){
    for(int k=0;k<nck;k++){
      for(int j=0;j<len;j++){
        for(int i=0;i<num;i++){
          ekinNHC += mNHC[j]*vNHC[k][i][j]*vNHC[k][i][j];
          potNHC  += xNHC[k][i][j];
        }}}//endfor
  }//endif

  for(int i=istrt0;i<istrt;i++){ekin += v[i].getMagSqr()*m[i];}       // g=0
  for(int i=istrt;i<iend;i++)  {ekin += v[i].getMagSqr()*(2.0*m[i]);} // gx=0
  for(int i=iend;i<n;i++)      {ekin += v[i].getMagSqr()*m[i];}       // gx!=0
  ekin     *=0.5;
  ekinNHC  *=0.5;

  ekin_ret[0]    = ekin;
  ekinNHC_ret[0] = ekinNHC;
  potNHC_ret[0]  = potNHC;

  //============================================================================
}//end routine
//============================================================================



