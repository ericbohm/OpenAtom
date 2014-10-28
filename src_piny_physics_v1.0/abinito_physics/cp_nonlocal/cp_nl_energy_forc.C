#include "standard_include.h"
#include "ckcomplex.h"

#include "../../../include/Atoms.h"
#include "../../../src_mathlib/mathlib.h"    

#include "../class_defs/allclass_gen.h"
#include "../class_defs/allclass_cp.h"

#include "../class_defs/CP_OPERATIONS/class_cpnonlocal.h"

//============================================================================
// This routine is depricated old junk
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//
// The following are the inputs into the function:
// - zMatrixRow: the portion of the reduced zmatrix corresponding to a 
//   particular state.
// - rowSize: the size of zMatrixRow
// - forces: a placeholder to put the computed forces, with 
// - gspace: the portion of gspace over which the forces should be computed
// - arrays k_x, k_y, k_z provide the g_x, g_y, g_z values at which the forces
//   should be computed.
//
//============================================================================

void CPNONLOCAL::CP_enl_force_calc(complex* zMatrixRow, int forcesSize, 
    int *k_x, int *k_y, int *k_z, 
    complex *StructFact,complex *forces,
    int state_ind, int mydoublePack,int numSfGrps, int indexSfGrp)

  //============================================================================
{/* Begin function */
  //============================================================================
  // Local variable declarations   

  GENERAL_DATA *general_data = GENERAL_DATA::get();
  CP           *cp           = CP::get();

#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_cp.h"

  double ecut   = cpcoeffs_info->ecut_psi; // its in Rydberg already
  double vol    = gencell->vol;  
  double *hmati = gencell->hmati;  
  double tpi    = 2.0*M_PI;
  int natm_nl   = cppseudo->nonlocal.natm;

#ifdef GJM_DEBUG_SIZE
  PRINTF(" %d : coefs in CP_enl_force_calc\n",forcesSize);
#endif

  if(indexSfGrp>=numSfGrps ||  indexSfGrp<0 || numSfGrps>natm_nl){
    PRINTF("Incorrect SF index %d %d %d\n",indexSfGrp,numSfGrps,natm_nl);
    EXIT(1);
  }//endif

  int natm_nl_grp,istrt,iend,ioff;
  get_grp_params(natm_nl,numSfGrps,indexSfGrp,&natm_nl_grp,&istrt,&iend,&ioff);

  //============================================================================
  // Compute force on the states the fast way

#define NEW_FORCE_JUNK

#ifdef  NEW_FORCE_JUNK

  char trans    = 'N';

  if(mydoublePack==0){

    int incx      = 1; 
    int incy      = 1;
    complex alpha = complex(-1.0,0.0);
    complex beta  = complex(1.0,0.0); 
    ZGEMV(&trans,&forcesSize,&natm_nl_grp,&alpha,&(StructFact[0]),
        &forcesSize,&(zMatrixRow[0]),&incx,&beta,&(forces[0]),&incy);
    PRINTF("Imaginary part of non-local forces is broken");
    PRINTF("when not double packing for atm groups > 1");
    for(int i=0;i<forcesSize;i++){forces[i].im=-forces[i].im;}
  }else{

    int incx     =  2; 
    int incy     =  1;
    double alpha = -1.0;
    double beta  =  1.0; //Forces from previous iter add in
    int nsize    =  2*forcesSize;
    double *dStructFact = reinterpret_cast<double *>(StructFact);
    double *dzMatrixRow = reinterpret_cast<double *>(zMatrixRow);
    double *dforces     = reinterpret_cast<double *>(forces);
    DGEMV(&trans,&nsize,&natm_nl_grp,&alpha,&(dStructFact[0]),
        &nsize,&(dzMatrixRow[0]),&incx,&beta,&(dforces[0]),&incy);

  }//endif

#ifdef GJM_DBG_NONLOCAL
  for(int i = 0; i < forcesSize; i++){
    //    if(k_x[i]==0 && k_y[i]==1 && k_z[i]==4){
    if(k_x[i]==0 && k_y[i]==0 && k_z[i]==0){
      CkPrintf("%d : %d %d %d : %g %g \n",state_ind,
          k_x[i],k_y[i],k_z[i],forces[i].re,forces[i].im);
    }
  }
#endif

#endif

  //============================================================================
  // Compute force on the states the fast way

#ifndef NEW_FORCE_JUNK

  for(int i = 0; i < forcesSize; i++){

    double gx,gy,gz,g2;
    gx = tpi*(k_x[i]*hmati[1] + k_y[i]*hmati[2] + k_z[i]*hmati[3]);
    gy = tpi*(k_x[i]*hmati[4] + k_y[i]*hmati[5] + k_z[i]*hmati[6]);
    gz = tpi*(k_x[i]*hmati[7] + k_y[i]*hmati[8] + k_z[i]*hmati[9]);
    g2 = gx*gx + gy*gy + gz*gz;

    //-----------------------------------------------------------------------
    //  For this state, compute the force using structure factor Zmatrix

    forces[i] = complex(0.0,0.0);
    if(g2<=ecut){


      for(int iNL = 0; iNL < natm_nl_grp; iNL++){
        int s_ind    = forcesSize*iNL + i;
        double tmp0  = zMatrixRow[iNL].re;
        complex tmp1 = (StructFact[s_ind]*tmp0).conj();
#ifndef GLENN_DEBUG_ON
        forces[i] -= tmp1; 
#endif
      }/* endfor atom loop */
    }//endif

    //-----------------------------------------------------------------------
  }/* endfor g-space loop */
#endif

  //============================================================================
  } /* End function */
  //============================================================================


  //============================================================================
  // This routine is is the ticket
  //============================================================================
  //cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  //============================================================================

  void CPNONLOCAL::CP_enl_atm_forc_calc(int numSfGrps, int indexSfGrp, FastAtoms *atoms,
      complex *zmatrixSum,complex *zmatrixSum_fx,complex *zmatrixSum_fy,
      complex *zmatrixSum_fz,double *enl_ret,int mydoublePack, int istate)

    //============================================================================
  {// begin routine
    //============================================================================

    GENERAL_DATA *general_data = GENERAL_DATA::get();
    CP           *cp           = CP::get();

    int natm_nl     = cp->cppseudo.nonlocal.natm;
    int *map_nl     = cp->cppseudo.nonlocal.map_nl;

    double *vnorm_0 = cp->cppseudo.nonlocal.vnorm_0;
    double *occ     = cp->cpcoeffs_info.occ_up;

    double occ_now  = occ[istate+1];
    double vol      = general_data->gencell.vol;
    double fpi      = 4.0*M_PI;
    double y00      = 1.0/sqrt(fpi);

    double *fx    = atoms->fx;
    double *fy    = atoms->fy;
    double *fz    = atoms->fz;


    if(indexSfGrp>=numSfGrps ||  indexSfGrp<0 || numSfGrps>natm_nl){
      PRINTF("Incorrect SF index %d %d %d\n",indexSfGrp,numSfGrps,natm_nl);
      EXIT(1);
    }//endif

    int natm_nl_grp,istrt,iend,ioff;
    get_grp_params(natm_nl,numSfGrps,indexSfGrp,&natm_nl_grp,&istrt,&iend,&ioff);

    //============================================================================

    double enl = 0.0;

    if(mydoublePack==0){

      for (int i = 0; i < natm_nl_grp; i++){
        int k   = i+istrt;
        double zmag = zmatrixSum[i].getMagSqr();
        double norm = vnorm_0[k]/vol;
        enl        += (zmag*norm);
        zmatrixSum[i].re *= (2.0*norm);
        zmatrixSum[i].im *= (2.0*norm);
      }//endfor

      for (int i = 0; i < natm_nl_grp; i++){
        int k   = i+istrt;
        int ind = map_nl[k]-1;
        fx[ind] -= (zmatrixSum[i]*zmatrixSum_fx[i].conj()).re;
        fy[ind] -= (zmatrixSum[i]*zmatrixSum_fy[i].conj()).re;
        fz[ind] -= (zmatrixSum[i]*zmatrixSum_fz[i].conj()).re;
      }//endfor

    }else{

      double *dzmatrixSum    = reinterpret_cast<double *>(zmatrixSum);
      double *dzmatrixSum_fx = reinterpret_cast<double *>(zmatrixSum_fx);
      double *dzmatrixSum_fy = reinterpret_cast<double *>(zmatrixSum_fy);
      double *dzmatrixSum_fz = reinterpret_cast<double *>(zmatrixSum_fz);

      for (int i = 0,j=0; i < natm_nl_grp; i++,j+=2){
        int k           = i+istrt;
        double zmag     = dzmatrixSum[j]*dzmatrixSum[j];
        double norm     = vnorm_0[k]/vol;
        enl            += (zmag*norm);
        dzmatrixSum[j] *= (2.0*norm);
      }//endfor

      for (int i = 0,j=0; i < natm_nl_grp; i++,j+=2){
        int k   = i+istrt;
        int ind = map_nl[k]-1;
        fx[ind] -= dzmatrixSum[j]*dzmatrixSum_fx[j];
        fy[ind] -= dzmatrixSum[j]*dzmatrixSum_fy[j];
        fz[ind] -= dzmatrixSum[j]*dzmatrixSum_fz[j];
      }//endfor

    }//endif

    (*enl_ret) = enl*occ_now;

    //============================================================================
  } /* End function */
  //============================================================================

