#include "standard_include.h"
#include "ckcomplex.h"

#include "../../../include/Atoms.h"    
#include "../../../src_mathlib/mathlib.h"    
#include "../class_defs/allclass_cp.h"

#include "../class_defs/CP_OPERATIONS/class_cpnonlocal.h"

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
//
//============================================================================

void CPNONLOCAL::CP_enl_matrix_calc(int gSpaceSize, complex *gspace, 
    int *k_x, int *k_y, int *k_z, 
    complex *StructFact,
    complex *StructFact_fx,complex *StructFact_fy,
    complex *StructFact_fz,
    complex *zmatrix, 
    complex *zmatrix_fx,complex *zmatrix_fy, 
    complex *zmatrix_fz, 
    int state_ind, int mydoublePack,
    int numSfGrps, int indexSfGrp)

  //============================================================================
{/* Begin function */
  //============================================================================

  CP *cp = CP::get();
#include "../class_defs/allclass_strip_cp.h"

  /*------------------*/
  /* Atom information */

  int natm_nl  = cppseudo->nonlocal.natm;
  double ecut  = 2.0*cpcoeffs_info->ecut; // in Rydbergs -> factor of 2
#ifdef GJM_DEBUG_SIZE
  PRINTF(" %d : coefs in CP_enl_matrix_calc : state %d\n",gSpaceSize,state_ind);
#endif

  //============================================================================

  if(indexSfGrp>=numSfGrps ||  indexSfGrp<0 || numSfGrps>natm_nl){
    PRINTF("Incorrect SF index %d %d %d\n",indexSfGrp,numSfGrps,natm_nl);
    EXIT(1);
  }//endif

  int natm_nl_grp,istrt,iend,ioff;
  get_grp_params(natm_nl,numSfGrps,indexSfGrp,&natm_nl_grp,&istrt,&iend,&ioff);

  //============================================================================
  // Compute non-local matrix using zgemvs

#define NEW_JUNK
#ifdef  NEW_JUNK

  char trans    ='t';
  if(mydoublePack==0){
    complex alpha = complex(1.0,0.0);
    complex beta  = complex(0.0,0.0);
    int incx      = 1; 
    int incy      = 1;

    ZGEMV(&trans,&gSpaceSize,&natm_nl_grp,&alpha,&(StructFact[0]),
        &gSpaceSize,&(gspace[0]),&incx,&beta,&(zmatrix[0]),&incy);

    ZGEMV(&trans,&gSpaceSize,&natm_nl_grp,&alpha,&(StructFact_fx[0]),
        &gSpaceSize,&(gspace[0]),&incx,&beta,&(zmatrix_fx[0]),&incy);

    ZGEMV(&trans,&gSpaceSize,&natm_nl_grp,&alpha,&(StructFact_fy[0]),
        &gSpaceSize,&(gspace[0]),&incx,&beta,&(zmatrix_fy[0]),&incy);

    ZGEMV(&trans,&gSpaceSize,&natm_nl_grp,&alpha,&(StructFact_fz[0]),
        &gSpaceSize,&(gspace[0]),&incx,&beta,&(zmatrix_fz[0]),&incy);

  }else{

    int incx      = 1; 
    int incy      = 2;
    int nsize     = 2*gSpaceSize;
    double alpha  = 1.0;
    double beta   = 0.0;
    memset(zmatrix,0,sizeof(complex)*natm_nl_grp);
    memset(zmatrix_fx,0,sizeof(complex)*natm_nl_grp);
    memset(zmatrix_fy,0,sizeof(complex)*natm_nl_grp);
    memset(zmatrix_fz,0,sizeof(complex)*natm_nl_grp);
    double *dgspace       = reinterpret_cast<double *>(gspace);
    double *dStructFact   = reinterpret_cast<double *>(StructFact);
    double *dStructFact_fx= reinterpret_cast<double *>(StructFact_fx);
    double *dStructFact_fy= reinterpret_cast<double *>(StructFact_fy);
    double *dStructFact_fz= reinterpret_cast<double *>(StructFact_fz);
    double *dzmatrix      = reinterpret_cast<double *>(zmatrix);
    double *dzmatrix_fx   = reinterpret_cast<double *>(zmatrix_fx);
    double *dzmatrix_fy   = reinterpret_cast<double *>(zmatrix_fy);
    double *dzmatrix_fz   = reinterpret_cast<double *>(zmatrix_fz);

    DGEMV(&trans,&nsize,&natm_nl_grp,&alpha,&(dStructFact[0]),
        &nsize,&(dgspace[0]),&incx,&beta,&(dzmatrix[0]),&incy);

    DGEMV(&trans,&nsize,&natm_nl_grp,&alpha,&(dStructFact_fx[0]),
        &nsize,&(dgspace[0]),&incx,&beta,&(dzmatrix_fx[0]),&incy);

    DGEMV(&trans,&nsize,&natm_nl_grp,&alpha,&(dStructFact_fy[0]),
        &nsize,&(dgspace[0]),&incx,&beta,&(dzmatrix_fy[0]),&incy);

    DGEMV(&trans,&nsize,&natm_nl_grp,&alpha,&(dStructFact_fz[0]),
        &nsize,&(dgspace[0]),&incx,&beta,&(dzmatrix_fz[0]),&incy);

  }//endif

#endif

  //============================================================================
  // Compute the non-local pseudo slowly

#ifndef  NEW_JUNK

  for(int i = 0; i < gSpaceSize; i++){

    //============================================================================
    // Only compute something if g^2 is nonzero

    gx = tpi*(k_x[i]*hmati[1] + k_y[i]*hmati[2] + k_z[i]*hmati[3]);
    gy = tpi*(k_x[i]*hmati[4] + k_y[i]*hmati[5] + k_z[i]*hmati[6]);
    gz = tpi*(k_x[i]*hmati[7] + k_y[i]*hmati[8] + k_z[i]*hmati[9]);
    g2 = gx*gx + gy*gy + gz*gz;

    if(g2<=ecut ){

      //------------------------------------------------------------------
      // Compute elements of the Z matrix using plane-wave coeffs and 
      // structure factor scaled by spherical harmonic and radial projector

      for(int iNL = 0; iNL < natm_nl_grp; iNL++){
        int s_ind     = gSpaceSize*iNL + i;
        complex CbyS  = gspace[i]*StructFact[s_ind];
        complex CbySx = gspace[i]*StructFact_fx[s_ind];
        complex CbySy = gspace[i]*StructFact_fy[s_ind];
        complex CbySz = gspace[i]*StructFact_fz[s_ind];
#ifdef CHECK_ZMATRIX
        if(state_ind == 0 && k_x[i] == 1 && k_y[i] == 0 && k_z[i] == 0){
          CkPrintf("zmat | %g %g : %g %g : %g %g\n",
              gspace[i].re,gspace[i].im,
              StructFact[s_ind].re,StructFact[s_ind].im
              CbyS.re,CbyS.im);
        }//endif
#endif
        zmatrix[iNL]    += CbyS;
        zmatrix_fx[iNL] += CbySx;
        zmatrix_fy[iNL] += CbySy;
        zmatrix_fz[iNL] += CbySz;
      }/* endfor atom loop */

    }// g2 within range
    //------------------------------------------------------------------
  }/* endfor g-space loop */
#endif

  //----------------------------------------------------------------------------
} /* End function */
//============================================================================

