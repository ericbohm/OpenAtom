#include "standard_include.h"
#include "../../../include/debug_flags.h"
#include "../../../include/CPcharmParaInfo.h"
#include "../class_defs/allclass_gen.h"
#include "../class_defs/allclass_cp.h"
#include "../include/class_defs/allclass_mdatoms.h"
#include "../class_defs/PINY_INIT/PhysicsParamTrans.h"

//========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//========================================================================
// This function is called once on the main processor. All initial
// values can be put in here. These values will be accessed from the
// CPcharmParaInfo class.
//========================================================================

void PhysicsParamTransfer::ParaInfoInit(CPcharmParaInfo *sim)

//========================================================================
  {//begin routine
//========================================================================
// Local variables 

  GENERAL_DATA *general_data = GENERAL_DATA::get();
  CP *cp                     = CP::get();
  MDATOMS      *mdatoms      = MDATOMS::get();
#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_cp.h"

  int ndump_frq       = genfilenames->iwrite_dump;
  int istart_typ_cp   = gensimopts->istart_cp;
  int cp_opt          = (gensimopts->cp+gensimopts->cp_wave);
  int cp_min_opt      = (gensimopts->cp_wave_min+gensimopts->cp_min);
  int cp_std          = gensimopts->cp;
  int cp_wave         = gensimopts->cp_wave;
  int cp_min_cg       = genminopts->cp_min_cg;
  int cp_min_std      = genminopts->cp_min_std;
  int *kmax           = cpewald->kmax_cp_dens_cp_box;
  int nkf1            = 4*(kmax[1]+1);
  int nkf2            = 4*(kmax[2]+1);
  int nkf3            = 4*(kmax[3]+1);
  int nstates         = cpcoeffs_info->nstate_up;
  int ntime           = gentimeinfo->ntime;
  int ibinary_opt     = cpopts->iread_coef_binary;
  int ibinary_write_opt= cpopts->iwrite_coef_binary;
  int natm_tot        = (mdatoms->mdclatoms_info.natm_tot);
  int natm_nl         = (cppseudo->nonlocal.natm);
  int natm_typ        = cppseudo->natm_typ;
  int cp_grad_corr_on = cpopts->cp_gga;

  double vol          = gencell->vol;
  double dt           = gentimeinfo->dt;
  double tol_norb     = cpconstrnt->c_tolnorb;
  double tol_cp_min   = genminopts->tol_coef;
  double tol_cp_dyn   = cpopts->tol_coef;

  int ees_eext_on     = cppseudo->nonlocal.ees_eext_on;
  int ees_nloc_on     = cppseudo->nonlocal.ees_on;
  int ngrid_eext_a    = cppseudo->ngrid_eext_a;
  int ngrid_eext_b    = cppseudo->ngrid_eext_b;
  int ngrid_eext_c    = cppseudo->ngrid_eext_c;
  int ngrid_nloc_a    = cppseudo->nonlocal.ngrid_a;
  int ngrid_nloc_b    = cppseudo->nonlocal.ngrid_b;
  int ngrid_nloc_c    = cppseudo->nonlocal.ngrid_c;

  int nlIters         = (cppseudo->nonlocal.nl_iter);
  int nmem_zmat_tot   = cppseudo->nonlocal.ntot_zmat;
  int *nmem_zmat      = cppseudo->nonlocal.n_zmat;
  int *ioff_zmat      = cppseudo->nonlocal.ioff_zmat;

//========================================================================

   if(cp_opt==0 && cp_min_opt==0){
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
     PRINTF("The present code only does cp applications\n");
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
     EXIT(1);
   }//endif

   if(cp_opt==1 && istart_typ_cp ==0){
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
     PRINTF("No gen-wave restarts for cp-dynamics\n");
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
     EXIT(1);
   }//endif

   if(istart_typ_cp ==0){
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
     PRINTF("No gen-wave restarts yet. Sorry\n");
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
     EXIT(1);
   }//endif

//========================================================================

   sim->nstates        = nstates;
   sim->natm_typ       = natm_typ;
   sim->natm_tot       = natm_tot;
   sim->natm_nl        = natm_nl;
   sim->nlIters        = nlIters;
   sim->ntime          = ntime;
   sim->dt             = dt;
   sim->vol            = vol;

   sim->cp_min_opt     = cp_min_opt;
   sim->cp_min_cg      = cp_min_cg;
   sim->cp_min_std     = cp_min_std;
   sim->cp_opt         = cp_opt;
   sim->cp_std         = cp_std;
   sim->cp_wave        = cp_wave;
   sim->cp_grad_corr_on= cp_grad_corr_on;

   sim->tol_norb       = tol_norb;
   sim->tol_cp_min     = tol_cp_min;
   sim->tol_cp_dyn     = tol_cp_dyn;

   sim->ndump_frq      = ndump_frq;
   sim->istart_typ_cp  = istart_typ_cp;
   sim->ibinary_opt    = ibinary_opt;
   sim->ibinary_write_opt= ibinary_write_opt;

   sim->sizeX          = nkf1;
   sim->sizeY          = nkf2;
   sim->sizeZ          = nkf3;

   sim->ees_eext_on    = ees_eext_on;
   sim->ees_nloc_on    = ees_nloc_on;
   sim->ngrid_nloc_a   = ngrid_nloc_a; 
   sim->ngrid_nloc_b   = ngrid_nloc_b; 
   sim->ngrid_nloc_c   = ngrid_nloc_c; 
   sim->ngrid_eext_a   = ngrid_eext_a; 
   sim->ngrid_eext_b   = ngrid_eext_b; 
   sim->ngrid_eext_c   = ngrid_eext_c; 

   sim->nmem_zmat_tot  = nmem_zmat_tot;
   sim->ioff_zmat      = new int [nlIters];
   sim->nmem_zmat      = new int [nlIters];
   for(int i=0;i<nlIters;i++){
     sim->ioff_zmat[i] = ioff_zmat[(i+1)];
     sim->nmem_zmat[i] = nmem_zmat[(i+1)];
   }//endfor

//-----------------------------------------------------------------------
  }//end routine
//========================================================================


//========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//========================================================================

void PhysicsParamTransfer::control_mapping_function(CPcharmParaInfo *sim,
                                                    int mydoublePack)

//========================================================================
  {//begin routine
//========================================================================
// Local variables 

  CP *cp                     = CP::get();
#include "../class_defs/allclass_strip_cp.h"

     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
     PRINTF("The present load balance function is not this one\n");
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
     EXIT(1);

  int *kmax      = cpewald->kmax_cp_dens_cp_box;
  int nkf1       = 4*(kmax[1]+1);
  int nkf2       = 4*(kmax[2]+1);
  int nkf3       = 4*(kmax[3]+1);
  int nktot      = cpewald->nktot_sm;
  int *ka        = cpewald->kastr_sm;
  int *kb        = cpewald->kbstr_sm;
  int *kc        = cpewald->kcstr_sm;

  int cp_para_typ=0;

  double *lines_per_plane;
  double *pts_per_plane;
  lines_per_plane = new double[nkf1];
  pts_per_plane   = new double[nkf1];
  
//========================================================================

   PRINTF("\n");PRINT_LINE_STAR;
   PRINTF("Creating G-space mapping functions\n");
   PRINT_LINE_DASH;PRINTF("\n");

   if(mydoublePack==0){

     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
     PRINTF("The present code only compute doublePack maps for CP\n");
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
     EXIT(1);

   }//endif

//========================================================================
   
  if(mydoublePack==1){

    if(cp_para_typ==1){
      compute_lines_per_plane_half_plane(nkf1,nkf2,nkf3,nktot,lines_per_plane,
                                         pts_per_plane,ka,kb,kc);
    }else{
      compute_lines_per_plane_half_sphere(nkf1,nkf2,nkf3,nktot,lines_per_plane,
                                        pts_per_plane,ka,kb,kc);
    }//endif

  }//endif

  sim->nplane_x        = nkf1;
  sim->lines_per_chareG = lines_per_plane;
  sim->pts_per_chareG   = pts_per_plane;

//========================================================================

  PRINTF("\n");PRINT_LINE_DASH;
  PRINTF("Completed g-space map set up\n");
  PRINT_LINE_STAR;PRINTF("\n");

//========================================================================
  }//end routine
//========================================================================


//========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//========================================================================

void PhysicsParamTransfer::compute_lines_per_plane_half_plane(
            int nkf1, int nkf2, int nkf3,int nktot, 
            double *lines, double *pts, int *ka, int *kb, int *kc
           )

//========================================================================
   {//begin routine
//========================================================================
//                 Cubic box salute
// Note Piny g-space isn't set up to match driver g-space
// but who cares because driver is a cubic box code. 04/15/2005
//      For PINY
// ka is the half space variable
// kc varies on full range.
// compute pts on kc-planes
//========================================================================
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
     PRINTF("The present load balance function is not this one\n");
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
     EXIT(1);

  if(nkf1!=nkf2 || nkf2!=nkf3){
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
     PRINTF("The present load balance function is for cubic boxes\n");
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
     EXIT(1);
  }//endif

  PRINTF("Computing mapping function for dble packed g-space planes.\n");
  PRINTF("The full number of planes are 1/2 filled.\n");

//========================================================================
// compute points on each kc plane

  for(int i=0;i<nkf3;i++){pts[i]=0.0;}
  for(int i=1;i<=nktot;i++){
    int idx = (kc[i] < 0 ? kc[i]+nkf3 : kc[i]);
    pts[idx]+=1.0;
  }//endfor
  pts[0]+=1.0; //g=0

  double rnorm_p = 1.0/((double)(nktot+1));
  for(int i=0;i<nkf3;i++){pts[i] *= rnorm_p;}
   
//========================================================================
// On each kc plane, compute the number of FFTs along kb direction.
// The number kb ffts is equal to the number of ka values in the plane.

  int *ka_on = (int *)cmalloc(nkf1*sizeof(int),"compute_lines_per_plane");
  int n2      = nkf3/2;

  for(int i=0;i<nkf3;i++){lines[i]=0.0;}
  double sum = 0.0;
  for(int j=0;j<nkf3;j++){
    int kcc = (j>n2 ? j-nkf3 : j);
    for(int i=0;i<nkf1;i++){ka_on[i]=0;}
    for(int i=1;i<=nktot;i++){
      if(kc[i]==kcc){
        int k = (ka[i] < 0 ? ka[i]+nkf1 : ka[i]);
        ka_on[k]=1;
      }//endif
    }//endfor
    for(int i=0;i<nkf1;i++){
      double add = ((double)ka_on[i]);
      lines[j] += add; sum += add;
    }//endfor
  }//endfor

  double rnorm_l = 1.0/sum;
  for(int i=0;i<nkf3;i++){lines[i]*=rnorm_l;}

  cfree(ka_on,"compute_lines_per_plane");

//========================================================================
// Just a little debug, early in the morning.

#ifdef _CP_DEBUG_LINES_PER_PLANE_
  FILE *fp = fopen("map_gjm_half_plane.out","w");
    double sum_p = 0.0;
    double sum_l = 0.0;
    for(int i=0;i<nkf3;i++){
      fprintf(fp,"line %d : pts %g %g: lines %g %g\n",
                 i,pts[i],pts[i]/rnorm_p,lines[i],lines[i]/rnorm_l);
      sum_p += pts[i];
      sum_l += lines[i];
    }//endfor
    fprintf(fp,"total : %g %g\n",sum_p,sum_l);
  fclose(fp);
#endif

//-----------------------------------------------------------------------
  }//end routine
//========================================================================


//========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//========================================================================

void PhysicsParamTransfer::compute_lines_per_plane_half_sphere(
            int nkf1, int nkf2, int nkf3,int nktot, 
            double *lines, double *pts, int *ka, int *kb, int *kc
           )

//========================================================================
   {//begin routine
//========================================================================
//                 Cubic box salute
// Note Piny g-space isn't set up to match driver g-space
// but who cares because driver is a cubic box code. 04/15/2005
//      For PINY
// ka is the half space variable
// kc varies on full range.
// compute pts on kc-planes
//========================================================================

     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
     PRINTF("The present load balance function is not this one\n");
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
     EXIT(1);
  if(nkf1!=nkf2 || nkf2!=nkf3){
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
     PRINTF("The present load balance function is for cubic boxes\n");
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
     EXIT(1);
  }//endif

  PRINTF("Computing mapping function for dble packed g-space planes.\n");
  PRINTF("Half the number of planes are filled.\n");
   
//========================================================================
// compute points on each ka plane : cpaimd keeps a little too much data

  double psum = 0.0;
  for(int i=0;i<nkf1;i++){pts[i]=0.0;}
  for(int i=1;i<=nktot;i++){
    int idx = (ka[i] < 0 ? ka[i]+nkf1 : ka[i]); // ka[i] >= 0
    if(ka[i]==0){
      psum += 2.0;
      pts[idx]+=2.0;
    }else{
      psum += 1.0;
      pts[idx]+=1.0;
    }
  }//endfor
  pts[0]+=1.0; //g=0
  psum += 1.0;

  double rnorm_p = 1.0/psum;
  for(int i=0;i<nkf1;i++){pts[i] *= rnorm_p;}
   
//========================================================================
// On each ka plane, compute the number of FFTs along kb direction.
// The number of kb ffts is equal to the number of kc values in the plane.
// Maximum number is nkf3 ffts of size nkf2

  int *kb_on = (int *)cmalloc(nkf3*sizeof(int),"compute_lines_per_plane");
  int n2      = nkf1/2;

  for(int i=0;i<nkf3;i++){lines[i]=0.0;}
  double sum = 0.0;
  for(int j=0;j<nkf1;j++){
    int kaa = (j>n2 ? j-nkf1 : j);
    for(int i=0;i<nkf3;i++){kb_on[i]=0;}
    for(int i=1;i<=nktot;i++){
      if(ka[i]==kaa){
        int k = (kc[i] < 0 ? kc[i]+nkf3 : kc[i]);
        kb_on[k]=1;
      }//endif
    }//endfor
    for(int i=0;i<nkf3;i++){
      double add = ((double)kb_on[i]);
      lines[j] += add; sum += add;
    }//endfor
  }//endfor

  double rnorm_l = 1.0/sum;
  for(int i=0;i<nkf1;i++){lines[i]*=rnorm_l;}

  cfree(kb_on,"compute_lines_per_plane");

//========================================================================
// Just a little debug, early in the morning.

#ifdef _CP_DEBUG_LINES_PER_PLANE_
  FILE *fp = fopen("map_gjm_half_sphere.out","w");
    double sum_p = 0.0;
    double sum_l = 0.0;
    for(int i=0;i<nkf3;i++){
      fprintf(fp,"line %d : pts %g %g: lines %g %g\n",
                 i,pts[i],pts[i]/rnorm_p,lines[i],lines[i]/rnorm_l);
      sum_p += pts[i];
      sum_l += lines[i];
    }//endfor
    fprintf(fp,"total : %g %g\n",sum_p,sum_l);
  fclose(fp);
#endif

//-----------------------------------------------------------------------
  }//end routine
//========================================================================



//========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//========================================================================

void PhysicsParamTransfer::get_Sfgrp_params(int natm_nl, int numSfGrps, 
       int indexSfGrp,int *n_ret, int *istrt_ret, int *iend_ret)
{
   int n     = (natm_nl/numSfGrps);
   int m     = (natm_nl % numSfGrps);

   int istrt = n*indexSfGrp + 1;
   if(indexSfGrp>=m){istrt += m;}
   if(indexSfGrp<m) {istrt += indexSfGrp;}
   if(indexSfGrp<m) {n++;}
   int iend  = n+istrt-1;

   (*n_ret)     = n;
   (*istrt_ret) = istrt;
   (*iend_ret)  = iend;
}
//========================================================================



//========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//========================================================================

void PhysicsParamTransfer::get_Sfgrp_max(int natm_nl, int numSfGrps, 
                                         int *natm_nl_grp_max)
{
     int n = (natm_nl / numSfGrps);
     int m = (natm_nl % numSfGrps);
     if(m!=0){n++;}

     (*natm_nl_grp_max) = n;
}
//========================================================================


//========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//========================================================================

void PhysicsParamTransfer::control_new_mapping_function(CPcharmParaInfo *sim,
                                                    int mydoublePack)

//========================================================================
  {//begin routine
//========================================================================
// Local variables 

   int *nlines_per_chareG = sim->nlines_per_chareG;
   int *npts_per_chareG   = sim->npts_per_chareG;
   int nchareG            = sim->nchareG;
   int sizeX              = sim->sizeX;
   int nlines_tot         = sim->nlines_tot;
   int npts_tot           = sim->npts_tot;

   int nsize = nchareG;
   double *pts_per_chareG   = new double [nsize];
   double *lines_per_chareG = new double [nsize];

   int npts   = 0;
   int nlines = 0;
   for(int i=0;i<nchareG;i++){
     pts_per_chareG[i]= ((double)npts_per_chareG[i])/((double)npts_tot);
     npts += npts_per_chareG[i];
     lines_per_chareG[i]= ((double)nlines_per_chareG[i])/((double)nlines_tot);
     nlines += nlines_per_chareG[i];
   }//endfor

   if(npts!=npts_tot || nlines != nlines_tot){
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
     PRINTF("Inconsistency between nlines and pts per plane and totals\n");
     PRINTF("%d %d : %d %d \n",npts,npts_tot,nlines,nlines_tot);
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
     EXIT(1);
   }//endif

   sim->pts_per_chareG   = pts_per_chareG;
   sim->lines_per_chareG = lines_per_chareG;

//========================================================================
  }
//========================================================================
