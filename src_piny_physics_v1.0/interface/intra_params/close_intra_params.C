/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: close_intra_parms.c                          */
/*                                                                          */
/* This subprogram tidies things up                                         */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../class_defs/allclass_gen.h"
#include "../class_defs/allclass_mdintegrate.h"
#include "../class_defs/allclass_mdatoms.h"
#include "../class_defs/allclass_mdinter.h"
#include "../class_defs/allclass_mdintra.h"
#include "../class_defs/allclass_cp.h"
#include "../class_defs/typedefs_par.h"
#include "../proto_defs/proto_intra_params_entry.h"
#include "../proto_defs/proto_intra_params_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void close_intra_params(MDCLATOMS_INFO *clatoms_info,
    MDCLATOMS_PIMD *clatoms_pimd,
    CPATOM_MAPS *cpatommaps,
    MDATOM_MAPS *atommaps,
    BUILD_INTRA *build_intra,
    MDINTRA *mdintra,
    NULL_INTER_PARSE *null_inter_parse,
    double *tot_memory,GENSIMOPTS *simopts)

  /*=======================================================================*/
  /*        Begin routine                                                 */
{/*begin routine*/
  /*=======================================================================*/
  /*           Local Variables                                             */

  int i,ntyp_pow,iii,inow,know;
  double now_memory;
  int natm_typ_mall;
  int natm_mall;
  int nghost_mall;
  int nfreeze_mall;
  int ncomp_mall;
  int nchrg_mall;
  int ngrp_21_mall;
  int ngrp_typ_21_mall;
  int ngrp_33_mall;
  int ngrp_typ_33_mall;
  int ngrp_watt_33_mall;
  int ngrp_typ_watt_33_mall;
  int ngrp_43_mall;
  int ngrp_typ_43_mall;
  int ngrp_23_mall;
  int ngrp_typ_23_mall;
  int ngrp_46_mall;
  int ngrp_typ_46_mall;
  int nbond_pow_mall;
  int nbond_typ_pow_mall;
  int nbond_con_mall;
  int nbond_typ_con_mall;
  int nbend_pow_mall;
  int nbend_typ_pow_mall;
  int nbend_con_mall;
  int nbend_typ_con_mall;
  int nbend_bnd_mall;
  int nbend_bnd_typ_mall;
  int ntors_pow_mall;
  int ntors_typ_pow_mall;
  int ntors_con_mall;
  int ntors_typ_con_mall;
  int nonfo_mall;
  int nonfo_typ_mall;
  int pimd_on;
  int nghost_old,nfreeze_old,ncomp_old;
  int pi_beads = clatoms_info->pi_beads; 
  int pimd_freez_typ = clatoms_pimd->pimd_freez_typ;

#include "../class_defs/allclass_strip_mdintra.h"


  /*======================================================================*/
  /* 0) Make the length of each vector an odd number                      */

  nghost_old          = build_intra->nghost_tot_max;
  nfreeze_old         = build_intra->nfreeze_max;
  ncomp_old           = NCOEF_GHOST_MAX;
  natm_mall           = clatoms_info->natm_tot;
  natm_typ_mall       = atommaps->natm_typ;
  nghost_mall         = mdghost_atoms->nghost_tot;
  nfreeze_mall        = mdconstrnt->nfreeze;
  ncomp_mall          = mdghost_atoms->natm_comp_max;
  ngrp_21_mall        = mdgrp_bond_con->num_21;
  ngrp_33_mall        = mdgrp_bond_con->num_33;
  ngrp_watt_33_mall   = mdgrp_bond_watts->num_33;
  ngrp_43_mall        = mdgrp_bond_con->num_43;
  ngrp_23_mall        = mdgrp_bond_con->num_23;
  ngrp_46_mall        = mdgrp_bond_con->num_46;
  ngrp_typ_21_mall    = mdgrp_bond_con->ntyp_21;
  ngrp_typ_33_mall    = mdgrp_bond_con->ntyp_33;
  ngrp_typ_watt_33_mall  = mdgrp_bond_watts->ntyp_33;
  ngrp_typ_43_mall    = mdgrp_bond_con->ntyp_43;
  ngrp_typ_23_mall    = mdgrp_bond_con->ntyp_23;
  ngrp_typ_46_mall    = mdgrp_bond_con->ntyp_46;
  nbond_pow_mall      = mdbond->npow;
  nbond_typ_pow_mall  = mdbond->ntyp_pow;
  nbond_con_mall      = mdbond->ncon;
  nbond_typ_con_mall  = mdbond->ntyp_con;
  nbend_pow_mall      = mdbend->npow;
  nbend_typ_pow_mall  = mdbend->ntyp_pow;
  nbend_bnd_mall      = mdbend_bnd->num;
  nbend_bnd_typ_mall  = mdbend_bnd->ntyp;
  ntors_pow_mall      = mdtors->npow;
  ntors_typ_pow_mall  = mdtors->ntyp_pow;
  nonfo_mall          = mdonfo->num;
  nonfo_typ_mall      = mdonfo->ntyp;

  /*======================================================================*/
  /* I) Atm stuff: Note constrain flag set here                           */
  /*               and charged atoms are counted                          */

  mdconstrnt->iconstrnt=0;
  if(mdbond->ncon > 0) mdconstrnt->iconstrnt=1;
  if(mdgrp_bond_con->num_21 > 0) mdconstrnt->iconstrnt=1;
  if(mdgrp_bond_con->num_23 > 0) mdconstrnt->iconstrnt=1;
  if(mdgrp_bond_con->num_33 > 0) mdconstrnt->iconstrnt=1;
  if(mdgrp_bond_watts->num_33 > 0) mdconstrnt->iconstrnt=1;
  if(mdgrp_bond_con->num_43 > 0) mdconstrnt->iconstrnt=1;
  if(mdgrp_bond_con->num_46 > 0) mdconstrnt->iconstrnt=1;

  clatoms_info->ichrg      = (int *)
    cmalloc(natm_mall*sizeof(int),"close_intra_params")-1;

  if(clatoms_info->cp_grimme==1){
    clatoms_info->c6Grimme     = (double *)
                cmalloc(natm_mall*sizeof(double),"close_intra_params")-1;
    clatoms_info->r0Grimme     = (double *)
                cmalloc(natm_mall*sizeof(double),"close_intra_params")-1;
  }//endif

  clatoms_info->nfree      = 0;
  clatoms_pimd->nfree_pimd = 0;
  for(i=1;i<=atommaps->nmol_typ;i++){
    inow                 = atommaps->nfree_1mol_jmol_typ[i]
      *atommaps->nmol_jmol_typ[i];
    know                 = 3*atommaps->natm_1mol_jmol_typ[i]
      *atommaps->nmol_jmol_typ[i];
    clatoms_info->nfree += inow;
    if(pimd_freez_typ==2){
      clatoms_pimd->nfree_pimd += inow*pi_beads;
    }else{
      clatoms_pimd->nfree_pimd += (inow + know*(pi_beads-1));
    }/*endif*/
  }/*endfor*/

  clatoms_info->nchrg = 0;
  for(i=1;i<=clatoms_info->natm_tot;i++){
    clatoms_info->mass[i] *= PROT_MASS;
    if(clatoms_info->q[i] != 0) {
      clatoms_info->nchrg++;
      clatoms_info->ichrg[(clatoms_info->nchrg)] = i;
    }
  }
  nchrg_mall = clatoms_info->nchrg;


  now_memory = ((double)(((natm_mall)*sizeof(double)*13+sizeof(int)*4)+
        natm_typ_mall*(sizeof(double)*0+sizeof(int)*25))
      *1.e-06);
  *tot_memory += now_memory;

  PRINTF("Atom allocation: %g Mbytes; Total memory: %g Mbytes\n",
      now_memory,*tot_memory);

  atommaps->atm_typ = (NAME *)crealloc(&(atommaps->atm_typ)[1],
      natm_typ_mall*sizeof(NAME),"close_intra_params")-1;
  clatoms_info->mass     = (double *)crealloc(&(clatoms_info->mass)[1],
      natm_mall*sizeof(double),"close_intra_params")-1;
  clatoms_info->q        = (double *)crealloc(&(clatoms_info->q)[1],
      natm_mall*sizeof(double),"close_intra_params")-1;
  cpatommaps->cp_vlnc_up  = (int *)crealloc(&(cpatommaps->cp_vlnc_up)[1],
      natm_mall*sizeof(int),"close_intra_params")-1;
  cpatommaps->cp_vlnc_dn  = (int *)crealloc(&(cpatommaps->cp_vlnc_dn)[1],
      natm_mall*sizeof(int),"close_intra_params")-1;
  cpatommaps->cp_vlnc_true_up  = (int *)crealloc(&(cpatommaps->cp_vlnc_true_up)[1],
      natm_mall*sizeof(int),"close_intra_params")-1;
  cpatommaps->cp_vlnc_true_dn  = (int *)crealloc(&(cpatommaps->cp_vlnc_true_dn)[1],
      natm_mall*sizeof(int),"close_intra_params")-1;
  cpatommaps->cp_atm_flag  = (int *)crealloc(&(cpatommaps->cp_atm_flag)[1],
      natm_mall*sizeof(int),"close_intra_params")-1;
  clatoms_info->ichrg  = (int *)crealloc(&(clatoms_info->ichrg)[1],
      nchrg_mall*sizeof(int),"close_intra_params")-1;
  clatoms_info->alp_pol  = (double *)crealloc(&(clatoms_info->alp_pol)[1],
      natm_mall*sizeof(double),"close_intra_params")-1;
  clatoms_info->b_neut   = (double *)crealloc(&(clatoms_info->b_neut)[1],
      natm_mall*sizeof(double),"close_intra_params")-1;
  clatoms_info->text_atm = (double *)crealloc(&(clatoms_info->text_atm)[1],
      natm_mall*sizeof(double),"close_intra_params")-1;
  atommaps->iatm_mol_typ  = (int *)crealloc(&(atommaps->iatm_mol_typ)[1],
      natm_mall*sizeof(int),"close_intra_params")-1;
  atommaps->iatm_atm_typ  = (int *)crealloc(&(atommaps->iatm_atm_typ)[1],
      natm_mall*sizeof(int),"close_intra_params")-1;
  atommaps->iatm_res_typ  = (int *)crealloc(&(atommaps->iatm_res_typ)[1],
      natm_mall*sizeof(int),"close_intra_params")-1;
  atommaps->iatm_mol_num  = (int *)crealloc(&(atommaps->iatm_mol_num)[1],
      natm_mall*sizeof(int),"close_intra_params")-1;
  atommaps->iatm_res_num  = (int *)crealloc(&(atommaps->iatm_res_num)[1],
      natm_mall*sizeof(int),"close_intra_params")-1;

  if(clatoms_info->pi_beads>1){
    clatoms_pimd->prekf = (double *)
      cmalloc(natm_mall*sizeof(double),"close_intra_params")-1;
  }

  /*========================================================================*/
  /* I) Ghost_atm stuff */

  mdghost_atoms->ighost_map  = (int *)crealloc(&(mdghost_atoms->ighost_map)[1],
      nghost_mall*sizeof(int),"close_intra_params")-1;
  mdghost_atoms->natm_comp   = (int *)crealloc(&(mdghost_atoms->natm_comp)[1],
      nghost_mall*sizeof(int),"close_intra_params")-1;
  mdghost_atoms->ighost_flag = (int *)crealloc(&(mdghost_atoms->ighost_flag)[1],
      natm_mall*sizeof(int),"close_intra_params")-1;

  /* need to set up reallocate matrix */
  mdghost_atoms->iatm_comp  = creall_int_mat(mdghost_atoms->iatm_comp,
      1,ncomp_old ,1,nghost_old,
      1,ncomp_mall,1,nghost_mall,"close_intra_params");
  mdghost_atoms->coef       = creall_mat(mdghost_atoms->coef,
      1,ncomp_old ,1,nghost_old,
      1,ncomp_mall,1,nghost_mall,"close_intra_params");

  /*========================================================================*/
  /* I) Freeze_atm stuff */

  mdconstrnt->freeze_map  = (int *) crealloc(&(mdconstrnt->freeze_map)[1],
      nfreeze_mall*sizeof(int),"close_intra_params")-1;
  mdconstrnt->freeze_flag = (int *)   crealloc(&(mdconstrnt->freeze_flag)[1],
      natm_mall*sizeof(int),"close_intra_params")-1;
  mdconstrnt->atom_label  = (int *)   crealloc(&(mdconstrnt->atom_label)[1],
      natm_mall*sizeof(int),"close_intra_params")-1;

  /*========================================================================*/
  /* II) Grp constrnts */



  mdgrp_bond_con->j1_21 = 
    (int *)crealloc(&(mdgrp_bond_con->j1_21)[1],
        ngrp_21_mall*sizeof(int),"close_intra_params")-1;
  mdgrp_bond_con->j2_21 = 
    (int *)   crealloc(&(mdgrp_bond_con->j2_21)[1],
        ngrp_21_mall*sizeof(int),"close_intra_params")-1;  
  mdgrp_bond_con->jtyp_21 = 
    (int *) crealloc(&(mdgrp_bond_con->jtyp_21)[1],
        ngrp_21_mall*sizeof(int),"close_intra_params")-1;
  mdgrp_bond_con->j1_33 = 
    (int *)crealloc(&(mdgrp_bond_con->j1_33)[1],
        ngrp_33_mall*sizeof(int),"close_intra_params")-1;
  mdgrp_bond_con->j2_33 = 
    (int *)   crealloc(&(mdgrp_bond_con->j2_33)[1],
        ngrp_33_mall*sizeof(int),"close_intra_params")-1;  
  mdgrp_bond_con->j3_33 = 
    (int *)   crealloc(&(mdgrp_bond_con->j3_33)[1],
        ngrp_33_mall*sizeof(int),"close_intra_params")-1;  
  mdgrp_bond_con->jtyp_33 = 
    (int *) crealloc(&(mdgrp_bond_con->jtyp_33)[1],
        ngrp_33_mall*sizeof(int),"close_intra_params")-1;
  mdgrp_bond_watts->j1_33 = 
    (int *)crealloc(&(mdgrp_bond_watts->j1_33)[1],
        ngrp_watt_33_mall*sizeof(int),"close_intra_params")-1;
  mdgrp_bond_watts->j2_33 = 
    (int *)   crealloc(&(mdgrp_bond_watts->j2_33)[1],
        ngrp_watt_33_mall*sizeof(int),"close_intra_params")-1;  
  mdgrp_bond_watts->j3_33 = 
    (int *)   crealloc(&(mdgrp_bond_watts->j3_33)[1],
        ngrp_watt_33_mall*sizeof(int),"close_intra_params")-1;  
  mdgrp_bond_watts->jtyp_33 = 
    (int *) crealloc(&(mdgrp_bond_watts->jtyp_33)[1],
        ngrp_watt_33_mall*sizeof(int),"close_intra_params")-1;
  mdgrp_bond_con->j1_23 = 
    (int *)   crealloc(&(mdgrp_bond_con->j1_23)[1],
        ngrp_23_mall*sizeof(int),"close_intra_params")-1;  
  mdgrp_bond_con->j2_23 = 
    (int *)   crealloc(&(mdgrp_bond_con->j2_23)[1],
        ngrp_23_mall*sizeof(int),"close_intra_params")-1;  
  mdgrp_bond_con->j3_23 = 
    (int *)   crealloc(&(mdgrp_bond_con->j3_23)[1],
        ngrp_23_mall*sizeof(int),"close_intra_params")-1;  
  mdgrp_bond_con->jtyp_23 = 
    (int *) crealloc(&(mdgrp_bond_con->jtyp_23)[1],
        ngrp_23_mall*sizeof(int),"close_intra_params")-1;
  mdgrp_bond_con->j1_43 = 
    (int *)crealloc(&(mdgrp_bond_con->j1_43)[1],
        ngrp_43_mall*sizeof(int),"close_intra_params")-1;
  mdgrp_bond_con->j2_43 = 
    (int *)   crealloc(&(mdgrp_bond_con->j2_43)[1],
        ngrp_43_mall*sizeof(int),"close_intra_params")-1;  
  mdgrp_bond_con->j3_43 = 
    (int *)   crealloc(&(mdgrp_bond_con->j3_43)[1],
        ngrp_43_mall*sizeof(int),"close_intra_params")-1;  
  mdgrp_bond_con->j4_43 = 
    (int *)   crealloc(&(mdgrp_bond_con->j4_43)[1],
        ngrp_43_mall*sizeof(int),"close_intra_params")-1;  
  mdgrp_bond_con->jtyp_43 = 
    (int *) crealloc(&(mdgrp_bond_con->jtyp_43)[1],
        ngrp_43_mall*sizeof(int),"close_intra_params")-1;
  mdgrp_bond_con->j1_46 = 
    (int *)   crealloc(&(mdgrp_bond_con->j1_46)[1],
        ngrp_46_mall*sizeof(int),"close_intra_params")-1;  
  mdgrp_bond_con->j2_46 = 
    (int *)   crealloc(&(mdgrp_bond_con->j2_46)[1],
        ngrp_46_mall*sizeof(int),"close_intra_params")-1;  
  mdgrp_bond_con->j3_46 = 
    (int *)   crealloc(&(mdgrp_bond_con->j3_46)[1],
        ngrp_46_mall*sizeof(int),"close_intra_params")-1;  
  mdgrp_bond_con->j4_46 = 
    (int *)   crealloc(&(mdgrp_bond_con->j4_46)[1],
        ngrp_46_mall*sizeof(int),"close_intra_params")-1;  
  mdgrp_bond_con->jtyp_46 = 
    (int *) crealloc(&(mdgrp_bond_con->jtyp_46)[1],
        ngrp_46_mall*sizeof(int),"close_intra_params")-1;
  mdgrp_bond_con->al_21 = cmall_mat(1,1,1,ngrp_21_mall,"close_intra_params");
  mdgrp_bond_con->al_33 = cmall_mat(1,3,1,ngrp_33_mall,"close_intra_params");
  mdgrp_bond_con->al_43 = cmall_mat(1,3,1,ngrp_43_mall,"close_intra_params");
  mdgrp_bond_con->al_23 = cmall_mat(1,2,1,ngrp_23_mall,"close_intra_params");
  mdgrp_bond_con->al_46 = cmall_mat(1,6,1,ngrp_46_mall,"close_intra_params");
  mdgrp_bond_con->eq_21 = 
    cmall_mat(1,1,1,ngrp_typ_21_mall,"close_intra_params");
  mdgrp_bond_con->eq_33 = 
    cmall_mat(1,3,1,ngrp_typ_33_mall,"close_intra_params");
  mdgrp_bond_watts->eq_33 = 
    cmall_mat(1,3,1,ngrp_typ_watt_33_mall,"close_intra_params");
  mdgrp_bond_con->eq_43   = 
    cmall_mat(1,3,1,ngrp_typ_43_mall,"close_intra_params");
  mdgrp_bond_con->eq_23   = 
    cmall_mat(1,2,1,ngrp_typ_23_mall,"close_intra_params");
  mdgrp_bond_con->eq_46   = 
    cmall_mat(1,6,1,ngrp_typ_46_mall,"close_intra_params");

  mdgrp_bond_watts->cos_thet0_2   = 
    (double *)cmalloc(ngrp_typ_watt_33_mall*sizeof(int),"close_intra_params")-1;
  mdgrp_bond_watts->sin_thet0_2   = 
    (double *)cmalloc(ngrp_typ_watt_33_mall*sizeof(int),"close_intra_params")-1;
  mdgrp_bond_watts->c_0_33    = 
    cmall_mat(1,3,1,ngrp_typ_watt_33_mall,"close_intra_params");
  mdgrp_bond_watts->c_1_33  = 
    cmall_mat(1,3,1,ngrp_typ_watt_33_mall,"close_intra_params");
  mdgrp_bond_watts->c_2_33  = 
    cmall_mat(1,3,1,ngrp_typ_watt_33_mall,"close_intra_params");
  mdgrp_bond_watts->c_3_33  = 
    cmall_mat(1,3,1,ngrp_typ_watt_33_mall,"close_intra_params");
  mdgrp_bond_watts->c_4_33  = 
    cmall_mat(1,3,1,ngrp_typ_watt_33_mall,"close_intra_params");
  mdgrp_bond_watts->c_5_33  = 
    cmall_mat(1,3,1,ngrp_typ_watt_33_mall,"close_intra_params");
  mdgrp_bond_watts->c_6_33  = 
    cmall_mat(1,3,1,ngrp_typ_watt_33_mall,"close_intra_params");
  mdgrp_bond_watts->dc_0_33 = 
    cmall_mat(1,3,1,ngrp_typ_watt_33_mall,"close_intra_params");
  mdgrp_bond_watts->dc_1_33 = 
    cmall_mat(1,3,1,ngrp_typ_watt_33_mall,"close_intra_params");
  mdgrp_bond_watts->dc_2_33 = 
    cmall_mat(1,3,1,ngrp_typ_watt_33_mall,"close_intra_params");
  mdgrp_bond_watts->dc_3_33 = 
    cmall_mat(1,3,1,ngrp_typ_watt_33_mall,"close_intra_params");
  mdgrp_bond_watts->dc_4_33 = 
    cmall_mat(1,3,1,ngrp_typ_watt_33_mall,"close_intra_params");
  mdgrp_bond_watts->dc_5_33 = 
    cmall_mat(1,3,1,ngrp_typ_watt_33_mall,"close_intra_params");
  mdgrp_bond_watts->dc_6_33 = 
    cmall_mat(1,3,1,ngrp_typ_watt_33_mall,"close_intra_params");

  /*======================================================================*/
  /*   II) Bonds  */

  now_memory = ((nbond_pow_mall)*(sizeof(double)*0+sizeof(int)*3)+
      (nbond_con_mall)*(sizeof(double)*1 +  sizeof(int)*3)+
      (nbond_typ_pow_mall)*(sizeof(double)*14+ sizeof(int)*0 )+
      (nbond_typ_con_mall)*(sizeof(double)*1 +  sizeof(int)*0)+
      (nbend_pow_mall)*(sizeof(double)*0 +  sizeof(int)*4)+
      (nbend_typ_pow_mall)*(sizeof(double)*29+  sizeof(int)*0)+
      (nbend_bnd_mall)*(sizeof(double)*0 +  sizeof(int)*4)+
      (nbend_bnd_typ_mall)*(sizeof(double)*44 +  sizeof(int)*0)+
      (ntors_pow_mall)*(sizeof(double)*0 +  sizeof(int)*5)+
      (ntors_typ_pow_mall)*(sizeof(double)*30 +  sizeof(int)*0)+
      (nonfo_mall)*(sizeof(double)*0 +  sizeof(int)*3)+
      (nonfo_typ_mall)*(sizeof(double)*3 +  sizeof(int)*0)
      )*1.e-06;

  *tot_memory += now_memory;
  PRINTF("Intramolecular allocation: %g Mbytes; Total memory: %g Mbytes\n",
      now_memory,*tot_memory);
  ntyp_pow = nbond_typ_pow_mall;


  mdbond->j1_pow = (int *)    crealloc(&(mdbond->j1_pow)[1],
      nbond_pow_mall*sizeof(int),"close_intra_params")-1;
  mdbond->j2_pow = (int *)    crealloc(&(mdbond->j2_pow)[1],
      nbond_pow_mall*sizeof(int),"close_intra_params")-1;
  mdbond->jtyp_pow   = (int *)crealloc(&(mdbond->jtyp_pow)[1],
      nbond_pow_mall*sizeof(int),"close_intra_params")-1;
  mdbond->j1_con = (int *)    crealloc(&(mdbond->j1_con)[1],
      nbond_con_mall*sizeof(int),"close_intra_params")-1;
  mdbond->j2_con = (int *)    crealloc(&(mdbond->j2_con)[1],
      nbond_con_mall*sizeof(int),"close_intra_params")-1;
  mdbond->jtyp_con   = (int *)crealloc(&(mdbond->jtyp_con)[1],
      nbond_con_mall*sizeof(int),"close_intra_params")-1;
  mdbond->al_con = (double *)
    cmalloc(nbond_con_mall*sizeof(double),"close_intra_params")-1;
  mdbond->eq_pow = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdbond->c_0    = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdbond->c_1    = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdbond->c_2    = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdbond->c_3    = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdbond->c_4    = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdbond->c_5    = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdbond->c_6    = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdbond->dc_0   = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdbond->dc_1   = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdbond->dc_2   = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdbond->dc_3   = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdbond->dc_4   = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdbond->dc_5   = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdbond->dc_6   = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdbond->eq_con = (double *)cmalloc(nbond_typ_con_mall
      *sizeof(double),"close_intra_params")-1;

  null_inter_parse->jbond1_nul= (int *)
    crealloc(&(null_inter_parse->jbond1_nul)[1],
        null_inter_parse->nbond_nul*sizeof(int),"close_intra_params")-1;
  null_inter_parse->jbond2_nul= (int *)
    crealloc(&(null_inter_parse->jbond2_nul)[1],
        null_inter_parse->nbond_nul*sizeof(int),"close_intra_params")-1;


  /*=======================================================================*/
  /*  III) Bends  */

  ntyp_pow = nbend_typ_pow_mall;
  mdbend->j1_pow = (int *)   crealloc(&(mdbend->j1_pow)[1],
      nbend_pow_mall*sizeof(int),"close_intra_params")-1;
  mdbend->j2_pow = (int *)   crealloc(&(mdbend->j2_pow)[1],
      nbend_pow_mall*sizeof(int),"close_intra_params")-1;
  mdbend->j3_pow = (int *)   crealloc(&(mdbend->j3_pow)[1],
      nbend_pow_mall*sizeof(int),"close_intra_params")-1;
  mdbend->jtyp_pow = (int *)crealloc(&(mdbend->jtyp_pow)[1],
      nbend_pow_mall*sizeof(int),"close_intra_params")-1;
  mdbend->eq_pow = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdbend->c_0    = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdbend->c_1    = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdbend->c_2    = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdbend->c_3    = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdbend->c_4    = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdbend->c_5    = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdbend->c_6    = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdbend->s_0    = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdbend->s_1    = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdbend->s_2    = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdbend->s_3    = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdbend->s_4    = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdbend->s_5    = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdbend->s_6    = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdbend->dc_0   = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdbend->dc_1   = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdbend->dc_2   = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdbend->dc_3   = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdbend->dc_4   = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdbend->dc_5   = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdbend->dc_6   = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdbend->ds_0   = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdbend->ds_1   = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdbend->ds_2   = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdbend->ds_3   = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdbend->ds_4   = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdbend->ds_5   = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdbend->ds_6   = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  null_inter_parse->jbend1_nul= (int *)  
    crealloc(&(null_inter_parse->jbend1_nul)[1],
        null_inter_parse->nbend_nul*sizeof(int),"close_intra_params")-1;
  null_inter_parse->jbend2_nul= (int *)  
    crealloc(&(null_inter_parse->jbend2_nul)[1],
        null_inter_parse->nbend_nul*sizeof(int),"close_intra_params")-1;
  null_inter_parse->jbend3_nul= (int *)  
    crealloc(&(null_inter_parse->jbend3_nul)[1],
        null_inter_parse->nbend_nul*sizeof(int),"close_intra_params")-1;
  /*=======================================================================*/
  /*  III) Bend_bnds  */

  mdbend_bnd->j1 = (int *)   crealloc(&(mdbend_bnd->j1)[1],
      nbend_bnd_mall*sizeof(int),"close_intra_params")-1;
  mdbend_bnd->j2 = (int *)   crealloc(&(mdbend_bnd->j2)[1],
      nbend_bnd_mall*sizeof(int),"close_intra_params")-1;
  mdbend_bnd->j3 = (int *)   crealloc(&(mdbend_bnd->j3)[1],
      nbend_bnd_mall*sizeof(int),"close_intra_params")-1;
  mdbend_bnd->jtyp = (int *)crealloc(&(mdbend_bnd->jtyp)[1],
      nbend_bnd_mall*sizeof(int),"close_intra_params")-1;
  mdbend_bnd->eq_bond = (double *)cmalloc(nbend_bnd_typ_mall*
      sizeof(double),"close_intra_params")-1;
  mdbend_bnd->eq_bend = (double *)cmalloc(nbend_bnd_typ_mall*
      sizeof(double),"close_intra_params")-1;
  mdbend_bnd->cbond_0    = (double *)cmalloc(nbend_bnd_typ_mall
      *sizeof(double),"close_intra_params")-1;
  mdbend_bnd->cbond_1    = (double *)cmalloc(nbend_bnd_typ_mall
      *sizeof(double),"close_intra_params")-1;
  mdbend_bnd->cbond_2    = (double *)cmalloc(nbend_bnd_typ_mall
      *sizeof(double),"close_intra_params")-1;
  mdbend_bnd->cbond_3    = (double *)cmalloc(nbend_bnd_typ_mall
      *sizeof(double),"close_intra_params")-1;
  mdbend_bnd->cbond_4    = (double *)cmalloc(nbend_bnd_typ_mall
      *sizeof(double),"close_intra_params")-1;
  mdbend_bnd->cbond_5    = (double *)cmalloc(nbend_bnd_typ_mall
      *sizeof(double),"close_intra_params")-1;
  mdbend_bnd->cbond_6    = (double *)cmalloc(nbend_bnd_typ_mall
      *sizeof(double),"close_intra_params")-1;
  mdbend_bnd->dcbond_0   = (double *)cmalloc(nbend_bnd_typ_mall
      *sizeof(double),"close_intra_params")-1;
  mdbend_bnd->dcbond_1   = (double *)cmalloc(nbend_bnd_typ_mall
      *sizeof(double),"close_intra_params")-1;
  mdbend_bnd->dcbond_2   = (double *)cmalloc(nbend_bnd_typ_mall
      *sizeof(double),"close_intra_params")-1;
  mdbend_bnd->dcbond_3   = (double *)cmalloc(nbend_bnd_typ_mall
      *sizeof(double),"close_intra_params")-1;
  mdbend_bnd->dcbond_4   = (double *)cmalloc(nbend_bnd_typ_mall
      *sizeof(double),"close_intra_params")-1;
  mdbend_bnd->dcbond_5   = (double *)cmalloc(nbend_bnd_typ_mall
      *sizeof(double),"close_intra_params")-1;
  mdbend_bnd->dcbond_6   = (double *)cmalloc(nbend_bnd_typ_mall
      *sizeof(double),"close_intra_params")-1;

  mdbend_bnd->cbend_0    = (double *)cmalloc(nbend_bnd_typ_mall
      *sizeof(double),"close_intra_params")-1;
  mdbend_bnd->cbend_1    = (double *)cmalloc(nbend_bnd_typ_mall
      *sizeof(double),"close_intra_params")-1;
  mdbend_bnd->cbend_2    = (double *)cmalloc(nbend_bnd_typ_mall
      *sizeof(double),"close_intra_params")-1;
  mdbend_bnd->cbend_3    = (double *)cmalloc(nbend_bnd_typ_mall
      *sizeof(double),"close_intra_params")-1;
  mdbend_bnd->cbend_4    = (double *)cmalloc(nbend_bnd_typ_mall
      *sizeof(double),"close_intra_params")-1;
  mdbend_bnd->cbend_5    = (double *)cmalloc(nbend_bnd_typ_mall
      *sizeof(double),"close_intra_params")-1;
  mdbend_bnd->cbend_6    = (double *)cmalloc(nbend_bnd_typ_mall
      *sizeof(double),"close_intra_params")-1;
  mdbend_bnd->sbend_0    = (double *)cmalloc(nbend_bnd_typ_mall
      *sizeof(double),"close_intra_params")-1;
  mdbend_bnd->sbend_1    = (double *)cmalloc(nbend_bnd_typ_mall
      *sizeof(double),"close_intra_params")-1;
  mdbend_bnd->sbend_2    = (double *)cmalloc(nbend_bnd_typ_mall
      *sizeof(double),"close_intra_params")-1;
  mdbend_bnd->sbend_3    = (double *)cmalloc(nbend_bnd_typ_mall
      *sizeof(double),"close_intra_params")-1;
  mdbend_bnd->sbend_4    = (double *)cmalloc(nbend_bnd_typ_mall
      *sizeof(double),"close_intra_params")-1;
  mdbend_bnd->sbend_5    = (double *)cmalloc(nbend_bnd_typ_mall
      *sizeof(double),"close_intra_params")-1;
  mdbend_bnd->sbend_6    = (double *)cmalloc(nbend_bnd_typ_mall
      *sizeof(double),"close_intra_params")-1;
  mdbend_bnd->dcbend_0   = (double *)cmalloc(nbend_bnd_typ_mall
      *sizeof(double),"close_intra_params")-1;
  mdbend_bnd->dcbend_1   = (double *)cmalloc(nbend_bnd_typ_mall
      *sizeof(double),"close_intra_params")-1;
  mdbend_bnd->dcbend_2   = (double *)cmalloc(nbend_bnd_typ_mall
      *sizeof(double),"close_intra_params")-1;
  mdbend_bnd->dcbend_3   = (double *)cmalloc(nbend_bnd_typ_mall
      *sizeof(double),"close_intra_params")-1;
  mdbend_bnd->dcbend_4   = (double *)cmalloc(nbend_bnd_typ_mall
      *sizeof(double),"close_intra_params")-1;
  mdbend_bnd->dcbend_5   = (double *)cmalloc(nbend_bnd_typ_mall
      *sizeof(double),"close_intra_params")-1;
  mdbend_bnd->dcbend_6   = (double *)cmalloc(nbend_bnd_typ_mall
      *sizeof(double),"close_intra_params")-1;
  mdbend_bnd->dsbend_0   = (double *)cmalloc(nbend_bnd_typ_mall
      *sizeof(double),"close_intra_params")-1;
  mdbend_bnd->dsbend_1   = (double *)cmalloc(nbend_bnd_typ_mall
      *sizeof(double),"close_intra_params")-1;
  mdbend_bnd->dsbend_2   = (double *)cmalloc(nbend_bnd_typ_mall
      *sizeof(double),"close_intra_params")-1;
  mdbend_bnd->dsbend_3   = (double *)cmalloc(nbend_bnd_typ_mall
      *sizeof(double),"close_intra_params")-1;
  mdbend_bnd->dsbend_4   = (double *)cmalloc(nbend_bnd_typ_mall
      *sizeof(double),"close_intra_params")-1;
  mdbend_bnd->dsbend_5   = (double *)cmalloc(nbend_bnd_typ_mall
      *sizeof(double),"close_intra_params")-1;
  mdbend_bnd->dsbend_6   = (double *)cmalloc(nbend_bnd_typ_mall
      *sizeof(double),"close_intra_params")-1;

  /*======================================================================*/
  /*  IV) Tors  */

  ntyp_pow = ntors_typ_pow_mall;
  mdtors->j1_pow = (int *)   crealloc(&(mdtors->j1_pow)[1],
      ntors_pow_mall*sizeof(int),"close_intra_params")-1;
  mdtors->j2_pow = (int *)   crealloc(&(mdtors->j2_pow)[1],
      ntors_pow_mall*sizeof(int),"close_intra_params")-1;
  mdtors->j3_pow = (int *)   crealloc(&(mdtors->j3_pow)[1],
      ntors_pow_mall*sizeof(int),"close_intra_params")-1;
  mdtors->j4_pow = (int *)   crealloc(&(mdtors->j4_pow)[1],
      ntors_pow_mall*sizeof(int),"close_intra_params")-1;
  mdtors->jtyp_pow = (int *)crealloc(&(mdtors->jtyp_pow)[1],
      ntors_pow_mall*sizeof(int),"close_intra_params")-1;
  mdtors->eq_pow = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdtors->c_0    = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdtors->c_1    = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdtors->c_2    = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdtors->c_3    = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdtors->c_4    = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdtors->c_5    = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdtors->c_6    = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdtors->s_0    = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdtors->s_1    = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdtors->s_2    = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdtors->s_3    = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdtors->s_4    = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdtors->s_5    = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdtors->s_6    = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdtors->dc_0   = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdtors->dc_1   = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdtors->dc_2   = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdtors->dc_3   = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdtors->dc_4   = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdtors->dc_5   = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdtors->dc_6   = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdtors->ds_0   = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdtors->ds_1   = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdtors->ds_2   = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdtors->ds_3   = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdtors->ds_4   = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdtors->ds_5   = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  mdtors->ds_6   = (double *)
    cmalloc(ntyp_pow*sizeof(double),"close_intra_params")-1;
  null_inter_parse->jtors1_nul= (int *) 
    crealloc(&(null_inter_parse->jtors1_nul)[1],
        null_inter_parse->ntors_nul*sizeof(int),"close_intra_params")-1;
  null_inter_parse->jtors2_nul= (int *) 
    crealloc(&(null_inter_parse->jtors2_nul)[1],
        null_inter_parse->ntors_nul*sizeof(int),"close_intra_params")-1;
  null_inter_parse->jtors3_nul= (int *) 
    crealloc(&(null_inter_parse->jtors3_nul)[1],
        null_inter_parse->ntors_nul*sizeof(int),"close_intra_params")-1;
  null_inter_parse->jtors4_nul= (int *) 
    crealloc(&(null_inter_parse->jtors4_nul)[1],
        null_inter_parse->ntors_nul*sizeof(int),"close_intra_params")-1;

  /*=======================================================================*/
  /*   V) Onfo  */

  mdonfo->j1     = (int *)   crealloc(&(mdonfo->j1)[1],
      nonfo_mall*sizeof(int),"close_intra_params")-1;
  mdonfo->j2     = (int *)   crealloc(&(mdonfo->j2)[1],
      nonfo_mall*sizeof(int),"close_intra_params")-1;
  mdonfo->jtyp   = (int *)   crealloc(&(mdonfo->jtyp)[1],
      nonfo_mall*sizeof(int),"close_intra_params")-1;
  mdonfo->feps   = (double *)
    cmalloc(nonfo_typ_mall*sizeof(double),"close_intra_params")-1;
  mdonfo->s6 = (double *)
    cmalloc(nonfo_typ_mall*sizeof(double),"close_intra_params")-1;
  mdonfo->sc = (double *)
    cmalloc(nonfo_typ_mall*sizeof(double),"close_intra_params")-1;
  null_inter_parse->jonfo1_nul= (int *) 
    crealloc(&(null_inter_parse->jonfo1_nul)[1],
        null_inter_parse->nonfo_nul*sizeof(int),"close_intra_params")-1;
  null_inter_parse->jonfo2_nul  = (int *) 
    crealloc(&(null_inter_parse->jonfo2_nul)[1],
        null_inter_parse->nonfo_nul*sizeof(int),"close_intra_params")-1;

  /*-----------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/












