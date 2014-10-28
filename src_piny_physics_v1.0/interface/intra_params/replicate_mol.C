/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: Build_map.c                                  */
/*                                                                          */
/* This subprogram builds some atoms maps                                   */
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
#include "../proto_defs/proto_intra_params_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void replicate_mol(MDCLATOMS_INFO *clatoms_info,
    CPATOM_MAPS *cpatom_maps,MDATOM_MAPS *atommaps,
    BUILD_INTRA *build_intra,MDINTRA *mdintra,
    NULL_INTER_PARSE *null_inter_parse,
    START_INDEX *start_index,int jmol_typ)

  /*========================================================================*/
  /*     Begin routine                                                      */
{/*begin routine*/

  /*========================================================================*/
  /*     Local Variables                                                    */

#include "../class_defs/allclass_strip_mdintra.h"   

  int imol,i,ishift,nmol,nmol1,k,iii,ig_shift,if_shift;
  int nbond_pow_add,nbond_con_add,nbond_null_add;
  int nbend_pow_add,nbend_null_add;
  int ntors_pow_add,ntors_null_add;
  int nonfo_add,nonfo_null_add,nbend_bnd_add;
  int natm_add;
  int nghost_add;
  int nfreeze_add;
  int ngrp_21_add; 
  int ngrp_33_add; 
  int ngrp_watt_33_add; 
  int ngrp_43_add; 
  int ngrp_23_add;
  int ngrp_46_add;
  int ibond_pow_off,ibond_con_off,ibond_null_off;
  int ibend_pow_off,ibend_null_off;
  int itors_pow_off,itors_null_off;
  int ionfo_off,ionfo_null_off,ibend_bnd_off;
  int iatm_off;
  int ighost_off;
  int ifreeze_off;
  int igrp_21_off; 
  int igrp_33_off; 
  int igrp_watt_33_off; 
  int igrp_43_off; 
  int igrp_23_off;
  int igrp_46_off;
  int ibond_pow_offn,ibond_con_offn,ibond_null_offn;
  int ibend_pow_offn,ibend_null_offn;
  int itors_pow_offn,itors_null_offn;
  int ionfo_offn,ionfo_null_offn,ibend_bnd_offn;
  int iatm_offn;
  int ighost_offn;
  int ifreeze_offn;
  int igrp_21_offn; 
  int igrp_33_offn; 
  int igrp_watt_33_offn; 
  int igrp_43_offn; 
  int igrp_23_offn;
  int igrp_46_offn;

  /*========================================================================*/
  /* I) Calculate how many of what to replicate                             */

  nbond_pow_add = mdbond->npow            - start_index->nbond_pow;
  nbond_con_add = mdbond->ncon            - start_index->nbond_con;
  nbond_null_add = null_inter_parse->nbond_nul - start_index->nbond_nul;
  nbend_pow_add = mdbend->npow            - start_index->nbend_pow;
  nbend_null_add = null_inter_parse->nbend_nul - start_index->nbend_nul;
  ntors_pow_add = mdtors->npow            - start_index->ntors_pow;
  ntors_null_add = null_inter_parse->ntors_nul - start_index->ntors_nul;
  nonfo_add     = mdonfo->num             - start_index->nonfo;
  nonfo_null_add = null_inter_parse->nonfo_nul - start_index->nonfo_nul;
  nbend_bnd_add = mdbend_bnd->num         - start_index->nbend_bnd;
  natm_add      = clatoms_info->natm_tot       - start_index->natm;
  nghost_add    = mdghost_atoms->nghost_tot      - start_index->nghost_tot;
  nfreeze_add   = mdconstrnt->nfreeze      - start_index->nfreeze;
  ngrp_21_add   = mdgrp_bond_con->num_21         - start_index->ngrp_21;
  ngrp_23_add   = mdgrp_bond_con->num_23         - start_index->ngrp_23;
  ngrp_43_add   = mdgrp_bond_con->num_43         - start_index->ngrp_43;
  ngrp_33_add   = mdgrp_bond_con->num_33         - start_index->ngrp_33;
  ngrp_watt_33_add   = mdgrp_bond_watts->num_33  
    - start_index->ngrp_watt_33;
  ngrp_46_add   = mdgrp_bond_con->num_46         - start_index->ngrp_46;


  /*========================================================================*/
  /* I) Rellocate                                                 */

  reallocate_intra_list(clatoms_info,cpatom_maps,
      atommaps,build_intra,mdintra,null_inter_parse,
      jmol_typ,nbond_pow_add, nbond_con_add, nbond_null_add,
      nbend_pow_add, nbend_null_add,
      ntors_pow_add, ntors_null_add,
      nonfo_add, nonfo_null_add, nbend_bnd_add,
      natm_add,nghost_add, nfreeze_add,ngrp_43_add,
      ngrp_33_add,ngrp_watt_33_add,
      ngrp_21_add,ngrp_23_add, ngrp_46_add);

  /*========================================================================*/
  /* II) Store the offsets                                                 */

  ibond_pow_off  = start_index->nbond_pow;
  ibond_con_off  = start_index->nbond_con;
  ibond_null_off = start_index->nbond_nul;
  ibend_pow_off  = start_index->nbend_pow;
  ibend_null_off = start_index->nbend_nul;
  itors_pow_off  = start_index->ntors_pow;
  itors_null_off = start_index->ntors_nul;
  ionfo_off      = start_index->nonfo;
  ionfo_null_off = start_index->nonfo_nul;
  ibend_bnd_off  = start_index->nbend_bnd;
  iatm_off       = start_index->natm;
  ighost_off     = start_index->nghost_tot;
  ifreeze_off    = start_index->nfreeze;
  igrp_21_off    = start_index->ngrp_21;
  igrp_23_off    = start_index->ngrp_23;
  igrp_33_off    = start_index->ngrp_33;
  igrp_watt_33_off    = start_index->ngrp_watt_33;
  igrp_43_off    = start_index->ngrp_43;
  igrp_46_off    = start_index->ngrp_46;

  ibond_pow_offn  = start_index->nbond_pow;
  ibond_con_offn  = start_index->nbond_con;
  ibond_null_offn = start_index->nbond_nul;
  ibend_pow_offn  = start_index->nbend_pow;
  ibend_null_offn = start_index->nbend_nul;
  itors_pow_offn  = start_index->ntors_pow;
  itors_null_offn = start_index->ntors_nul;
  ionfo_offn      = start_index->nonfo;
  ionfo_null_offn = start_index->nonfo_nul;
  ibend_bnd_offn  = start_index->nbend_bnd;
  iatm_offn       = start_index->natm;
  ighost_offn     = start_index->nghost_tot;
  ifreeze_offn     = start_index->nfreeze;
  igrp_21_offn    = start_index->ngrp_21;
  igrp_23_offn    = start_index->ngrp_23;
  igrp_33_offn    = start_index->ngrp_33;
  igrp_watt_33_offn    = start_index->ngrp_watt_33;
  igrp_43_offn    = start_index->ngrp_43;
  igrp_46_offn    = start_index->ngrp_46;

  /*========================================================================*/
  /* III) Replicate the bonds                                               */

  nmol = atommaps->nmol_jmol_typ[jmol_typ];
  for(imol=1;imol<=(nmol-1);imol++){
    ishift = imol*natm_add;
    ibond_pow_offn += nbond_pow_add;
    for(i=1;i<=nbond_pow_add;i++){
      mdbond->j1_pow[(i+ibond_pow_offn)] = 
        mdbond->j1_pow[(i+ibond_pow_off)]+ishift;
      mdbond->j2_pow[(i+ibond_pow_offn)] = 
        mdbond->j2_pow[(i+ibond_pow_off)]+ishift;
      mdbond->jtyp_pow[(i+ibond_pow_offn)] = 
        mdbond->jtyp_pow[(i+ibond_pow_off)];
    }/*endfor*/
    ibond_con_offn += nbond_con_add;
    for(i=1;i<=nbond_con_add;i++){
      mdbond->j1_con[(i+ibond_con_offn)] = 
        mdbond->j1_con[(i+ibond_con_off)]+ishift;
      mdbond->j2_con[(i+ibond_con_offn)] = 
        mdbond->j2_con[(i+ibond_con_off)]+ishift;
      mdbond->jtyp_con[(i+ibond_con_offn)] = 
        mdbond->jtyp_con[(i+ibond_con_off)];
    }/*endfor*/
    ibond_null_offn += nbond_null_add;
    for(i=1;i<=nbond_null_add;i++){
      null_inter_parse->jbond1_nul[(i+ibond_null_offn)] = 
        null_inter_parse->jbond1_nul[(i+ibond_null_off)]+ishift;
      null_inter_parse->jbond2_nul[(i+ibond_null_offn)] = 
        null_inter_parse->jbond2_nul[(i+ibond_null_off)]+ishift;
    }/*endfor*/
  }/*endfor*/

  /*========================================================================*/
  /* IV) Replicate the bends                                                */

  nmol = atommaps->nmol_jmol_typ[jmol_typ];
  for(imol=1;imol<=(nmol-1);imol++){
    ishift = imol*natm_add;
    ibend_pow_offn += nbend_pow_add;
    for(i=1;i<=nbend_pow_add;i++){
      mdbend->j1_pow[(i+ibend_pow_offn)] = 
        mdbend->j1_pow[(i+ibend_pow_off)]+ishift;
      mdbend->j2_pow[(i+ibend_pow_offn)] = 
        mdbend->j2_pow[(i+ibend_pow_off)]+ishift;
      mdbend->j3_pow[(i+ibend_pow_offn)] = 
        mdbend->j3_pow[(i+ibend_pow_off)]+ishift;
      mdbend->jtyp_pow[(i+ibend_pow_offn)] = 
        mdbend->jtyp_pow[(i+ibend_pow_off)];
    }/*endfor*/
    ibend_null_offn += nbend_null_add;
    for(i=1;i<=nbend_null_add;i++){
      null_inter_parse->jbend1_nul[(i+ibend_null_offn)] = 
        null_inter_parse->jbend1_nul[(i+ibend_null_off)]+ishift;
      null_inter_parse->jbend2_nul[(i+ibend_null_offn)] = 
        null_inter_parse->jbend2_nul[(i+ibend_null_off)]+ishift;
      null_inter_parse->jbend3_nul[(i+ibend_null_offn)] = 
        null_inter_parse->jbend3_nul[(i+ibend_null_off)]+ishift;
    }/*endfor*/
  }/*endfor*/

  /*========================================================================*/
  /* V) Replicate the torsions                                              */

  nmol = atommaps->nmol_jmol_typ[jmol_typ];
  for(imol=1;imol<=(nmol-1);imol++){
    ishift = imol*natm_add;
    itors_pow_offn += ntors_pow_add;
    for(i=1;i<=ntors_pow_add;i++){
      mdtors->j1_pow[(i+itors_pow_offn)] = 
        mdtors->j1_pow[(i+itors_pow_off)]+ishift;
      mdtors->j2_pow[(i+itors_pow_offn)] = 
        mdtors->j2_pow[(i+itors_pow_off)]+ishift;
      mdtors->j3_pow[(i+itors_pow_offn)] = 
        mdtors->j3_pow[(i+itors_pow_off)]+ishift;
      mdtors->j4_pow[(i+itors_pow_offn)] = 
        mdtors->j4_pow[(i+itors_pow_off)]+ishift;
      mdtors->jtyp_pow[(i+itors_pow_offn)] = 
        mdtors->jtyp_pow[(i+itors_pow_off)];
    }/*endfor*/
    itors_null_offn += ntors_null_add;
    for(i=1;i<=ntors_null_add;i++){
      null_inter_parse->jtors1_nul[(i+itors_null_offn)] = 
        null_inter_parse->jtors1_nul[(i+itors_null_off)]+ishift;
      null_inter_parse->jtors2_nul[(i+itors_null_offn)] = 
        null_inter_parse->jtors2_nul[(i+itors_null_off)]+ishift;
      null_inter_parse->jtors3_nul[(i+itors_null_offn)] = 
        null_inter_parse->jtors3_nul[(i+itors_null_off)]+ishift;
      null_inter_parse->jtors4_nul[(i+itors_null_offn)] = 
        null_inter_parse->jtors4_nul[(i+itors_null_off)]+ishift;
    }/*endfor*/
  }/*endfor*/

  /*========================================================================*/
  /* VI) Replicate the onfos                                                */

  nmol = atommaps->nmol_jmol_typ[jmol_typ];
  for(imol=1;imol<=(nmol-1);imol++){
    ishift = imol*natm_add;
    ionfo_offn += nonfo_add;
    for(i=1;i<=nonfo_add    ;i++){
      mdonfo->j1[(i+ionfo_offn)] = 
        mdonfo->j1[(i+ionfo_off)]+ishift;
      mdonfo->j2[(i+ionfo_offn)] = 
        mdonfo->j2[(i+ionfo_off)]+ishift;
      mdonfo->jtyp[(i+ionfo_offn)] = 
        mdonfo->jtyp[(i+ionfo_off)];
    }/*endfor*/
    ionfo_null_offn += nonfo_null_add;
    for(i=1;i<=nonfo_null_add;i++){
      null_inter_parse->jonfo1_nul[(i+ionfo_null_offn)] = 
        null_inter_parse->jonfo1_nul[(i+ionfo_null_off)]+ishift;
      null_inter_parse->jonfo2_nul[(i+ionfo_null_offn)] = 
        null_inter_parse->jonfo2_nul[(i+ionfo_null_off)]+ishift;
    }/*endfor*/
  }/*endfor*/

  /*========================================================================*/
  /* VII) Replicate the bend_bnds                                           */

  nmol = atommaps->nmol_jmol_typ[jmol_typ];
  for(imol=1;imol<=(nmol-1);imol++){
    ishift = imol*natm_add;
    ibend_bnd_offn += nbend_bnd_add;
    for(i=1;i<=nbend_bnd_add;i++){
      mdbend_bnd->j1[(i+ibend_bnd_offn)] = 
        mdbend_bnd->j1[(i+ibend_bnd_off)]+ishift;
      mdbend_bnd->j2[(i+ibend_bnd_offn)] = 
        mdbend_bnd->j2[(i+ibend_bnd_off)]+ishift;
      mdbend_bnd->j3[(i+ibend_bnd_offn)] = 
        mdbend_bnd->j3[(i+ibend_bnd_off)]+ishift;
      mdbend_bnd->jtyp[(i+ibend_bnd_offn)] = 
        mdbend_bnd->jtyp[(i+ibend_bnd_off)];
    }/*endfor*/
  }/*endfor*/

  /*========================================================================*/
  /* VIII) Replicate the atoms                                              */

  nmol = atommaps->nmol_jmol_typ[jmol_typ];
  for(imol=1;imol<=(nmol-1);imol++){
    iatm_offn += natm_add;
    ishift = imol*natm_add;
    ig_shift = (nghost_add)*imol;
    if_shift = (nfreeze_add)*imol;
    for(i=1;i<=natm_add     ;i++){
      clatoms_info->mass[(i+iatm_offn)]   = clatoms_info->mass[(i+iatm_off)];
      clatoms_info->q[(i+iatm_offn)]      = clatoms_info->q[(i+iatm_off)];
      cpatom_maps->cp_vlnc_up[(i+iatm_offn)]= 
        cpatom_maps->cp_vlnc_up[(i+iatm_off)];
      cpatom_maps->cp_vlnc_dn[(i+iatm_offn)]= 
        cpatom_maps->cp_vlnc_dn[(i+iatm_off)];
      cpatom_maps->cp_vlnc_true_up[(i+iatm_offn)]= 
        cpatom_maps->cp_vlnc_true_up[(i+iatm_off)];
      cpatom_maps->cp_vlnc_true_dn[(i+iatm_offn)]= 
        cpatom_maps->cp_vlnc_true_dn[(i+iatm_off)];
      cpatom_maps->cp_atm_flag[(i+iatm_offn)]= 
        cpatom_maps->cp_atm_flag[(i+iatm_off)];
      clatoms_info->alp_pol[(i+iatm_offn)]= clatoms_info->alp_pol[(i+iatm_off)];
      clatoms_info->b_neut[(i+iatm_offn)] = clatoms_info->b_neut[(i+iatm_off)];
      clatoms_info->text_atm[(i+iatm_offn)] = 
        clatoms_info->text_atm[(i+iatm_off)];
      atommaps->iatm_atm_typ[(i+iatm_offn)] = 
        atommaps->iatm_atm_typ[(i+iatm_off)];
      atommaps->iatm_mol_typ[(i+iatm_offn)] = jmol_typ;
      atommaps->iatm_mol_num[(i+iatm_offn)] = imol+1;
      atommaps->iatm_res_num[(i+iatm_offn)] = 
        atommaps->iatm_res_num[(i+iatm_off)];
      atommaps->iatm_res_typ[(i+iatm_offn)] = 
        atommaps->iatm_res_typ[(i+iatm_off)];
      mdconstrnt->atom_label[(i+iatm_offn)] = 
        mdconstrnt->atom_label[(i+iatm_off)];
      mdconstrnt->freeze_flag[(i+iatm_offn)]  =  0;
      if(mdconstrnt->freeze_flag[(i+iatm_off)] != 0) {
        mdconstrnt->freeze_flag[(i+iatm_offn)]  =  
          mdconstrnt->freeze_flag[(i+iatm_off)]+if_shift;}
      mdghost_atoms->ighost_flag[(i+iatm_offn)]  =  0;
      if(mdghost_atoms->ighost_flag[(i+iatm_off)] != 0) {
        mdghost_atoms->ighost_flag[(i+iatm_offn)]  = 
          mdghost_atoms->ighost_flag[(i+iatm_off)]+ig_shift;}
    }/*endfor*/
  }/*endfor*/

  /*========================================================================*/
  /* IX) Replicate the ghosts                                              */

  nmol = atommaps->nmol_jmol_typ[jmol_typ];
  for(imol=1;imol<=(nmol-1);imol++){
    ighost_offn += nghost_add;
    ishift = imol*natm_add;
    for(i=1;i<=nghost_add     ;i++){
      mdghost_atoms->ighost_map[(i+ighost_offn)] = 
        mdghost_atoms->ighost_map[(i+ighost_off)] + ishift;
      mdghost_atoms->natm_comp[(i+ighost_offn)] = 
        mdghost_atoms->natm_comp[(i+ighost_off)];
      for(k=1;k<=mdghost_atoms->natm_comp[(i+ighost_off)];k++){
        mdghost_atoms->iatm_comp[k][(i+ighost_offn)] = 
          mdghost_atoms->iatm_comp[k][(i+ighost_off)]+ishift;
        mdghost_atoms->coef[k][(i+ighost_offn)] = 
          mdghost_atoms->coef[k][(i+ighost_off)];
      }/*endfor*/
    }/*endfor*/
  }/*endfor*/

  /*========================================================================*/
  /* IX) Replicate the frozen atoms                                         */

  nmol = atommaps->nmol_jmol_typ[jmol_typ];
  for(imol=1;imol<=(nmol-1);imol++){
    ifreeze_offn += nfreeze_add;
    ishift = imol*natm_add;
    for(i=1;i<=nfreeze_add     ;i++){
      mdconstrnt->freeze_map[(i+ifreeze_offn)] = 
        mdconstrnt->freeze_map[(i+ifreeze_off)] + ishift;
    }/*endfor*/
  }/*endfor*/

  /*========================================================================*/
  /* X) Replicate the grp_cons                                              */

  nmol = atommaps->nmol_jmol_typ[jmol_typ];
  for(imol=1;imol<=(nmol-1);imol++){
    ishift = imol*natm_add;
    igrp_21_offn += ngrp_21_add;
    for(i=1;i<=ngrp_21_add;i++){
      mdgrp_bond_con->j1_21[(i+igrp_21_offn)] = 
        mdgrp_bond_con->j1_21[(i+igrp_21_off)]+ishift;
      mdgrp_bond_con->j2_21[(i+igrp_21_offn)] = 
        mdgrp_bond_con->j2_21[(i+igrp_21_off)]+ishift;
      mdgrp_bond_con->jtyp_21[(i+igrp_21_offn)] = 
        mdgrp_bond_con->jtyp_21[(i+igrp_21_off)];
    }/*endfor*/
    igrp_23_offn += ngrp_23_add;
    for(i=1;i<=ngrp_23_add;i++){
      mdgrp_bond_con->j1_23[(i+igrp_23_offn)] = 
        mdgrp_bond_con->j1_23[(i+igrp_23_off)]+ishift;
      mdgrp_bond_con->j2_23[(i+igrp_23_offn)] = 
        mdgrp_bond_con->j2_23[(i+igrp_23_off)]+ishift;
      mdgrp_bond_con->j3_23[(i+igrp_23_offn)] = 
        mdgrp_bond_con->j3_23[(i+igrp_23_off)]+ishift;
      mdgrp_bond_con->jtyp_23[(i+igrp_23_offn)] = 
        mdgrp_bond_con->jtyp_23[(i+igrp_23_off)];
    }/*endfor*/
    igrp_33_offn += ngrp_33_add;
    for(i=1;i<=ngrp_33_add;i++){
      mdgrp_bond_con->j1_33[(i+igrp_33_offn)] = 
        mdgrp_bond_con->j1_33[(i+igrp_33_off)]+ishift;
      mdgrp_bond_con->j2_33[(i+igrp_33_offn)] = 
        mdgrp_bond_con->j2_33[(i+igrp_33_off)]+ishift;
      mdgrp_bond_con->j3_33[(i+igrp_33_offn)] = 
        mdgrp_bond_con->j3_33[(i+igrp_33_off)]+ishift;
      mdgrp_bond_con->jtyp_33[(i+igrp_33_offn)] = 
        mdgrp_bond_con->jtyp_33[(i+igrp_33_off)];
    }/*endfor*/
    igrp_watt_33_offn += ngrp_watt_33_add;
    for(i=1;i<=ngrp_watt_33_add;i++){
      mdgrp_bond_watts->j1_33[(i+igrp_watt_33_offn)] = 
        mdgrp_bond_watts->j1_33[(i+igrp_watt_33_off)]+ishift;
      mdgrp_bond_watts->j2_33[(i+igrp_watt_33_offn)] = 
        mdgrp_bond_watts->j2_33[(i+igrp_watt_33_off)]+ishift;
      mdgrp_bond_watts->j3_33[(i+igrp_watt_33_offn)] = 
        mdgrp_bond_watts->j3_33[(i+igrp_watt_33_off)]+ishift;
      mdgrp_bond_watts->jtyp_33[(i+igrp_watt_33_offn)] = 
        mdgrp_bond_watts->jtyp_33[(i+igrp_watt_33_off)];
    }/*endfor*/
    igrp_43_offn += ngrp_43_add;
    for(i=1;i<=ngrp_43_add;i++){
      mdgrp_bond_con->j1_43[(i+igrp_43_offn)] = 
        mdgrp_bond_con->j1_43[(i+igrp_43_off)]+ishift;
      mdgrp_bond_con->j2_43[(i+igrp_43_offn)] = 
        mdgrp_bond_con->j2_43[(i+igrp_43_off)]+ishift;
      mdgrp_bond_con->j3_43[(i+igrp_43_offn)] = 
        mdgrp_bond_con->j3_43[(i+igrp_43_off)]+ishift;
      mdgrp_bond_con->j4_43[(i+igrp_43_offn)] = 
        mdgrp_bond_con->j4_43[(i+igrp_43_off)]+ishift;
      mdgrp_bond_con->jtyp_43[(i+igrp_43_offn)] = 
        mdgrp_bond_con->jtyp_43[(i+igrp_43_off)];
    }/*endfor*/
    igrp_46_offn += ngrp_46_add;
    for(i=1;i<=ngrp_46_add;i++){
      mdgrp_bond_con->j1_46[(i+igrp_46_offn)] = 
        mdgrp_bond_con->j1_46[(i+igrp_46_off)]+ishift;
      mdgrp_bond_con->j2_46[(i+igrp_46_offn)] = 
        mdgrp_bond_con->j2_46[(i+igrp_46_off)]+ishift;
      mdgrp_bond_con->j3_46[(i+igrp_46_offn)] = 
        mdgrp_bond_con->j3_46[(i+igrp_46_off)]+ishift;
      mdgrp_bond_con->j4_46[(i+igrp_46_offn)] = 
        mdgrp_bond_con->j4_46[(i+igrp_46_off)]+ishift;
      mdgrp_bond_con->jtyp_46[(i+igrp_46_offn)] = 
        mdgrp_bond_con->jtyp_46[(i+igrp_46_off)];
    }/*endfor*/
  }/*endfor*/

  /*========================================================================*/
  /* XI) Increase the list counters                                         */

  nmol1 = nmol-1;


  mdbond->npow += nbond_pow_add*nmol1;
  mdbond->ncon += nbond_con_add*nmol1;
  null_inter_parse->nbond_nul += nbond_null_add*nmol1;

  mdbend->npow += nbend_pow_add*nmol1;
  null_inter_parse->nbend_nul += nbend_null_add*nmol1;

  mdbend_bnd->num += nbend_bnd_add*nmol1;

  mdtors->npow += ntors_pow_add*nmol1;
  null_inter_parse->ntors_nul += ntors_null_add*nmol1;

  mdonfo->num += nonfo_add*nmol1;
  null_inter_parse->nonfo_nul += nonfo_null_add*nmol1;

  clatoms_info->natm_tot += natm_add*nmol1;
  mdghost_atoms->nghost_tot += nghost_add*nmol1;
  mdconstrnt->nfreeze += nfreeze_add*nmol1;

  mdgrp_bond_con->num_21 += ngrp_21_add*nmol1;
  mdgrp_bond_con->num_23 += ngrp_23_add*nmol1;
  mdgrp_bond_con->num_33 += ngrp_33_add*nmol1;
  mdgrp_bond_watts->num_33 += ngrp_watt_33_add*nmol1;
  mdgrp_bond_con->num_43 += ngrp_43_add*nmol1;   
  mdgrp_bond_con->num_46 += ngrp_46_add*nmol1;   


  /*==========================================================================*/
}/*end routine*/ 
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void reallocate_intra_list(MDCLATOMS_INFO *clatoms_info,
    CPATOM_MAPS *cpatom_maps,MDATOM_MAPS *atommaps,
    BUILD_INTRA *build_intra,MDINTRA *mdintra,
    NULL_INTER_PARSE *null_inter_parse,
    int jmol_typ,
    int nbond_pow_add,int nbond_con_add,int nbond_null_add,
    int nbend_pow_add,int nbend_null_add,
    int ntors_pow_add,int ntors_null_add,
    int nonfo_add,int nonfo_null_add,int nbend_bnd_add,
    int natm_add,int nghost_add,int nfreeze_add,int ngrp_43_add,
    int ngrp_33_add,int ngrp_watt_33_add,int ngrp_21_add,
    int ngrp_23_add,int ngrp_46_add)

  /*==========================================================================*/
{ /*begin routine */
  /*==========================================================================*/

#include "../class_defs/allclass_strip_mdintra.h"   

  int nmol,nmol1,mem_add_now,iii;
  int nfreeze_old,nfreeze_new;
  int nghost_new,nghost_old,ncomp_new,ncomp_old;

  /*==========================================================================*/
  /* I) Reallocate the atoms                                                  */

  nmol = atommaps->nmol_jmol_typ[jmol_typ];
  nmol1 = (nmol-1)*natm_add;
  if(clatoms_info->natm_tot+nmol1 > build_intra->natm_tot_max){
    mem_add_now = clatoms_info->natm_tot+nmol1 - build_intra->natm_tot_max;
    build_intra->natm_tot_max += MAX(mem_add_now,NMEM_MIN);

    clatoms_info->mass = (double *)crealloc(&((clatoms_info->mass[1])),
        build_intra->natm_tot_max*sizeof(double),"realloc_inta_list")-1;
    clatoms_info->q        = (double *)crealloc(&((clatoms_info->q[1])),
        build_intra->natm_tot_max*sizeof(double),"realloc_inta_list")-1;
    cpatom_maps->cp_vlnc_up  = (int *)crealloc(&((cpatom_maps->cp_vlnc_up[1])),
        build_intra->natm_tot_max*sizeof(int),"realloc_intra_list")-1;
    cpatom_maps->cp_vlnc_dn  = (int *)crealloc(&((cpatom_maps->cp_vlnc_dn[1])),
        build_intra->natm_tot_max*sizeof(int),"realloc_intra_list")-1;
    cpatom_maps->cp_vlnc_true_up  = (int *)crealloc(&((cpatom_maps->cp_vlnc_true_up[1])),
        build_intra->natm_tot_max*sizeof(int),"realloc_intra_list")-1;
    cpatom_maps->cp_vlnc_true_dn  = (int *)crealloc(&((cpatom_maps->cp_vlnc_true_dn[1])),
        build_intra->natm_tot_max*sizeof(int),"realloc_intra_list")-1;
    cpatom_maps->cp_atm_flag =(int *)crealloc(&((cpatom_maps->cp_atm_flag[1])),
        build_intra->natm_tot_max*sizeof(int),"realloc_intra_list")-1;
    clatoms_info->alp_pol  = (double *)crealloc(&((clatoms_info->alp_pol[1])),
        build_intra->natm_tot_max*sizeof(double),"realloc_inta_list")-1;
    clatoms_info->b_neut   = (double *)crealloc(&((clatoms_info->b_neut[1])),
        build_intra->natm_tot_max*sizeof(double),"realloc_inta_list")-1;
    clatoms_info->text_atm =(double *)crealloc(&((clatoms_info->text_atm[1])),
        build_intra->natm_tot_max*sizeof(double),"realloc_inta_list")-1;
    atommaps->iatm_mol_typ  = (int *) crealloc(&(atommaps->iatm_mol_typ[1]),
        build_intra->natm_tot_max*sizeof(int),"realloc_intra_list")-1;
    atommaps->iatm_atm_typ  = (int *) crealloc(&(atommaps->iatm_atm_typ[1]),
        build_intra->natm_tot_max*sizeof(int),"realloc_intra_list")-1;
    atommaps->iatm_res_typ  = (int *) crealloc(&(atommaps->iatm_res_typ[1]),
        build_intra->natm_tot_max*sizeof(int),"realloc_intra_list")-1;
    atommaps->iatm_mol_num  = (int *) crealloc(&(atommaps->iatm_mol_num[1]),
        build_intra->natm_tot_max*sizeof(int),"realloc_intra_list")-1;
    atommaps->iatm_res_num  = (int *) crealloc(&(atommaps->iatm_res_num[1]),
        build_intra->natm_tot_max*sizeof(int),"realloc_intra_list")-1;
    mdghost_atoms->ighost_flag   = (int *)crealloc(&(mdghost_atoms->ighost_flag[1]),
        build_intra->natm_tot_max*sizeof(int),"realloc_intra_list")-1;
    mdconstrnt->freeze_flag   = (int *)crealloc(&(mdconstrnt->freeze_flag[1]),
        build_intra->natm_tot_max*sizeof(int),"realloc_intra_list")-1;
    mdconstrnt->atom_label   = (int *)crealloc(&(mdconstrnt->atom_label[1]),
        build_intra->natm_tot_max*sizeof(int),"realloc_intra_list")-1;
  }/*endif*/

  /*==========================================================================*/
  /* I) Reallocate the ghosts                                                 */

  nmol = atommaps->nmol_jmol_typ[jmol_typ];
  nmol1 = (nmol-1)*nghost_add;
  if(mdghost_atoms->nghost_tot+nmol1 > build_intra->nghost_tot_max){
    nghost_old  = build_intra->nghost_tot_max;
    ncomp_old   = NCOEF_GHOST_MAX;
    build_intra->nghost_tot_max+=
      MAX((mdghost_atoms->nghost_tot+nmol1-build_intra->nghost_tot_max),NMEM_MIN);
    nghost_new  = build_intra->nghost_tot_max;
    ncomp_new   = NCOEF_GHOST_MAX;
    mdghost_atoms->ighost_map = (int *) crealloc(
        &(mdghost_atoms->ighost_map)[1],nghost_new*sizeof(int),
        "realloc_intra_list")-1;
    mdghost_atoms->natm_comp  = (int *) crealloc(  
        &(mdghost_atoms->natm_comp)[1],nghost_new*sizeof(int),
        "realloc_intra_list")-1;
    mdghost_atoms->iatm_comp  = creall_int_mat(mdghost_atoms->iatm_comp,
        1,ncomp_old,1,nghost_old,
        1,ncomp_new,1,nghost_new,
        "reallocate_intra_list");
    mdghost_atoms->coef       = creall_mat(mdghost_atoms->coef,
        1,ncomp_old,1,nghost_old,
        1,ncomp_new,1,nghost_new,
        "reallocate_intra_list");
  }/*endif*/

  /*==========================================================================*/
  /* I) Reallocate the freeze                                                 */

  nmol = atommaps->nmol_jmol_typ[jmol_typ];
  nmol1 = (nmol-1)*nfreeze_add;
  if(mdconstrnt->nfreeze+nmol1 > build_intra->nfreeze_max){
    nfreeze_old  = build_intra->nfreeze_max;
    build_intra->nfreeze_max+=
      MAX((mdconstrnt->nfreeze+nmol1-build_intra->nfreeze_max),NMEM_MIN);
    nfreeze_new  = build_intra->nfreeze_max;
    mdconstrnt->freeze_map = (int *) crealloc(
        &(mdconstrnt->freeze_map)[1],
        nfreeze_new*sizeof(int),"realloc_intra_list")-1;
  }/*endif*/

  /*==========================================================================*/
  /* I) Reallocate the grpcons                                                */


  nmol = atommaps->nmol_jmol_typ[jmol_typ];
  nmol1 = (nmol-1)*ngrp_21_add;

  if(mdgrp_bond_con->num_21+nmol1 > build_intra->ngrp_21_max){
    build_intra->ngrp_21_max += 
      MAX((mdgrp_bond_con->num_21+nmol1-build_intra->ngrp_21_max),NMEM_MIN);
    mdgrp_bond_con->j1_21     = 
      (int *) crealloc(&(mdgrp_bond_con->j1_21)[1],
          build_intra->ngrp_21_max*sizeof(int),"realloc_intra_list")-1;
    mdgrp_bond_con->j2_21     = 
      (int *) crealloc(&(mdgrp_bond_con->j2_21)[1],
          build_intra->ngrp_21_max*sizeof(int),"realloc_intra_list")-1;
    mdgrp_bond_con->jtyp_21   = 
      (int *) crealloc(&(mdgrp_bond_con->jtyp_21)[1],
          build_intra->ngrp_21_max*sizeof(int),"realloc_intra_list")-1;
  }/*endif*/
  nmol = atommaps->nmol_jmol_typ[jmol_typ];
  nmol1 = (nmol-1)*ngrp_23_add;
  if(mdgrp_bond_con->num_23+nmol1 > build_intra->ngrp_23_max){
    build_intra->ngrp_23_max += 
      MAX((mdgrp_bond_con->num_23+nmol1-build_intra->ngrp_23_max),NMEM_MIN);
    mdgrp_bond_con->j1_23     = 
      (int *) crealloc(&(mdgrp_bond_con->j1_23)[1],
          build_intra->ngrp_23_max*sizeof(int),"realloc_intra_list")-1;
    mdgrp_bond_con->j2_23     = 
      (int *) crealloc(&(mdgrp_bond_con->j2_23)[1],
          build_intra->ngrp_23_max*sizeof(int),"realloc_intra_list")-1;
    mdgrp_bond_con->j3_23     = 
      (int *) crealloc(&(mdgrp_bond_con->j3_23)[1],
          build_intra->ngrp_23_max*sizeof(int),"realloc_intra_list")-1;
    mdgrp_bond_con->jtyp_23   = 
      (int *) crealloc(&(mdgrp_bond_con->jtyp_23)[1],
          build_intra->ngrp_23_max*sizeof(int),"realloc_intra_list")-1;
  }/*endif*/
  nmol = atommaps->nmol_jmol_typ[jmol_typ];
  nmol1 = (nmol-1)*ngrp_33_add;
  if(mdgrp_bond_con->num_33+nmol1 > build_intra->ngrp_33_max){
    build_intra->ngrp_33_max += 
      MAX((mdgrp_bond_con->num_33+nmol1-build_intra->ngrp_33_max),NMEM_MIN);
    mdgrp_bond_con->j1_33     = 
      (int *) crealloc(&(mdgrp_bond_con->j1_33)[1],
          build_intra->ngrp_33_max*sizeof(int),"realloc_intra_list")-1;
    mdgrp_bond_con->j2_33     = 
      (int *) crealloc(&(mdgrp_bond_con->j2_33)[1],
          build_intra->ngrp_33_max*sizeof(int),"realloc_intra_list")-1;
    mdgrp_bond_con->j3_33     = 
      (int *) crealloc(&(mdgrp_bond_con->j3_33)[1],
          build_intra->ngrp_33_max*sizeof(int),"realloc_intra_list")-1;
    mdgrp_bond_con->jtyp_33   = 
      (int *) crealloc(&(mdgrp_bond_con->jtyp_33)[1],
          build_intra->ngrp_33_max*sizeof(int),"realloc_intra_list")-1;
  }/*endif*/
  nmol = atommaps->nmol_jmol_typ[jmol_typ];
  nmol1 = (nmol-1)*ngrp_watt_33_add;
  if(mdgrp_bond_watts->num_33+nmol1 > build_intra->ngrp_watt_33_max){
    build_intra->ngrp_watt_33_max += 
      MAX((mdgrp_bond_watts->num_33+nmol1-build_intra->ngrp_watt_33_max),
          NMEM_MIN);
    mdgrp_bond_watts->j1_33     = 
      (int *) crealloc(&(mdgrp_bond_watts->j1_33)[1],
          build_intra->ngrp_watt_33_max*sizeof(int),"realloc_intra_list")-1;
    mdgrp_bond_watts->j2_33     = 
      (int *) crealloc(&(mdgrp_bond_watts->j2_33)[1],
          build_intra->ngrp_watt_33_max*sizeof(int),"realloc_intra_list")-1;
    mdgrp_bond_watts->j3_33     = 
      (int *) crealloc(&(mdgrp_bond_watts->j3_33)[1],
          build_intra->ngrp_watt_33_max*sizeof(int),"realloc_intra_list")-1;
    mdgrp_bond_watts->jtyp_33   = 
      (int *) crealloc(&(mdgrp_bond_watts->jtyp_33)[1],
          build_intra->ngrp_watt_33_max*sizeof(int),"realloc_intra_list")-1;
  }/*endif*/
  nmol = atommaps->nmol_jmol_typ[jmol_typ];
  nmol1 = (nmol-1)*ngrp_43_add;
  if(mdgrp_bond_con->num_43+nmol1 > build_intra->ngrp_43_max){
    build_intra->ngrp_43_max += 
      MAX((mdgrp_bond_con->num_43+nmol1-build_intra->ngrp_43_max),NMEM_MIN);
    mdgrp_bond_con->j1_43     = 
      (int *) crealloc(&(mdgrp_bond_con->j1_43)[1],
          build_intra->ngrp_43_max*sizeof(int),"realloc_intra_list")-1;
    mdgrp_bond_con->j2_43     = 
      (int *) crealloc(&(mdgrp_bond_con->j2_43)[1],
          build_intra->ngrp_43_max*sizeof(int),"realloc_intra_list")-1;
    mdgrp_bond_con->j3_43     = 
      (int *) crealloc(&(mdgrp_bond_con->j3_43)[1],
          build_intra->ngrp_43_max*sizeof(int),"realloc_intra_list")-1;
    mdgrp_bond_con->j4_43     =  
      (int *) crealloc(&(mdgrp_bond_con->j4_43)[1],
          build_intra->ngrp_43_max*sizeof(int),"realloc_intra_list")-1;
    mdgrp_bond_con->jtyp_43   = 
      (int *) crealloc(&(mdgrp_bond_con->jtyp_43)[1],
          build_intra->ngrp_43_max*sizeof(int),"realloc_intra_list")-1;
  }/*endif*/
  nmol = atommaps->nmol_jmol_typ[jmol_typ];
  nmol1 = (nmol-1)*ngrp_46_add;
  if(mdgrp_bond_con->num_46+nmol1 > build_intra->ngrp_46_max){
    build_intra->ngrp_46_max += 
      MAX((mdgrp_bond_con->num_46+nmol1-build_intra->ngrp_46_max),NMEM_MIN);
    mdgrp_bond_con->j1_46     = 
      (int *) crealloc(&(mdgrp_bond_con->j1_46)[1],
          build_intra->ngrp_46_max*sizeof(int),"realloc_intra_list")-1;
    mdgrp_bond_con->j2_46     = 
      (int *) crealloc(&(mdgrp_bond_con->j2_46)[1],
          build_intra->ngrp_46_max*sizeof(int),"realloc_intra_list")-1;
    mdgrp_bond_con->j3_46     = 
      (int *) crealloc(&(mdgrp_bond_con->j3_46)[1],
          build_intra->ngrp_46_max*sizeof(int),"realloc_intra_list")-1;
    mdgrp_bond_con->j4_46     =  
      (int *) crealloc(&(mdgrp_bond_con->j4_46)[1],
          build_intra->ngrp_46_max*sizeof(int),"realloc_intra_list")-1;
    mdgrp_bond_con->jtyp_46   = 
      (int *) crealloc(&(mdgrp_bond_con->jtyp_46)[1],
          build_intra->ngrp_46_max*sizeof(int),"realloc_intra_list")-1;
  }/*endif*/

  /*===========================================================================*/
  /* I) Reallocate the bonds */

  nmol = atommaps->nmol_jmol_typ[jmol_typ];
  nmol1 = (nmol-1)*nbond_pow_add;
  if(mdbond->npow+nmol1 > build_intra->nbond_pow_max){
    build_intra->nbond_pow_max += 
      MAX((mdbond->npow+nmol1-build_intra->nbond_pow_max),NMEM_MIN);
    mdbond->j1_pow = (int *)crealloc(&(mdbond->j1_pow)[1],
        build_intra->nbond_pow_max*sizeof(int),"realloc_intra_list")-1;
    mdbond->j2_pow = (int *)crealloc(&(mdbond->j2_pow)[1],
        build_intra->nbond_pow_max*sizeof(int),"realloc_intra_list")-1;
    mdbond->jtyp_pow=(int *)crealloc(&(mdbond->jtyp_pow)[1],
        build_intra->nbond_pow_max*sizeof(int),"realloc_intra_list")-1;
  }/*endif*/

  nmol1 = (nmol-1)*nbond_con_add;
  if(mdbond->ncon+nmol1 > build_intra->nbond_con_max){
    build_intra->nbond_con_max += 
      MAX((mdbond->ncon+nmol1-build_intra->nbond_con_max),NMEM_MIN);
    mdbond->j1_con = (int *)crealloc(&(mdbond->j1_con)[1],
        build_intra->nbond_con_max*sizeof(int),"realloc_intra_list")-1;
    mdbond->j2_con = (int *)crealloc(&(mdbond->j2_con)[1],
        build_intra->nbond_con_max*sizeof(int),"realloc_intra_list")-1;
    mdbond->jtyp_con=(int *)crealloc(&(mdbond->jtyp_con)[1],
        build_intra->nbond_con_max*sizeof(int),"realloc_intra_list")-1;
  }/*endif*/

  nmol1 = (nmol-1)*nbond_null_add;
  if((null_inter_parse->nbond_nul)+nmol1 > build_intra->nbond_nul_max) {
    build_intra->nbond_nul_max   += 
      MAX((null_inter_parse->nbond_nul+nmol1-build_intra->nbond_nul_max),NMEM_MIN);
    null_inter_parse->jbond1_nul  = 
      (int *) crealloc(&(null_inter_parse->jbond1_nul)[1],
          build_intra->nbond_nul_max*sizeof(int),"realloc_intra_list")-1;
    null_inter_parse->jbond2_nul  = 
      (int *) crealloc(&(null_inter_parse->jbond2_nul)[1],
          build_intra->nbond_nul_max*sizeof(int),"realloc_intra_list")-1;
  }  /*endif*/

  /*====================================================================*/
  /* II) Reallocate the bends                                           */

  nmol1 = (nmol-1)*nbend_pow_add;
  if(mdbend->npow+nmol1 > build_intra->nbend_pow_max){
    build_intra->nbend_pow_max += 
      MAX((mdbend->npow+nmol1-build_intra->nbend_pow_max),NMEM_MIN);
    mdbend->j1_pow =(int*)crealloc(&(mdbend->j1_pow)[1], 
        build_intra->nbend_pow_max*sizeof(int),"realloc_intra_list")-1;
    mdbend->j2_pow =(int*)crealloc(&(mdbend->j2_pow)[1],
        build_intra->nbend_pow_max*sizeof(int),"realloc_intra_list")-1;
    mdbend->j3_pow =(int*)crealloc(&(mdbend->j3_pow)[1],
        build_intra->nbend_pow_max*sizeof(int),"realloc_intra_list")-1;
    mdbend->jtyp_pow=(int*)crealloc(&(mdbend->jtyp_pow)[1],
        build_intra->nbend_pow_max*sizeof(int),"realloc_intra_list")-1;
  }/*endif*/

  nmol1 = (nmol-1)*nbend_null_add;
  if(null_inter_parse->nbend_nul+nmol1 > build_intra->nbend_nul_max){
    build_intra->nbend_nul_max += 
      MAX((null_inter_parse->nbend_nul+nmol1-build_intra->nbend_nul_max),NMEM_MIN);
    null_inter_parse->jbend1_nul     = 
      (int *) crealloc(&(null_inter_parse->jbend1_nul)[1],
          build_intra->nbend_nul_max*sizeof(int),"realloc_intra_list")-1;
    null_inter_parse->jbend2_nul     = 
      (int *) crealloc(&(null_inter_parse->jbend2_nul)[1],
          build_intra->nbend_nul_max*sizeof(int),"realloc_intra_list")-1;
    null_inter_parse->jbend3_nul     = 
      (int *) crealloc(&(null_inter_parse->jbend3_nul)[1],
          build_intra->nbend_nul_max*sizeof(int),"realloc_intra_list")-1;
  }     /*endif*/

  /*====================================================================*/
  /* IV) Reallocate the torsions                                        */

  nmol1 = (nmol-1)*ntors_pow_add;
  if((mdtors->npow)+nmol1 > build_intra->ntors_pow_max){
    build_intra->ntors_pow_max += 
      MAX((mdtors->npow+nmol1-build_intra->ntors_pow_max),NMEM_MIN);
    mdtors->j1_pow =(int*)crealloc(&(mdtors->j1_pow)[1],
        build_intra->ntors_pow_max*sizeof(int),"realloc_intra_list")-1;
    mdtors->j2_pow =(int*)crealloc(&(mdtors->j2_pow)[1],
        build_intra->ntors_pow_max*sizeof(int),"realloc_intra_list")-1;
    mdtors->j3_pow =(int*)crealloc(&(mdtors->j3_pow)[1],
        build_intra->ntors_pow_max*sizeof(int),"realloc_intra_list")-1;
    mdtors->j4_pow =(int*)crealloc(&(mdtors->j4_pow)[1],
        build_intra->ntors_pow_max*sizeof(int),"realloc_intra_list")-1;
    mdtors->jtyp_pow=(int*)crealloc(&(mdtors->jtyp_pow)[1],
        build_intra->ntors_pow_max*sizeof(int),"realloc_intra_list")-1;
  }/*endif*/

  nmol1 = (nmol-1)*ntors_null_add;
  if(null_inter_parse->ntors_nul+nmol1 > build_intra->ntors_nul_max){
    build_intra->ntors_nul_max += 
      MAX((null_inter_parse->ntors_nul+nmol1-build_intra->ntors_nul_max),NMEM_MIN);
    null_inter_parse->jtors1_nul = 
      (int *) crealloc(&(null_inter_parse->jtors1_nul[1]),
          build_intra->ntors_nul_max*sizeof(int),"realloc_intra_list")-1;
    null_inter_parse->jtors2_nul     = 
      (int *) crealloc(&(null_inter_parse->jtors2_nul[1]),
          build_intra->ntors_nul_max*sizeof(int),"realloc_intra_list")-1;
    null_inter_parse->jtors3_nul     = 
      (int *) crealloc(&(null_inter_parse->jtors3_nul[1]),
          build_intra->ntors_nul_max*sizeof(int),"realloc_intra_list")-1;
    null_inter_parse->jtors4_nul     = 
      (int *) crealloc(&(null_inter_parse->jtors4_nul[1]),
          build_intra->ntors_nul_max*sizeof(int),"realloc_intra_list")-1;
  }      /*endif*/

  /*====================================================================*/
  /* V) Reallocate the onfos                                         */

  nmol1 = (nmol-1)*nonfo_add;
  if(mdonfo->num+nmol1 > build_intra->nonfo_max){
    build_intra->nonfo_max += 
      MAX((mdonfo->num+nmol1-build_intra->nonfo_max),NMEM_MIN);
    mdonfo->j1 =(int *) crealloc(&(mdonfo->j1)[1],
        build_intra->nonfo_max*sizeof(int),"realloc_intra_list")-1;
    mdonfo->j2 =(int *) crealloc(&(mdonfo->j2)[1],
        build_intra->nonfo_max*sizeof(int),"realloc_intra_list")-1;
    mdonfo->jtyp =(int *) crealloc(&(mdonfo->jtyp)[1],
        build_intra->nonfo_max*sizeof(int),"realloc_intra_list")-1;
  }/*endif*/

  nmol1 = (nmol-1)*nonfo_null_add;
  if(null_inter_parse->nonfo_nul+nmol1 > build_intra->nonfo_nul_max){
    build_intra->nonfo_nul_max += 
      MAX((null_inter_parse->nonfo_nul+nmol1-build_intra->nonfo_nul_max),NMEM_MIN);
    null_inter_parse->jonfo1_nul = 
      (int *) crealloc(&(null_inter_parse->jonfo1_nul)[1],
          build_intra->nonfo_nul_max*sizeof(int),"realloc_intra_list")-1;
    null_inter_parse->jonfo2_nul     = 
      (int *) crealloc(&(null_inter_parse->jonfo2_nul)[1],
          build_intra->nonfo_nul_max*sizeof(int),"realloc_intra_list")-1;
  }/*endif*/

  /*====================================================================*/
  /* VI) Reallocate the Uri-Bradleys                                    */

  nmol1 = (nmol-1)*nbend_bnd_add;
  if(mdbend_bnd->num+nmol1 > build_intra->nbend_bnd_max){
    mem_add_now = mdbend_bnd->num+nmol1-build_intra->nbend_bnd_max;
    build_intra->nbend_bnd_max += MAX(mem_add_now,NMEM_MIN);
    mdbend_bnd->j1 =(int*)crealloc(&(mdbend_bnd->j1)[1], 
        build_intra->nbend_bnd_max*sizeof(int),"realloc_intra_list")-1;
    mdbend_bnd->j2 =(int*)crealloc(&(mdbend_bnd->j2)[1],
        build_intra->nbend_bnd_max*sizeof(int),"realloc_intra_list")-1;
    mdbend_bnd->j3 =(int*)crealloc(&(mdbend_bnd->j3)[1],
        build_intra->nbend_bnd_max*sizeof(int),"realloc_intra_list")-1;
    mdbend_bnd->jtyp=(int*)crealloc(&(mdbend_bnd->jtyp)[1],
        build_intra->nbend_bnd_max*sizeof(int),"realloc_intra_list")-1;
  }/*endif*/

  /*==========================================================================*/
} /*end routine*/
/*==========================================================================*/





