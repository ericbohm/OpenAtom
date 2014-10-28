/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: fetch_freeze.c                               */
/*                                                                          */
/* This subprogram sets up frozen atoms                                     */
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

#define DEBUG_OFF

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void fetch_freeze(CLASS_PARSE *class_parse, MDATOM_MAPS *atommaps, 
    MDCONSTRNT *constrnt,BUILD_INTRA *build_intra, 
    START_INDEX *start_index, MDCLATOMS_INFO *clatoms_info, 
    int jmol_typ)

  /*==========================================================================*/
  /*begin routine*/{

    /*==========================================================================*/
    /* Define Local variable */

    int i,iii,nfreeze_now,iresidue_off,nresidue,count;
    int count_atm,ires;
    double *mass = clatoms_info->mass;


    /*==========================================================================*/
    /* I) No atoms frozen */ 

    if(class_parse->mol_freeze_opt[jmol_typ]==0){
      for(i=start_index->natm+1;i<=clatoms_info->natm_tot;i++){
        constrnt->freeze_flag[i] = 0;
      }/*endfor*/
    }/*endif*/

    /*==========================================================================*/
    /* II) All atoms frozen */ 
    if(class_parse->mol_freeze_opt[jmol_typ]==1){
      for(i=start_index->natm+1;i<=clatoms_info->natm_tot;i++){
        constrnt->freeze_flag[i] = constrnt->nfreeze+i-start_index->natm;
      }/*endfor*/
      nfreeze_now = clatoms_info->natm_tot - start_index->natm;
      constrnt->nfreeze += nfreeze_now;

      if((constrnt->nfreeze) > build_intra->nfreeze_max){
        build_intra->nfreeze_max += MAX(NMEM_MIN,nfreeze_now);
        constrnt->freeze_map  = (int *) crealloc(&(constrnt->freeze_map[1]),
            build_intra->nfreeze_max*
            sizeof(int),"fetch_freeze")-1;
      }/*endif*/
      for(i=start_index->nfreeze+1;i<=constrnt->nfreeze;i++){
        constrnt->freeze_map[i] = start_index->natm+i-start_index->nfreeze;
      }/*endfor*/

      /* Do the degrees of freedom for each molecule             */
      iresidue_off = atommaps->jres_jmol_typ_strt[jmol_typ]-1;
      atommaps->nfree_1mol_jmol_typ[jmol_typ] = 0; 
      constrnt->icons_jmol_typ[jmol_typ] = 0;
      nresidue = MAX(atommaps->nres_1mol_jmol_typ[jmol_typ],1);
      for(i=1;i<=nresidue;i++){
        atommaps->nfree_jres_jmol_typ[(iresidue_off+i)] = 0;
        /*endfor*/}

    }/*endif*/

    /*==========================================================================*/
    /* III) Backbone atoms frozen */ 

    if(class_parse->mol_freeze_opt[jmol_typ]==2){

      /* Count the freeze map and fill the freeze_flag list */
      nfreeze_now = constrnt->nfreeze;
      for(i=start_index->natm+1;i<=clatoms_info->natm_tot;i++){
        if(constrnt->atom_label[i] == 1){
          nfreeze_now++;
          constrnt->freeze_flag[i] = i;
        }/*endif*/
      }/*endfor*/

      /* Realloc the freeze map  */
      if(nfreeze_now > build_intra->nfreeze_max){
        build_intra->nfreeze_max += MAX(NMEM_MIN,nfreeze_now-constrnt->nfreeze);
        constrnt->freeze_map  = (int *) crealloc(&(constrnt->freeze_map[1]),
            build_intra->nfreeze_max*
            sizeof(int),"fetch_freeze")-1;
      }/*endif*/

      /* Fill the freeze map and do the degrees of freedom for each molecule */
      count = 0;
      count_atm = 0;
      iresidue_off = atommaps->jres_jmol_typ_strt[jmol_typ]-1;
      ires = iresidue_off+1;
      for(i=start_index->natm+1;i<=clatoms_info->natm_tot;i++){
        count_atm++;
        if(count_atm > atommaps->natm_jres_jmol_typ[ires]){
          count_atm = 1;
          ires++;
        }/*endif*/
        if(constrnt->atom_label[i] == 1){
          count++;
          constrnt->freeze_map[(constrnt->nfreeze+count)] = i;
          atommaps->nfree_1mol_jmol_typ[jmol_typ] -= 3; 
          atommaps->nfree_jres_jmol_typ[ires] -= 3;
        }/*endif*/
      }/*endfor*/
      constrnt->nfreeze = nfreeze_now;
      constrnt->nfreeze = nfreeze_now;

#ifdef DEBUG
      for(i=1;i<=constrnt->nfreeze;i++){
        PRINTF("This is present freeze map %d\n",constrnt->freeze_map[i]);
        scanf("%d",&iii); 
      }

      PRINTF("This is molecule degrees of freedom %d\n",
          atommaps->nfree_1mol_jmol_typ[jmol_typ]);
      scanf("%d",&iii); 

      nresidue = MAX(atommaps->nres_1mol_jmol_typ[jmol_typ],1);
      for(i=1;i<=nresidue;i++){
        PRINTF("This is residue degrees of freedom %d %d\n",
            atommaps->nfree_jres_jmol_typ[i+iresidue_off],i);
        scanf("%d",&iii); 
      }
#endif

    }/*endif*/

    /*==========================================================================*/
    /* IV) All Heavy atoms frozen */ 

    if(class_parse->mol_freeze_opt[jmol_typ]==3){

      /* Count the freeze map and fill the freeze_flag list */
      nfreeze_now = constrnt->nfreeze;

      for(i=start_index->natm+1;i<=clatoms_info->natm_tot;i++){
        if(mass[i] > 1.5){
          nfreeze_now++;
          constrnt->freeze_flag[i] = i;
        }/*endif*/
      }/*endfor*/

      /* Realloc the freeze map  */

      if(nfreeze_now > build_intra->nfreeze_max){
        PRINTF(" reallocing freeze_map \n");
        build_intra->nfreeze_max += MAX(NMEM_MIN,nfreeze_now-constrnt->nfreeze);
        constrnt->freeze_map  = (int *) crealloc(&(constrnt->freeze_map[1]),
            build_intra->nfreeze_max*
            sizeof(int),"fetch_freeze")-1;
      }/*endif*/

      /* Fill the freeze map and do the degrees of freedom for each molecule */
      count = 0;
      count_atm = 0;
      iresidue_off = atommaps->jres_jmol_typ_strt[jmol_typ]-1;
      ires = iresidue_off+1;

      for(i=start_index->natm+1;i<=clatoms_info->natm_tot;i++){
        count_atm++;
        if(count_atm > atommaps->natm_jres_jmol_typ[ires]){
          count_atm = 1;
          ires++;
        }/*endif*/
        if(mass[i] > 1.5){
          count++;
          constrnt->freeze_map[(constrnt->nfreeze+count)] = i;
          atommaps->nfree_1mol_jmol_typ[jmol_typ] -= 3; 
          atommaps->nfree_jres_jmol_typ[ires] -= 3;
        }/*endif*/
      }/*endfor*/
      constrnt->nfreeze = nfreeze_now;

      /*    PRINTF("nfreeze_now %d \n",nfreeze_now);
            PRINTF("EXITing from fetch_freeze \n"); EXIT(1); */


    }/*endif*/
    /*--------------------------------------------------------------------------*/
  }/*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void freeze_con_check(MDATOM_MAPS *atommaps,MDBOND *bond, MDBEND *bend, 
    MDTORS *tors,MDGRP_BOND_CON *grp_bond_con,
    MDCONSTRNT *constrnt)

  /*==========================================================================*/
  /*begin routine*/{
    /*==========================================================================*/
    /* Define Local variable */
    int k1,k2,k3,k4,i;
    /*==========================================================================*/
    /* I) Check all constrained bonds */ 
    for(i=1;i<=bond->ncon;i++){
      k1 = bond->j1_con[i];
      k2 = bond->j2_con[i];
      if((constrnt->freeze_flag[k1] > 0) || (constrnt->freeze_flag[k2] > 0)){
        PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("You are freezing atoms that are involved in a  \n");
        PRINTF("constrained bond.                              \n");
        PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        FFLUSH(stdout);
        EXIT(1);
      }/*endif*/
    }/*endfor*/      
    /*==========================================================================*/
    /* IV) Check all 2-1 group constraints */ 
    for(i=1;i<=grp_bond_con->num_21;i++){
      k1 = grp_bond_con->j1_21[i];
      k2 = grp_bond_con->j2_21[i];
      if((constrnt->freeze_flag[k1] > 0) || 
          (constrnt->freeze_flag[k2] > 0)) {
        PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("You are freezing atoms that are involved in a  \n");
        PRINTF("constrained group bond.                        \n");
        PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        FFLUSH(stdout);
        EXIT(1);
      }/*endif*/
    }/*endfor*/      
    /*==========================================================================*/
    /* IV.V) Check all 2-3 group constraints */ 
    for(i=1;i<=grp_bond_con->num_23;i++){
      k1 = grp_bond_con->j1_23[i];
      k2 = grp_bond_con->j2_23[i];
      k3 = grp_bond_con->j3_23[i];
      if((constrnt->freeze_flag[k1] > 0) || 
          (constrnt->freeze_flag[k2] > 0) || 
          (constrnt->freeze_flag[k3] > 0)){
        PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("You are freezing atoms that are involved in a  \n");
        PRINTF("constrained group bond.                        \n");
        PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        FFLUSH(stdout);
        EXIT(1);
      }/*endif*/
    }/*endfor*/      
    /*==========================================================================*/
    /* V) Check all 3-3 group constraints */ 
    for(i=1;i<=grp_bond_con->num_33;i++){
      k1 = grp_bond_con->j1_33[i];
      k2 = grp_bond_con->j2_33[i];
      k3 = grp_bond_con->j3_33[i];
      if((constrnt->freeze_flag[k1] > 0) || 
          (constrnt->freeze_flag[k2] > 0) || 
          (constrnt->freeze_flag[k3] > 0)){
        PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("You are freezing atoms that are involved in a  \n");
        PRINTF("constrained group bond.                        \n");
        PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        FFLUSH(stdout);
        EXIT(1);
      }/*endif*/
    }/*endfor*/      
    /*==========================================================================*/
    /* V.V) Check all 4-3 group constraints */ 
    for(i=1;i<=grp_bond_con->num_43;i++){
      k1 = grp_bond_con->j1_43[i];
      k2 = grp_bond_con->j2_43[i];
      k3 = grp_bond_con->j3_43[i];
      k4 = grp_bond_con->j4_43[i];
      if((constrnt->freeze_flag[k1] > 0)|| 
          (constrnt->freeze_flag[k2] > 0) || 
          (constrnt->freeze_flag[k3] > 0) ||
          (constrnt->freeze_flag[k4] > 0)){
        PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("You are freezing atoms that are involved in a  \n");
        PRINTF("constrained group bond.                        \n");
        PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        FFLUSH(stdout);
        EXIT(1);
      }/*endif*/
    }/*endfor*/      
    /*==========================================================================*/
    /* V) Check all 4-6 group constraints */ 
    for(i=1;i<=grp_bond_con->num_46;i++){
      k1 = grp_bond_con->j1_46[i];
      k2 = grp_bond_con->j2_46[i];
      k3 = grp_bond_con->j3_46[i];
      k4 = grp_bond_con->j4_46[i];
      if((constrnt->freeze_flag[k1] > 0) || 
          (constrnt->freeze_flag[k2] > 0) || 
          (constrnt->freeze_flag[k3] > 0)||(constrnt->freeze_flag[k4] > 0)){
        PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("You are freezing atoms that are involved in a  \n");
        PRINTF("constrained group bond.                        \n");
        PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        FFLUSH(stdout);
        EXIT(1);
      }/*endif*/
    }/*endfor*/      
    /*--------------------------------------------------------------------------*/
  }/*end routine*/ 
/*==========================================================================*/
