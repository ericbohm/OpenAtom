/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                  Mask routines                                           */
/*                                                                          */
/*==========================================================================*/


#include "standard_include.h"
#include "../class_defs/allclass_mdintegrate.h"
#include "../class_defs/allclass_mdatoms.h"
#include "../class_defs/allclass_mdinter.h"
#include "../class_defs/allclass_mdintra.h"
#include "../class_defs/allclass_cp.h"
#include "../class_defs/typedefs_par.h"
#include "../proto_defs/proto_intra_params_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_handle_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void check_atm_mask(BUILD_INTRA *build_intra,DICT_INTRA *dict_intra,
    MDATOM_MAPS *atommaps, char fun_key[], char filename[],
    int iresidue, int iresidue_off)

  /*========================================================================*/
  /*     Begin routine                                                      */
{/*begin routine*/
  int index,i;
  /*========================================================================*/
  /* Make sure you got the number of atoms in pure residue correct          */

  build_intra->natmind_1res_now = build_intra->natm_1res_now;
  if(build_intra->natm_1res_now!=
      atommaps->natm_jres_jmol_typ[iresidue_off+iresidue]){
    strcpy(fun_key,"res_name_def");
    index = 3;
    keyarg_barf(dict_intra->mol_name_dict,filename,fun_key,index);
  }/*endif*/

  /*========================================================================*/
  /* Make sure you got the atom indices of pure residue correct             */

  for(i=1;i<=build_intra->natm_1res_now;i++){
    if(build_intra->mask_atm[i]!=1){
      strcpy(fun_key,"res_name_def");
      index = 3;
      keyarg_barf(dict_intra->atm_dict,filename,fun_key,index);
    }/*endif*/
  }/*endfor*/

  /*-----------------------------------------------------------------------*/
}/*end routine*/
/*========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void create_atm_ind(MDCLATOMS_INFO *clatoms_info,CPATOM_MAPS *cpatom_maps,
    MDATOM_MAPS *atommaps,MDCONSTRNT *mdconstrnt,
    MDGHOST_ATOMS *mdghost_atoms, BUILD_INTRA *build_intra,
    DICT_INTRA *dict_intra,char fun_key[],char filename[])

  /*========================================================================*/
  /*     Begin routine                                                      */
{/*begin routine*/
  int natm_now,i,index;
  /*========================================================================*/
  /* Check masks                                                            */

  natm_now = 0;
  for(i=1;i<=build_intra->natmind_1res_now;i++){
    if((build_intra->mask_atm[i]>1)||(build_intra->mask_atm[i]<0)){
      strcpy(fun_key,"atom_create or atom_destroy");
      index = 2;
      keyarg_barf(dict_intra->atm_dict,filename,fun_key,index);
    }/*endif*/
    natm_now += (build_intra->mask_atm)[i] ;
  }/*endfor*/

  if(natm_now!=build_intra->natm_1res_now){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("INTERNAL ERROR: atom masks incorrectly set ");
    PRINTF("in control_res_params.c\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);  
    EXIT(1);
  }/*endfor*/

  /*========================================================================*/
  /* Build indicies                                                         */

  build_intra->index_atm[1] = (build_intra->mask_atm)[1];
  for(i=2;i<=build_intra->natmind_1res_now;i++){
    build_intra->index_atm[i] = (build_intra->mask_atm[i] + 
        build_intra->index_atm[(i-1)]);
  }/*endfor*/

  /*========================================================================*/
  /* Reallocate atom memory                                                 */

  if((clatoms_info->natm_tot+build_intra->natm_1res_now) > 
      build_intra->natm_tot_max){
    build_intra->natm_tot_max += MAX(NMEM_MIN,natm_now);
    clatoms_info->mass = (double *)crealloc(&((clatoms_info->mass[1])),
        build_intra->natm_tot_max*sizeof(double),"create_atm_ind")-1;
    clatoms_info->q         = (double *)crealloc(&((clatoms_info->q[1])),
        build_intra->natm_tot_max*sizeof(double),"create_atm_ind")-1;
    cpatom_maps->cp_vlnc_up =(int *)crealloc(&((cpatom_maps->cp_vlnc_up[1])),
        build_intra->natm_tot_max*sizeof(int),"create_atm_ind")-1;
    cpatom_maps->cp_vlnc_dn =(int *)crealloc(&((cpatom_maps->cp_vlnc_dn[1])),
        build_intra->natm_tot_max*sizeof(int),"create_atm_ind")-1;
    cpatom_maps->cp_vlnc_true_up =(int *)crealloc(&((cpatom_maps->cp_vlnc_true_up[1])),
        build_intra->natm_tot_max*sizeof(int),"create_atm_ind")-1;
    cpatom_maps->cp_vlnc_true_dn =(int *)crealloc(&((cpatom_maps->cp_vlnc_true_dn[1])),
        build_intra->natm_tot_max*sizeof(int),"create_atm_ind")-1;
    cpatom_maps->cp_atm_flag=
      (int *)crealloc(&((cpatom_maps->cp_atm_flag[1])),
          build_intra->natm_tot_max*sizeof(int),"create_atm_ind")-1;
    clatoms_info->alp_pol =(double *)crealloc(&((clatoms_info->alp_pol[1])),
        build_intra->natm_tot_max*sizeof(double),"create_atm_ind")-1;
    clatoms_info->b_neut  = (double *)crealloc(&((clatoms_info->b_neut[1])),
        build_intra->natm_tot_max*sizeof(double),"create_atm_ind")-1;
    clatoms_info->text_atm=(double *)crealloc(&((clatoms_info->text_atm[1])),
        build_intra->natm_tot_max*sizeof(double),"create_atm_ind")-1;
    atommaps->iatm_mol_typ  = (int *) crealloc(&(atommaps->iatm_mol_typ[1]),
        build_intra->natm_tot_max*sizeof(int),"create_atm_ind")-1;
    atommaps->iatm_atm_typ  = (int *) crealloc(&(atommaps->iatm_atm_typ[1]),
        build_intra->natm_tot_max*sizeof(int),"create_atm_ind")-1;
    atommaps->iatm_res_typ  = (int *) crealloc(&(atommaps->iatm_res_typ[1]),
        build_intra->natm_tot_max*sizeof(int),"create_atm_ind")-1;
    atommaps->iatm_mol_num  = (int *) crealloc(&(atommaps->iatm_mol_num[1]),
        build_intra->natm_tot_max*sizeof(int),"create_atm_ind")-1;
    atommaps->iatm_res_num  = (int *) crealloc(&(atommaps->iatm_res_num[1]),
        build_intra->natm_tot_max*sizeof(int),"create_atm_ind")-1;
    mdghost_atoms->ighost_flag  = 
      (int *) crealloc(&(mdghost_atoms->ighost_flag[1]),
          build_intra->natm_tot_max*sizeof(int),"create_atm_ind")-1;
    mdconstrnt->freeze_flag  = (int *)crealloc(&(mdconstrnt->freeze_flag[1]),
        build_intra->natm_tot_max*sizeof(int),"create_atm_ind")-1;
    mdconstrnt->atom_label  = (int *)crealloc(&(mdconstrnt->atom_label[1]),
        build_intra->natm_tot_max*sizeof(int),"create_atm_ind")-1;
  }/*endif*/

  /*========================================================================*/
  /* Reallocate bond_site memory                                            */

  if(build_intra->natm_1res_now>build_intra->natm_1res_max){
    build_intra->natm_1res_max += 
      MAX(NMEM_MIN,build_intra->natm_1res_now-build_intra->natm_1res_max);
    build_intra->bond_site  = (BOND_SITE *)
      crealloc(&(build_intra->bond_site[1]),build_intra->natm_1res_max*
          sizeof(BOND_SITE),"create_atm_ind")-1;
    build_intra->iatm_ind_chk = (int *)
      crealloc(&(build_intra->iatm_ind_chk[1]),build_intra->natm_1res_max*
          sizeof(int),"create_atm_ind")-1;
  }

  /*-----------------------------------------------------------------------*/
}/*end routine*/
/*========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void check_atm_ind(BUILD_INTRA *build_intra)

  /*========================================================================*/
  /*     Begin routine                                                      */
{/*begin routine*/
  int i,iii;
  /*========================================================================*/

  for(i=1;i<=build_intra->natm_1res_now;i++){
    if(build_intra->iatm_ind_chk[i]!=1){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("INTERNAL ERROR: atom indices off ");
      PRINTF("in control_res_params.c\n");
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);  
      EXIT(1);
    }/*endif*/
  }/*endfor*/

  /*-----------------------------------------------------------------------*/
}/*end routine*/
/*========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void init_build_intra(BUILD_INTRA *build_intra, MDATOM_MAPS *atommaps,
    int iresidue, int iresidue_off)

  /*========================================================================*/
  /*     Begin routine                                                      */
{/*begin routine*/
  int i;
  /*========================================================================*/

  build_intra->natm_1res_now      = 0;
  build_intra->natmind_1res_now   = 0;
  build_intra->natm_1res_pure_now = 
    atommaps->natm_jres_jmol_typ[iresidue_off+iresidue];

  if(build_intra->natm_1res_pure_now > 
      build_intra->natmind_1res_max){
    build_intra->natmind_1res_max += 
      MAX(NMEM_MIN,(atommaps->natm_jres_jmol_typ)[iresidue_off+iresidue]
          - (build_intra->natmind_1res_max));
    build_intra->mask_atm = (int *)
      crealloc(&(build_intra->mask_atm[1]),build_intra->natmind_1res_max*
          sizeof(int),"init_build_intra")-1;
    build_intra->index_atm = (int *)
      crealloc(&(build_intra->index_atm[1]),build_intra->natmind_1res_max*
          sizeof(int),"init_build_intra")-1;
  }/*endif*/

  for(i=1;i<=build_intra->natmind_1res_max;i++){build_intra->mask_atm[i]=0;}

  /*-----------------------------------------------------------------------*/
}/*end routine*/
/*========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_atm_mask(DICT_WORD atm_dict[],int num_atm_dict,
    char fun_key[],char filename[],
    BUILD_INTRA *build_intra,int num, int jres,int jmol_typ)

  /*=======================================================================*/
{/*begin routine*/
  /*=======================================================================*/
  /*          Local Variables                                             */
  int index;
  int iatm_ind;
  /*=======================================================================*/
  /* Error check */

  if(num==9 || num==10){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Atom creates and atom destroys not permitted \n");
    PRINTF("in the residue parameter file %s\n",filename);
    PRINTF("of the %d residue of the %d molecule\n",jres,jmol_typ);
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }/*endif*/

  /*=======================================================================*/
  /* Set the mask increment atoms counters                                 */

  sscanf(atm_dict[2].keyarg,"%d",&iatm_ind);
  index = 2;
  if((iatm_ind<0)||(iatm_ind>build_intra->natm_1res_pure_now)){
    keyarg_barf(atm_dict,filename,fun_key,index);
  }
  build_intra->mask_atm[iatm_ind]++;
  build_intra->natm_1res_now    ++;

  /*=======================================================================*/
}  /*end routine*/
/*==========================================================================*/




/*==========================================================================*/

void set_atm_mask_rb(DICT_WORD atm_dict[],int num_atm_dict,
    char *fun_key,char *filename,
    BUILD_INTRA *build_intra, int num, int jres,int jmol_typ)

  /*==========================================================================*/
{/*begin routine*/
  /*==========================================================================*/
  /*              Local Variables                                             */
  int i,index;
  int iatm_ind,itemp;
  /*==========================================================================*/

  if(num==2){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Atom definitions not permitted \n");
    PRINTF("in the residue morphing file %s\n",filename);
    PRINTF("of the %d residue of the %d molecule\n",jres,jmol_typ);
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }

  /*==========================================================================*/
  /* Get atom index                                                           */

  sscanf(atm_dict[2].keyarg,"%d",&iatm_ind);
  index = 2;

  /*==========================================================================*/
  /* Create new atom or destroy existing atoms                                */

  if(num==9){
    if((iatm_ind<0)||(iatm_ind<=build_intra->natm_1res_pure_now)){
      keyarg_barf(atm_dict,filename,fun_key,index);}
    build_intra->natm_1res_now += 1;
  }/*endif*/

  if(num==8){
    if((iatm_ind<0)||(iatm_ind>build_intra->natm_1res_pure_now)){
      keyarg_barf(atm_dict,filename,fun_key,index);}
    build_intra->natm_1res_now     -= 1;
    /*endif*/}

  /*==========================================================================*/
  /* Check array sizes                                                        */

  itemp = MAX(build_intra->natmind_1res_now,iatm_ind);
  build_intra->natmind_1res_now=itemp;
  if(build_intra->natmind_1res_now>build_intra->natmind_1res_max){
    build_intra->natmind_1res_max += 
      MAX(NMEM_MIN,build_intra->natmind_1res_now-build_intra->natmind_1res_max);
    build_intra->mask_atm = (int *)
      crealloc(&(build_intra->mask_atm[1]),
          build_intra->natmind_1res_max*sizeof(int),"set_atm_mask_rb")-1;
    build_intra->index_atm = (int *)
      crealloc(&(build_intra->index_atm[1]),
          build_intra->natmind_1res_max*sizeof(int),"set_atm_mask_rb")-1;
    build_intra->bond_site  = (BOND_SITE *)
      crealloc(&(build_intra->bond_site[1]),
          build_intra->natmind_1res_max*sizeof(BOND_SITE),"set_atm_mask_rb")-1;
    build_intra->iatm_ind_chk= (int *)
      crealloc(&(build_intra->bond_site[1]),
          build_intra->natmind_1res_max*sizeof(int),"set_atm_mask_rb")-1;
    for(i=build_intra->natmind_1res_now+1;
        i<=build_intra->natmind_1res_max;i++){build_intra->mask_atm[i]=0;}
    /*endif*/}

  /*==========================================================================*/
  /* Set atom masks                                                           */

  if(strcasecmp(fun_key,"atom_create")==0){
    build_intra->mask_atm[iatm_ind]++;
  }/*endif*/
  if(strcasecmp(fun_key,"atom_destroy")==0){
    build_intra->mask_atm[iatm_ind]--;
    /*endif*/}

  /*==========================================================================*/
  /*end routine*/}
  /*==========================================================================*/





