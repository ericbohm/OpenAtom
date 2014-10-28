/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: control_set_mol_params.c                     */
/*                                                                          */
/* This subprogram sets molecule/cp setup parameters                        */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../class_defs/allclass_gen.h"
#include "../class_defs/allclass_cp.h"
#include "../class_defs/allclass_mdintegrate.h"
#include "../class_defs/allclass_mdatoms.h"
#include "../class_defs/allclass_mdinter.h"
#include "../class_defs/allclass_mdintra.h"
#include "../class_defs/typedefs_par.h"
#include "../proto_defs/proto_mol_params_local.h"
#include "../proto_defs/proto_handle_entry.h"
#include "../proto_defs/proto_friend_lib_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void control_set_mol_params(MDATOM_MAPS *mdatom_maps,CPOPTS *cpopts,
    CPCOEFFS_INFO *cpcoeffs_info,CP_PARSE *cp_parse,
    CLASS_PARSE *class_parse,
    MDINTRA *mdintra,MDSURFACE *mdsurface,
    FILENAME_PARSE *filename_parse,
    FREE_PARSE *free_parse,
    DICT_MOL *dict_mol,DICT_WORD *word,
    char *fun_key,int *nfun_key,
    int iextend,double t_ext,int ifirst,int pi_beads)
  /*=======================================================================*/
{/*begin routine*/
  /*=======================================================================*/
  /*             Local variable declarations                               */

  int   nline,nkey,i,num,ierr;
  FILE *fp;
  int  *mol_ind_chk;                                 /* Index check      */
  int  *ifound;

  int nmol_typ            = mdatom_maps->nmol_typ;
  int *nmol_jmol_typ      = mdatom_maps->nmol_jmol_typ;
  int *nres_1mol_jmol_typ = mdatom_maps->nres_1mol_jmol_typ;

  int num_fun_dict        = dict_mol->num_fun_dict;
  char *molsetname        = filename_parse->molsetname;

  int bond_free_num       = mdintra->mdbond_free.num;
  int bend_free_num       = mdintra->mdbend_free.num;
  int tors_free_num       = mdintra->mdtors_free.num;
  int rbar_sig_free_iopt  = mdintra->mdrbar_sig_free.iopt;

  int *imoltyp_rbar1_free;
  int *imoltyp_rbar2_free;
  int *imol_rbar1_free;
  int *imol_rbar2_free;
  int *ires_rbar1_free;
  int *ires_rbar2_free;
  int nbar_bond;

  int *imoltyp_bond_free = free_parse->imoltyp_bond_free;
  int *imol_bond_free    = free_parse->imol_bond_free;
  int *ires_bond_free    = free_parse->ires_bond_free;

  int *imoltyp_bend_free = free_parse->imoltyp_bend_free;
  int *imol_bend_free    = free_parse->imol_bend_free;
  int *ires_bend_free    = free_parse->ires_bend_free;

  int *imoltyp_tors_free = free_parse->imoltyp_tors_free;
  int *imol_tors_free    = free_parse->imol_tors_free;
  int *ires_tors_free    = free_parse->ires_tors_free; 

  /*=======================================================================*/
  /* 0) Set up molecular index checking memory                             */

  mol_ind_chk     = (int *) cmalloc(nmol_typ*sizeof(int),"control_set_mol_params")-1;
  ifound          = (int *) cmalloc(num_fun_dict*sizeof(int),"control_set_mol_params")-1;
  for(i=1;i<=nmol_typ;i++){ mol_ind_chk[i]=0;}
  for(i=1;i<=num_fun_dict;i++){ ifound[i]=0;}

  /*=======================================================================*/
  /* I) Fetch a valid functional key word from molset file                 */

  fp = cfopen((const char *) molsetname,"r");

  nline          = 1;
  *nfun_key      = 0;
  word->iuset    = 0;
  word->key_type = 0;
  while(get_fun_key(fp,fun_key,&nline,nfun_key,molsetname)){
    get_fun_key_index(fun_key,dict_mol->num_fun_dict,dict_mol->fun_dict,
        nline,*nfun_key,molsetname,&num);

    /*=======================================================================*/
    /* II) Fetch the key words and key args of the functional key word       */
    /*     and stick them into the correct dictionary                        */

    if(num==1){
      set_mol_dict(&(dict_mol->mol_dict),&(dict_mol->num_mol_dict),
          iextend,class_parse->tau_nhc_def,t_ext,ifirst);
    }
    nkey=0;
    while(get_word(fp,word,&nline,&nkey,*nfun_key,molsetname)){
      switch(num){
        case 1:
          put_word_dict(word,dict_mol->mol_dict,dict_mol->num_mol_dict,
              fun_key,nline,nkey,*nfun_key,molsetname);
          break;
        case 2:
          put_word_dict(word,dict_mol->wave_dict,
              dict_mol->num_wave_dict,fun_key,nline,nkey,
              *nfun_key,molsetname);
          break;
        case 3:
          put_word_dict(word,dict_mol->bond_free_dict,
              dict_mol->num_bond_free_dict,fun_key,nline,
              nkey,*nfun_key,molsetname);
          break;
        case 4:
          put_word_dict(word,dict_mol->bend_free_dict,
              dict_mol->num_bend_free_dict,fun_key,nline,
              nkey,*nfun_key,molsetname);
          break;
        case 5:
          put_word_dict(word,dict_mol->tors_free_dict,
              dict_mol->num_tors_free_dict,fun_key,nline,
              nkey,*nfun_key,molsetname);
          break;
        case 6:
          put_word_dict(word,dict_mol->def_base_dict,
              dict_mol->num_def_base_dict,fun_key,nline,
              nkey,*nfun_key,molsetname);
          break;
        case 7:
          put_word_dict(word,dict_mol->user_base_dict,
              dict_mol->num_user_base_dict,fun_key,nline,
              nkey,*nfun_key,molsetname);
          break;
        case 8:
          put_word_dict(word,dict_mol->rbar_free_dict,
              dict_mol->num_rbar_free_dict,fun_key,nline,
              nkey,*nfun_key,molsetname);
          break;
        case 9:
          put_word_dict(word,dict_mol->surface_dict,
              dict_mol->num_surface_dict,fun_key,nline,
              nkey,*nfun_key,molsetname);
          break;
      } /* end switch dictionary fills*/
    } /* endwhile fetching keywords and keyargs*/

    /*=====================================================================*/
    /* III) Assign the key args of the key words to the appropriate        */
    /*      program variables                                              */

    ifound[num]=1;
    switch(num){
      case 1:
        set_mol_params(filename_parse,fun_key,dict_mol->mol_dict,
            dict_mol->num_mol_dict,
            class_parse,mdatom_maps,mol_ind_chk,pi_beads);
        break;

      case 2:
        set_wave_params(molsetname,fun_key,
            dict_mol->wave_dict,dict_mol->num_wave_dict,
            cpopts,cpcoeffs_info,cp_parse);
        break;

      case 3:
        set_bond_free_params(molsetname,fun_key,
            dict_mol->bond_free_dict,
            dict_mol->num_bond_free_dict,
            &(mdintra->mdbond_free),free_parse,
            mdatom_maps->nmol_typ);
        break;
      case 4:
        set_bend_free_params(molsetname,fun_key,
            dict_mol->bend_free_dict,
            dict_mol->num_bend_free_dict,
            &(mdintra->mdbend_free),free_parse,
            mdatom_maps->nmol_typ);
        break;

      case 5:
        set_tors_free_params(molsetname,fun_key,
            dict_mol->tors_free_dict,
            dict_mol->num_tors_free_dict,
            &(mdintra->mdtors_free),free_parse,
            mdatom_maps->nmol_typ);
        break;

      case 6:
        set_def_base_params(filename_parse,dict_mol->def_base_dict,
            dict_mol->num_def_base_dict);
        break;

      case 7:
        set_user_base_params(filename_parse,dict_mol->user_base_dict,
            dict_mol->num_user_base_dict);
        break;

      case 8:
        set_rbar_free_params(molsetname,fun_key,
            dict_mol->rbar_free_dict,
            dict_mol->num_rbar_free_dict,
            &(mdintra->mdrbar_sig_free),free_parse,
            mdatom_maps->nmol_typ);
        break;

      case 9:
        set_surf_params(molsetname,fun_key,
            dict_mol->surface_dict,
            dict_mol->num_surface_dict,
            mdsurface); 
        break;
    } /* end switch assigning dict stuff to variables */

  } /*endwhile getting functional keywords data*/

  /*=====================================================================*/
  /* IV) Make sure everything that is required has been found            */ 
  /*     if it hasn't either die or define defaults                      */ 

  if(ifound[1]==0){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("  Required functional keyword %s not found     \n",
        dict_mol->fun_dict[1].keyword);
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }/*endif*/
  if(ifound[6]==0){
    set_def_base_params(filename_parse,dict_mol->def_base_dict,
        dict_mol->num_def_base_dict);
  }/*endif*/
  if(ifound[7]==0){
    set_user_base_params(filename_parse,dict_mol->user_base_dict,
        dict_mol->num_user_base_dict);
  }/*endif*/

  /*=======================================================================*/
  /* IV) Check  indices                                                    */

  for(i = 1;i<=nmol_typ;i++){
    if(mol_ind_chk[i]!=1){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@  \n");
      PRINTF("Molecule number %d specified %d times in file %s \n",
          i,mol_ind_chk[i],molsetname);
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@  \n");
      FFLUSH(stdout);
      EXIT(1);
    }/*endif*/
  } /*endfor*/

  /*=======================================================================*/
  /* V) Check  free energy indices                                         */

  ierr = 0;

  if(bond_free_num>0){
    for(i=1;i<=2;i++){
      if(imol_bond_free[i]>nmol_jmol_typ[imoltyp_bond_free[i]]){ierr=i;}
      if(ires_bond_free[i]>MAX(nres_1mol_jmol_typ[imoltyp_bond_free[i]],1))
      {ierr=i;}
    }/*endfor*/
  }/*endif*/

  if(bend_free_num>0){
    for(i=1;i<=3;i++){
      if(imol_bend_free[i]>nmol_jmol_typ[imoltyp_bend_free[i]]){ierr=i;}
      if(ires_bend_free[i]>MAX(nres_1mol_jmol_typ[imoltyp_bend_free[i]],1))
      {ierr=i;}
    }/*endfor*/
  }/*endif*/

  if(tors_free_num>0){tors_free_num=mdintra->mdtors_free.num;}
  if(tors_free_num==1){
    for(i=1;i<=4;i++){
      if(imol_tors_free[i]>nmol_jmol_typ[imoltyp_tors_free[i]]){ierr=i;}
      if(ires_tors_free[i]>MAX(nres_1mol_jmol_typ[imoltyp_tors_free[i]],1))
      {ierr=i;}
    }/*endfor*/
  }/*endif*/
  if(tors_free_num==2){
    for(i=1;i<=8;i++){
      if(imol_tors_free[i]>nmol_jmol_typ[imoltyp_tors_free[i]]){ierr=i;}
      if(ires_tors_free[i]>MAX(nres_1mol_jmol_typ[imoltyp_tors_free[i]],1))
      {ierr=i;}
    }/*endfor*/
  }/*endif*/

  if(rbar_sig_free_iopt>0){

    imoltyp_rbar1_free = free_parse->imoltyp_rbar1_free;
    imoltyp_rbar2_free = free_parse->imoltyp_rbar2_free;
    imol_rbar1_free    = free_parse->imol_rbar1_free;
    imol_rbar2_free    = free_parse->imol_rbar2_free;
    ires_rbar1_free    = free_parse->ires_rbar1_free;
    ires_rbar2_free    = free_parse->ires_rbar2_free;
    nbar_bond          = free_parse->nbar_bond;

    for(i=1;i<=nbar_bond;i++){
      if(imol_rbar1_free[i]>nmol_jmol_typ[imoltyp_rbar1_free[i]]){ierr=1;}
      if(imol_rbar2_free[i]>nmol_jmol_typ[imoltyp_rbar2_free[i]]){ierr=2;}
      if(ires_rbar1_free[i]>
          MAX(nres_1mol_jmol_typ[imoltyp_rbar1_free[i]],1)){ierr=1;}
      if(ires_rbar2_free[i]>
          MAX(nres_1mol_jmol_typ[imoltyp_rbar2_free[i]],1)){ierr=2;}
    }/*endfor*/
  }/*endif*/

  if(ierr>0){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Error specifing the %dth atom in a free energy def \n",ierr);
    PRINTF("in set up file %s \n",molsetname);
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }/*endif*/ 

  /*=======================================================================*/
  /* VI) Free memory                                                       */

  cfree(&mol_ind_chk[1],"control_set_mol_params");
  cfree(&ifound[1],"control_set_mol_params");

  /*-----------------------------------------------------------------------*/
}/*end routine*/
/*=======================================================================*/


