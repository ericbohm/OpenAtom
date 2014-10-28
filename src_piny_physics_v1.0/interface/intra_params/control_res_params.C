/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: control_res_parms.c                          */
/*                                                                          */
/* This subprogram reads in molecular parameter files and sets              */
/* molecular and intramolecular data sets                                   */
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

void control_res_params(double *tot_memory,MDCLATOMS_INFO *clatoms_info,
    CPATOM_MAPS *cpatom_maps,MDATOM_MAPS *atommaps,
    MDINTRA *mdintra,
    RESBOND_PARSE *resbond_parse,
    BUILD_INTRA *build_intra,
    FILENAME_PARSE *filename_parse,
    FREE_PARSE *free_parse,CLASS_PARSE *class_parse,
    NULL_INTER_PARSE *null_inter_parse,
    char filename[],DICT_INTRA *dict_intra,
    char fun_key[],int jmol_typ)

  /*========================================================================*/
  /*     Begin routine                                                      */
{/*begin routine*/
  /*========================================================================*/
  /*             Local variable declarations                                */

  int iresidue,nresidue,iparm,iresidue_off;
  int i,istart,iend,jatm,iii;
  int mol_only,natm_mol,mol_or_res,nres_bond;

#include "../class_defs/allclass_strip_mdintra.h"

  /*=======================================================================*/
  /* 0) Initialize offsets and flags                                       */

  iresidue_off = atommaps->jres_jmol_typ_strt[jmol_typ]-1;
  nresidue     = MAX(atommaps->nres_1mol_jmol_typ[jmol_typ],1);
  mol_or_res   = 2;
  mol_only     = 0;
  if(atommaps->nres_1mol_jmol_typ[jmol_typ]==0){mol_only=1;}
  atommaps->jatm_jmol_typ_strt[jmol_typ] = clatoms_info->natm_tot+1;

  /*=======================================================================*/
  /* I) Setup  residues of molecule jmol_typ                          */

  for(iresidue=1;iresidue<=nresidue;iresidue++){
    atommaps->jatm_jres_1mol_jmol_typ_strt[(iresidue+iresidue_off)]=
      clatoms_info->natm_tot;
    build_intra->nghost_now = 0;
    /*--------------------------------------------------------------------------*/
    /*  0) Output                                                               */

#ifdef DEBUG
    strcpy(filename, filename_parse->res_param_name[iresidue]);
    PRINTF("**************************************************************\n");
    PRINTF("Reading from residue parameter file %s\n",filename);
    PRINTF("--------------------------------------------------------------\n");
    PRINTF("\n");
#endif

    /*--------------------------------------------------------------------------*/
    /*  A) Get residue name and number of atoms if not a single molecule        */
    if(mol_only != 1){
      /*    i) Read from residue parm file                                 */
      iparm = 1;
      strcpy(filename, filename_parse->res_param_name[iresidue]);
      check_parmfile(filename,&(dict_intra->num_fun_dict),
          &(dict_intra->fun_dict),fun_key,
          nresidue,&nres_bond,mol_or_res);
      fetch_residue_name(filename,dict_intra,atommaps,fun_key,build_intra,
          jmol_typ,iresidue,iresidue_off,iparm,mol_only);
      /*   ii) Read from residue bond file                                 */
      iparm = 0;
      istart = resbond_parse->res_bond1_off[iresidue];
      for(i=1;i<=resbond_parse->nres_bond1_jres[iresidue];i++){
        strcpy(filename,resbond_parse->resbond_prm[(istart+i)].res1_file);
        if(strlen(filename)!=0){
          check_parmfile(filename,&(dict_intra->num_fun_dict),
              &(dict_intra->fun_dict),fun_key,
              nresidue,&nres_bond,mol_or_res);
          fetch_residue_name(filename,dict_intra,atommaps,fun_key,
              build_intra,jmol_typ,iresidue,iresidue_off,iparm,mol_only);
        }/*endif*/
      }/*endfor*/
      iparm = 0;
      istart = resbond_parse->res_bond2_off[iresidue];
      for(i=1;i<=resbond_parse->nres_bond2_jres[iresidue];i++){
        strcpy(filename,resbond_parse->resbond_prm[(istart+i)].res2_file);
        if(strlen(filename)!=0){
          check_parmfile(filename,&(dict_intra->num_fun_dict),
              &(dict_intra->fun_dict),fun_key,
              nresidue,&nres_bond,mol_or_res);
          fetch_residue_name(filename,dict_intra,atommaps,fun_key,
              build_intra,jmol_typ,iresidue,iresidue_off,iparm,mol_only);
        }/*endif*/
      }/*endfor*/
      /*   iii) Read from residue fix file                                 */
      iparm = 0;
      strcpy(filename,filename_parse->res_fix_name[iresidue]);
      if(strlen(filename)!=0){
        check_parmfile(filename,&(dict_intra->num_fun_dict),
            &(dict_intra->fun_dict),fun_key,
            nresidue,&nres_bond,mol_or_res);
        fetch_residue_name(filename,dict_intra,atommaps,fun_key,
            build_intra,jmol_typ,iresidue,iresidue_off,iparm,mol_only);
      }/*endif*/
    }else{
      atommaps->natm_jres_jmol_typ[(iresidue_off+1)] = 
        atommaps->natm_1mol_jmol_typ[jmol_typ];
    }/*endif*/
    /*-------------------------------------------------------------------------*/
    /*  B) Get atm indicies: set_params/set_atm_mask.c and fetch_residue.c     */
    /*    i) Initialize build intra                                            */
    init_build_intra(build_intra,atommaps,iresidue,iresidue_off);
    /*   ii) Read from residue parm file :                                     */
    iparm = 1;
    strcpy(filename, filename_parse->res_param_name[iresidue]);
    fetch_residue_atm_masks(filename,dict_intra,atommaps,
        fun_key,build_intra,jmol_typ,iresidue, 
        iresidue_off,iparm);
    /*   iii) Check atom masks                                                 */
    check_atm_mask(build_intra,dict_intra,atommaps,fun_key,filename,
        iresidue,iresidue_off);
    /*   iv) Read from residue bond files                                      */
    iparm = 0;
    istart = resbond_parse->res_bond1_off[iresidue];
    for(i=1;i<=resbond_parse->nres_bond1_jres[iresidue];i++){
      strcpy(filename,resbond_parse->resbond_prm[(istart+i)].res1_file);
      if(strlen(filename)!=0){
        fetch_residue_atm_masks(filename,dict_intra,atommaps,
            fun_key,build_intra,jmol_typ,iresidue, 
            iresidue_off,iparm);
      }/*endif*/
    }/*endfor*/
    iparm = 0;
    istart = resbond_parse->res_bond2_off[iresidue];
    for(i=1;i<=resbond_parse->nres_bond2_jres[iresidue];i++){
      strcpy(filename,resbond_parse->resbond_prm[(istart+i)].res2_file);
      if(strlen(filename)!=0){
        fetch_residue_atm_masks(filename,dict_intra,atommaps,
            fun_key,build_intra,jmol_typ,iresidue, 
            iresidue_off,iparm);
      }/*endif*/
    }/*endfor*/
    /*  v) Read from the fix file                                           */
    iparm = 0;
    strcpy(filename, filename_parse->res_fix_name[iresidue]);
    if(strlen(filename)!=0){
      fetch_residue_atm_masks(filename,dict_intra,atommaps,
          fun_key,build_intra,jmol_typ,iresidue, 
          iresidue_off,iparm);
    }/*endif*/
    /*  vi) Make the atom indices                                               */
    create_atm_ind(clatoms_info,cpatom_maps,atommaps,mdconstrnt,
        mdghost_atoms,build_intra,dict_intra,fun_key,filename);
    /*------------------------------------------------------------------------*/
    /*  C) Get atms parameters                                                */

    /*    i) Initialize the atom index checker                               */
    for(i=1;i<=build_intra->natm_1res_now;i++){
      (build_intra->iatm_ind_chk)[i]=0;
    }/*endfor*/
    /*   ii) Read from residue parm file                                      */
    iparm = 1;
    strcpy(filename, filename_parse->res_param_name[iresidue]);
    fetch_residue_atm_prms(filename,dict_intra,cpatom_maps,atommaps,
        clatoms_info,mdghost_atoms,mdconstrnt,fun_key,build_intra,
        jmol_typ,iresidue,iresidue_off,iparm);
    /*  iii) Read residue bond files                                          */
    iparm = 0;
    istart = resbond_parse->res_bond1_off[iresidue];
    for(i=1;i<=resbond_parse->nres_bond1_jres[iresidue];i++){
      strcpy(filename,resbond_parse->resbond_prm[(istart+i)].res1_file);
      if(strlen(filename)!=0){
        fetch_residue_atm_prms(filename,dict_intra,cpatom_maps,atommaps,
            clatoms_info,mdghost_atoms,mdconstrnt,fun_key,build_intra,
            jmol_typ,iresidue,iresidue_off,iparm);
      }/*endif*/
    }/*endfor*/
    iparm = 0;
    istart = resbond_parse->res_bond2_off[iresidue];
    for(i=1;i<=resbond_parse->nres_bond2_jres[iresidue];i++){
      strcpy(filename,resbond_parse->resbond_prm[(istart+i)].res2_file);
      if(strlen(filename)!=0){
        fetch_residue_atm_prms(filename,dict_intra,cpatom_maps,atommaps,
            clatoms_info,mdghost_atoms,mdconstrnt,fun_key,build_intra,
            jmol_typ,iresidue,iresidue_off,iparm);
      }/*endif*/
    }/*endfor*/
    /*  iv) Read the residue fix file                                    */
    iparm = 0;
    strcpy(filename, filename_parse->res_fix_name[iresidue]);
    if(strlen(filename)!=0){
      fetch_residue_atm_prms(filename,dict_intra,cpatom_maps,atommaps,
          clatoms_info,mdghost_atoms,mdconstrnt,fun_key,build_intra,
          jmol_typ,iresidue,iresidue_off,iparm);
    }/*endif*/
    /*    v) Check atm indicies and set atom numbers */
    check_atm_ind(build_intra);
    atommaps->natm_jres_jmol_typ[(iresidue_off+iresidue)]=
      build_intra->natm_1res_now;
    /*------------------------------------------------------------------------*/
    /*  D) Get bond/bend/tors/onfo/grp_con params                             */

    /*  i) Initial degrees of freedom and constraint options                  */
    atommaps->nfree_jres_jmol_typ[(iresidue_off+iresidue)] = 
      3*(atommaps->natm_jres_jmol_typ[(iresidue_off+iresidue)]
          - build_intra->nghost_now);
    mdconstrnt->icons_jres_jmol_typ[(iresidue_off+iresidue)] = 0;
    /*  ii) Read from residue parm file                                      */
    iparm = 1;
    strcpy(filename, filename_parse->res_param_name[iresidue]);
    fetch_residue_connectivity(filename,dict_intra,
        atommaps,fun_key,clatoms_info,
        build_intra,mdintra,
        null_inter_parse,jmol_typ,iresidue,iresidue_off,iparm,
        class_parse->mol_hydrog_con_opt[jmol_typ]);
    /*  iii) Read residue bond files                                          */
    iparm = 0;
    istart = resbond_parse->res_bond1_off[iresidue];
    for(i=1;i<=resbond_parse->nres_bond1_jres[iresidue];i++){
      strcpy(filename,resbond_parse->resbond_prm[(istart+i)].res1_file);
      if(strlen(filename)!=0){
        fetch_residue_connectivity(filename,dict_intra,
            atommaps,fun_key,clatoms_info,build_intra,mdintra,
            null_inter_parse,jmol_typ,iresidue,iresidue_off,iparm,
            class_parse->mol_hydrog_con_opt[jmol_typ]);
      }/*endif*/
    }/*endfor*/
    iparm = 0;
    istart = resbond_parse->res_bond2_off[iresidue];
    for(i=1;i<=resbond_parse->nres_bond2_jres[iresidue];i++){
      strcpy(filename,resbond_parse->resbond_prm[(istart+i)].res2_file);
      if(strlen(filename)!=0){
        fetch_residue_connectivity(filename,dict_intra,
            atommaps,fun_key,clatoms_info, build_intra,mdintra,
            null_inter_parse,jmol_typ,iresidue,iresidue_off,iparm,
            class_parse->mol_hydrog_con_opt[jmol_typ]);
      }/*endif*/
    }/*endfor*/
    /*  iv) Read residue fix files                                          */
    iparm = 0;
    strcpy(filename, filename_parse->res_fix_name[iresidue]);
    if(strlen(filename)!=0){
      fetch_residue_connectivity(filename,dict_intra,
          atommaps,fun_key,clatoms_info,
          build_intra,mdintra,
          null_inter_parse,jmol_typ,iresidue,iresidue_off,iparm,
          class_parse->mol_hydrog_con_opt[jmol_typ]);
    }/*endif*/
    /*------------------------------------------------------------------------*/
    /*  E) Use connectivity information stored in the bondsite structure      */
    /*     to fill the res_bondrpm structure for each residue bond to which   */
    /*     the current residue is a member                                    */
    if(resbond_parse->nres_bond1_jres[iresidue]>0 ||
        resbond_parse->nres_bond2_jres[iresidue]>0){
      fetch_resbond_prm(resbond_parse,build_intra->bond_site,
          jmol_typ,iresidue,build_intra->natm_1res_now,
          clatoms_info->natm_tot);
    }/*endif*/
    /*------------------------------------------------------------------------*/
    /*  F) Set free energy indicies                                           */
    natm_mol = atommaps->natm_1mol_jmol_typ[jmol_typ];
    fetch_free_energy_index(build_intra,free_parse,mdintra,
        jmol_typ,iresidue,clatoms_info->natm_tot,natm_mol);
    /*------------------------------------------------------------------------*/
    /*  G) Increment total number of atoms in system                        */
    clatoms_info->natm_tot += build_intra->natm_1res_now;
    /*------------------------------------------------------------------------*/
    /*  H) Output                                                             */
#ifdef DEBUG
    PRINTF("\n");
    PRINTF("--------------------------------------------------------------\n");
    strcpy(filename, filename_parse->res_param_name[iresidue]);
    PRINTF("Completed read from residue parameter file %s\n",filename);
    PRINTF("**************************************************************\n");
    PRINTF("\n");     
#endif

  }/*endfor:iresidue*/

  /*=======================================================================*/
  /* II) Get the number of degrees of freedom of each molecule             */

  iresidue_off = atommaps->jres_jmol_typ_strt[jmol_typ]-1;
  atommaps->nfree_1mol_jmol_typ[jmol_typ] = 0; 
  atommaps->natm_1mol_jmol_typ[jmol_typ]  = 0;
  mdconstrnt->icons_jmol_typ[jmol_typ] = 0;
  nresidue = MAX(atommaps->nres_1mol_jmol_typ[jmol_typ],1);
  for(i=1;i<=nresidue;i++){
    atommaps->nfree_1mol_jmol_typ[jmol_typ] += 
      atommaps->nfree_jres_jmol_typ[(iresidue_off+i)];
    atommaps->natm_1mol_jmol_typ[jmol_typ] += 
      atommaps->natm_jres_jmol_typ[(iresidue_off+i)];
    if(mdconstrnt->icons_jres_jmol_typ[(iresidue_off+i)]==1){
      mdconstrnt->icons_jmol_typ[jmol_typ] = 1;
    }
#ifdef DEBUG
    PRINTF("jmol,ires,mol_free,res_free,natm_1mol,natm_res,icons_mol,icons_res\n");
    PRINTF(" %d   %d     %d        %d       %d          %d        %d      %d\n",
        jmol_typ,i,
        atommaps->nfree_1mol_jmol_typ[jmol_typ],
        atommaps->nfree_jres_jmol_typ[(iresidue_off+i)],
        atommaps->natm_1mol_jmol_typ[jmol_typ],
        atommaps->natm_jres_jmol_typ[(iresidue_off+i)],
        mdconstrnt->icons_jmol_typ[jmol_typ],
        mdconstrnt->icons_jres_jmol_typ[(iresidue_off+i)]);
    mal_verify(1);
    scanf("%d",&iii);
#endif
  }/*endfor*/
  if(natm_mol!=atommaps->natm_1mol_jmol_typ[jmol_typ]){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Error in number of atoms in molecule number %d\n",jmol_typ);
    PRINTF("%d expected %d found\n",natm_mol,
        atommaps->natm_1mol_jmol_typ[jmol_typ]);
    PRINTF("in file %s \n",filename_parse->mol_param_name[jmol_typ]);
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }/*endif*/

  /*=======================================================================*/
  /* III) Build temperature stuff                                          */

  istart = atommaps->jatm_jmol_typ_strt[jmol_typ];
  iend = atommaps->natm_1mol_jmol_typ[jmol_typ]+istart-1;
  for(jatm=istart;jatm<=iend;jatm++){
    clatoms_info->text_atm[jatm] =  class_parse->text_nhc_mol[jmol_typ];
  }/*endfor*/

  /*=======================================================================*/
  /* IV) Debug                                                             */
#ifdef DEBUG
  vomit_intra_list(clatoms_info,atommaps,mdintra,null_inter_parse);
#endif
  /*-----------------------------------------------------------------------*/
} /*end routine */
/*=======================================================================*/


#ifdef DEBUG


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void vomit_intra_list (MDCLATOMS_INFO *clatoms_info,  
    MDATOM_MAPS *atommaps, MDINTRA *mdintra,
    NULL_INTER_PARSE *null_inter_parse)

  /*==========================================================================*/
{ /*begin routine */
  /*==========================================================================*/

  int i,iii,j;

#include "../class_defs/allclass_strip_mdintra.h"

  /*==========================================================================*/
  /* i) Atoms */


  PRINTF("------\n");mal_verify(1);
  PRINTF("Atoms: %d\n",clatoms_info->natm_tot);
  PRINTF("------\n");
  for(i=1;i<=clatoms_info->natm_tot;i++){
    PRINTF("%d mass    %g\n",i,clatoms_info->mass[i]);
    PRINTF("%d q       %g\n",i,clatoms_info->q[i]);            
    PRINTF("%d alp_pol %g\n",i,clatoms_info->alp_pol[i]);     
    PRINTF("%d b_neut  %g\n",i,clatoms_info->b_neut[i]);       
    PRINTF("%d text    %g\n",i,clatoms_info->text_atm[i]);
    PRINTF("%d mol num     %d\n",i,atommaps->iatm_mol_num[i]);
    PRINTF("%d res num     %d\n",i,atommaps->iatm_res_num[i]);
    PRINTF("%d ind mol_typ %d\n",i,atommaps->iatm_mol_typ[i]);
    PRINTF("%d ind res_typ %d\n",i,atommaps->iatm_res_typ[i]);
    PRINTF("%d ind atm_typ %d\n",i,atommaps->iatm_atm_typ[i]);
    PRINTF("%d ind ghost_flag %d\n",i,atommaps->ighost_flag[i]);
    PRINTF("%d ind mol_typ %s\n",i,
        atommaps->mol_typ[atommaps->iatm_mol_typ[i]]);
    PRINTF("%d ind res_typ %s\n",i,
        atommaps->res_typ[atommaps->iatm_res_typ[i]]);
    PRINTF("%d ind atm_typ %s\n",i,
        atommaps->atm_typ[atommaps->iatm_atm_typ[i]]);
    if((i%10)==0){scanf("%d",&iii);}
  }/*endfor*/
  mal_verify(1);PRINTF("Enter an integer ");SCANF("%d",&iii);

  /*==========================================================================*/
  /* ii) Bonds */

  PRINTF("------\n");
  PRINTF("pow bonds: %d\n ",bond->npow);
  PRINTF("------\n");
  for(i=1;i<=bond->npow;i++){
    PRINTF("%d atom1 %d atom2 %d type %d \n",i,
        bond->j1_pow[i],bond->j2_pow[i],
        bond->jtyp_pow[i]);
  }/*endfor*/
  mal_verify(1);PRINTF("Enter an integer ");SCANF("%d",&iii);

  PRINTF("------\n");
  PRINTF("con bonds: %d\n ",bond->ncon);
  PRINTF("------\n");
  for(i=1;i<=bond->ncon;i++){
    PRINTF("%d atom1 %d atom2 %d type %d \n",i,
        bond->j1_con[i],bond->j2_con[i],
        bond->jtyp_con[i]);
  }/*endfor*/
  mal_verify(1);PRINTF("Enter an integer ");SCANF("%d",&iii);

  PRINTF("------\n");
  PRINTF("null bonds: %d\n ",null_inter_parse->nbond_nul);
  PRINTF("------\n");
  for(i=1;i<=null_inter_parse->nbond_nul;i++){
    PRINTF("%d atom1 %d atom2 %d \n",i,
        null_inter_parse->jbond1_nul[i],
        null_inter_parse->jbond2_nul[i]);
  }/*endfor*/
  mal_verify(1);PRINTF("Enter an integer ");SCANF("%d",&iii);

  /*==========================================================================*/
  /* ii) Bends */

  PRINTF("------\n");
  PRINTF("pow bends: %d\n",bend->npow);
  PRINTF("------\n");
  for(i=1;i<=bend->npow;i++){
    PRINTF("%d atom1 %d atom2 %d atom3 %d type %d \n",i,
        bend->j1_pow[i],bend->j2_pow[i],
        bend->j3_pow[i],bend->jtyp_pow[i]);
  }/*endfor*/
  mal_verify(1);PRINTF("Enter an integer ");SCANF("%d",&iii);

  PRINTF("------\n");
  PRINTF("con bends: %d\n",bend->ncon);
  PRINTF("------\n");
  for(i=1;i<=bend->ncon;i++){
    PRINTF("%d atom1 %d atom2 %d atom 3 %d type %d \n",i,
        bend->j1_con[i],bend->j2_con[i],
        bend->j3_con[i],bend->jtyp_con[i]);
  }/*endfor*/
  mal_verify(1);PRINTF("Enter an integer ");SCANF("%d",&iii);

  PRINTF("------\n");
  PRINTF("null bends: %d\n ",null_inter_parse->nbend_nul);
  PRINTF("------\n");
  for(i=1;i<=null_inter_parse->nbend_nul;i++){
    PRINTF("%d atom1 %d atom2 %d atom3 %d \n",i,
        null_inter_parse->jbend1_nul[i],
        null_inter_parse->jbend2_nul[i],
        null_inter_parse->jbend3_nul[i]);
  }/*endfor*/
  mal_verify(1);PRINTF("Enter an integer ");SCANF("%d",&iii);

  /*==========================================================================*/
  /* ii) Bend_bnds */

  PRINTF("------\n");
  PRINTF("bend_bnd: %d\n",bend_bnd->num);
  PRINTF("------\n");
  for(i=1;i<=bend_bnd->num;i++){
    PRINTF("%d atom1 %d atom2 %d atom3 %d type %d \n",i,
        bend_bnd->j1[i],bend_bnd->j2[i],
        bend_bnd->j3[i],bend_bnd->jtyp[i]);
  }/*endfor*/
  mal_verify(1);PRINTF("Enter an integer ");SCANF("%d",&iii);

  /*==========================================================================*/
  /* ii) Tors */
  PRINTF("------\n");
  PRINTF("pow torsions: %d\n",tors->npow);
  PRINTF("------\n");
  for(i=1;i<=tors->npow;i++){
    PRINTF("%d atom1 %d atom2 %d atom3 %d atom 4 %d type %d \n",i,
        tors->j1_pow[i],tors->j2_pow[i],
        tors->j3_pow[i],tors->j4_pow[i],
        tors->jtyp_pow[i]);
  }/*endfor*/
  mal_verify(1);PRINTF("Enter an integer ");SCANF("%d",&iii);

  PRINTF("------\n");
  PRINTF("con torsions: %d\n",tors->ncon);
  PRINTF("------\n");
  for(i=1;i<=tors->ncon;i++){
    PRINTF("%d atom1 %d atom2 %d atom3 %d atom 4 %d type %d \n",i,
        tors->j1_con[i],tors->j2_con[i],
        tors->j3_con[i],tors->j4_con[i],
        tors->jtyp_con[i]);
  }/*endfor*/
  mal_verify(1);PRINTF("Enter an integer ");SCANF("%d",&iii);

  PRINTF("------\n");
  PRINTF("null torsions: %d\n ",null_inter_parse->ntors_nul);
  PRINTF("------\n");
  for(i=1;i<=null_inter_parse->ntors_nul;i++){
    PRINTF("%d atom1 %d atom2 %d atom3 %d atom4 %d\n",i,
        null_inter_parse->jtors1_nul[i],
        null_inter_parse->jtors2_nul[i],
        null_inter_parse->jtors3_nul[i],
        null_inter_parse->jtors4_nul[i]);
  }/*endfor*/
  mal_verify(1);PRINTF("Enter an integer ");SCANF("%d",&iii);

  /*==========================================================================*/
  /* ii) Onfos */

  PRINTF("------\n");
  PRINTF("onfos: %d\n",onfo->num);
  PRINTF("------\n");
  for(i=1;i<=onfo->num;i++){
    PRINTF("%d atom1 %d atom2 %d type %d \n",i,
        onfo->j1[i],onfo->j2[i],
        onfo->jtyp[i]);
  }/*endfor*/
  mal_verify(1);PRINTF("Enter an integer ");SCANF("%d",&iii);

  PRINTF("------\n");
  PRINTF("null onfos: %d\n ",null_inter_parse->nbond_nul);
  PRINTF("------\n");
  for(i=1;i<=null_inter_parse->nbond_nul;i++){
    PRINTF("%d atom1 %d atom2 %d  \n",i,
        null_inter_parse->jbond1_nul[i],
        null_inter_parse->jbond2_nul[i]);
  }/*endfor*/
  mal_verify(1);PRINTF("Enter an integer ");SCANF("%d",&iii);

  /*==========================================================================*/
  /* ii) Ghosts */

  PRINTF("------\n");
  PRINTF("ghost atoms: %d\n",ghost_atoms->nghost_tot);
  PRINTF("------\n");

  for(i=1;i<=ghost_atoms->nghost_tot;i++){
    PRINTF("atm %d is ghost number %d and is formed from %d atms\n",
        ghost_atoms->ighost_map[i],i,ghost_atoms->natm_comp[i]);
    for(j=1;j<=ghost_atoms->natm_comp[i];j++){
      PRINTF("atm %d atm_ind %d coef %g\n",
          j,ghost_atoms->iatm_comp[j][i],ghost_atoms->coef[j][i]);
    }/*endfor*/
  }/*endfor*/
  mal_verify(1);PRINTF("Enter an integer ");SCANF("%d",&iii);

  /*==========================================================================*/
  /* i) Grp con 21 */

  PRINTF("------\n");mal_verify(1);
  PRINTF("Grp con 21: %d \n",grp_bond_con->num_21);
  PRINTF("------\n");

  for(i=1;i<=grp_bond_con->num_21;i++){
    PRINTF("%d atom1 %d atom2 %d type %d \n",i,
        grp_bond_con->j1_21[i],grp_bond_con->j2_21[i],
        grp_bond_con->jtyp_21[i]);
  }/*endfor*/
  mal_verify(1);PRINTF("Enter an integer ");SCANF("%d",&iii);

  /*==========================================================================*/
  /* i) Grp con 23 */

  PRINTF("------\n");mal_verify(1);
  PRINTF("Grp con 23: %d \n",grp_bond_con->num_23);
  PRINTF("------\n");

  for(i=1;i<=grp_bond_con->num_23;i++){
    PRINTF("%d atom1 %d atom2 %d atom3 %d type %d \n",i,
        grp_bond_con->j1_23[i],grp_bond_con->j2_23[i],
        grp_bond_con->j3_23[i],
        grp_bond_con->jtyp_23[i]);
  }/*endfor*/
  mal_verify(1);PRINTF("Enter an integer ");SCANF("%d",&iii);

  /*==========================================================================*/
  /* i) Grp con 33 */

  PRINTF("------\n");mal_verify(1);
  PRINTF("Grp con 33: %d \n",grp_bond_con->num_33);
  PRINTF("------\n");

  for(i=1;i<=grp_bond_con->num_33;i++){
    PRINTF("%d atom1 %d atom2 %d atom3 %d type %d \n",i,
        grp_bond_con->j1_33[i],grp_bond_con->j2_33[i],
        grp_bond_con->j3_33[i],
        grp_bond_con->jtyp_33[i]);
  }/*endfor*/
  mal_verify(1);PRINTF("Enter an integer ");SCANF("%d",&iii);

  /*==========================================================================*/
  /* i) Grp Watts 33 */

  PRINTF("------\n");mal_verify(1);
  PRINTF("Grp Watts 33: %d \n",grp_bond_watts->num_33);
  PRINTF("------\n");

  for(i=1;i<=grp_bond_watts->num_33;i++){
    PRINTF("%d atom1 %d atom2 %d atom3 %d type %d \n",i,
        grp_bond_watts->j1_33[i],grp_bond_watts->j2_33[i],
        grp_bond_watts->j3_33[i],
        grp_bond_watts->jtyp_33[i]);
  }/*endfor*/
  mal_verify(1);PRINTF("Enter an integer ");SCANF("%d",&iii);

  /*==========================================================================*/
  /* i) Grp con 43 */

  PRINTF("------\n");mal_verify(1);
  PRINTF("Grp con 43: %d \n",grp_bond_con->num_43);
  PRINTF("------\n");

  for(i=1;i<=grp_bond_con->num_43;i++){
    PRINTF("%d atom1 %d atom2 %d atom3 %d atom4 %d type %d \n",i,
        grp_bond_con->j1_43[i],grp_bond_con->j2_43[i],
        grp_bond_con->j3_43[i],
        grp_bond_con->j4_43[i],
        grp_bond_con->jtyp_43[i]);
  }/*endfor*/
  mal_verify(1);PRINTF("Enter an integer ");SCANF("%d",&iii);

  /*==========================================================================*/
  /* i) Grp con 46 */

  PRINTF("------\n");mal_verify(1);
  PRINTF("Grp con 46: %d \n",grp_bond_con->num_46);
  PRINTF("------\n");

  for(i=1;i<=grp_bond_con->num_46;i++){
    PRINTF("%d atom1 %d atom2 %d atom3 %d atom4 %d type %d \n",i,
        grp_bond_con->j1_46[i],grp_bond_con->j2_46[i],
        grp_bond_con->j3_46[i],grp_bond_con->j4_46[i],
        grp_bond_con->jtyp_46[i]);
  }/*endfor*/

  mal_verify(1);PRINTF("Enter an integer ");SCANF("%d",&iii);

  /*==========================================================================*/
  /* ii) Free energy */

  PRINTF("------\n");mal_verify(1);
  PRINTF("Bond free energy: %d \n",bond_free->num);
  PRINTF("------\n");
  if(bond_free->num>0){
    PRINTF("file name %s atom1 %d atom2 %d\n",
        bond_free->file,bond_free->j1,bond_free->j2);
    PRINTF("fk %g eq %g npow %d \n",
        bond_free->fk,bond_free->eq,bond_free->npow);
    PRINTF("nhist %d del %g rmin %g rmax %g\n", 
        bond_free->nhist,bond_free->del,
        bond_free->rmin,bond_free->rmax);
  }/*endif*/

  PRINTF("------\n");mal_verify(1);
  PRINTF("Bend free energy: %d \n",bend_free->num);
  PRINTF("------\n");
  if(bend_free->num>0){
    PRINTF("file name %s atom1 %d atom2 %d atom3 %d\n",
        bend_free->file,bend_free->j1,bend_free->j2,
        bend_free->j3);
    PRINTF("fk %g eq %g npow %d \n",
        bend_free->fk,bend_free->eq,bend_free->npow);
    PRINTF("nhist %d del %g \n", 
        bend_free->nhist,bend_free->del);
  }/*endif*/

  PRINTF("------\n");mal_verify(1);
  PRINTF("Tors free energy: %d \n",tors_free->num);
  PRINTF("------\n");
  if(tors_free->num>0){
    PRINTF("file name %s atom1 %d atom2 %d atom3 %d atom4 %d\n",
        tors_free->file,tors_free->j1,tors_free->j2,
        tors_free->j3,tors_free->j4);
    PRINTF("fk %g eq %g npow %d \n",
        tors_free->fk,tors_free->eq,tors_free->npow);
    PRINTF("nhist %d del %g \n", 
        tors_free->nhist,tors_free->del);
  }/*endif*/
  mal_verify(1);PRINTF("Enter an integer ");SCANF("%d",&iii);

  /*=======================================================================*/
}/*end routine */
/*=======================================================================*/

void mal_verify(int i){}

#endif

