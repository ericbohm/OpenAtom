/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: set_sim_parms.c                              */
/*                                                                          */
/* This subprogram sets the simulation inputs                               */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../class_defs/allclass_gen.h"
#include "../class_defs/allclass_mdintegrate.h"
#include "../class_defs/allclass_mdatoms.h"
#include "../class_defs/allclass_cp.h"
#include "../class_defs/allclass_mdinter.h"
#include "../class_defs/allclass_mdintra.h"
#include "../class_defs/typedefs_par.h"
#include "../proto_defs/proto_sim_params_local.h"
#include "../proto_defs/proto_handle_entry.h"
#include "../proto_defs/proto_friend_lib_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void set_sim_params_temper(GENTEMPERING_CTRL *tempering_ctrl,
    GENSIMOPTS *gensimopts,DICT_WORD *dict,char *fun_key,
    char *input_name)
  /*==========================================================================*/
{/* begin routine */
  /*==========================================================================*/

  int ifound,index,i;
  int int_key_arg,iopt;
  int npara_temps;
  double real_key_arg;
  double *t_ext,*p_ext,*s_ext;
  FILE *fp;
  char *fname,*directory;

  /*==========================================================================*/
  /* Read in the 1st key argument and malloc filename space if necessary      */
  /*--------------------------------------------------------------------------*/
  /*  1)\num{#} */
  index=1;
  sscanf(dict[index].keyarg,"%d",&int_key_arg);
  tempering_ctrl->npara_temps = int_key_arg;
  gensimopts->ntemper         = int_key_arg;
  npara_temps                 = int_key_arg;
  if(int_key_arg <= 0){keyarg_barf(dict,input_name,fun_key,index);}

  if(npara_temps>1){
    tempering_ctrl->history_name     = (char *)cmalloc(MAXWORD*sizeof(char),"set_sim_prms_tmpr");
    tempering_ctrl->wgt_name         = (char *)cmalloc(MAXWORD*sizeof(char),"set_sim_prms_tmpr");
    tempering_ctrl->troyer_name      = (char *)cmalloc(MAXWORD*sizeof(char),"set_sim_prms_tmpr");
    tempering_ctrl->output_directory = (char *)cmalloc(MAXWORD*sizeof(char),"set_sim_prms_tmpr");
  }/*endif*/

  /*==========================================================================*/
  /* Read in the rest of the key arguments */
  /*-----------------------------------------------------------------------*/
  /*  2)\ens_opt{#} */
  index=2;
  ifound = 0;
  iopt = 0;
  tempering_ctrl->nvt = 0;  tempering_ctrl->npt = 0; tempering_ctrl->nst = 0;
  if(strcasecmp(dict[index].keyarg,"nvt")==0){tempering_ctrl->nvt = 1; iopt=1; ifound++;}
  if(strcasecmp(dict[index].keyarg,"npt")==0){tempering_ctrl->npt = 1; iopt=2; ifound++;}
  if(strcasecmp(dict[index].keyarg,"nst")==0){tempering_ctrl->nst = 1; iopt=3; ifound++;}
  if(ifound!=1){keyarg_barf(dict,input_name,fun_key,index);}
  if(iopt!=1){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Only simple NVT tempering supported\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }/*endif*/
  /*-----------------------------------------------------------------------*/
  /*  5)\rsmpl_opt{} */
  index=5;
  ifound = 0;
  tempering_ctrl->rsmpl_opt = 0;
  if(strcasecmp(dict[index].keyarg,"yes")==0){tempering_ctrl->rsmpl_opt = 1; ifound++;}
  if(strcasecmp(dict[index].keyarg,"no")==0) {tempering_ctrl->rsmpl_opt = 0; ifound++;}
  if(ifound!=1){keyarg_barf(dict,input_name,fun_key,index);}
  /*-----------------------------------------------------------------------*/
  /*  6)\switch_steps{#} */
  index=6;
  sscanf(dict[index].keyarg,"%d",&tempering_ctrl->switch_steps);
  if(tempering_ctrl->switch_steps<1){keyarg_barf(dict,input_name,fun_key,index);}
  /*-----------------------------------------------------------------------*/
  /*  7)\history_fname{#} */
  if(npara_temps>1){
    index=7; ifound = 0;
    check_for_slash(dict[index].keyarg,dict[index].keyword,&ifound);
    if(ifound!=0){keyarg_barf(dict,input_name,fun_key,index);}
    strcpy(tempering_ctrl->history_name,dict[index].keyarg);
  }/*endif*/
  /*-----------------------------------------------------------------------*/
  /*  8)\history_frq{#} */
  index=8;
  sscanf(dict[index].keyarg,"%d",&tempering_ctrl->history_frq);
  if(tempering_ctrl->history_frq<1){keyarg_barf(dict,input_name,fun_key,index);}
  /*-----------------------------------------------------------------------*/
  /*  9)\wgt_fname{#} */
  if(npara_temps>1){
    index=9; ifound = 0;
    check_for_slash(dict[index].keyarg,dict[index].keyword,&ifound);
    if(ifound!=0){keyarg_barf(dict,input_name,fun_key,index);}
    strcpy(tempering_ctrl->wgt_name,dict[index].keyarg);
  }/*endif*/
  /*-----------------------------------------------------------------------*/
  /*  10)\troyer_fname{#} */
  if(npara_temps>1){
    index=10; ifound = 0;
    check_for_slash(dict[index].keyarg,dict[index].keyword,&ifound);
    if(ifound!=0){keyarg_barf(dict,input_name,fun_key,index);}
    strcpy(tempering_ctrl->troyer_name,dict[index].keyarg);
  }/*endif*/
  /*-----------------------------------------------------------------------*/
  /*  11)\restart_par_temp{} */
  index=11;
  ifound = 0;
  if(strcasecmp(dict[index].keyarg,"yes")==0){tempering_ctrl->ipt_restart = 1; ifound++;}
  if(strcasecmp(dict[index].keyarg,"no")==0){tempering_ctrl->ipt_restart = 0; ifound++;}
  if(ifound!=1){keyarg_barf(dict,input_name,fun_key,index);}
  /*-----------------------------------------------------------------------*/
  /*  12)\output_directory{} */
  if(npara_temps>1){
    index=12;
    strcpy(tempering_ctrl->output_directory,dict[index].keyarg);
  }/*endif*/
  /*==========================================================================*/
  /*  If we have more than one temperer:  read in the state point info :      */
  /*                                      malloc some memory                  */
  /*                                      check the directory structure       */

  if(npara_temps>1){
    tempering_ctrl->t_ext = (double *)cmalloc(npara_temps*sizeof(double),"set_sim_prms_tmpr")-1;
    tempering_ctrl->p_ext = (double *)cmalloc(npara_temps*sizeof(double),"set_sim_prms_tmpr")-1;
    tempering_ctrl->s_ext = (double *)cmalloc(npara_temps*sizeof(double),"set_sim_prms_tmpr")-1;
    fname                 = (char   *)cmalloc(MAXWORD*sizeof(char),"set_sim_prms_tmpr");

    strcpy(fname,dict[3].keyarg);
    t_ext     = tempering_ctrl->t_ext;
    p_ext     = tempering_ctrl->p_ext;
    s_ext     = tempering_ctrl->s_ext;
    directory = tempering_ctrl->output_directory;

    PRINTF("     Reading tempering state points from : %s\n",fname);
    PRINTF("      Tempering state points are:\n");

    fp = cfopen((const char *)fname,"r");
    for(i=1;i<=npara_temps;i++){
      switch(iopt){
        case 1: fscanf(fp,"%lg",&t_ext[i]);
                PRINTF("       %g\n",t_ext[i]);
                break;
        case 2: fscanf(fp,"%lg %lg",&t_ext[i],&p_ext[i]);
                PRINTF("       %g %g\n",t_ext[i],p_ext[i]);
                p_ext[i] *= PCONV;
                break;
        case 3: fscanf(fp,"%lg %lg %lg",&t_ext[i],&p_ext[i],&s_ext[i]);
                PRINTF("       %g %g %g\n",t_ext[i],p_ext[i],s_ext[i]);
                p_ext[i] *= PCONV;  s_ext[i] *= STENS_CONV;
                break;
      }/*end switch*/
      readtoendofline_check(fp,fname,i,npara_temps);
    }/*endfor*/
    fclose(fp);
    PRINTF("     Done reading from file: %s\n",fname);
    fflush(stdout);

    for(i=1;i<=npara_temps;i++){
      sprintf (fname, "%s/Temper.%d/ChkDirExistOAtom",directory,i);
      FILE *fpck = fopen(fname,"w");
      if(fpck==NULL){
        sprintf (fname, "%s/Temper.%d",directory,i);
        PRINTF("    @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("    Output directory, %s , is not present\n",fname);
        PRINTF("    @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        EXIT(1);
      }//endif
      fclose(fpck);
    }/*endfor*/

    cfree(fname,"set_sim_params_ctrl_temper"); /* file not needed again */

  }/*endif : We have tempering*/

  /*==========================================================================*/
}/* end routine */
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_sim_params_list(MDINTEGRATE *mdintegrate, MDATOMS *mdatoms,
    MDINTER *mdinter,GENERAL_DATA *general_data,
    MDINTRA *mdintra,
    CP *cp,CLASS_PARSE *class_parse,CP_PARSE *cp_parse,
    FILENAME_PARSE *filename_parse,
    DICT_WORD *dict,char *fun_key)

  /*=======================================================================*/
  /*             Begin routine                                              */
{/*begin routine */
  /*=======================================================================*/
  /*             Local variable declarations                                */

#include "../class_defs/allclass_strip_mdatoms.h"
#include "../class_defs/allclass_strip_mdintegrate.h"
#include "../class_defs/allclass_strip_mdinter.h"
#include "../class_defs/allclass_strip_mdintra.h"

  int ifound = 0;
  int index;
  int ilimit = 50; 
  double real_key_arg;

  /*========================================================================*/
  /* 1)\verlist_pad{#} */
  sscanf(dict[1].keyarg,"%lg",&real_key_arg);
  mdverlist->jver_pad = (int)(real_key_arg);
  index = 1;
  if(mdverlist->jver_pad <= 0)
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if((mdverlist->jver_pad < 10)||  (mdverlist->jver_pad > ilimit) ){ 
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("You have requested a verlist_padding < 10 or > %d \n",ilimit);
    PRINTF("Are you certain this is what you would like to do?   \n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 2)\verlist_mem_safe{#} */
  sscanf(dict[2].keyarg,"%lg",&real_key_arg);
  mdverlist->mem_safe = (real_key_arg); 
  index=2;
  if(mdverlist->mem_safe <= 0) 
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if(mdverlist->mem_safe > 1.5){ 
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("You have requested a verlist_mem_safe > 1.5      \n");
    PRINTF("Are you certain this is what you would like to do?\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 3)\verlist_mem_min{#} */
  sscanf(dict[3].keyarg,"%lg",&real_key_arg);
  mdverlist->nmem_min_lst = (int)(real_key_arg);
  index=3;
  if(mdverlist->nmem_min_lst <= 0) 
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 4)\verlist_skin{#} */
  sscanf(dict[4].keyarg,"%lg",&(mdinteract->skin));
  index=4;
  if(mdinteract->skin < 0.){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  }
  if(mdinteract->skin>2.5){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("You have requested a verlist skin > 2.5A      \n");
    PRINTF("Are you certain this is what you would like to do?\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/
  mdinteract->skin /= BOHR;
  /*-----------------------------------------------------------------------*/ 
  /* 8)\lnkcell_cell_divs{#} */
  sscanf(dict[8].keyarg,"%lg",&real_key_arg);
  mdlnklist->ncell_div_avg = (int)(real_key_arg);
  index=8;
  if(mdlnklist->ncell_div_avg <= 0) 
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if(mdlnklist->ncell_div_avg>40){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("You have requested a lnk cell with > 40s per side \n");
    PRINTF("Are you certain this is what you would like to do?\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 9)\lnk_cell_excl_safe{#} */
  sscanf(dict[9].keyarg,"%lg",&real_key_arg);
  mdlnklist->lnk_excl_safe = (int)(real_key_arg);
  index=9;
  if(mdlnklist->lnk_excl_safe <= 0) 
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 10)\lnk_cell_vol_safe{#} */
  sscanf(dict[10].keyarg,"%lg",&real_key_arg);
  mdlnklist->lnk_vol_safe = (int)(real_key_arg);
  index=10;
  if(mdlnklist->lnk_vol_safe <= 0) 
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if(mdlnklist->lnk_vol_safe > 1.1){ 
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("You have requested a lnk_cell_vol_safe > 1.1  \n");
    PRINTF("Are you certain this is what you would like to do?\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 11)\lnkcell_force_odd{on,off} */
  ifound      = 0;
  if(strcasecmp(dict[11].keyarg,"on")==0) {
    mdlnklist->lnk_for_odd = 1; ifound++;}
  if(strcasecmp(dict[11].keyarg,"off")==0){
    mdlnklist->lnk_for_odd = 0; ifound++;}
  index=11;
  if(ifound == 0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/
  /* 12)\shave_skin_opt{#} */
  ifound      = 0; 
  if(strcasecmp(dict[12].keyarg,"on")==0) {
    mdinteract->ishave_opt = 1; ifound++;}
  if(strcasecmp(dict[12].keyarg,"off")==0){
    mdinteract->ishave_opt = 0; ifound++;}
  index=12;
  if(ifound == 0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /*  13)\neighbor_list{ver_list,lnk_list,no_list} */
  ifound = 0;
  mdnbr_list->iver = 0;
  mdnbr_list->ilnk = 0;
  mdnbr_list->nolst = 0;
  if(strcasecmp(dict[13].keyarg,"no_list")==0)  {
    mdnbr_list->nolst = 1;ifound++;}
  if(strcasecmp(dict[13].keyarg,"lnk_list")==0) {
    mdnbr_list->ilnk = 1;ifound++;}
  if(strcasecmp(dict[13].keyarg,"ver_list")==0) {
    mdnbr_list->iver = 1;ifound++;}
  mdnbr_list->mdlnklist.ilnk = mdnbr_list->ilnk;
  mdnbr_list->mdverlist.iver = mdnbr_list->iver;
  index=13;
  if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/  
  /*  14)\update_type{no_list,lnk_list} */
  ifound  = 0;
  if(strcasecmp(dict[14].keyarg,"no_list")==0){
    mdverlist->nolst_ver_update = 1;
    mdverlist->lnk_ver_update = 0;
    ifound++;
  }
  if(strcasecmp(dict[14].keyarg,"lnk_list")==0) {
    mdverlist->nolst_ver_update = 0;
    mdverlist->lnk_ver_update = 1;
    ifound++;
  }
  mdnbr_list->mdlnklist.lnk_ver_update = mdverlist->lnk_ver_update;
  index=14;
  if(ifound == 0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /*-----------------------------------------------------------------------*/ 
  /* 5)\brnch_root_list_opt{#} */
  ifound      = 0;
  if(strcasecmp(dict[5].keyarg,"on")==0) {
    mdnbr_list->brnch_root_list_opt = 1; ifound++;}
  if(strcasecmp(dict[5].keyarg,"off")==0){
    mdnbr_list->brnch_root_list_opt = 0; ifound++;}
  if(strcasecmp(dict[5].keyarg,"pare_down")==0){
    mdnbr_list->brnch_root_list_opt = 2; ifound++;}
  index=5;
  mdexcl->brnch_root_list_opt = mdnbr_list->brnch_root_list_opt;
  mdnbr_list->mdbrnch_root.brnch_root_list_opt = 
    mdnbr_list->brnch_root_list_opt;
  if(ifound == 0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if(mdnbr_list->iver ==0 && mdnbr_list->brnch_root_list_opt > 0){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Brnch-root list option only implemented for verlet lists.\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 6)\brnch_root_list_skin{#} */
  sscanf(dict[6].keyarg,"%lg",&(mdinteract->brnch_root_skin));
  index=6;
  if(mdinteract->brnch_root_skin < 0.){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  }
  mdinteract->brnch_root_skin /= BOHR;
  if(mdnbr_list->brnch_root_list_opt==0){
    mdinteract->brnch_root_skin = 0.0;
  }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 7)\brnch_root_cutoff{#} */
  sscanf(dict[7].keyarg,"%lg",&(mdinteract->brnch_root_cut));
  index=7;
  if(mdinteract->brnch_root_cut < 0.){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  }
  if(mdinteract->brnch_root_cut > 1.112){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("You have requested a brnch root cutoff > 1.112A \n");
    PRINTF("Are you certain this is what you would like to do?\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/
  mdinteract->brnch_root_cut /= BOHR;

  /*========================================================================*/
}/*end routine*/ 
/*========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_sim_params_cp(MDINTEGRATE *mdintegrate, MDATOMS *mdatoms,
    MDINTER *mdinter,GENERAL_DATA *general_data,
    MDINTRA *mdintra,
    CP *cp,CLASS_PARSE *class_parse,CP_PARSE *cp_parse,
    FILENAME_PARSE *filename_parse,
    DICT_WORD *dict,char *fun_key)

  /*=======================================================================*/
  /*          Begin routine                            */
{/*begin routine */
  /*=======================================================================*/
  /*          Local variable declarations                    */

#include "../class_defs/allclass_strip_mdatoms.h"
#include "../class_defs/allclass_strip_mdintegrate.h"
#include "../class_defs/allclass_strip_mdinter.h"
#include "../class_defs/allclass_strip_mdintra.h"
#include "../class_defs/allclass_strip_cp.h"
#include "../class_defs/allclass_strip_gen.h"

  int iii;
  int ifound,index,itemp,cp_orth_sum,int_key_arg;
  double real_key_arg;

  /*========================================================================*/

  cpcoeffs_info->ncoef_l   = 0;
  cpcoeffs_info->ncoef     = 0;
  cpcoeffs_info->nstate_up = 0;
  cpcoeffs_info->nstate_dn = 0;
  cptherm_info->num_c_nhc  = 0;
  cppseudo->n_ang_max      = 0;

  /*========================================================================*/
  /* Determine dft_typ before you do anything */
  /*-----------------------------------------------------------------------*/  
  /*  7)\cp_dft_typ{lsda,lda,gga_lda,gga_lsda} */
  ifound  = 0;
  cpopts->cp_lda  = 0;
  cpopts->cp_lsda = 0;
  cpopts->cp_gga  = 0;
  if(strcasecmp(dict[7].keyarg,"lsda")==0) {
    cpopts->cp_lsda = 1;ifound++;}
  if(strcasecmp(dict[7].keyarg,"lda")==0)  {
    cpopts->cp_lda  = 1;ifound++;}
  if(strcasecmp(dict[7].keyarg,"gga_lsda")==0){
    cpopts->cp_lsda = 1;cpopts->cp_gga=1; ifound++;}
  if(strcasecmp(dict[7].keyarg,"gga_lda")==0) {
    cpopts->cp_lda  = 1;cpopts->cp_gga=1; ifound++;}
  index=7;
  if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 1)\cp_vxc_typ{pz_lda,pw_lda,pz_lsda}                     */
  ifound = 0;
  index=1;
  if(strcasecmp(dict[1].keyarg,"pz_lda")==0) {ifound++;}
  if(strcasecmp(dict[1].keyarg,"pw_lda")==0) {ifound++;}
  if(strcasecmp(dict[1].keyarg,"pz_lsda")==0){ifound++;}
  sscanf(dict[1].keyarg,"%s",cppseudo->vxc_typ);
  if(ifound == 0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 2)\cp_ggax_typ{becke,pw91x,debug97x,fila_1x,fila_2x,pbe_x, */
  /*            brx89,brx2k,off}                */
  ifound = 0;
  itemp = cpopts->cp_gga;
  cpopts->cp_gga = 0;
  cpopts->cp_laplacian_on = 0;
  cpopts->cp_becke = 0;
  if(strcasecmp(dict[2].keyarg,"off")==0)    {
    ifound++;cpopts->cp_gga=0;}
  if(strcasecmp(dict[2].keyarg,"becke")==0)  {
    ifound++;cpopts->cp_gga=1; cpopts->cp_becke=1;}
  if(strcasecmp(dict[2].keyarg,"pw91x")==0)  {
    ifound++;cpopts->cp_gga=1;}
  if(strcasecmp(dict[2].keyarg,"fila_1x")==0)  {
    ifound++;cpopts->cp_gga=1;}
  if(strcasecmp(dict[2].keyarg,"fila_2x")==0)  {
    ifound++;cpopts->cp_gga=1;}
  if(strcasecmp(dict[2].keyarg,"pbe_x")==0)  {
    ifound++;cpopts->cp_gga=1;}
  if(strcasecmp(dict[2].keyarg,"brx89")==0 ||
      strcasecmp(dict[2].keyarg,"brx2k")==0) {
    ifound++;
    cpopts->cp_gga=1;
    cpopts->cp_ke_dens_on = 1;
    cpopts->cp_tau_functional = 1;
    cpopts->cp_laplacian_on   = 1;
  }/* endif brx89 (a tau functional) */
  if(strcasecmp(dict[2].keyarg,"debug97x")==0)  {
    if(gensimopts->debug_cp != 1 || 
        gensimopts->debug_cp_pimd != 1){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("You may not choose the debug97x functional unless you \n");
      PRINTF("actually doing a debug_cp or debug_cp_pimd\n");
      PRINTF("simulation. It isn't a real functional\n");
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }/* endif */
    ifound++;cpopts->cp_gga=1;
  }/* endif */
  index=2;
  if(ifound == 0){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  }
  sscanf(dict[2].keyarg,"%s",cppseudo->ggax_typ);
  if(itemp!=cpopts->cp_gga){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  }

  if(strcasecmp(cppseudo->vxc_typ,"pz_lda")==0 &&
      strcasecmp(cppseudo->ggax_typ,"pbe_x") == 0) {
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("In order to use the PBE GGA functional, \n");
    PRINTF("You must choose pw_lda as the cp_vxc_typ\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }/* endif */
  /*-----------------------------------------------------------------------*/ 
  /* 3)\cp_ggac_typ{lyp,lypm1,pw91c,pbe_c,tau1_c,off}             */
  ifound = 0;
  itemp = cpopts->cp_gga;
  cpopts->cp_gga = 0;
  cpopts->cp_lyp = 0;
  index=3;
  if(strcasecmp(dict[3].keyarg,"off")==0)    {
    ifound++;cpopts->cp_gga=1;}
  if(strcasecmp(dict[3].keyarg,"pw91c")==0)  {
    ifound++;cpopts->cp_gga=1;}
  if(strcasecmp(dict[3].keyarg,"lyp")==0)    {
    ifound++;cpopts->cp_gga=1; cpopts->cp_lyp=1;}
  if(strcasecmp(dict[3].keyarg,"pbe_c")==0)    {
    ifound++;cpopts->cp_gga=1;}
  if(strcasecmp(dict[3].keyarg,"lypm1")==0)    {
    ifound++;cpopts->cp_gga=1;cpopts->cp_laplacian_on=1;}
  if(strcasecmp(dict[3].keyarg,"tau1_c")==0 ){
    ifound++;
    cpopts->cp_gga=1;
    cpopts->cp_ke_dens_on     = 1;
    cpopts->cp_tau_functional = 1;
    cpopts->cp_laplacian_on   = 1;
  }/* endif tau1 (a tau- and Laplacian-dependent  functional) */
  if(ifound == 0){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  }
  sscanf(dict[3].keyarg,"%s",cppseudo->ggac_typ);
  if(strcasecmp(dict[3].keyarg,"off")==0 && 
      strcasecmp(dict[2].keyarg,"off")==0){
    cpopts->cp_gga=0;
  }
  if(itemp!=cpopts->cp_gga){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  }
  if(strcasecmp(cppseudo->vxc_typ,"pz_lda" )==0 &&
      strcasecmp(cppseudo->ggac_typ,"pbe_c") == 0) {
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("In order to use the PBE GGA functional, \n");
    PRINTF("You must choose pw_lda as the cp_vxc_typ\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }/* endif */
  /*-----------------------------------------------------------------------*/  
  /*  4)\cp_sic{on,off} */
  ifound  = 0;
  if(strcasecmp(dict[4].keyarg,"on")==0)  {cpopts->cp_sic = 1;ifound++;}
  if(strcasecmp(dict[4].keyarg,"off")==0) {cpopts->cp_sic = 0;ifound++;}
  index=4;
  if(ifound == 0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/  
  /*  5)\cp_e_e_interact{on,off} */
  ifound      = 0;
  if(strcasecmp(dict[5].keyarg,"on")==0) {cpopts->cp_nonint=0;ifound++;}
  if(strcasecmp(dict[5].keyarg,"off")==0){cpopts->cp_nonint=1;ifound++;}
  index=5;
  if(ifound == 0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/  
  /*  6)\cp_norb{full_ortho,norm_only,no_constrnt,off} */
  ifound      = 0;
  if(strcasecmp(dict[6].keyarg,"full_ortho")==0) 
  {cpopts->cp_norb = 1;ifound++;}
  if(strcasecmp(dict[6].keyarg,"norm_only")==0) 
  {cpopts->cp_norb = 2;ifound++;}
  if(strcasecmp(dict[6].keyarg,"no_constrnt")==0)
  {cpopts->cp_norb = 3;ifound++;}
  if(strcasecmp(dict[6].keyarg,"off")==0){cpopts->cp_norb = 0;ifound++;}
  index=6;
  if(ifound == 0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 8)\cp_nl_list{on,off} */
  ifound      = 0;
  if(strcasecmp(dict[8].keyarg,"on")==0) {
    cppseudo->nl_cut_on = 1;ifound++;}
  if(strcasecmp(dict[8].keyarg,"off")==0){
    cppseudo->nl_cut_on = 0;ifound++;}
  index=8;
  if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 9)\cp_mass_tau_def{#} */
  sscanf(dict[9].keyarg,"%lg",&(cp_parse->cp_mass_tau_def));
  cpcoeffs_info->tau_mass =  cp_parse->cp_mass_tau_def;
  index=9;
  if(cp_parse->cp_mass_tau_def<=0.)
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 10)\cp_mass_cut_def{#} */
  sscanf(dict[10].keyarg,"%lg",&(cp_parse->cp_mass_cut_def));
  cpcoeffs_info->ecut_mass =  cp_parse->cp_mass_cut_def;
  index=10;
  if(cp_parse->cp_mass_cut_def<=0.)
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 11)\cp_energy_cut_def{#} */
  sscanf(dict[11].keyarg,"%lg",&(cp_parse->cp_ecut_def));
  index=11;
  if(cp_parse->cp_ecut_def <= 0. ) 
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 12)\cp_fict_KE{#} */
  sscanf(dict[12].keyarg,"%lg",&(cpopts->te_ext));
  index=12;
  if(cpopts->te_ext <= 0. ) 
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if((cpopts->te_ext > 3000.0) || (cpopts->te_ext <100.0)) {
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("The recommended fictious cp_fict_KE range is 100K - 3000K\n"); 
    PRINTF("per 4.1x10^6 coefs\n"); 
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
  }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /*  13)\cp_ptens{on,off} */
  ifound      = 0;
  if(strcasecmp(dict[13].keyarg,"on")==0)
  {cpopts->cp_ptens_calc=1;ifound++;}
  if(strcasecmp(dict[13].keyarg,"off")==0)
  {cpopts->cp_ptens_calc=0;ifound++;}
  cppseudo->cp_ptens_calc = cpopts->cp_ptens_calc;
  index=13;
  if(ifound == 0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/  
  /* 14)\cp_init_orthog{on,off} */
  ifound  = 0;
  if(strcasecmp(dict[14].keyarg,"on")==0)  
  {cpopts->cp_init_orthog = 1;ifound++;}
  if(strcasecmp(dict[14].keyarg,"off")==0) 
  {cpopts->cp_init_orthog = 0;ifound++;}
  index=14;
  if(ifound == 0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 15)\cp_cg_line_min_len{#} */
  sscanf(dict[15].keyarg,"%lg",&real_key_arg);
  genminopts->cp_cg_line_min_len = (int)(real_key_arg);
  index=15;
  if(genminopts->cp_cg_line_min_len<0)
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 16)\cp_minimize_typ{min_std,min_cg,min_diis} */
  ifound = 0;
  cpopts->cp_hess_calc=0;
  genminopts->cp_min_std = 0;
  genminopts->cp_min_cg = 0;
  genminopts->cp_min_diis = 0;
  if(strcasecmp(dict[16].keyarg,"min_std")==0) {
    genminopts->cp_min_std  = 1; ifound++;}
  if(strcasecmp(dict[16].keyarg,"min_cg")==0) {
    genminopts->cp_min_cg   = 1; ifound++;}
  if(strcasecmp(dict[16].keyarg,"min_diis")==0){
    genminopts->cp_min_diis = 1;
    cpopts->cp_hess_calc=1; ifound++;}
  index=16;
  if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 17)\cp_diis_hist_len{#} */
  sscanf(dict[17].keyarg,"%lg",&real_key_arg);
  genminopts->cp_diis_hist_len = (int)(real_key_arg);
  index=17;
  if(genminopts->cp_diis_hist_len<0)
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 18)\cp_orth_meth{gram_schmidt,lowdin,normalize,none} */
  ifound  = 0;
  cpopts->cp_gs = 0;
  cpopts->cp_low = 0;
  cpopts->cp_normalize = 0;
  if(strcasecmp(dict[18].keyarg,"lowdin")==0) {
    cpopts->cp_low = 1;ifound++;}
  if(strcasecmp(dict[18].keyarg,"gram_schmidt")==0) {
    cpopts->cp_gs=1;ifound++;}
  if(strcasecmp(dict[18].keyarg,"normalize")==0) {
    cpopts->cp_normalize=1;ifound++;}
  if(strcasecmp(dict[18].keyarg,"none")==0) {
    ifound++;}
  index=18;
  if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if(gensimopts->cp_min == 1 || 
      gensimopts->cp_wave_min == 1 ||
      gensimopts->cp_wave_min_pimd == 1 || 
      cpopts->cp_init_orthog == 1){
    if(cpopts->cp_normalize == 1 && cpopts->cp_norb !=2){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("You may only use the normalization option under\n");
      PRINTF("norb minimization or md with the norm_only option\n");
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }/* endif */
    cp_orth_sum = cpopts->cp_gs + cpopts->cp_low + cpopts->cp_normalize;
    if(cp_orth_sum == 0 && cpopts->cp_norb !=3){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("You may only use the none option under\n");
      PRINTF("norb minimization or md with the none option\n");
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }/* endif */
  }/* endif general_data */
  /*-----------------------------------------------------------------------*/ 
  /* 19)\cp_restart_type{initial,restart_pos,restart_posvel,restart_all}*/
  ifound = 0;
  if(strcasecmp(dict[19].keyarg,"gen_wave")==0)    {
    cp_parse->istart_cp = 0; ifound++;}
  if(strcasecmp(dict[19].keyarg,"initial")==0)    {
    cp_parse->istart_cp = 1; ifound++;}
  if(strcasecmp(dict[19].keyarg,"restart_pos")==0)  {
    cp_parse->istart_cp = 2; ifound++;}
  if(strcasecmp(dict[19].keyarg,"restart_posvel")==0){
    cp_parse->istart_cp= 3; ifound++;}
  if(strcasecmp(dict[19].keyarg,"restart_all")==0)  {
    cp_parse->istart_cp = 4; ifound++;}
  index=19;
  if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  gensimopts->istart_cp = cp_parse->istart_cp;
  /*-----------------------------------------------------------------------*/ 
  /* 20)\diis_hist_len{#} */
  sscanf(dict[20].keyarg,"%lg",&real_key_arg);
  genminopts->diis_hist_len = (int)(real_key_arg);
  index=20;
  if(genminopts->diis_hist_len<0)
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 21)\nlvps_list_skin{#}   */
  sscanf(dict[21].keyarg,"%lg",&real_key_arg);
  cppseudo->nlvps_skin = (real_key_arg);
  index=21;
  if(cppseudo->nlvps_skin<0){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);}
  cppseudo->nlvps_skin /= BOHR;
  /*-----------------------------------------------------------------------*/ 
  /* 22)\gradient_cutoff{#}   */
  sscanf(dict[22].keyarg,"%lg",&real_key_arg);
  cppseudo->gga_cut = (real_key_arg);
  index=22;
  if(cppseudo->gga_cut<0){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);}
  /*-----------------------------------------------------------------------*/
  /* 23)\zero_cp_vel{initial,periodic,no}*/
  ifound = 0;
  if(strcasecmp(dict[23].keyarg,"initial")==0)    {
    cpopts->zero_cp_vel = 1; ifound++;}
  if(strcasecmp(dict[23].keyarg,"periodic")==0)    {
    cpopts->zero_cp_vel = 2; ifound++;}
  if(strcasecmp(dict[23].keyarg,"no")==0)  {
    cpopts->zero_cp_vel = 0; ifound++;}
  index=23;
  if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 24)\cp_check_perd_size{#}   */
  ifound = 0;
  if(strcasecmp(dict[24].keyarg,"on")==0)    {
    cpopts->icheck_perd_size = 1; ifound++;}
  if(strcasecmp(dict[24].keyarg,"off")==0)    {
    cpopts->icheck_perd_size = 0; ifound++;}
  index=24;
  if(ifound != 1){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);}
  /*-----------------------------------------------------------------------*/ 
  /* 25)\cp_tol_edge_dist{#}   */
  sscanf(dict[25].keyarg,"%lg",&real_key_arg);
  cpopts->tol_edge_dist = (real_key_arg);
  index=25;
  if(cpopts->tol_edge_dist <= 0){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);}

  if( (cpopts->tol_edge_dist < 3.0) || 
      (cpopts->tol_edge_dist > 6.0) ) {
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("The recommended value of the CP cluster size tolerance\n");
    PRINTF("is 3-6 A\n"); 
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
  }/*endif*/
  cpopts->tol_edge_dist /= BOHR;
  /*-----------------------------------------------------------------------*/ 
  /* 26)\cp_dual_grid_opt{#} */
  ifound = 0;
  if(strcasecmp(dict[26].keyarg,"off")==0)    {
    cpopts->cp_dual_grid_opt = 0; ifound++;}
  if(strcasecmp(dict[26].keyarg,"prop")==0)    {
    cpopts->cp_dual_grid_opt = 1; ifound++;}
  if(strcasecmp(dict[26].keyarg,"not_prop")==0)    {
    cpopts->cp_dual_grid_opt = 2; ifound++;}
  index=26;
  if(ifound != 1){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);}
  /*-----------------------------------------------------------------------*/ 
  /* 27)\cp_check_perd_size{#}   */
  ifound = 0;
  if(strcasecmp(dict[27].keyarg,"on")==0)    {
    cpopts->icheck_dual_size = 1; ifound++;}
  if(strcasecmp(dict[27].keyarg,"off")==0)    {
    cpopts->icheck_dual_size = 0; ifound++;}
  index=27;
  if(ifound != 1){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);}
  /*-----------------------------------------------------------------------*/ 
  /* 28)\cp_move_dual_box_opt{#}   */
  ifound = 0;
  if(strcasecmp(dict[28].keyarg,"on")==0)    {
    gencell->imov_cp_box = 1; ifound++;}
  if(strcasecmp(dict[28].keyarg,"off")==0)    {
    gencell->imov_cp_box = 0; ifound++;}
  index=28;
  if(ifound != 1){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);}
  /*-----------------------------------------------------------------------*/ 
  /* 29)\cp_energy_cut_dual_grid_def{#} */
  sscanf(dict[29].keyarg,"%lg",&(cp_parse->cp_ecut_dual_grid_def));
  index=29;
  if(cp_parse->cp_ecut_dual_grid_def <= 0. ) 
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 30)\cp_alpha_conv_dual{#} */
  sscanf(dict[30].keyarg,"%lg",&(cppseudo->alpha_conv_dual));
  index=30;
  if(cppseudo->alpha_conv_dual <= 0. ) 
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if( cppseudo->alpha_conv_dual < 7.0 ||
      cppseudo->alpha_conv_dual > 10.0){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("The recommended range of cp_alpha_conv_dual is 7-10 \n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
  }
  /*-----------------------------------------------------------------------*/ 
  /* 31)\interp_pme_dual{#} */
  sscanf(dict[31].keyarg,"%lg",&real_key_arg);
  cppseudo->n_interp_pme_dual= (int)(real_key_arg);
  index=31;
  if((cppseudo->n_interp_pme_dual  <4 )||
      ((cppseudo->n_interp_pme_dual  %2)!=0))
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);

  /*-----------------------------------------------------------------------*/ 
  /* 32)\cp_elf_calc_frq{#} */
  sscanf(dict[32].keyarg,"%lg",&real_key_arg);
  cpopts->cp_elf_calc_frq =  (int)(real_key_arg);
  index=32;
  if((cpopts->cp_elf_calc_frq  <0 ))
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if(cpopts->cp_elf_calc_frq > 0) {cpopts->cp_ke_dens_on = 1;}

  /*-----------------------------------------------------------------------*/ 
  /* 33)\cp_ngrid_skip{#} */
  sscanf(dict[33].keyarg,"%lg",&real_key_arg);
  cpopts->cp_ngrid_skip =  (int)(real_key_arg);
  index=33;
  if((cpopts->cp_ngrid_skip  < 0 ))
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);

  /*-----------------------------------------------------------------------*/ 
  /* 34)\cp_isok_opt{#}   */
  ifound = 0;
  if(strcasecmp(dict[34].keyarg,"on")==0)    {
    cpopts->cp_isok_opt = 1; ifound++;}
  if(strcasecmp(dict[34].keyarg,"off")==0)    {
    cpopts->cp_isok_opt = 0; ifound++;}
  index=34;
  if(ifound != 1){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);}

  /*-----------------------------------------------------------------------*/ 
  /* 35)\cp_nonloc_ees_opt{on,off} */
  index=35;
  ifound = 0;
  if(strcasecmp(dict[index].keyarg,"on")==0)    {
    cppseudo->nonlocal.ees_on = 1; ifound++;}
  if(strcasecmp(dict[index].keyarg,"off")==0)    {
    cppseudo->nonlocal.ees_on = 0; ifound++;}
  if(ifound != 1){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);}
  cppseudo->ees_nonloc_on = cppseudo->nonlocal.ees_on;
  /*-----------------------------------------------------------------------*/ 
  /* 36)\cp_pseudo_ees_order{4} */
  index=36;
  sscanf(dict[index].keyarg,"%lg",&real_key_arg);
  int_key_arg = (int) real_key_arg;
  cppseudo->n_interp_ps          = int_key_arg;
  cppseudo->nonlocal.n_interp  = int_key_arg;
  cppseudo->nonlocal.n_interp2 = (int_key_arg*int_key_arg);
  if( ((int_key_arg % 2)!=0) || (int_key_arg < 4) ){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  }/*endif*/
  if(int_key_arg>8){
    printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$$$$$$\n");
    printf("The recommended range of cp_nonloc_pme_order is 4-8 \n");
    printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
  }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 37)\cp_pseudo_pme_scale{1.2} */
  index=37;
  sscanf(dict[index].keyarg,"%lg",&real_key_arg);
  cppseudo->fft_size_scale_ps         = real_key_arg;
  cppseudo->nonlocal.fft_size_scale = real_key_arg;
  if(real_key_arg <= 1.1){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  }/*endif*/
  if(real_key_arg>1.4){
    printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$$$$$$\n");
    printf("The recommended range of cp_nonloc_pme_scale is 1.1-1.4 \n");
    printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 38)\cp_eext_ees_opt{on,off} */
  index=38;
  ifound = 0;
  if(strcasecmp(dict[index].keyarg,"on")==0)    {
    cppseudo->nonlocal.ees_eext_on = 1; ifound++;}
  if(strcasecmp(dict[index].keyarg,"off")==0)    {
    cppseudo->nonlocal.ees_eext_on = 0; ifound++;}
  if(ifound != 1){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);}
  cppseudo->ees_eext_on = cppseudo->nonlocal.ees_eext_on;
  /*-----------------------------------------------------------------------*/ 
  /* 39)\cp_min_update{on,off} */
  index=39;
  ifound = 0;
  if(strcasecmp(dict[index].keyarg,"on")==0)    {
    gensimopts->cp_min_update = 1; ifound++;}
  if(strcasecmp(dict[index].keyarg,"off")==0)    {
    gensimopts->cp_min_update = 0; ifound++;}
  if(ifound != 1){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);}
  if(gensimopts->cp_min_update==0){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
    PRINTF("You are in danger! Your cp_min_update is off! \n");
    PRINTF("Are you debugging some CS code? \n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n\n");
  }//

  /*-----------------------------------------------------------------------*/ 
  /* 40)\cp_force_complex_psi{on,off} */
  index=40;
  ifound = 0;
  if(strcasecmp(dict[index].keyarg,"on")==0)    {
    cpopts->cp_force_complex_psi = 1; ifound++;}
  if(strcasecmp(dict[index].keyarg,"off")==0)    {
    cpopts->cp_force_complex_psi = 0; ifound++;}
  if(ifound != 1){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);}
  if(cpopts->cp_force_complex_psi==1){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
    PRINTF("You have forced the KS states to be complex.\n");
    PRINTF("If you are at the Gamma point, the code will be slow.\n");
    PRINTF("Are you debugging?\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n\n");
  }//endif
  /*-----------------------------------------------------------------------*/ 
  /* 41)\cp_allow_duplicate_kpts{on,off} */
  index=41;
  ifound = 0;
  if(strcasecmp(dict[index].keyarg,"on")==0)    {
    cpopts->cp_allow_dup_kpts = 1; ifound++;}
  if(strcasecmp(dict[index].keyarg,"off")==0)    {
    cpopts->cp_allow_dup_kpts = 0; ifound++;}
  if(ifound != 1){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);}
  if(cpopts->cp_allow_dup_kpts==1){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
    PRINTF("You have permitted the use of duplicate kpoints.\n");
    PRINTF("Are you debugging?\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n\n");
  }//endif
  /*-----------------------------------------------------------------------*/ 
  /* 42)\cp_bomd_max_minimization_steps{#} */
  index=42;
  sscanf(dict[index].keyarg,"%lg",&real_key_arg);
  gentimeinfo->btime = (int)real_key_arg;
  // TODO: Sanity checks here for <= 0.0 > 1000000
  /*-----------------------------------------------------------------------*/ 
  /* 43)\cp_bomd_timestep_scale{#} */
  index=43;
  sscanf(dict[index].keyarg,"%lg",&real_key_arg);
  gentimeinfo->bomd_scale = real_key_arg;
  // TODO: Sanity checks here for <= 0.0 > 1.0
  /*========================================================================*/
  // CP mass warning messages for the uninitiated

  PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  double cmass_tau_def = cp_parse->cp_mass_tau_def;
  double cmass_cut_def = cp_parse->cp_mass_cut_def;
  double tau_true      = cmass_tau_def*sqrt(cmass_cut_def/CP_EMAGIC);
  PRINTF("The true cp_mass_tau is %g. The value set %g\n",tau_true,cmass_tau_def);
  PRINTF("is scaled by sqrt(cmass_cut_def/%g)=%g.\n",1.0/(2.0*CP_EMAGIC),
      sqrt(cmass_cut_def/(2.0*CP_EMAGIC)));
  PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n\n");

  /*========================================================================*/
}/*end routine*/ 
/*========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_sim_params_gen(MDINTEGRATE *mdintegrate, MDATOMS *mdatoms,
    MDINTER *mdinter,GENERAL_DATA *general_data,
    MDINTRA *mdintra,
    CP *cp,CLASS_PARSE *class_parse,CP_PARSE *cp_parse,
    FILENAME_PARSE *filename_parse,
    DICT_WORD *dict,char *fun_key)

  /*=======================================================================*/
  /*          Begin routine                            */
{/*begin routine */
  /*=======================================================================*/
  /*          Local variable declarations                    */

#include "../class_defs/allclass_strip_mdatoms.h"
#include "../class_defs/allclass_strip_mdintegrate.h"
#include "../class_defs/allclass_strip_mdinter.h"
#include "../class_defs/allclass_strip_mdintra.h"
#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_cp.h"

  int ifound,index;
  double real_key_arg;
  int   int_key_arg;

  /*========================================================================*/
  /*-----------------------------------------------------------------------*/ 
  /*  1)\simulation_typ{md,minimize,pimd,cp,cp_wave,cp_min,cp_wave_min,
      cp_wave_min_pimd,cp_pimd,cp_wave_pimd,cp_bomd
      debug,debug_pimd,debug_cp} */

  ifound =0;

  gensimopts->minimize = 0;
  gensimopts->md = 0;
  gensimopts->pimd = 0;
  gensimopts->cp = 0;
  gensimopts->cp_wave = 0;
  gensimopts->cp_wave_pimd = 0;
  gensimopts->cp_pimd = 0; 
  gensimopts->cp_min = 0; 
  gensimopts->cp_bomd = 0;
  gensimopts->cp_wave_min = 0;
  gensimopts->cp_wave_min_pimd = 0;
  gensimopts->debug = 0;
  gensimopts->debug_pimd = 0;
  gensimopts->debug_cp = 0;
  gensimopts->debug_cp_pimd = 0; 

  if(strcasecmp(dict[1].keyarg,"minimize")==0)  {
    gensimopts->minimize = 1; ifound++;}
  if(strcasecmp(dict[1].keyarg,"md")==0){
    gensimopts->md = 1; ifound++;}
  if(strcasecmp(dict[1].keyarg,"pimd")==0){
    gensimopts->pimd = 1; ifound++;}
  if(strcasecmp(dict[1].keyarg,"cp")==0){
    gensimopts->cp = 1; ifound++;}
  if(strcasecmp(dict[1].keyarg,"cp_pimd")==0){
    gensimopts->cp_pimd = 1; ifound++;}
  if(strcasecmp(dict[1].keyarg,"cp_wave")==0){
    gensimopts->cp_wave = 1; ifound++;}
  if(strcasecmp(dict[1].keyarg,"cp_wave_pimd")==0){
    gensimopts->cp_wave_pimd = 1; ifound++;}
  if(strcasecmp(dict[1].keyarg,"cp_bomd")==0){
    gensimopts->cp_bomd = 1; ifound++;}
  if(strcasecmp(dict[1].keyarg,"cp_min")==0){
    gensimopts->cp_min = 1; ifound++;}
  if(strcasecmp(dict[1].keyarg,"cp_wave_min")==0){
    gensimopts->cp_wave_min=1;ifound++;}
  if(strcasecmp(dict[1].keyarg,"cp_wave_min_pimd")==0){
    gensimopts->cp_wave_min_pimd=1;ifound++;}
  if(strcasecmp(dict[1].keyarg,"debug")==0) {
    gensimopts->debug = 1; ifound++;}
  if(strcasecmp(dict[1].keyarg,"debug_cp")==0) {
    gensimopts->debug_cp = 1; ifound++;}
  if(strcasecmp(dict[1].keyarg,"debug_cp_pimd")==0) {
    gensimopts->debug_cp_pimd = 1; ifound++;}
  if(strcasecmp(dict[1].keyarg,"debug_pimd")==0) {
    gensimopts->debug_pimd = 1; ifound++;}
  index=1;
  if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if(gensimopts->cp_bomd == 1) gensimopts->cp_wave_min = 1;

  /*-----------------------------------------------------------------------*/ 
  /*  2)\ensemble_typ{nve,nvt,npt_i,npt_f,nst} */
  ifound = 0;
  genensopts->nve = 0;
  genensopts->nvt = 0;
  genensopts->npt_i = 0;
  genensopts->npt_f = 0;
  genensopts->nst = 0;

  if(strcasecmp(dict[2].keyarg,"nve")==0)  {
    genensopts->nve = 1;ifound++;}
  if(strcasecmp(dict[2].keyarg,"nvt")==0)  {
    genensopts->nvt = 1;ifound++;}
  if(strcasecmp(dict[2].keyarg,"npt_i")==0){
    genensopts->npt_i = 1;
    mdintegrate->mdbaro.iopt = 1; ifound++;}
  if(strcasecmp(dict[2].keyarg,"npt_f")==0){
    genensopts->npt_f = 1;ifound++;}

  mdintegrate->mdpar_rahman.npt_f = 
    genensopts->npt_f;
  index=2;
  if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 3)\num_time_step{#} */
  sscanf(dict[3].keyarg,"%lg",&real_key_arg);
  gentimeinfo->ntime = (int)(real_key_arg);
  index=3;
  if(gentimeinfo->ntime<0)
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if(gentimeinfo->ntime>10000000){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("You have requested over 10 million time steps    \n");
    PRINTF("Are you certain this is what you would like to do?\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/

  /*-----------------------------------------------------------------------*/ 
  /* 4)\time_step{#} */
  sscanf(dict[4].keyarg,"%lg",&real_key_arg);
  gentimeinfo->dt = real_key_arg;
  index=4;
  if(gentimeinfo->dt==0.0)
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if(gentimeinfo->dt>50.0){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("You have requested a time step greater than 50fs\n");
    PRINTF("Are you certain this is what you would like to do?\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/
  gentimeinfo->dt /= TIME_CONV;

  /*-----------------------------------------------------------------------*/ 
  /* 5)\temperature{#} */
  sscanf(dict[5].keyarg,"%lg",&real_key_arg);
  genstatepoint->t_ext = real_key_arg;
  index=5;
  if(genstatepoint->t_ext<0.0)
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if(genstatepoint->t_ext>1000.0){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("You have requested a temperature greater than 1000K\n");
    PRINTF("Are you certain this is what you would like to do?\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/

  /*-----------------------------------------------------------------------*/ 
  /* 6)\pressure{#} */
  sscanf(dict[6].keyarg,"%lg",&real_key_arg);
  genstatepoint->pext = real_key_arg;
  index=6;
  if(genstatepoint->pext<0.0)
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if(genstatepoint->pext>100000.0){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("You have requested a pressure greater than 0.1Mbar\n");
    PRINTF("Are you certain this is what you would like to do?\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/
  genstatepoint->pext *= PCONV;
  /*-----------------------------------------------------------------------*/  
  /* 7)\restart_type{initial,restart_pos,restart_posvel,restart_all}*/
  ifound = 0;
  if(strcasecmp(dict[7].keyarg,"initial")==0)    {
    class_parse->istart = 1; ifound++;}
  if(strcasecmp(dict[7].keyarg,"restart_pos")==0)  {
    class_parse->istart = 2; ifound++;}
  if(strcasecmp(dict[7].keyarg,"restart_posvel")==0){
    class_parse->istart= 3; ifound++;}
  if(strcasecmp(dict[7].keyarg,"restart_all")==0)  {
    class_parse->istart = 4; ifound++;}
  gensimopts->istart = class_parse->istart;
  index=7;
  if(ifound != 1)keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/  
  /* 8)\minimize_typ{min_std,min_cg,min_diis} */
  ifound = 0;
  genminopts->min_std = 0;
  genminopts->min_cg = 0;
  genminopts->min_diis = 0;
  if(strcasecmp(dict[8].keyarg,"min_std")==0) {
    genminopts->min_std  = 1; ifound++;}
  if(strcasecmp(dict[8].keyarg,"min_cg")==0) {
    genminopts->min_cg   = 1; ifound++;}
  if(strcasecmp(dict[8].keyarg,"min_diis")==0){
    genminopts->min_diis = 1; ifound++;}
  index=8;
  if(ifound != 1)keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 9)\annealing_rate{#} */
  sscanf(dict[9].keyarg,"%lg",&real_key_arg);
  gensimopts->ann_rate = (double)(real_key_arg);
  index=9;
  if(gensimopts->ann_rate<=0)
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 14)\rndm_seed{#} */
  sscanf(dict[14].keyarg,"%d",&(int_key_arg));
  mdvel_samp->iseed = (long) int_key_arg;
  cpvel_samp->iseed = (long) int_key_arg;
  mdvel_samp->qseed = (double) mdvel_samp->iseed;
  cpvel_samp->qseed = (double) cpvel_samp->iseed;
  /*-----------------------------------------------------------------------*/ 
  /* 15)\rndm_seed2{#} */
  sscanf(dict[15].keyarg,"%ld",&(mdvel_samp->iseed2));
  sscanf(dict[15].keyarg,"%ld",&(cpvel_samp->iseed2));
  /*-----------------------------------------------------------------------*/ 
  /* 16)\generic_fft_opt{on,off} */
  ifound = 0;
  if(strcasecmp(dict[16].keyarg,"fftw")==0) {
    gensimopts->fftopt  = 0; ifound++;}
  if(strcasecmp(dict[16].keyarg,"essl")==0) {
    gensimopts->fftopt  = 1; ifound++;}
  index=16;
  if(ifound != 1)keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 17)\alpha_clus{#}   */
  sscanf(dict[17].keyarg,"%lg",&real_key_arg);
  genewald->alp_clus = (real_key_arg);
  index=17;
  if(genewald->alp_clus<0){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);}
  if( (genewald->alp_clus < 7.0) || 
      (genewald->alp_clus > 22.0) ) {
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("The recommended value of the dimensionless parameter\n");
    PRINTF("alp_clus is 8-22\n"); 
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
  }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 18)\ecut_clus{#}   */
  sscanf(dict[18].keyarg,"%lg",&real_key_arg);
  genewald->ecut_clus = (real_key_arg);
  index=18;
  if(genewald->ecut_clus < 0){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);}
  /*-----------------------------------------------------------------------*/
  /* 19)\surf_tens{#} */
  sscanf(dict[19].keyarg,"%lg",&real_key_arg);
  genstatepoint->stens_ext= real_key_arg;
  index=19;
  if(genstatepoint->stens_ext>1000.0){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
    PRINTF("You have requested a surface tension greater than 1000 N/m\n");
    PRINTF("Are you certain this is what you would like to do?\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/
  genstatepoint->stens_ext*= STENS_CONV;

  /*-----------------------------------------------------------------------*/
  /* 20)\min_num_atoms_per_proc{#}   */

  /*-----------------------------------------------------------------------*/
  /* 23)\annealing_opt{on,off} */

  ifound = 0;
  if(strcasecmp(dict[23].keyarg,"on")==0)  {
    gensimopts->anneal_opt = 1;        
    ifound++;
  }
  if(strcasecmp(dict[23].keyarg,"off")==0)  {
    gensimopts->anneal_opt = 0;        
    ifound++;
  }
  index=23;
  if(ifound != 1){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  }

  if((gensimopts->anneal_opt == 1) && 
      (gensimopts->ann_rate == 1.0)){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Dude, you selected the annealing option but have an \n");
    PRINTF("annealing rate of 1.0.  Maybe you want to reconsider \n");
    PRINTF("something more sensible for the latter\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }/*endif*/

  /*-----------------------------------------------------------------------*/
  /* 24)\ann_start_temperature{#} */

  sscanf(dict[24].keyarg,"%lg",&real_key_arg);
  gensimopts->ann_start_temp = real_key_arg;
  index=24;
  if(gensimopts->ann_start_temp < 0.0){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);}


  /*-----------------------------------------------------------------------*/
  /* 25)\ann_target_temperature{#} */

  sscanf(dict[25].keyarg,"%lg",&real_key_arg);
  gensimopts->ann_target_temp = real_key_arg;
  index=25;
  if(gensimopts->ann_start_temp < 0.0){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);}

  /*========================================================================*/
}/*end routine*/ 
/*========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_sim_params_vpot(MDINTEGRATE *mdintegrate, MDATOMS *mdatoms,
    MDINTER *mdinter,GENERAL_DATA *general_data,
    MDINTRA *mdintra,
    CP *cp,CLASS_PARSE *class_parse,CP_PARSE *cp_parse,
    FILENAME_PARSE *filename_parse,
    DICT_WORD *dict,char *fun_key)

  /*=======================================================================*/
  /*          Begin routine                            */
{/*begin routine */
  /*=======================================================================*/
  /*          Local variable declarations                    */

#include "../class_defs/allclass_strip_mdatoms.h"
#include "../class_defs/allclass_strip_mdintegrate.h"
#include "../class_defs/allclass_strip_mdinter.h"
#include "../class_defs/allclass_strip_mdintra.h"
#include "../class_defs/allclass_strip_cp.h"
#include "../class_defs/allclass_strip_gen.h"

  int ifound,index;
  int    int_key_arg;
  double real_key_arg;

  /*========================================================================*/
  /* 1)\shift_inter_pe{on,off} */
  ifound      = 0;
  mdenergy_ctrl->iswit_vdw=0;
  if(strcasecmp(dict[1].keyarg,"swit")==0) {
    class_parse->ishift_pot = 2;mdenergy_ctrl->iswit_vdw=1;ifound++;}
  if(strcasecmp(dict[1].keyarg,"on")==0) {
    class_parse->ishift_pot = 1;ifound++;}
  if(strcasecmp(dict[1].keyarg,"off")==0){
    class_parse->ishift_pot = 0;ifound++;}
  mdinteract->iswit_vdw=mdenergy_ctrl->iswit_vdw;
  genmdstat_avg->iswit_vdw = mdenergy_ctrl->iswit_vdw;
  index=1;
  if(ifound != 1){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  }
  if(class_parse->ishift_pot==0){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("Potential shift is off.  Energy conservation will degrade.\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/
  if(class_parse->ishift_pot==1){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("Potential shift is on.  Potential energy will change.\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 2)\inter_spline_pts{#} */
  sscanf(dict[2].keyarg,"%lg",&real_key_arg);
  mdinteract->nsplin  = (int)(real_key_arg);
  mdecor->nsplin     = (int)(real_key_arg);
  mdecor->nsplin_m2  = mdecor->nsplin-2;
  index=2;
  if(mdinteract->nsplin <= 0) 
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if(mdinteract->nsplin > 2000 ){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("You have requested more than 1000 spline points   \n");
    PRINTF("Are you certain this is what you would like to do?\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 3)\intra_block_min{} */
  /*-----------------------------------------------------------------------*/ 
  /* 4)\pten_inter_respa{#} */
  index=4;
  sscanf(dict[4].keyarg,"%lg",&real_key_arg);
  mdinteract->pten_inter_guess = real_key_arg*PCONV;
  /*-----------------------------------------------------------------------*/ 
  /* 5)\pten_kin_respa{#} */
  index=5;
  sscanf(dict[5].keyarg,"%lg",&real_key_arg);
  mdinteract->pten_kin_guess = real_key_arg*PCONV;
  /*-----------------------------------------------------------------------*/ 
  /* 6)\pseud_spline_pts{#} */
  sscanf(dict[6].keyarg,"%lg",&real_key_arg);
  cppseudo->nsplin_g = (int)(real_key_arg);
  index=6;
  if(cppseudo->nsplin_g <= 0) 
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 13)\ewald_interp_pme{#} */
  sscanf(dict[13].keyarg,"%lg",&real_key_arg);
  mdpart_mesh->n_interp= (int)(real_key_arg);
  cpatom_pme->n_interp = mdpart_mesh->n_interp;
  index=13;
  if((mdpart_mesh->n_interp  <4 )||
      ((mdpart_mesh->n_interp  %2)!=0))
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 7)\scratch_length{#} */
  /*-----------------------------------------------------------------------*/ 
  /* 8)\ewald_alpha{#} */
  sscanf(dict[8].keyarg,"%lg",&(genewald->alp_ewd));
  index=8;
  if(genewald->alp_ewd <= 0.)
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if(genewald->alp_ewd>40){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("You have requested an alpha ewald parameter > 40\n");
    PRINTF("Are you certain this is what you would like to do?\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 9)\ewald_kmax{#} */
  sscanf(dict[9].keyarg,"%lg",&real_key_arg);
  class_parse->kmax_ewd = (int)(real_key_arg);
  index=9;
  if(class_parse->kmax_ewd <0) 
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if(class_parse->kmax_ewd>40){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("You have requested a kmax ewald parameter > 40\n");
    PRINTF("Are you certain this is what you would like to do?\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/
  if((double)(class_parse->kmax_ewd)<genewald->alp_ewd){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("You have requested a kmax ewald parameter < alpha ewald\n");
    PRINTF("Are you certain this is what you would like to do?\n");
    PRINTF("For well converged forces kmax ewald >= alpha ewald \n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 10)\ewald_respa_kmax{#} */
  sscanf(dict[10].keyarg,"%lg",&real_key_arg);
  class_parse->kmax_res = (int)(real_key_arg);
  index=10;
  if((class_parse->kmax_res < 0)||
      (class_parse->kmax_res>class_parse->kmax_ewd)){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);}
  /*-----------------------------------------------------------------------*/ 
  /* 11)\ewald_pme_opt{#} */
  ifound      = 0;
  if(strcasecmp(dict[11].keyarg,"on")==0) {
    mdpart_mesh->pme_on = 1; ifound++;}
  if(strcasecmp(dict[11].keyarg,"off")==0){
    mdpart_mesh->pme_on = 0; ifound++;}
  cpatom_pme->pme_on = mdpart_mesh->pme_on;
  index=11;
  if(ifound == 0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 12)\ewald_kmax_pme{#} */
  sscanf(dict[12].keyarg,"%lg",&real_key_arg);
  mdpart_mesh->kmax_pme = (int)(real_key_arg);
  index=12;
  if(mdpart_mesh->kmax_pme  <0) 
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if(mdpart_mesh->pme_on == 1){
    if(mdpart_mesh->kmax_pme >40){
      PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
      PRINTF("You have requested a PME kmax ewald parameter > 40\n");
      PRINTF("Are you certain this is what you would like to do?\n");
      PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
    }/*endif*/
    if(class_parse->kmax_ewd>mdpart_mesh->kmax_pme){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("The ewald PME kmax parameter < ewald kmax\n");
      PRINTF("This is not allowed\n");
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    }/*endif*/
  }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 14)\ewald_respa_pme_opt{#} */
  ifound      = 0;
  if(strcasecmp(dict[14].keyarg,"on")==0) {
    mdpart_mesh->pme_res_on = 1; ifound++;}
  if(strcasecmp(dict[14].keyarg,"off")==0){
    mdpart_mesh->pme_res_on = 0; ifound++;}
  index=14;
  if(ifound == 0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 15)\ewald_respa_kmax_pme{#} */
  sscanf(dict[15].keyarg,"%lg",&real_key_arg);
  mdpart_mesh->kmax_pme_res = (int)(real_key_arg);
  index=15;
  if(mdpart_mesh->kmax_pme_res  <0) 
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if(mdpart_mesh->pme_res_on == 1){
    if(mdpart_mesh->kmax_pme_res >40){
      PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
      PRINTF("You have requested a RESPA PME kmax ewald parameter > 40\n");
      PRINTF("Are you certain this is what you would like to do?\n");
      PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
    }/*endif*/
    if(class_parse->kmax_res>mdpart_mesh->kmax_pme_res){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("The RESPA PME kmax parameter assigned < ewald respa kmax\n");
      PRINTF("This is not allowed\n");      
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    }/*endif*/
  }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 16)\ewald_respa_interp_pme{#} */
  sscanf(dict[16].keyarg,"%lg",&real_key_arg);
  mdpart_mesh->n_interp_res= (int)(real_key_arg);
  index=16;
  if((mdpart_mesh->n_interp_res  <4 )||
      ((mdpart_mesh->n_interp_res  %2)!=0))
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);

  /*-----------------------------------------------------------------------*/ 
  /* 17)\sep_VanderWaals{on,off} */
  ifound      = 0;
  if(strcasecmp(dict[17].keyarg,"on")==0) {
    mdenergy_ctrl->isep_vvdw = 1;ifound++;}
  if(strcasecmp(dict[17].keyarg,"off")==0){
    mdenergy_ctrl->isep_vvdw = 0;ifound++;}
  index=17;
  if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 18)\dielectric_opt{on,off} */
  ifound      = 0;
  if(strcasecmp(dict[18].keyarg,"on")==0) {
    mdinteract->dielectric_opt = 1;ifound++;}
  if(strcasecmp(dict[18].keyarg,"off")==0){
    mdinteract->dielectric_opt = 0;ifound++;}
  index=18;
  if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/
  /* 19)\dielectric_rheal{#} */
  sscanf(dict[19].keyarg,"%lg",&real_key_arg);
  mdinteract->dielectric_rheal = (double)(real_key_arg)/0.529177;
  index=19;
  if(mdinteract->dielectric_rheal<0)
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/
  /* 20)\dielectric_cut{#} */
  sscanf(dict[20].keyarg,"%lg",&real_key_arg);
  mdinteract->dielectric_cut = (double)(real_key_arg)/0.529177;
  index=20;
  if(mdinteract->dielectric_cut<0)
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/
  /* 21)\dielectric_eps{#} */
  sscanf(dict[21].keyarg,"%lg",&real_key_arg);
  mdinteract->dielectric_eps = (double)(real_key_arg);
  index=21;
  if(mdinteract->dielectric_eps<1)
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/
  /* 22)\std_intra_block{on,off} */
  /*-----------------------------------------------------------------------*/
  /* 23)\con_intra_block{on,off} */
  /*-----------------------------------------------------------------------*/
  /* 24)\inter_PE_calc_freq */
  sscanf(dict[24].keyarg,"%lg",&real_key_arg); 
  gentimeinfo->iget_pe_real_inter_freq = (int)(real_key_arg);
  mdenergy_ctrl->iget_pe_real_inter_freq     = (int)(real_key_arg);
  mdenergy_ctrl->itime        = 0;
  mdenergy_ctrl->iget_pe_real_inter = 1;
  mdenergy_ctrl->iget_pv_real_inter = 1;
  index=24;
  if(gentimeinfo->iget_pe_real_inter_freq < 1) 
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if(gentimeinfo->iget_pe_real_inter_freq > 10){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("You have requested an inter PE calc freq > 10 steps \n");
    PRINTF("Are you certain this is what you would like to do?\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/
  /*-----------------------------------------------------------------------*/
  /*  25)\inter_calc_opt{#} */
  ifound = 0;
  if(strcasecmp(dict[25].keyarg,"spline")==0)    {
    mdenergy_ctrl->lj_coul_explicit_opt = 0; ifound++;}
  if(strcasecmp(dict[25].keyarg,"grind")==0)    {
    mdenergy_ctrl->lj_coul_explicit_opt = 1; ifound++;}
  mdinteract->lj_coul_explicit_opt = mdenergy_ctrl->lj_coul_explicit_opt;

  index=25;
  if(ifound != 1){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);}
  /*-----------------------------------------------------------------------*/ 
  /*  26)\nkvec_perd_fix{#} */
  index=26;
  sscanf(dict[index].keyarg,"%d",&int_key_arg); 
  genewald->nfix_para = int_key_arg;
  if(int_key_arg<=0){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);}
  /*-----------------------------------------------------------------------*/ 
  /*  27)\nkvec_perd_expnd{#} */
  index=27;
  sscanf(dict[index].keyarg,"%d",&int_key_arg); 
  genewald->nfix_perp = int_key_arg;
  if(int_key_arg<=0){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);}

  /*========================================================================*/
}/*end routine*/ 
/*========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_sim_params_run(MDINTEGRATE *mdintegrate, MDATOMS *mdatoms,
    MDINTER *mdinter,GENERAL_DATA *general_data,
    MDINTRA *mdintra,
    CP *cp,CLASS_PARSE *class_parse,CP_PARSE *cp_parse,
    FILENAME_PARSE *filename_parse,
    DICT_WORD *dict,char *fun_key)

  /*=======================================================================*/
  /*          Begin routine                            */
{/*begin routine */
  /*=======================================================================*/
  /*          Local variable declarations                    */

#include "../class_defs/allclass_strip_mdatoms.h"
#include "../class_defs/allclass_strip_mdintegrate.h"
#include "../class_defs/allclass_strip_mdinter.h"
#include "../class_defs/allclass_strip_mdintra.h"
#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_cp.h"

  int ifound,index,pimd_on;
  int iii,int_key_arg;
  double real_key_arg;

  pimd_on = gensimopts->pimd + gensimopts->cp_pimd 
    + gensimopts->cp_wave_pimd 
    + gensimopts->debug_pimd 
    + gensimopts->debug_cp_pimd
    + gensimopts->cp_wave_min_pimd;


  /*========================================================================*/
  /* 1)\init_resmpl_atm_vel{on,off} */
  ifound      = 0;
  if(strcasecmp(dict[1].keyarg,"on")==0) {
    class_parse->ivx_smpl = 1;ifound++;}
  if(strcasecmp(dict[1].keyarg,"off")==0){
    class_parse->ivx_smpl = 0;ifound++;}
  index=1;
  if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 2)\init_resmpl_cp_vel{on,off} */
  ifound      = 0;
  if(strcasecmp(dict[2].keyarg,"on")==0)   {
    cp_parse->ivc_smpl = 1;ifound++;}
  if(strcasecmp(dict[2].keyarg,"off")==0)  {
    cp_parse->ivc_smpl = 0;ifound++;}
  index=2;
  if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 3)\init_resmpl_cp_nhc{on,off} */
  ifound      = 0;
  if(strcasecmp(dict[3].keyarg,"on")==0)  {
    cp_parse->ivcnhc_smpl = 1;ifound++;}
  if(strcasecmp(dict[3].keyarg,"off")==0) {
    cp_parse->ivcnhc_smpl = 0;ifound++;}
  index=3;
  if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 4a)\resmpl_frq_atm_vel{#} */
  sscanf(dict[4].keyarg,"%lg",&real_key_arg);
  mdvel_samp->nvx_smpl = (int)(real_key_arg);
  index=4;
  if(mdvel_samp->nvx_smpl<0)
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 4b)\rescale_frq_atm_vel{#} */
  sscanf(dict[25].keyarg,"%lg",&real_key_arg);
  mdvel_samp->nvx_scale = (int)(real_key_arg);
  index=25;
  if(mdvel_samp->nvx_scale<0)
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);


  /*-----------------------------------------------------------------------*/ 
  /* 5)\respa_steps_lrf{#} */
  sscanf(dict[5].keyarg,"%lg",&real_key_arg);
  gentimeinfo->nres_ter = (int)(real_key_arg);
  index=5;
  if(gentimeinfo->nres_ter<0){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);}

  if(gentimeinfo->nres_ter==1){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("You have requested one long-range forces RESPA step.  \n");
    PRINTF("This will result in additional overhead\n");
    PRINTF("with no gain in efficiency!! \n");
    PRINTF("Are you certain this is what you would like to do?\n");
    PRINTF("It would be better not to specify the \n");
    PRINTF("key word or give it the value zero.   \n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }
  gentimeinfo->int_res_ter   = 1;
  genewald->int_res_ter      = 1;
  mdenergy_ctrl->int_res_ter = 1;
  if(gentimeinfo->nres_ter==0){
    gentimeinfo->int_res_ter   = 0;
    genewald->int_res_ter      = 0;
    mdenergy_ctrl->int_res_ter = 0;
    gentimeinfo->nres_ter      = 1;
  }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 6)\respa_steps_torsion{#} */
  sscanf(dict[6].keyarg,"%lg",&real_key_arg);
  gentimeinfo->nres_tor = (int)(real_key_arg);
  index=6;
  if(gentimeinfo->nres_tor<0)
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  gentimeinfo->int_res_tor=1;
  if(gentimeinfo->nres_tor==0){
    gentimeinfo->int_res_tor=0;
    gentimeinfo->nres_tor=1;
  }
  /*-----------------------------------------------------------------------*/ 
  /* 7)\respa_steps_intra{#} */
  sscanf(dict[7].keyarg,"%lg",&real_key_arg);
  gentimeinfo->nres_tra = (int)(real_key_arg);
  index=7;
  if(gentimeinfo->nres_tra<0)
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  gentimeinfo->int_res_tra=1;
  mdenergy_ctrl->int_res_tra=1;
  if(gentimeinfo->nres_tra==0){
    gentimeinfo->int_res_tra=0;
    mdenergy_ctrl->int_res_tra=0;
    gentimeinfo->nres_tra=1;
  }
  /*-----------------------------------------------------------------------*/ 
  /* 8)\respa_rheal{#} */
  sscanf(dict[8].keyarg,"%lg",&(mdinteract->rheal_res));
  index=8;
  if(mdinteract->rheal_res<=0)
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if(mdinteract->rheal_res>2.5){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("You have requested RESPA healing length > 2.5A   \n");
    PRINTF("Are you certain this is what you would like to do?\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/
  mdinteract->rheal_res  /= BOHR;
  /*-----------------------------------------------------------------------*/ 
  /* 9)\shake_tol{#} */
  sscanf(dict[9].keyarg,"%lg",&(mdconstrnt->tolshake));
  index=9;
  if(mdconstrnt->tolshake<=0)
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if(mdconstrnt->tolshake<1.0e-08) {
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
    PRINTF("You have requested a SHAKE tolerence < 1.0e-08    \n");
    PRINTF("Are you certain this is what you would like to do?\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }
  /*-----------------------------------------------------------------------*/ 
  /* 10)\rattle_tol{#} */
  sscanf(dict[10].keyarg,"%lg",&(mdconstrnt->tolratl));
  index=10;
  if(mdconstrnt->tolratl <= 0)
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if(mdconstrnt->tolratl<1.0e-08) {
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
    PRINTF("You have requested a RATTLE tolerence < 1.0e-08    \n");
    PRINTF("Are you certain this is what you would like to do?\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }
  /*-----------------------------------------------------------------------*/ 
  /* 11)\max_constrnt_iter{#} */
  sscanf(dict[11].keyarg,"%lg",&real_key_arg);
  mdconstrnt->max_iter = (int)(real_key_arg);
  index=11;
  if(mdconstrnt->max_iter <= 0) 
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  mdgrp_bond_con->max_iter = mdconstrnt->max_iter;
  if(mdconstrnt->max_iter<100){ 
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("You have requested a max_constrnt_iter < 100  \n");
    PRINTF("Are you certain this is what you would like to do?\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 12)\init_rescale_atm_vel{on,off} */
  ifound      = 0;
  if(strcasecmp(dict[12].keyarg,"on")==0) {
    class_parse->ivx_scale = 1;ifound++;}
  if(strcasecmp(dict[12].keyarg,"off")==0){
    class_parse->ivx_scale = 0;ifound++;}
  index=12;
  if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 13)\init_rescale_atm_nhc{on,off} */
  ifound      = 0;
  if(strcasecmp(dict[13].keyarg,"on")==0) {
    class_parse->ivnhc_scale = 1;ifound++;}
  if(strcasecmp(dict[13].keyarg,"off")==0){
    class_parse->ivnhc_scale = 0;ifound++;}
  index=13;
  if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 14)\init_rescale_cp_vel{on,off} */
  ifound      = 0;
  if(strcasecmp(dict[14].keyarg,"on")==0)   {
    cp_parse->ivc_scale = 1;ifound++;}
  if(strcasecmp(dict[14].keyarg,"off")==0)  {
    cp_parse->ivc_scale = 0;ifound++;}
  index=14;
  if((cpopts->cp_isok_opt == 1) && (cp_parse->ivc_scale != 1)){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("When using the CP isokinetic method, it is advisable\n");
    PRINTF("to do an initial rescale of your CP velocities\n");
    PRINTF("Therefore, I am setting this flag for you.\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
  }/* endif */
  if(cpopts->cp_isok_opt == 1) cp_parse->ivc_scale = 1;
  if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 15)\resmpl_frq_cp_vel{ # } */
  sscanf(dict[15].keyarg,"%lg",&real_key_arg);
  cpvel_samp->nvc_smpl = (int)(real_key_arg);
  index=15;
  if(cpvel_samp->nvc_smpl<0)
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 16)\group_con_tol{#} */
  sscanf(dict[16].keyarg,"%lg",&(mdgrp_bond_con->tol));
  index=16;
  if(mdgrp_bond_con->tol <= 0)
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if(mdgrp_bond_con->tol<1.0e-08) {
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
    PRINTF("You have requested a GROUP CONSTRAINT tolerence < 1.0e-08\n");
    PRINTF("Are you certain this is what you would like to do?\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }
  /*-----------------------------------------------------------------------*/ 
  /* 17)\cp_norb_tol{#} */
  sscanf(dict[17].keyarg,"%lg",&(cpconstrnt->c_tolnorb));
  index=17;
  if(cpconstrnt->c_tolnorb<=0)
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);

  if(cpconstrnt->c_tolnorb>0.01&&cpopts->cp_norb==2){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("The tolerence must be quite low with the    \n");
    PRINTF("norm_only norb option. The rotation does not\n");
    PRINTF("preserve the individual norms only the sum.  \n");
    PRINTF("Work is in progress to overcome this.     \n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }/*endif*/

  /*-----------------------------------------------------------------------*/ 
  /* 18)\cp_ks_rot_frq{#} */
  cpopts->ks_rot_on = 1;      
  sscanf(dict[18].keyarg,"%lg",&real_key_arg);
  cpopts->n_ks_rot = (int)(real_key_arg);
  index=18;
  if(cpopts->n_ks_rot<0)
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if(cpopts->n_ks_rot==0){
    cpopts->ks_rot_on=0;
    cpopts->n_ks_rot=1;
  }
  if(cpopts->ks_rot_on == 1 && cpopts->cp_norb > 0){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");    
    PRINTF("You are studying the Kohn-Sham eigenvalues under cp_norb \n");
    PRINTF("Eigenvalues are not guaranteed to be correct\n");
    PRINTF("Turn norb off to insure the correct eigenvalues\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }
  /*-----------------------------------------------------------------------*/ 
  /* 19)\cp_shake_tol{#} */
  sscanf(dict[19].keyarg,"%lg",&(cpconstrnt->c_tolshake));
  index=19;
  if(cpconstrnt->c_tolshake<=0)
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 20)\cp_rattle_tol{#} */
  sscanf(dict[20].keyarg,"%lg",&(cpconstrnt->c_tolratl));
  index=20;
  if(cpconstrnt->c_tolratl <= 0) 
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/
  /* 24)\cp_min_tol{#} */
  sscanf(dict[24].keyarg,"%lg",&(genminopts->tol_coef));
  index=24;
  if(genminopts->tol_coef<=0)
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/

  /*-----------------------------------------------------------------------*/ 
  /* 21)\cp_run_tol{#} */
  sscanf(dict[21].keyarg,"%lg",&(cpopts->tol_coef));
  index=21;
  if(cpopts->tol_coef<=0)
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if(cpopts->tol_coef < genminopts->tol_coef ){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("You have requested a cp run tolerence less than your \n");
    PRINTF("minimization tolerence %g %g \n",cpopts->tol_coef,
        genminopts->tol_coef);
    PRINTF("Are you certain this is what you would like to do?\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/
  if(cpopts->tol_coef > 10000.0*genminopts->tol_coef ){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("You have requested a cp run tolerence greater than  \n");
    PRINTF("10000 times the minimization tolerence %g %g \n",
        cpopts->tol_coef,
        genminopts->tol_coef);
    PRINTF("Are you certain this is what you would like to do?\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/
  /*-----------------------------------------------------------------------*/
  /* 22)\zero_com_vel{yes,no}*/
  ifound = 0;
  if(strcasecmp(dict[22].keyarg,"yes")==0){
    class_parse->zero_com_vel=1;ifound++;}
  if(strcasecmp(dict[22].keyarg,"no" )==0){
    class_parse->zero_com_vel=0;ifound++;}
  index = 22;
  if(ifound != 1){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);  
  }
  /*-----------------------------------------------------------------------*/
  /* 23)\min_tol{#} */
  sscanf(dict[23].keyarg,"%lg",&(genminopts->tol_atom));
  index=23;
  if(genminopts->tol_atom<=0)
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 26)hess_opt{full,unit} */
  ifound      = 0;
  if(strcasecmp(dict[26].keyarg,"full")==0)   {
    cpopts->cp_hess_calc = 1;ifound++;}
  if(strcasecmp(dict[26].keyarg,"unit")==0)  {
    cpopts->cp_hess_calc = 0;ifound++;}
  index=26;
  if(ifound != 1){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  }
  if(gensimopts->cp_min != 1 || 
      gensimopts->debug_cp != 1){
    if((genminopts->min_cg == 1) &&
        cpopts->cp_hess_calc > 0 ){
      PRINTF("$$$$$$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$$$$\n");
      PRINTF("You have chosen to do a conjugate gradient minimization \n");
      PRINTF("with a nonunit hessian preconditioner.            \n");
      PRINTF("There is no rigorous basis for this approach, and it    \n");
      PRINTF("is likely not to work.  I will let you try as an     \n");
      PRINTF("experiment, however, if you find that the minimizer is  \n");
      PRINTF("is not behaving well, consider the \"unit\" option\n");
      PRINTF("for the \"hess_opt\" keyword\n");
      PRINTF("$$$$$$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$$$$\n");
    }/*endif*/
  }/*endif*/
  /*-----------------------------------------------------------------------*/
  /* 27)\class_mass_scale_fact{#}  (Masses divided by this factor) */
  sscanf(dict[27].keyarg,"%lg",&(mdclatoms_info->mass_sc_fact));
  index=27;
  if(mdclatoms_info->mass_sc_fact<=0)
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if((gensimopts->cp_min + 
        gensimopts->minimize == 0) &&
      mdclatoms_info->mass_sc_fact != 1.0 ){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("Dude, if you are not doing a minimization run     \n");
    PRINTF("then, like, why are you scaling your atom masses??\n");
    PRINTF("It's cool!  I'll reset the scaling factor for you!\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
    mdclatoms_info->mass_sc_fact = 1.0;
  }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 28)\hmat_int_typ{normal,upper_triangle} */
  ifound      = 0;
  if(strcasecmp(dict[28].keyarg,"normal")==0)   {
    gencell->hmat_int_typ = 0;ifound++;}
  if(strcasecmp(dict[28].keyarg,"upper_triangle")==0)  {
    gencell->hmat_int_typ = 1;ifound++;}
  index=28;

  if(genensopts->npt_f==1){
    if((pimd_on>0)&&(gencell->hmat_int_typ==0)){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("You are running a path integral simulation under\n");
      PRINTF("NPTF.  You must use the upper diagonal form of the\n");
      PRINTF("cell matrix.\n");
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }/*endif*/
  }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 29)\hmat_cons_typ{none,ortho_rhom,mono_clin} */
  ifound      = 0;
  if(strcasecmp(dict[29].keyarg,"none")==0)   {
    gencell->hmat_cons_typ = 0;ifound++;}
  if(strcasecmp(dict[29].keyarg,"ortho_rhom")==0)  {
    gencell->hmat_cons_typ = 1;ifound++;}
  if(strcasecmp(dict[29].keyarg,"mono_clin")==0)  {
    gencell->hmat_cons_typ = 2;ifound++;}

  index=29;
  if(ifound == 0){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  }/*endif*/

  /*-----------------------------------------------------------------------*/ 
  /* 30)\hmat_cons_typ{none,ortho_rhom,mono_clin} */
  ifound      = 0;
  if(strcasecmp(dict[30].keyarg,"yes")==0)   {
    genminopts->min_atm_com_fix_opt = 1;ifound++;}
  if(strcasecmp(dict[30].keyarg,"no")==0)  {
    genminopts->min_atm_com_fix_opt = 0;ifound++;}
  index=30;
  if(ifound == 0){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  }/*endif*/

  /*-----------------------------------------------------------------------*/ 
  /* 31)\rescale_frq_cp_vel{ # } */
  sscanf(dict[31].keyarg,"%lg",&real_key_arg);
  cpvel_samp->nvc_scal = (int)(real_key_arg);
  index=31;
  if(cpvel_samp->nvc_scal<0)
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 32)\auto_rescale_cp_vel{ on/off } */
  ifound = 0;
  if(strcasecmp(dict[32].keyarg,"on")==0)   {
    cpvel_samp->iauto_vc_scal_opt = 1;ifound++;}
  if(strcasecmp(dict[32].keyarg,"off")==0)  {
    cpvel_samp->iauto_vc_scal_opt = 0;ifound++;}
  index=32;
  if(ifound == 0){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 33)\auto_rescale_cp_vel_tol{ # } */
  sscanf(dict[33].keyarg,"%lg",&real_key_arg);
  cpvel_samp->vc_scal_tol = real_key_arg;
  index=33;
  if(cpvel_samp->vc_scal_tol<=1.0){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);}
  if(cpvel_samp->vc_scal_tol<3.0){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("Coef velocity Rescale tolerence set less than 3kT\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
  }/*endif*/
  if(cpvel_samp->vc_scal_tol>10.0){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("Coef velocity Rescale tolerence set greater than 10kT\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
  }/*endif*/

  /*-----------------------------------------------------------------------*/ 
  /* 34)\cp_norb_rot_kescal */
  ifound = 0;
  if(strcasecmp(dict[34].keyarg,"on")==0)   {
    cpopts->cp_norb_rot_kescal = 1; ifound++;}
  if(strcasecmp(dict[34].keyarg,"off")==0)  {
    cpopts->cp_norb_rot_kescal = 0; ifound++;}
  index=34;
  if(ifound == 0){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 35)\norb_fin_diff_order */
  ifound = 0;
  sscanf(dict[35].keyarg,"%d",&int_key_arg);
  cpopts->norb_fin_diff_order = int_key_arg;
  index=35;
  if(int_key_arg<0 || int_key_arg>2){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  }/*endif*/

  /*-----------------------------------------------------------------------*/ 
  /* 36)\cp_min_diagonalize */
  ifound = 0;
  if(strcasecmp(dict[36].keyarg,"on")==0)   {
    cpopts->cp_min_diagonalize = 1; ifound++;}
  if(strcasecmp(dict[36].keyarg,"off")==0)  {
    cpopts->cp_min_diagonalize = 0; ifound++;}
  index=36;
  if(ifound == 0){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  }/*endif*/

  /*========================================================================*/
}/*end routine*/ 
/*========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_sim_params_nhc(MDINTEGRATE *mdintegrate, MDATOMS *mdatoms,
    MDINTER *mdinter,GENERAL_DATA *general_data,
    MDINTRA *mdintra,
    CP *cp,CLASS_PARSE *class_parse,CP_PARSE *cp_parse,
    FILENAME_PARSE *filename_parse,
    DICT_WORD *dict,char *fun_key)

  /*=======================================================================*/
  /*          Begin routine                            */
{/*begin routine */
  /*=======================================================================*/
  /*          Local variable declarations                    */

#include "../class_defs/allclass_strip_mdatoms.h"
#include "../class_defs/allclass_strip_mdintegrate.h"
#include "../class_defs/allclass_strip_mdinter.h"
#include "../class_defs/allclass_strip_mdintra.h"
#include "../class_defs/allclass_strip_cp.h"
#include "../class_defs/allclass_strip_gen.h"

  int ifound,index,int_key_arg;
  double real_key_arg;

  /*========================================================================*/
  /* 1)\respa_steps_nhc{#} */
  sscanf(dict[1].keyarg,"%lg",&real_key_arg);
  mdtherm_info_bead->nres_nhc = (int)(real_key_arg);
  mdtherm_info->nres_nhc = (int)(real_key_arg);
  index=1;
  if(mdtherm_info->nres_nhc < 1) 
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if(mdtherm_info->nres_nhc>5){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("You have requested more than 5 NHC RESPA steps    \n");
    PRINTF("Are you certain this is what you would like to do?\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 2)\yosh_steps_nhc{#} */
  sscanf(dict[2].keyarg,"%lg",&real_key_arg);
  mdtherm_info_bead->nyosh_nhc = (int)(real_key_arg);
  mdtherm_info->nyosh_nhc = (int)(real_key_arg);
  index=2;
  if(mdtherm_info->nyosh_nhc != 1 && 
      mdtherm_info->nyosh_nhc != 3 && 
      mdtherm_info->nyosh_nhc != 5 && 
      mdtherm_info->nyosh_nhc != 7 &&
      mdtherm_info->nyosh_nhc != 9 )
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 3)\atm_nhc_tau_def{#} */
  sscanf(dict[3].keyarg,"%lg",&(class_parse->tau_nhc_def));
  index=3;
  if(class_parse->tau_nhc_def <= 0) 
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if(class_parse->tau_nhc_def < 10.0 ){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("You have requested a NHC particle default tau < 10fs\n");
    PRINTF("Are you certain this is what you would like to do?\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 4)\cp_nhc_len{#} */
  sscanf(dict[4].keyarg,"%lg",&real_key_arg);
  cptherm_info->len_c_nhc = (int)(real_key_arg);
  index=4;
  if(cptherm_info->len_c_nhc != 2) {
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /* 5)\resmpl_frq_cp_nhc{ #} */
  sscanf(dict[5].keyarg,"%lg",&real_key_arg);
  cpvel_samp->nvcnhc_smpl = (int)(real_key_arg);
  index=5;
  if(cpvel_samp->nvcnhc_smpl<0)
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/  
  /* 6)\init_resmpl_atm_nhc{on,off} */
  ifound      = 0;
  if(strcasecmp(dict[6].keyarg,"on")==0) {
    class_parse->ivnhc_smpl = 1;ifound++;}
  if(strcasecmp(dict[6].keyarg,"off")==0){
    class_parse->ivnhc_smpl = 0;ifound++;}
  index=6;
  if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 7)\init_resmpl_cp_nhc{on,off} */
  ifound      = 0;
  if(strcasecmp(dict[7].keyarg,"on")==0)  {
    cp_parse->ivcnhc_scale = 1;ifound++;}
  if(strcasecmp(dict[7].keyarg,"off")==0) {
    cp_parse->ivcnhc_scale = 0;ifound++;}
  index=7;
  if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 8)\resmpl_atm_nhc{#} */
  sscanf(dict[8].keyarg,"%lg",&real_key_arg);
  mdvel_samp->nvnhc_smpl = (int)(real_key_arg);
  index=8;
  if(mdvel_samp->nvnhc_smpl<0)
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 9)\atm_nhc_len{#} */              
  sscanf(dict[9].keyarg,"%lg",&real_key_arg);
  mdtherm_info_bead->len_nhc  = (int)(real_key_arg);
  mdtherm_info->len_nhc = (int)(real_key_arg);
  mdintegrate->mdbaro.len_nhc      = (int)(real_key_arg);
  index=9;
  if(mdtherm_info_bead->len_nhc < 0){ 
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);}
  if(mdtherm_info_bead->therm_typ == 1
      && mdtherm_info_bead->len_nhc > 3 ){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
    PRINTF("You have requested a NHC length > 3      \n");
    PRINTF("Are you certain this is what you would like to do?\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/
  /*-----------------------------------------------------------------------*/
  /* 10)\cp_nhc_tau_def{#} */
  sscanf(dict[10].keyarg,"%lg",&(cp_parse->cp_tau_nhc_def));
  index=10;
  if(cp_parse->cp_tau_nhc_def<=0) 
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  cpopts->cp_tau_nhc  = cp_parse->cp_tau_nhc_def;
  /*-----------------------------------------------------------------------*/ 
  /* 11)\cp_respa_steps_nhc{#} */
  sscanf(dict[11].keyarg,"%lg",&real_key_arg);
  cptherm_info->nres_c_nhc = (int)(real_key_arg);
  index=11;
  if(cptherm_info->nres_c_nhc < 0) 
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 12)\cp_yosh_steps_nhc{#} */
  sscanf(dict[12].keyarg,"%lg",&real_key_arg);
  cptherm_info->nyosh_c_nhc = (int)(real_key_arg);
  index=12;
  if(cptherm_info->nyosh_c_nhc != 1 && cptherm_info->nyosh_c_nhc != 3 && 
      cptherm_info->nyosh_c_nhc != 5 && cptherm_info->nyosh_c_nhc != 7
      && cptherm_info->nyosh_c_nhc != 9)
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 13)\respa_xi_opt{#} */
  sscanf(dict[13].keyarg,"%lg",&real_key_arg);
  gentimeinfo->ix_respa = (int)(real_key_arg);
  index=13;
  if(gentimeinfo->ix_respa <= 0 
      || gentimeinfo->ix_respa > 5) 
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if(gentimeinfo->ix_respa>1) {
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
    PRINTF("WARNING: xi respa option chosen!\n");
    PRINTF("energy drifts may occur if the number \n");
    PRINTF("of respa steps is taken too large\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }
  /*-----------------------------------------------------------------------*/
  /* 14)\thermostat_type{#} */

  if(strcasecmp(dict[14].keyarg,"NHC")==0)  {
    mdtherm_info_bead->therm_typ = 1;
    mdtherm_info->therm_typ = 1;
  }/*endif*/

  if(strcasecmp(dict[14].keyarg,"GGMT")==0)  {
    mdtherm_info_bead->therm_typ = 2;
    mdtherm_info->therm_typ = 2;
  }/*endif*/

  index=14;
  if(mdtherm_info->therm_typ <=0
      || mdtherm_info->therm_typ >= 3) {
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  }

  if(mdtherm_info_bead->therm_typ == 2
      && mdtherm_info_bead->len_nhc > 3 ){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("You have requested a GGMT length > 3      \n");
    PRINTF("Available GGMT lengths are 2 or 3    \n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    EXIT(1);
  }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 15)\cp_therm_heat_fact{#} */
  sscanf(dict[15].keyarg,"%lg",&real_key_arg);
  index=15;
  if(real_key_arg<1.0){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  }/*endif*/

  /*-----------------------------------------------------------------------*/ 
  /* 16)\cp_num_nhc_iso{#} */
  sscanf(dict[16].keyarg,"%d",&int_key_arg);
  cptherm_info->num_c_nhc_iso = int_key_arg;
  index=16;
  if(int_key_arg<1){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 17)\atm_isokin_opt{on/off} */
  ifound      = 0;
  if(strcasecmp(dict[17].keyarg,"on")==0)  {
    mdtherm_info->isokin_opt = 1; ;ifound++;}
  if(strcasecmp(dict[17].keyarg,"off")==0) {
    mdtherm_info->isokin_opt = 0; ;ifound++;}
  index=17;
  if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if(mdtherm_info->isokin_opt==1 && mdtherm_info->len_nhc<3){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("With isokin_opt you need len_nhc >=3 \n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    EXIT(1);
  }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 18)\atm_num_nhc_iso{#} */
  sscanf(dict[18].keyarg,"%d",&int_key_arg);
  mdtherm_info->num_nhc_iso = int_key_arg;
  index=18;
  if(int_key_arg<1){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  }/*endif*/
  /*---------------------------------------------------------------------*/
  /*  19)\cp_nhc_chunk{#} */
  index=19;
  sscanf(dict[19].keyarg,"%d",&int_key_arg);
  cptherm_info->nck_c_nhc = int_key_arg;
  if(int_key_arg<1){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  }/*endif*/

  /*========================================================================*/
}/*end routine*/ 
/*========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_sim_params_vol(MDINTEGRATE *mdintegrate, MDATOMS *mdatoms,
    MDINTER *mdinter,GENERAL_DATA *general_data,
    MDINTRA *mdintra,
    CP *cp,CLASS_PARSE *class_parse,CP_PARSE *cp_parse,
    FILENAME_PARSE *filename_parse,
    DICT_WORD *dict,char *fun_key)

  /*=======================================================================*/
  /*          Begin routine                            */
{/*begin routine */
  /*=======================================================================*/
  /*          Local variable declarations                    */

#include "../class_defs/allclass_strip_mdatoms.h"
#include "../class_defs/allclass_strip_mdintegrate.h"
#include "../class_defs/allclass_strip_mdinter.h"
#include "../class_defs/allclass_strip_mdintra.h"
#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_cp.h"

  int ifound,index;
  double real_key_arg;
  int cp_dual_grid_opt_on = cpopts->cp_dual_grid_opt;

  /*========================================================================*/
  /* 1)\volume_tau{#} */
  sscanf(dict[1].keyarg,"%lg",&(class_parse->tau_vol));
  index=1;
  if(class_parse->tau_vol <= 0) 
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if(class_parse->tau_vol > 10000.0){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("You have requested a volume tau > 10ps        \n");
    PRINTF("Are you certain this is what you would like to do?\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/
  if(class_parse->tau_vol < 100.0 ){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("You have requested a volume tau < 100fs        \n");
    PRINTF("Are you certain this is what you would like to do?\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 2)\volume_nhc_tau{#} */
  sscanf(dict[2].keyarg,"%lg",&(class_parse->tau_vol_nhc));
  index=2;
  if(class_parse->tau_vol_nhc <= 0) 
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if(class_parse->tau_vol_nhc > 10000.0){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("You have requested a NHC volume tau  > 10ps      \n");
    PRINTF("Are you certain this is what you would like to do?\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/
  if(class_parse->tau_vol_nhc < 100.0 ){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("You have requested a NHC volume tau < 100fs      \n");
    PRINTF("Are you certain this is what you would like to do?\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 3)\periodicity{0,1,2,3} */
  if(strcasecmp(dict[3].keyarg,"0_ewald")==0)  {
    gencell->iperd = 4;
  }else{
    sscanf(dict[3].keyarg,"%lg",&real_key_arg);
    gencell->iperd = (int)(real_key_arg);
  }/*endif*/
  genewald->iperd = gencell->iperd;
  index=3;
  if(gencell->iperd == 4 && cp_dual_grid_opt_on >= 1){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("cluster ewald  not supported for dualed systems yet. \n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }/*endif*/
  if((gencell->iperd<0)||(gencell->iperd>4))
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 4)\intra_perds{on,off} */
  ifound      = 0;
  if(strcasecmp(dict[4].keyarg,"on")==0) {
    gencell->intra_perds = 1;ifound++;}
  if(strcasecmp(dict[4].keyarg,"off")==0){
    gencell->intra_perds = 0;ifound++;}
  index=4;
  if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/
  /*========================================================================*/
}/*end routine*/ 
/*========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_sim_params_write(MDINTEGRATE *mdintegrate, MDATOMS *mdatoms,
    MDINTER *mdinter,GENERAL_DATA *general_data,
    MDINTRA *mdintra,
    CP *cp,CLASS_PARSE *class_parse,CP_PARSE *cp_parse,
    FILENAME_PARSE *filename_parse,
    DICT_WORD *dict,char *fun_key)

  /*=======================================================================*/
  /*          Begin routine                            */
{/*begin routine */
  /*=======================================================================*/
  /*          Local variable declarations                    */

#include "../class_defs/allclass_strip_mdatoms.h"
#include "../class_defs/allclass_strip_mdintegrate.h"
#include "../class_defs/allclass_strip_mdinter.h"
#include "../class_defs/allclass_strip_mdintra.h"
#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_cp.h"

  int ifound,index;
  double real_key_arg;
  char *strip,*strip2;

  /*========================================================================*/

  strip = (char *)cmalloc(MAXWORD*sizeof(char),"set_sim_params_write");
  strip2 = (char *)cmalloc(MAXWORD*sizeof(char),"set_sim_params_write");

  /*========================================================================*/
  /* 1)\write_binary_cp_coef{on,off} */
  ifound  = 0;
  if(strcasecmp(dict[1].keyarg,"off")==0) 
  {cpopts->iwrite_coef_binary = 0;ifound++;}
  if(strcasecmp(dict[1].keyarg,"on")==0)  
  {cpopts->iwrite_coef_binary = 1;ifound++;}
  if(strcasecmp(dict[1].keyarg,"off_gzip")==0) 
  {cpopts->iwrite_coef_binary = 2;ifound++;}
  if(strcasecmp(dict[1].keyarg,"on_gzip")==0) 
  {cpopts->iwrite_coef_binary = 3;ifound++;}
  index=1;
  if(ifound == 0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 2)\write_force_freq{#} */
  sscanf(dict[2].keyarg,"%lg",&real_key_arg);
  genfilenames->iwrite_atm_for = (int)(real_key_arg);
  index=2;
  if(genfilenames->iwrite_atm_for < 1) 
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 3)\write_screen_freq */
  sscanf(dict[3].keyarg,"%lg",&real_key_arg);
  genfilenames->iwrite_screen = (int)(real_key_arg);
  index=3;
  if(genfilenames->iwrite_screen < 1) 
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if(genfilenames->iwrite_screen > 1000){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("You have requested a write screen freq > 1000 steps \n");
    PRINTF("Are you certain this is what you would like to do?\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 4)\write_dump_freq */
  sscanf(dict[4].keyarg,"%lg",&real_key_arg);
  genfilenames->iwrite_dump = (int)(real_key_arg);
  index=4;
  if(genfilenames->iwrite_dump < 1) 
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if(genfilenames->iwrite_dump > 2000){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("You have requested a write dump freq > 2000 steps \n");
    PRINTF("Are you certain this is what you would like to do?\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 5)\write_inst_freq */
  sscanf(dict[5].keyarg,"%lg",&real_key_arg);
  genfilenames->iwrite_inst = (int)(real_key_arg);
  index=5;
  if(genfilenames->iwrite_inst < 1) 
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if(genfilenames->iwrite_inst > 2000){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("You have requested a write inst freq > 2000 steps \n");
    PRINTF("Are you certain this is what you would like to do?\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 6)\write_pos_freq */
  sscanf(dict[6].keyarg,"%lg",&real_key_arg);
  genfilenames->iwrite_confp = (int)(real_key_arg);
  index=6;
  if(genfilenames->iwrite_confp < 1) 
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if(genfilenames->iwrite_confp > 2000){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("You have requested a write position freq > 2000 steps \n");
    PRINTF("Are you certain this is what you would like to do?\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 7)\write_vel_freq */
  sscanf(dict[7].keyarg,"%lg",&real_key_arg);
  genfilenames->iwrite_confv = (int)(real_key_arg);
  index=7;
  if(genfilenames->iwrite_confv < 1) 
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if(genfilenames->iwrite_confv > 2000){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("You have requested a write velocity freq > 2000 steps \n");
    PRINTF("Are you certain this is what you would like to do?\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 8)\write_cp_c_freq */
  sscanf(dict[8].keyarg,"%lg",&real_key_arg);
  genfilenames->iwrite_confc = (int)(real_key_arg);
  index=8;
  if(genfilenames->iwrite_confc < 1) 
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/  
  /* 9)\screen_output_units */
  ifound = 0;
  if(strcasecmp(dict[9].keyarg,"au")==0)    {
    genfilenames->iwrite_units=0; ifound++;}
  if(strcasecmp(dict[9].keyarg,"kcal_mol")==0)    {
    genfilenames->iwrite_units=1; ifound++;}
  if(strcasecmp(dict[9].keyarg,"kelvin")==0)    {
    genfilenames->iwrite_units=2; ifound++;}
  index=9;
  if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/  
  /* 10)\conf_file_format */
  ifound = 0;
  if(strcasecmp(dict[10].keyarg,"binary")==0)    {
    genfilenames->iwrite_conf_binary=1; ifound++;}
  if(strcasecmp(dict[10].keyarg,"formatted")==0)  {
    genfilenames->iwrite_conf_binary=0; ifound++;}
  index=10;
  if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/  
  /* 11)\path_cent_file */
  sscanf(dict[11].keyarg,"%s",genfilenames->centname);
  /*-----------------------------------------------------------------------*/ 
  /* 12)\atm_force_file */
  sscanf(dict[12].keyarg,"%s",genfilenames->forcename);
  /*-----------------------------------------------------------------------*/ 
  /* 13)\cp_coef_file */
  sscanf(dict[13].keyarg,"%s",genfilenames->ccname);
  /*-----------------------------------------------------------------------*/ 
  /* 14)\conf_partial_freq */
  sscanf(dict[14].keyarg,"%lg",&real_key_arg);
  genfilenames->iwrite_par_confp = (int)(real_key_arg);
  index=14;
  if(genfilenames->iwrite_par_confp < 1) 
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if(genfilenames->iwrite_par_confp > 2000){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("You have requested a partial write pos freq > 2000 steps \n");
    PRINTF("Are you certain this is what you would like to do?\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 15)\path_cent_freq */
  sscanf(dict[15].keyarg,"%lg",&real_key_arg);
  genfilenames->iwrite_path_cent = (int)(real_key_arg);
  index=15;
  if(genfilenames->iwrite_path_cent < 1) 
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if(genfilenames->iwrite_path_cent > 2000){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("You have requested a write centroid freq > 2000 steps \n");
    PRINTF("Are you certain this is what you would like to do?\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 16)\conf_partial_limits */
  strcpy(strip,dict[16].keyarg);
  parse_part_lim(strip,strip2,&(genfilenames->low_lim_par),
      &(genfilenames->high_lim_par));     
  index=16;
  if((genfilenames->low_lim_par <= 0) || 
      (genfilenames->high_lim_par < 0)) 
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if((genfilenames->low_lim_par > 
        genfilenames->high_lim_par)&&
      (dict[16].iuset==1)){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("You have requested a lower limit higher than your \n");
    PRINTF("upper limit in the conf_partial_limits keyarg     \n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 17)\read_binary_cp_coef{on,off} */
  ifound  = 0;
  if(strcasecmp(dict[17].keyarg,"off")==0) 
  {cpopts->iread_coef_binary = 0;ifound++;}
  if(strcasecmp(dict[17].keyarg,"on")==0)  
  {cpopts->iread_coef_binary = 1;ifound++;}
  if(strcasecmp(dict[17].keyarg,"off_gzip")==0) 
  {cpopts->iread_coef_binary = 2;ifound++;}
  if(strcasecmp(dict[17].keyarg,"on_gzip")==0)  
  {cpopts->iread_coef_binary = 3;ifound++;}
  index=17;
  if(ifound == 0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/  
  /* 18)\sim_name */
  sscanf(dict[18].keyarg,"%s",filename_parse->simname);
  /*-----------------------------------------------------------------------*/ 
  /* 19)\out_restart_file */
  index = 19;
  check_for_slash(dict[index].keyarg,dict[index].keyword,&ifound);
  if(ifound!=0){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  }/*endif*/
  sscanf(dict[index].keyarg,"%s",genfilenames->dname);
  /*-----------------------------------------------------------------------*/ 
  /* 20)\in_restart_file */
  index = 20;
  check_for_slash(dict[index].keyarg,dict[index].keyword,&ifound);
  if(ifound!=0){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  }/*endif*/
  sscanf(dict[20].keyarg,"%s",filename_parse->dnamei);
  sscanf(dict[20].keyarg,"%s",genfilenames->dnamei);
  /*-----------------------------------------------------------------------*/ 
  /* 21)\instant_file */
  sscanf(dict[21].keyarg,"%s",genfilenames->iname);
  /*-----------------------------------------------------------------------*/ 
  /* 22)\atm_pos_file */
  index = 22;
  check_for_slash(dict[index].keyarg,dict[index].keyword,&ifound);
  if(ifound!=0){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  }/*endif*/
  sscanf(dict[22].keyarg,"%s",genfilenames->cpname);
  /*-----------------------------------------------------------------------*/ 
  /* 23)\atm_vel_file */
  index = 23;
  check_for_slash(dict[index].keyarg,dict[index].keyword,&ifound);
  if(ifound!=0){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  }/*endif*/
  sscanf(dict[23].keyarg,"%s",genfilenames->cvname);
  /*-----------------------------------------------------------------------*/ 
  /* 24)\conf_partial_file */
  index = 24;
  check_for_slash(dict[index].keyarg,dict[index].keyword,&ifound);
  if(ifound!=0){
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  }/*endif*/
  sscanf(dict[24].keyarg,"%s",genfilenames->cpparname);
  /*-----------------------------------------------------------------------*/ 
  /* 25)\mol_set_file */
  sscanf(dict[25].keyarg,"%s",filename_parse->molsetname);
  /*-----------------------------------------------------------------------*/ 
  /* 26)\cp_restart_out_file */
  sscanf(dict[26].keyarg,"%s",genfilenames->dnamec);
  /*-----------------------------------------------------------------------*/ 
  /* 27)\cp_restart_in_file */
  sscanf(dict[27].keyarg,"%s",filename_parse->dnameci);
  /*-----------------------------------------------------------------------*/ 
  /* 28)\cp_kseigs_file */
  sscanf(dict[28].keyarg,"%s",genfilenames->ksname);
  /*-----------------------------------------------------------------------*/ 
  /* 29)\cp_elf_file */
  sscanf(dict[29].keyarg,"%s",genfilenames->elfname);
  /*-----------------------------------------------------------------------*/ 
  /*  30)\atm_coord_dir_in */
  sscanf(dict[30].keyarg,"%s",genfilenames->atm_crd_dir_in);
  /*-----------------------------------------------------------------------*/ 
  /*  31)\atm_coord_dir_out */
  sscanf(dict[31].keyarg,"%s",genfilenames->atm_crd_dir_out);
  /*========================================================================*/

  cfree(strip,"set_sim_params_write");
  cfree(strip2,"set_sim_params_write");

  /*========================================================================*/
}/*end routine*/ 
/*========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_sim_params_pimd(MDINTEGRATE *mdintegrate, MDATOMS *mdatoms,
    MDINTER *mdinter,GENERAL_DATA *general_data,
    MDINTRA *mdintra,
    CP *cp,CLASS_PARSE *class_parse,CP_PARSE *cp_parse,
    FILENAME_PARSE *filename_parse,
    DICT_WORD *dict,char *fun_key)

  /*=======================================================================*/
  /*          Begin routine                            */
{/*begin routine */
  /*=======================================================================*/
  /*          Local variable declarations                    */

#include "../class_defs/allclass_strip_mdatoms.h"
#include "../class_defs/allclass_strip_mdintegrate.h"
#include "../class_defs/allclass_strip_mdinter.h"
#include "../class_defs/allclass_strip_mdintra.h"
#include "../class_defs/allclass_strip_cp.h"
#include "../class_defs/allclass_strip_gen.h"

  int ifound,index;
  double real_key_arg;

  /*========================================================================*/
  /* 1)\path_int_beads{#} */
  sscanf(dict[1].keyarg,"%lg",&real_key_arg);
  mdclatoms_info->pi_beads = (int)(real_key_arg);
  mdclatoms_pimd->pi_beads = (int)(real_key_arg);
  cpcoeffs_info->pi_beads = (int)(real_key_arg);
  gensimopts->pi_beads = (int)(real_key_arg);
  index=1;
  if((mdclatoms_info->pi_beads<0))
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);

  /*-----------------------------------------------------------------------*/ 
  /* 2)\path_int_gamma_adb{#} */
  sscanf(dict[2].keyarg,"%lg",&real_key_arg);
  mdclatoms_pimd->gamma_adb = (real_key_arg);
  index=2;
  if((mdclatoms_pimd->gamma_adb<0.0))
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if(mdclatoms_pimd->gamma_adb<1.0){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("You have requested a path integral          \n");
    PRINTF("adiabaticity parameter less than one.          \n");
    PRINTF("Are you certain this is what you would like to do?\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 3)\path_int_md_typ{staging,centroid} */
  ifound      = 0;
  if(strcasecmp(dict[3].keyarg,"staging")==0) {
    gensimopts->pi_md_typ = 1;ifound++;}
  if(strcasecmp(dict[3].keyarg,"centroid")==0){
    gensimopts->pi_md_typ = 2;ifound++;}
  index=3;
  if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if(gensimopts->pi_md_typ!=1){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("You have requested the centroid PIMD method      \n");
    PRINTF("We are not ready for that yet but if you         \n");
    PRINTF("implement it, please check it to the OpenSource  \n");
    PRINTF("repostitory (once it works, of course)           \n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
    EXIT(1);
  }/*endif*/
  if(mdclatoms_pimd->gamma_adb>1.0 && gensimopts->pi_md_typ!=2){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("You have requested the staging PIMD method       \n");
    PRINTF("with an adiabaticity parameter greater than one. \n");
    PRINTF("Staging is for sampling\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    EXIT(1);
  }//endif
  /*-----------------------------------------------------------------------*/
  /* 4)\pi_beads_level_full{#} */
  sscanf(dict[4].keyarg,"%lg",&real_key_arg);
  mdclatoms_pimd->pi_beads_full_ter = (int)(real_key_arg);
  index=4;
  if(mdclatoms_pimd->pi_beads_full_ter<1||
      mdclatoms_pimd->pi_beads_full_ter>mdclatoms_info->pi_beads)
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if(dict[4].iuset==0){
    mdclatoms_pimd->pi_beads_full_ter=mdclatoms_info->pi_beads;
  }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 5)\pi_beads_level_inter_short{#} */
  sscanf(dict[5].keyarg,"%lg",&real_key_arg);
  mdclatoms_pimd->pi_beads_res_ter = (int)(real_key_arg);
  index=5;
  if(mdclatoms_pimd->pi_beads_res_ter<0||
      mdclatoms_pimd->pi_beads_res_ter>mdclatoms_info->pi_beads)
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 6)\pi_beads_level_intra_res{#} */
  sscanf(dict[6].keyarg,"%lg",&real_key_arg);
  mdclatoms_pimd->pi_beads_res_tra = (int)(real_key_arg);
  index=6;
  if(mdclatoms_pimd->pi_beads_res_tra<0||
      mdclatoms_pimd->pi_beads_res_tra>mdclatoms_info->pi_beads)
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 7)\pi_beads_level_intra{#} */
  sscanf(dict[7].keyarg,"%lg",&real_key_arg);
  mdclatoms_pimd->pi_beads_full_tra = (int)(real_key_arg);
  index=7;
  if(mdclatoms_pimd->pi_beads_full_tra<0||
      mdclatoms_pimd->pi_beads_full_tra>mdclatoms_info->pi_beads)
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 8)\respa_steps_pimd{#} */
  sscanf(dict[8].keyarg,"%lg",&real_key_arg);
  gentimeinfo->nres_pimd = (int)(real_key_arg);
  index=8;
  if(gentimeinfo->nres_pimd<0)
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if(gentimeinfo->nres_pimd==0){gentimeinfo->nres_pimd=1;}
  if(gentimeinfo->nres_pimd>1){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("You have requested nres_pimd=%d>1 \n",gentimeinfo->nres_pimd);
    PRINTF("This will just be ingnored for now\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }//endif
  gentimeinfo->nres_pimd = 1;
  /*-----------------------------------------------------------------------*/ 
  /* 9)\initial_spread_size{} */
  sscanf(dict[9].keyarg,"%lg",&real_key_arg);
  mdclatoms_pimd->rcut_spread = (real_key_arg)/BOHR;
  index=9;
  if(mdclatoms_pimd->rcut_spread<0)
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if(mdclatoms_pimd->rcut_spread<1.0){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("You have requested an initial spread size < 1 Bohr \n");
    PRINTF("Are you certain this is what you would like to do?\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/
  /*-----------------------------------------------------------------------*/
  /* 10)\initial_spread_opt{on,off} */
  ifound      = 0;
  if(strcasecmp(dict[10].keyarg,"on")==0) {
    gensimopts->initial_spread_opt = 1; ifound++;}
  if(strcasecmp(dict[10].keyarg,"off")==0){
    gensimopts->initial_spread_opt = 0; ifound++;}
  index=10;
  if(ifound == 0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if(gensimopts->initial_spread_opt == 1 &&   
      gensimopts->pimd == 0){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("You have requested the initial spread option \n");
    PRINTF("\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/
  /*-----------------------------------------------------------------------*/
  /*  11)\pimd_freeze_type{centroid,all_mode} */
  ifound      = 0;
  if(strcasecmp(dict[11].keyarg,"centroid")==0) {
    mdclatoms_pimd->pimd_freez_typ = 1; ifound++;}
  if(strcasecmp(dict[11].keyarg,"all_mode")==0){
    mdclatoms_pimd->pimd_freez_typ = 2; ifound++;}
  mdconstrnt->pimd_freez_typ = mdclatoms_pimd->pimd_freez_typ;
  index=11;
  if(ifound == 0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 12)\pi_temperature{#} */
  sscanf(dict[12].keyarg,"%lg",&real_key_arg);
  mdclatoms_pimd->pi_temperature = real_key_arg;
  index=12;
  if(mdclatoms_pimd->pi_temperature < 0) 
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if(mdclatoms_pimd->pi_temperature == 0.0)
    mdclatoms_pimd->pi_temperature = genstatepoint->t_ext;
  if(mdclatoms_pimd->pi_temperature != genstatepoint->t_ext){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("You have requested a path integral temperature   \n");
    PRINTF("not equal to the external temperature!!          \n");
    PRINTF("Are you CERTAIN this is what you would like to do?\n");
    PRINTF("Path Int temperature %.6g  External temperature %.6g\n",
        mdclatoms_pimd->pi_temperature,genstatepoint->t_ext);
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/
  /*========================================================================*/
}/*end routine*/ 
/*========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  set_sim_params_finale: Consistency checks                   */
/*==========================================================================*/

void set_sim_params_finale(MDINTEGRATE *mdintegrate, MDATOMS *mdatoms,
    MDINTER *mdinter,GENERAL_DATA *general_data,
    MDINTRA *mdintra,
    CP *cp,CLASS_PARSE *class_parse,CP_PARSE *cp_parse,
    FILENAME_PARSE *filename_parse)

  /*=======================================================================*/
  /*          Begin routine                            */
{/*begin routine */
  /*=======================================================================*/
  /*          Local variable declarations                    */

#include "../class_defs/allclass_strip_mdatoms.h"
#include "../class_defs/allclass_strip_mdintegrate.h"
#include "../class_defs/allclass_strip_mdinter.h"
#include "../class_defs/allclass_strip_mdintra.h"
#include "../class_defs/allclass_strip_cp.h"
#include "../class_defs/allclass_strip_gen.h"

  int p1,p2,p3,p4,i1,i2,i3,i4,psum,ip,iii;

  int cp_any_on,cp_on,cp_min_on;
  int pimd_on,cppimd_on,iextended_on;

  int iperd    = gencell->iperd;

  iextended_on = genensopts->nvt+genensopts->npt_i
    + genensopts->npt_f;
  cp_min_on = gensimopts->cp_min 
    + gensimopts->cp_wave_min 
    + gensimopts->cp_wave_min_pimd;
  cp_on     = gensimopts->cp
    + gensimopts->cp_wave
    + gensimopts->cp_pimd 
    + gensimopts->cp_wave_pimd
    + gensimopts->debug_cp
    + gensimopts->debug_cp_pimd;
  pimd_on  = gensimopts->pimd + gensimopts->cp_pimd 
    + gensimopts->cp_wave_pimd 
    + gensimopts->debug_pimd 
    + gensimopts->debug_cp_pimd
    + gensimopts->cp_wave_min_pimd;
  cppimd_on    = gensimopts->cp_pimd
    + gensimopts->cp_wave_pimd
    + gensimopts->cp_wave_min_pimd;
  cp_any_on    = cp_min_on + cp_on;

  /*========================================================================*/
  /* 0) Store the useful flags in many places                   */

  cpopts->cp_any_on        = cp_any_on;
  cpcoeffs_info->cp_any_on    = cp_any_on;
  cptherm_info->cp_any_on     = cp_any_on;
  cpewald->cp_any_on       = cp_any_on;
  cppseudo->cp_any_on      = cp_any_on;
  gensimopts->cp_any_on    = cp_any_on;
  gencell->cp_any_on       = cp_any_on;
  genewald->cp_any_on      = cp_any_on;

  mdclatoms_pimd->pimd_on     = pimd_on;
  mdclatoms_info->pimd_on     = pimd_on;

  mdtherm_info->iextended_on  = iextended_on;

  mdtherm_info_bead->my_type  = 1;
  if(pimd_on==1){
    mdtherm_info_bead->iextended_on  = iextended_on;
  }else{
    mdtherm_info_bead->iextended_on  = 0;
  }/*endif*/

  /*========================================================================*/
  /* I) General Consistency Checks                        */

  if(genensopts->nve==1){
    mdtherm_info_bead->len_nhc = 0;
    mdtherm_info->len_nhc = 0;
    mdtherm_info_bead->num_nhc = 0;
    mdtherm_info->num_nhc = 0;
  }

  if((mdclatoms_info->pi_beads > 1)&&(pimd_on == 0)){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("The number of path integral beads is > 1, but \n");
    PRINTF("you have not turned on a path integral simulation\n");
    PRINTF("type.  Therefore, this run will terminate, now.\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }/*endif*/

  if((genensopts->nve==1)&&(pimd_on == 1)){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("You are doing a path integral calculation .  The \n");
    PRINTF("microcanonical ensemble (NVE) is meaningless\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }/*endif*/

  if(pimd_on==0&&gentimeinfo->ix_respa==5){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("The xi_respa option must be 1,2,3, or 4 under \n");
    PRINTF("classical MD \n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }

  if(class_parse->istart>2&&gensimopts->initial_spread_opt==1){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Spread option is supported only for restart_pos and\n");
    PRINTF("initial restart types\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }

  if(mdinteract->dielectric_opt==1 && gencell->iperd != 0){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Position dependent dielectric constant only    \n");
    PRINTF("supported under cluster boundary conditions    \n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }


  /*========================================================================*/
  /* II) Respa consistency checks                         */

  if(gentimeinfo->int_res_ter==1 && 
      gentimeinfo->int_res_tra == 0){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Having inter molecular respa on and intra molecular respa off\n");
    PRINTF("doesn't makes sense.\n");
    PRINTF("The program therefore insists that if int_res_ter=1\n");
    PRINTF("then int_res_tra=1\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }
  if(gentimeinfo->int_res_ter==1 &&
      gentimeinfo->int_res_tor == 0){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Having inter molecular respa on and torsion respa off\n");
    PRINTF("doesn'tmakes sense. \n");
    PRINTF("The program therefore insists that if int_res_ter=1\n");
    PRINTF("then int_res_tor=1\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }
  if(gentimeinfo->int_res_tor==1 && 
      gentimeinfo->int_res_tra == 0){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Having torsional respa on and intra molecular respa off\n");
    PRINTF("doesn't makes sense.\n");
    PRINTF("The program therefore insists that if int_res_tor=1\n");
    PRINTF("then int_res_tra=1\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }

  /*========================================================================*/
  /* II.V) Thermostat Consistency Checks                     */

  if((mdtherm_info->therm_typ == 2) && 
      gensimopts->anneal_opt == 1){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Annealing with GGMT thermostats not available\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }

  if((mdtherm_info_bead->therm_typ == 2) && 
      gensimopts->anneal_opt == 1){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Annealing with GGMT thermostats not available\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }


  /*========================================================================*/
  /* III) CP consistency checks                           */

  if(gensimopts->cp==1 || gensimopts->cp_wave==1 || 
      gensimopts->cp_min == 1 || 
      gensimopts->cp_wave_min==1){
    if(cpopts->cp_lda == 1){
      if(strcasecmp(cppseudo->vxc_typ,"pz_lsda")==0){
        PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("If you are doing a cp calculation using lda you\n");
        PRINTF("cannot use an lsda vxc.\n"); 
        PRINTF("you have selected vxc type %s\n",cppseudo->vxc_typ);
        PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        FFLUSH(stdout);
        EXIT(1);
      } /*endif*/
    }/*endif*/

    if(gensimopts->cp_min==1 || 
        gensimopts->cp_wave_min==1 || 
        gensimopts->cp_wave_min_pimd==1){
      if(cpopts->cp_ptens_calc == 1){
        PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
        PRINTF("There really is not much use in computing the\n");
        PRINTF("pressure tensor under CP minimization.  I will\n"); 
        PRINTF("allow you to continue, but you might want to turn off\n");
        PRINTF("the cp_ptens option.\n");
        PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
        FFLUSH(stdout);
      } /*endif*/
    }/*endif*/

    if(cpopts->cp_lsda==1){
      if(strcasecmp(cppseudo->vxc_typ,"pz_lda")==0 ||
          strcasecmp(cppseudo->vxc_typ,"pw_lda")==0){
        PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("If you are doing a cp calculation using lsda you\n");
        PRINTF("cannot use an lda vxc, you have selected vxc type %s\n",
            cppseudo->vxc_typ);
        PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        FFLUSH(stdout);
        EXIT(1);
      }/*endif*/
    }/*endif*/

    if(cpopts->cp_sic==1 && cpopts->cp_norb > 0){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("If you are doing a CP calculation using sic\n");
      PRINTF("you cannot use norb integration\n");
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }/*endif*/
  }/*endif*/

  if(cp_min_on == 1 && cpopts->cp_norb > 0) {
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Norb minimization not supported yet.  Sorry, but\n");
    PRINTF("technical support is overworked and underpaid.  \n");
    PRINTF("See Dawn, she'll give you the references!!\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }/*endif*/

  if( (cp_parse->istart_cp == 0) && (cp_on==1)){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Dude, you are doing a full CP calculation with a \n");
    PRINTF("GEN_WAVE restart?????? NO WAY!!!\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }/*enddif*/

  if(cp_on ==1 && cpopts->cp_ptens_calc != 0){
    if(gencell->iperd == 0 || gencell->iperd == 4){ 
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("Dude, you are doing a CP calculation with cluster \n");
      PRINTF("boundary conditions. Why are you calculating the\n");
      PRINTF("the pressure tensor as well?\n");
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }/*endif*/
  }/*endif*/

  if(cpopts->ks_rot_on == 1){
    genfilenames->iwrite_kseigs=cpopts->n_ks_rot;
  } else {
    genfilenames->iwrite_kseigs=gentimeinfo->ntime+1;
  }/* endif */

  if( cpopts->cp_elf_calc_frq > 0){
    genfilenames->iwrite_elf = cpopts->cp_elf_calc_frq;
  } else {
    genfilenames->iwrite_elf = gentimeinfo->ntime+1;
  }/* endif */

  if( (cpopts->cp_elf_calc_frq > 0)&&(cp_min_on>0)){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("NO ELFing during cp-minimization \n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }/*endif*/

  if( (cpopts->cp_elf_calc_frq > 0)&& 
      (cpopts->cp_dual_grid_opt >= 1)){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("NO ELFing with dual gridding without debugging! \n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }/*endif*/

  if( (cpopts->icheck_perd_size != 1)  && (cp_on + cp_min_on > 0)) {
    if(gencell->iperd <= 2 || gencell->iperd == 4){ 
      PRINTF("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");    
      PRINTF("You are doing a CP calculation with reduced boundary\n");
      PRINTF("conditions (< 3) and have elected NOT to check the\n"); 
      PRINTF("inter-particle distances for instances of \n"); 
      PRINTF("particles escaping the box. \n"); 
      PRINTF("ARE YOU SURE THIS IS WHAT YOU WANT TO DO?????? \n"); 
      PRINTF("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");    
    }/*endif*/
  }/*endif*/

  if( (gencell->imov_cp_box == 1)  && (cp_on + cp_min_on > 0)) {
    if(cpopts->cp_dual_grid_opt < 2){ 
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("If you are doing a CP calculation with the small box \n");
      PRINTF("moving, you must use the not_prop option for the  \n");
      PRINTF("cp_dual_grid_opt. Good bye, and have a nice day \n");
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }/*endif*/
  }/*endif*/

  if( (cp_on + cp_min_on > 0)) {
    if(cpopts->cp_dual_grid_opt > 2 && cpopts->cp_ptens_calc == 1){ 
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("If you would like to do a dual grid NPT CP calculation \n");
      PRINTF("you will have to code up the pressure tensor \n");
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }/*endif*/
  }/*endif*/

  if((cp_on || cp_min_on) && (mdenergy_ctrl->iget_pe_real_inter_freq != 1)){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("For CP calculations, interatomic PE should be \n");
    PRINTF("Calculated every step.  Setting the frequency to 1\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n"); 
    mdenergy_ctrl->iget_pe_real_inter_freq = 1;   
  }

  if(cpopts->cp_dual_grid_opt >= 1 && cpopts->cp_hess_calc != 0){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Dude, if you would like to do a dual grid with the hessian \n");
    PRINTF("you will have to code it up \n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }/*endif*/

  if(cpopts->cp_hess_calc == 1  && cpopts->cp_ptens_calc != 0){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Simulataneous calculation of pressure tensor and \n");
    PRINTF("atomic Hessian under CP not currently available\n");
    PRINTF("If you like REALLY need this option, see support staff\n");
    PRINTF("and they will allocate the separate memory space\n");
    PRINTF("required for this.\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }/*endif*/

  if( (cp_on + cp_min_on > 1) && (mdinteract->dielectric_opt == 1)){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Dude, the dielectric option must be off for all  \n");
    PRINTF("CP runs \n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }/*endif*/

  if( (cp_on + cp_min_on > 1)){
    if(cpatom_pme->pme_on  == 1 && 
        cpopts->cp_dual_grid_opt != 2){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("Dude, you may only pme when you dual grid the \n");
      PRINTF("not proportional way \n");
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }/*endif*/
    if(cpatom_pme->pme_on  == 0 && 
        cpopts->cp_dual_grid_opt == 2){
      PRINTF("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");    
      PRINTF("Dude, why aren't you using pme when you are dual gridding\n");
      PRINTF("the not proportional way. It really is a time saver. \n");
      PRINTF("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");    
      FFLUSH(stdout);
    }/*endif*/

  }/*endif*/

  if( (cpopts->cp_dual_grid_opt >= 1) && 
      (gensimopts->cp_min== 1) &&
      (genminopts->min_atm_com_fix_opt==0)){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");
    PRINTF("Dual gridding is great BUT you might want to fix\n");    
    PRINTF("the com of the system using the min_atm_com_fix option\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");   
  }/*endif*/

  if((cpopts->cp_isok_opt == 1) && (cpopts->cp_norb == 0)){
    PRINTF("@@@@@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Sorry, but if you want to use the CP isokinetic option,\n");
    PRINTF("you need to turn on one of the nonorthogonal (cp_norb)\n");
    PRINTF("options.  Try `full_ortho' if you are not sure what to do\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }/* endif */


  /*========================================================================*/
  /* Reduced Periodicity Warnings */

  if( (iperd!=3) && (cp_on + cp_min_on > 0) ){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("Danger! Periodicity is not 3! The convergence of the\n");
    PRINTF("electronic structure with box size in the non-periodic\n");
    PRINTF("directions must be tested!\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");    
  }/*endif*/

  if( (iperd>0) && (iperd!=3) ){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("Danger! Periodicity is not 3 or 0-> The box side, L,\n");
    PRINTF("in the non-periodic direction(s) must be in the range\n");
    PRINTF("L>D where D is the relevent distance parameter, \n");
    PRINTF("clusters, D=4R, surfaces D=2T, wires, D=2l. Here,\n");
    PRINTF("R=cluster radius, T=surface thickness, l=wire length\n");
    PRINTF("formed by the atomic positions \n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");    
  }/*endif*/

  /*========================================================================*/
  /* IV) CP Ensemble checks                            */

  if(gensimopts->cp_wave==1||
      gensimopts->cp_wave_pimd == 1){
    if((genensopts->npt_i + genensopts->npt_f == 1)){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("If you are doing a CP calculation moving only the \n");
      PRINTF("wave functions, constant pressure cannnot be employed \n");
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }
  }

  if(cpopts->cp_ptens_calc==0 && cp_on ==1){
    if((genensopts->npt_i + genensopts->npt_f == 1)){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("If you are doing a CP calculation under constant pressure \n");
      PRINTF("You must turn on the CP pressure tensor calculation \n");
      PRINTF("(cp_ptens{on} \n");
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }
  }


  /*========================================================================*/
  /* V) Velocity sampling checks   */

  mdvel_samp->ivel_smpl_on=0;
  if(mdvel_samp->nvx_smpl>0){mdvel_samp->ivel_smpl_on=1;}
  if(mdvel_samp->nvnhc_smpl>0){mdvel_samp->ivel_smpl_on=1;}

  mdvel_samp->ivel_scale_on=0;
  if(mdvel_samp->nvx_scale>0){mdvel_samp->ivel_scale_on=1;}

  cpvel_samp->ivelc_smpl_on=0;
  if(cpvel_samp->nvc_smpl>0){cpvel_samp->ivelc_smpl_on=1;}
  if(cpvel_samp->nvcnhc_smpl>0){cpvel_samp->ivelc_smpl_on=1;}

  cpvel_samp->ivelc_scal_on=0;
  if(cpvel_samp->nvc_scal>0){cpvel_samp->ivelc_scal_on=1;}

  /*========================================================================*/
  /* V) Ensemble checks */

  if((gencell->iperd<2) || (gencell->iperd==4) ){
    if((genensopts->npt_f+genensopts->npt_i)==1){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("Constant Pressure is not implemented for\n");
      PRINTF("periodicity in less than 2 dimensions   \n");
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }/*endif*/
  }/*endif*/

  if((gensimopts->ann_rate!=1.0)&&
      ((genensopts->npt_f+genensopts->npt_i)==1)){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Constant Pressure is not implemented using\n");
    PRINTF("the temperature annealing option\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }/*endif*/

  if(genensopts->npt_f == 1){
    if(gencell->hmat_cons_typ==0){
      PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
      PRINTF("You have requested flexible constant pressure. \n");
      PRINTF("This only is appropriate for 3D solid systems.\n");
      PRINTF("For example, membranes won't work but protein crystals will.\n");
      PRINTF("Are you certain this is what you would like to do?\n");
      PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    }/*endif*/
    if(gencell->hmat_cons_typ==2){
      PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
      PRINTF("You have requested flexible constant pressure\n");
      PRINTF("with the mono_clinic constraint. For example,\n");
      PRINTF("this ensemble will be stable for some membranes systems\n");
      PRINTF("but not others. It is not a good idea for liquids.\n");
      PRINTF("Are you certain this is what you would like to do?\n");
      PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    }/*endif*/
    if(gencell->hmat_cons_typ==1){
      PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
      PRINTF("You have requested flexible constant pressure\n");
      PRINTF("with the orthorhombic constraint. This ensemble will be\n");
      PRINTF("stable for membranes but not liquids.\n");
      PRINTF("Are you certain this is what you would like to do?\n");
      PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    }/*endif*/
  }/*endif*/

  /*========================================================================*/
  /* VI) Path integral consistency checks */

  if((pimd_on==1)&&(gensimopts->initial_spread_opt==1)){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Initial spread option to be performed off line now\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }/*endif*/

  if((pimd_on==1)&&(gensimopts->pi_beads<=1)){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("You have chosen the Path Integral option with  \n");
    PRINTF("only one bead.                     \n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
    FFLUSH(stdout);
  }/*endif*/

#ifdef NO_PATH_INT_ANNEALING
  if((pimd_on==1)&&(gensimopts->ann_rate!=1.0)){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Temperature annealing with path integrals not allowed.\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }/*endif*/
#endif

  /*========================================================================*/
  /* VII) Bead RESPA checks*/

  p4 = mdclatoms_pimd->pi_beads_full_ter;
  p3 = mdclatoms_pimd->pi_beads_res_ter+1;
  p2 = mdclatoms_pimd->pi_beads_full_tra+1;
  p1 = mdclatoms_pimd->pi_beads_res_tra+1;
  psum = p4*p3*p2*p1;
  if(cppimd_on == 0){
    if(psum != mdclatoms_info->pi_beads){ 
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("Your have improperly set up your RESPA beads   \n");
      PRINTF("for a Path Integral MD simulation.  The product of \n");
      PRINTF("bead_level parameters %d  must equal pi_beads %d. \n",psum,
          mdclatoms_info->pi_beads);
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }/*endif*/
  }/*endif*/
  mdclatoms_pimd->ip_lab = (int *)
    cmalloc(mdclatoms_info->pi_beads*sizeof(int),"set_sim_params_finale")-1;
  ip=0;
  mdclatoms_pimd->pi_beads_full_ter_wght=0.0;
  mdclatoms_pimd->pi_beads_res_ter_wght=0.0;
  mdclatoms_pimd->pi_beads_full_tra_wght=0.0;
  for(i4=1;i4<=p4;i4++){
    for(i3=1;i3<=p3;i3++){
      for(i2=1;i2<=p2;i2++){
        for(i1=1;i1<=p1;i1++){
          ip++;
          if((i3==1)&&(i2==1)&&(i1==1)){mdclatoms_pimd->ip_lab[ip]=4;}
          if((i3!=1)&&(i2==1)&&(i1==1)){mdclatoms_pimd->ip_lab[ip]=3;}
          if((i3!=1)&&(i2!=1)&&(i1==1)){mdclatoms_pimd->ip_lab[ip]=2;}
          if((i3!=1)&&(i2!=1)&&(i1!=1)){mdclatoms_pimd->ip_lab[ip]=1;}
          if(mdclatoms_pimd->ip_lab[ip]==4){
            mdclatoms_pimd->pi_beads_full_ter_wght+=1.0;
            mdclatoms_pimd->pi_beads_res_ter_wght+=1.0;
            mdclatoms_pimd->pi_beads_full_tra_wght+=1.0;
          }/*endif*/
          if(mdclatoms_pimd->ip_lab[ip]==3){
            mdclatoms_pimd->pi_beads_res_ter_wght+=1.0;
            mdclatoms_pimd->pi_beads_full_tra_wght+=1.0;
          }/*endif*/
          if(mdclatoms_pimd->ip_lab[ip]==2){
            mdclatoms_pimd->pi_beads_full_tra_wght+=1.0;
          }/*endif*/
        }/*endfor*/
      }/*endfor*/
    }/*endfor*/
  }/*endfor*/
  mdclatoms_pimd->pi_beads_full_ter_wght/=
    ((double)mdclatoms_info->pi_beads);
  mdclatoms_pimd->pi_beads_res_ter_wght/=
    ((double)mdclatoms_info->pi_beads);
  mdclatoms_pimd->pi_beads_full_tra_wght/=
    ((double)mdclatoms_info->pi_beads);
  /*========================================================================*/
  /* VIII) Parallel checks */

  if(class_parse->istart!=4){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("Restart types other than restart_all employed\n");
    PRINTF("in parallel may effect comparisons with simulations \n");
    PRINTF("performed in scalar or with different numbers of procs\n");
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
  }/*endif*/
  if((cp_on+cp_min_on)==1){
    if(cp_parse->istart_cp!=4){
      PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
      PRINTF("CP-Restart types other than restart_all employed\n");
      PRINTF("in parallel may effect comparisons with simulations \n");
      PRINTF("performed in scalar or with different numbers of procs\n");
      PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    }/*endif*/
  }/*endif*/

  /*========================================================================*/
  /* IX) PME checks                                 */

  if(gencell->iperd > 0 && cp_on == 0){

    if((gentimeinfo->int_res_ter==0) &&
        (mdpart_mesh->pme_res_on == 1)){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("PME respa option requires long range forces respa \n");
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }/*endif*/

    if(mdpart_mesh->pme_res_on == 1 && mdpart_mesh->pme_on == 0){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("No PME respa option without the PME option \n");
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }/*endif*/
    if((gentimeinfo->int_res_ter==1) &&
        (mdpart_mesh->pme_res_on == 1)){ 
      if((mdpart_mesh->n_interp_res  > mdpart_mesh->n_interp)){
        PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("The RESPA PME n_interp parameter > PME n_interp \n");
        PRINTF("This is not allowed\n");      
        PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
        FFLUSH(stdout);
        EXIT(1);      
      }/*endif*/
      if((mdpart_mesh->kmax_pme  < mdpart_mesh->kmax_pme_res)){
        PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("The RESPA PME kmax parameter > PME kmax parameter\n");
        PRINTF("This is not allowed\n");      
        PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
        FFLUSH(stdout);
        EXIT(1);
      }/*endif*/
    }/*endif*/

    if((gentimeinfo->int_res_ter==1) &&
        (mdpart_mesh->pme_res_on == 0) && (mdpart_mesh->pme_on == 1) ){
      if(class_parse->kmax_res==class_parse->kmax_ewd){
        PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
        PRINTF("Inter respa with PME option on but PME respa opt is off.\n");
        PRINTF("Placing full reciprocal space calculation \n");
        PRINTF("inside the short range respa loop.\n");
        PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
      }else{
        if(class_parse->kmax_res!=0){
          PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
          PRINTF("PME option on but PME respa opt is off under LRF respa.\n");
          PRINTF("This is only OK if kmax_res=kmax or kmax_res=0\n");
          PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
          FFLUSH(stdout);
          EXIT(1);
        }/*endif*/
      }/*endif*/
    }/*endif*/

    if(mdpart_mesh->pme_res_on == 1 && mdpart_mesh->pme_on == 1){
      if(class_parse->kmax_res==class_parse->kmax_ewd){
        if((mdpart_mesh->kmax_pme_res == mdpart_mesh->kmax_pme) &&
            (mdpart_mesh->n_interp_res == mdpart_mesh->n_interp)){
          PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
          PRINTF("PME Respa identical to full PME\n");
          PRINTF("Placing full reciprocal space calculation \n");
          PRINTF("inside the short range respa loop.\n");
          PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
          mdpart_mesh->pme_res_on = 0; 
        }/*endif*/
      }/*endif*/
    }/*endif*/

  }/*endif:iperd > 0*/

  /*========================================================================*/
  /* Reset flags because charmness of being is easier that way     */

  if(pimd_on==1){
    gensimopts->cp_wave_min   = gensimopts->cp_wave_min_pimd; 
    gensimopts->cp            = gensimopts->cp_pimd;
    gensimopts->cp_wave       = gensimopts->cp_wave_pimd;
    gensimopts->cp_wave_min_pimd = 0; 
    gensimopts->cp_pimd          = 0;
    gensimopts->cp_wave_pimd     = 0;
  }/*endif*/

  /*========================================================================*/
}/*end routine*/ 
/*========================================================================*/
