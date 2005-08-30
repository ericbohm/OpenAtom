/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: control_intra_parms.c                        */
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
#include "../proto_defs/proto_intra_params_entry.h"
#include "../proto_defs/proto_intra_params_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_math.h"

#define DEBUG_OFF

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void control_intra_params(double *tot_memory,MDCLATOMS_INFO *clatoms_info,
                          MDCLATOMS_PIMD *clatoms_pimd,
                          CPATOM_MAPS *cpatom_maps,
                          MDATOM_MAPS *atommaps,MDINTRA *mdintra,
                          FILENAME_PARSE *filename_parse,
                          FREE_PARSE *free_parse,CLASS_PARSE *class_parse,
                          NULL_INTER_PARSE *null_inter_parse,
                          GENSIMOPTS *simopts, int isurf_on)

/*========================================================================*/
/*     Begin routine                                                      */
{/*begin routine*/

/*========================================================================*/
/*     Local Variables                                                    */

  DICT_INTRA dict_intra;
  BUILD_INTRA build_intra;
  RESBOND_PARSE resbond_parse;
  START_INDEX start_index;

#include "../class_defs/allclass_strip_mdintra.h"


  int i,iii;
  int jmol_typ;
  int ifirst,nresidue,mol_or_res;
  int nmol_tot,nres_bond;
  int nres_tot,ncon_tot,nmass_unphys;
  double cpu1,cpu2;

  char *filename,*fun_key;

  int pi_beads = clatoms_info->pi_beads;

/*=======================================================================*/
/*   I) Start the routine                                                */

  PRINTF("\n");
  PRINT_LINE_STAR;
  PRINTF("Reading molecular parameter files\n");
  PRINT_LINE_DASH;PRINTF("\n");
  cputime(&cpu1);


/*=======================================================================*/
/* II) Initialize/malloc intra stuff                                     */
/*      (init_intra_params.c)                                            */

  filename = (char *)cmalloc(MAXWORD*sizeof(char),"control_intr_params");  
  fun_key  = (char *)cmalloc(MAXWORD*sizeof(char),"control_intr_params");  

  init_intra_params(clatoms_info,cpatom_maps,atommaps,&build_intra,
                    mdintra,null_inter_parse,&resbond_parse,
                    filename_parse);

/*=======================================================================*/
/* III) Set up the dictionaries                                          */
/*      (set_intra_dict.c)                                               */

  ifirst = 1;
  dict_intra.word = (DICT_WORD *)cmalloc(sizeof(DICT_WORD),"control_intr_params");
  set_intra_fun_dict(&dict_intra.fun_dict,&dict_intra.num_fun_dict,ifirst);
  set_atm_dict(&dict_intra.atm_dict,&dict_intra.num_atm_dict,ifirst);
  set_intra_dict(&dict_intra.intra_dict,&dict_intra.num_intra_dict,ifirst);
  set_mol_name_dict(&dict_intra.mol_name_dict,&dict_intra.num_mol_name_dict,
                ifirst);
  set_res_name_dict(&dict_intra.res_name_dict,&dict_intra.num_res_name_dict,
                ifirst);
  set_res_def_dict(&dict_intra.res_def_dict,&dict_intra.num_res_def_dict,
                ifirst);
  set_res_bond_dict(&dict_intra.res_bond_dict,&dict_intra.num_res_bond_dict,
                ifirst);

  strcpy(atommaps->atm_typ[1],"HELP");

/*=======================================================================*/
/* V) Zero the list counters                                             */
 
   mdbond->npow             = 0;
   mdbond->ncon             = 0;
   null_inter_parse->nbond_nul   = 0;
   mdbend->npow             = 0;
   null_inter_parse->nbend_nul   = 0;
   mdtors->npow             = 0;
   mdtors->nimpr            = 0;
   null_inter_parse->ntors_nul   = 0;
   mdonfo->num              = 0;
   null_inter_parse->nonfo_nul   = 0;
   mdbend_bnd->num          = 0;
   clatoms_info->natm_tot        = 0;
   atommaps->nres_typ            = 0;
   mdconstrnt->nfreeze             = 0;
   atommaps->natm_typ            = 0;
   mdghost_atoms->nghost_tot       = 0;
   mdghost_atoms->natm_comp_max    = 0;
   mdgrp_bond_con->num_21   = 0;
   mdgrp_bond_con->num_33   = 0;
   mdgrp_bond_con->num_43   = 0;
   mdgrp_bond_con->num_23   = 0;
   mdgrp_bond_con->num_46   = 0;
   mdgrp_bond_con->ntyp_21  = 0;
   mdgrp_bond_con->ntyp_33  = 0;
   mdgrp_bond_con->ntyp_43  = 0;
   mdgrp_bond_con->ntyp_23  = 0;
   mdgrp_bond_con->ntyp_46  = 0;
   mdgrp_bond_watts->num_33 = 0;
   mdgrp_bond_watts->ntyp_33= 0;

/*=======================================================================*/
/* IV) Loop over molecular parameter files                               */

  for(jmol_typ=1;jmol_typ<=atommaps->nmol_typ;jmol_typ++){

/*-----------------------------------------------------------------------*/
/*  0) Write to the screen                                               */

    PRINTF("**************************************************************\n");
    PRINTF("Reading from molecular parameter file %s\n",
          filename_parse->mol_param_name[jmol_typ]);
    PRINTF("--------------------------------------------------------------\n");
    PRINTF("\n");

/*-----------------------------------------------------------------------*/
/* 1) Store the present list counter values                              */

    start_index.nbond_pow    = mdbond->npow;
    start_index.nbond_con    = mdbond->ncon;
    start_index.nbond_nul    = null_inter_parse->nbond_nul;
    start_index.nbend_pow    = mdbend->npow;
    start_index.nbend_nul    = null_inter_parse->nbend_nul;
    start_index.ntors_pow    = mdtors->npow;
    start_index.ntors_nul    = null_inter_parse->ntors_nul;
    start_index.nonfo        = mdonfo->num;
    start_index.nonfo_nul    = null_inter_parse->nonfo_nul;
    start_index.nbend_bnd    = mdbend_bnd->num;
    start_index.natm         = clatoms_info->natm_tot;
    start_index.nfreeze      = mdconstrnt->nfreeze;
    start_index.nghost_tot   = mdghost_atoms->nghost_tot;
    start_index.ngrp_21      = mdgrp_bond_con->num_21;
    start_index.ngrp_33      = mdgrp_bond_con->num_33;
    start_index.ngrp_43      = mdgrp_bond_con->num_43;
    start_index.ngrp_23      = mdgrp_bond_con->num_23;
    start_index.ngrp_46      = mdgrp_bond_con->num_46;
    start_index.ngrp_watt_33 = mdgrp_bond_watts->num_33;
   
/*------------------------------------------------------------------------*/
/*  2) Count and error check the molecular parm file:                     */
/*     (fetch_residue.c)                                                  */

    strcpy(filename,filename_parse->mol_param_name[jmol_typ]);
    nresidue = atommaps->nres_1mol_jmol_typ[jmol_typ];/* spec in mol_set_file*/

    mol_or_res = 1;
    check_parmfile(filename,&(dict_intra.num_fun_dict),&(dict_intra.fun_dict),
                    fun_key,nresidue,&nres_bond,mol_or_res);
                                          /* in fetch_residue.c */
    resbond_parse.nres_bond = nres_bond; 
    resbond_parse.nresidue  = nresidue;
    resbond_parse_realloc(&resbond_parse);
                                          /* in manipulate_res_bonds.c */

/*-----------------------------------------------------------------------*/
/* 3) Read in the molecule name functional keyword: molecule type        */
/*     number of atoms and/or number of residues                         */
/*     (fetch_residue.c)                                                 */


    fetch_molname(filename,&dict_intra,atommaps,
                  fun_key,jmol_typ,nresidue);
/*-----------------------------------------------------------------------*/
/* 4) Read in the molecule residue definitions:                          */
/*       index,natom,parm_file,name                                      */
/*     (fetch_residue.c)                                                 */

    if(nresidue>0){
       fetch_residue_defs(filename,&dict_intra,
                          atommaps,filename_parse,
                          fun_key,nresidue,
                          jmol_typ,&build_intra);

    }else{
      fetch_residue_def0(atommaps,&build_intra,jmol_typ);  
      strcpy(filename_parse->res_param_name[1],
             filename_parse->mol_param_name[jmol_typ]);
        /* assigns the molecule type to the residue type 
           when there are no explicit residues           */
    }/*endif*/

/*-----------------------------------------------------------------------*/
/* 5) Read in the molecule resbonds then map them to the residues:       */
/*     Residue types(pairs),residue indices (pairs), bond sites(pairs),  */
/*     bond files(pairs), bond modifier and bond labels read in          */
/*     (fetch_residue.c and manipulate_res_bonds.c)                      */

    if(nres_bond>0){
       fetch_residue_bonds(filename,&dict_intra,fun_key,
                           resbond_parse.resbond_prm,
                           atommaps,nresidue,
                           nres_bond,jmol_typ,pi_beads); 
                    /* in fetch_residue.c */
       map_residue_bonds(&resbond_parse); /* in manipulate_res_bonds.c */
    }else{
       for(i=1;i<=MAX(nresidue,1);i++){
         resbond_parse.nres_bond1_jres[i]=0;resbond_parse.res_bond1_off[i]=0;
         resbond_parse.nres_bond2_jres[i]=0;resbond_parse.res_bond2_off[i]=0;
       }/*endfor*/  
    }/*endif*/
/*-----------------------------------------------------------------------*/
/* 6) Read in the residues: Set the atoms,bonds,bends,torsions,etc.      */
/*                          of each residue. Modificiations of residues  */
/*                          involved in bonds performed.                 */
/*                          Molecules with no residues treated here also */
/*     (control_res_params.c)                                           */


    control_res_params(tot_memory,clatoms_info,cpatom_maps,
                       atommaps,mdintra,&resbond_parse,&build_intra,
                       filename_parse,free_parse,class_parse,
                       null_inter_parse,filename,&dict_intra,
                       fun_key,jmol_typ);

/*------------------------------------------------------------------------*/
/* 7) Bond the residues together and create intramolecular interactions   */
/*    based on the molecular connectivity                                 */
/*    (residue_bond.c)                                                    */

    if(nres_bond>0){
       resbond_parse.ionfo  = class_parse->ionfo_opt[jmol_typ];
       resbond_parse.iconv  = class_parse->ires_bond_conv[jmol_typ];
       residue_bond(clatoms_info,
                    atommaps,mdintra,&resbond_parse,
                    &build_intra,jmol_typ,
                    class_parse->mol_hydrog_con_opt[jmol_typ]);
    }/*endif*/

/*------------------------------------------------------------------------*/
/* 7.5) Print out progress in intramolecular connectivity lists           */

#ifdef DEBUG
    PRINTF("End Fetch: nbonds=%d,nbends=%d,ntors=%d,nimpr=%d,nonfo=%d\n",
     mdbond->npow,mdbend_bnd->num,mdtors->npow,mdtors->nimpr,mdonfo->num);
#endif

/*------------------------------------------------------------------------*/
/* 7.7) Implement the freeze option                                       */

 fetch_freeze(class_parse,atommaps,mdconstrnt,&build_intra,
              &start_index,clatoms_info,jmol_typ);

/*------------------------------------------------------------------------*/
/* 7.7.5) Check for consistency with zero_com_vel                          */

  if(mdconstrnt->nfreeze > 0 && class_parse->zero_com_vel==1){
    PRINTF("@@@@@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Frozen atoms not compatible with zeroing the center of mass\n");
    PRINTF("If you want to freeze atoms, set the zero_com_vel to `no'\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }/*endif*/

/*------------------------------------------------------------------------*/
/* 7.8) Implement the hydrog_mass option                                  */

  fetch_hydrog_mass(class_parse,atommaps,mdconstrnt,&build_intra,
                    &start_index,clatoms_info,jmol_typ);

/*-----------------------------------------------------------------------*/
/* 8) Replicate the molecule now that it has been constructed            */
/*    (replicate_mol.c)                                                    */

  if(atommaps->nmol_jmol_typ[jmol_typ] > 1){
     replicate_mol(clatoms_info,cpatom_maps,atommaps,&build_intra,mdintra,
                   null_inter_parse,&start_index,jmol_typ);
  }/*endif*/

/*-----------------------------------------------------------------------*/

 }/*endfor:jmol_typ*/

/*=======================================================================*/
/*  Check for frozen atoms involved in constraints */

  freeze_con_check(atommaps,mdbond,mdbend,mdtors,mdgrp_bond_con,mdconstrnt);

/*=======================================================================*/
/*  DEBUG */

#ifdef DEBUG
  PRINTF("Vomitting \n");mal_verify(1);
  vomit_intra_list(clatoms_info,atommaps,mdintra,null_inter_parse);
#endif

/*=======================================================================*/
/*  VI)Tidy up                                                           */

 close_intra_params(clatoms_info,clatoms_pimd,cpatom_maps,atommaps,
                    &build_intra,mdintra,null_inter_parse,tot_memory,
                    simopts);

/*=======================================================================*/
/*  VIII) Output to screen                                                */

  cputime(&cpu2);
  PRINTF("\n"); 
  PRINT_LINE_DASH;
  PRINTF("Completed reading molecular parameter files %g\n",cpu2-cpu1);
  PRINT_LINE_STAR;PRINTF("\n");

/*========================================================================*/
/* Check for physical masses */

  nmass_unphys = 0;
  for(i=1;i<=clatoms_info->natm_tot;i++){
    if(clatoms_info->mass[i]<900.0){nmass_unphys++;}
  }/*endfor*/

#define NONEXPERT
#ifdef NONEXPERT
  if(nmass_unphys>0){
    PRINTF("@@@@@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("There were %d atoms with masses less than 1/2 AMU\n",nmass_unphys);
    PRINTF("This might be OK for experts, but not in general.\n");
    PRINTF("If you are performing cp minimization, use the option\n");
    PRINTF("class_mass_scale_fact in sim_run_def, instead.\n");
    PRINTF("Expert users can disable this error by undefining the\n");
    PRINTF("NONEXPERT ifdef in control_intra_params.c on line 358\n");
    PRINTF("However, this is not recommended as the experts pass\n");
    PRINTF("on damaged parm files to students and postdocs which\n");
    PRINTF("causes collatoral damage. Better to clean your parms \n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }/*endif*/
#endif

/*========================================================================*/
/*  IX) Set the intramolecular potential params                          */

  set_intra_potent(mdintra,&build_intra,
               filename_parse->def_intra_name,
               filename_parse->user_intra_name);
  
#ifdef DEBUG
   PRINTF("Vomitting \n");mal_verify(1);
   vomit_intra_potent(mdintra,&build_intra);
#endif

/*========================================================================*/
/*  X) Free memory                */                                     


  cfree(fun_key,"control_intra_params");
  cfree(filename,"control_intra_params");

  cfree(&dict_intra.atm_dict[1],"control_intra_params");
  cfree(&dict_intra.intra_dict[1],"control_intra_params");
  cfree(&dict_intra.mol_name_dict[1],"control_intra_params");
  cfree(&dict_intra.word[0],"control_intra_params");
  cfree(&dict_intra.fun_dict[1],"control_intra_params");

  cfree(&build_intra.mask_atm[1],"control_intra_params");
  cfree(&build_intra.bond_site[1],"control_intra_params");
  cfree(&build_intra.index_atm[1],"control_intra_params");
  cfree(&build_intra.iatm_ind_chk[1],"control_intra_params");

  cfree(&build_intra.cbond_typ_pow[1],"control_intra_params");
  cfree(&build_intra.cbond_typ_con[1],"control_intra_params");
  cfree(build_intra.cbond_typ_now,"control_intra_params");

  cfree(&build_intra.cbend_typ_pow[1],"control_intra_params");
  cfree(build_intra.cbend_typ_now,"control_intra_params");

  cfree(&build_intra.ctors_typ_pow[1],"control_intra_params");
  cfree(build_intra.ctors_typ_now,"control_intra_params");

  cfree(&build_intra.confo_typ[1],"control_intra_params");
  cfree(build_intra.confo_typ_now,"control_intra_params"); 

  cfree(&build_intra.cbend_bnd_typ[1],"control_intra_params");
  cfree(build_intra.cbend_bnd_typ_now,"control_intra_params");

/*=======================================================================*/
/*  X) Write out synopsis                                               */
  
  PRINT_LINE_STAR;
  PRINTF("System Summary  \n");
  PRINT_LINE_DASH;PRINTF("\n");

  nmol_tot=0;
  nres_tot=0;
  for(i=1;i<=atommaps->nmol_typ;i++){
    nmol_tot+=(atommaps->nmol_jmol_typ)[i];
    nres_tot+=(atommaps->nres_1mol_jmol_typ)[i]*
              (atommaps->nmol_jmol_typ)[i];
  }/*endfor*/
  atommaps->nres_tot = nres_tot;   

  PRINTF("There are %d molecular units      \n",nmol_tot);
  PRINTF("There are %d molecular unit types \n",atommaps->nmol_typ);
  PRINTF("There are %d residue units        \n",nres_tot);
  PRINTF("There are %d residue unit types   \n",atommaps->nres_typ);
  PRINTF("There are %d atoms                \n",clatoms_info->natm_tot);
  PRINTF("There are %d degrees of freedom   \n",clatoms_info->nfree);
  PRINTF("There are %d ghost atoms          \n",mdghost_atoms->nghost_tot);
  PRINTF("There are %d atom types           \n",atommaps->natm_typ);
  PRINTF("There are %d charged atoms        \n",clatoms_info->nchrg);
  PRINTF("There are %d power series bonds   \n",mdbond->npow);
  PRINTF("There are %d constrained bonds    \n",mdbond->ncon);
  PRINTF("There are %d null bonds           \n",null_inter_parse->nbond_nul);
  PRINTF("There are %d power series bends   \n",mdbend->npow);
  PRINTF("There are %d null bends           \n",null_inter_parse->nbend_nul);
  PRINTF("There are %d Urey-Bradley bends   \n",mdbend_bnd->num);
  PRINTF("There are %d power series torsions\n",mdtors->npow
                                               -mdtors->nimpr);
  PRINTF("There are %d improper torsions    \n",mdtors->nimpr);
  PRINTF("There are %d null torsions        \n",null_inter_parse->ntors_nul);
  PRINTF("There are %d lj onefours          \n",mdonfo->num);
  PRINTF("There are %d null onefours        \n",null_inter_parse->nonfo_nul);
  PRINTF("There are %d free energy bonds    \n",mdbond_free->num);
  PRINTF("There are %d free energy bends    \n",mdbend_free->num);
  PRINTF("There are %d free energy torsions \n",mdtors_free->num);
  PRINTF("There are %d free energy rbar-sig \n",mdrbar_sig_free->nfree);
  PRINTF("There are %d 21 group constraints \n",mdgrp_bond_con->num_21);
  PRINTF("There are %d 23 group constraints \n",mdgrp_bond_con->num_23);
  PRINTF("There are %d 33 group constraints \n",mdgrp_bond_con->num_33);
  PRINTF("There are %d 43 group constraints \n",mdgrp_bond_con->num_43);
  PRINTF("There are %d 46 group constraints \n",mdgrp_bond_con->num_46);
  PRINTF("There are %d 33 group Watts       \n",mdgrp_bond_watts->num_33);
  PRINTF("There is  %d surface              \n",isurf_on);
  PRINTF("\n");
  
  PRINT_LINE_DASH;
  PRINTF("End System Summary\n");
  PRINT_LINE_STAR;PRINTF("\n");
  if((simopts->debug+simopts->debug_cp+simopts->debug_pimd)==1){
    PRINTF("Enter an integer ");SCANF("%d",&iii);
  }/*endif*/

/*=======================================================================*/
/*  XI) No atm minimization with constraints                             */

  ncon_tot = mdgrp_bond_con->num_21
           + mdgrp_bond_con->num_23
           + mdgrp_bond_con->num_33
           + mdgrp_bond_con->num_43
           + mdgrp_bond_con->num_46
           + mdbond->ncon;
  if(((simopts->cp_min)==1)&&(ncon_tot>0)){
    PRINTF("@@@@@@@@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Atomic position minimization with constraints under CP \n");
    PRINTF("not implemented\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }/*endif*/ 

/*--------------------------------------------------------------------------*/
  }/*end routine*/ 
/*==========================================================================*/


#ifdef DEBUG

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void  vomit_intra_potent(MDINTRA *mdintra, BUILD_INTRA *build_intra)

/*==========================================================================*/
{ /*begin routine */
/*==========================================================================*/

int i,iii,j;
#include "../class_defs/allclass_strip_mdintra.h"

/*==========================================================================*/
/* ii) Bonds */

PRINTF("------\n");
PRINTF("pow bonds typ: %d\n ",bonded->mdbond.ntyp_pow);
PRINTF("------\n");
   for(i=1;i<=bonded->mdbond.ntyp_pow;i++){
     PRINTF("%d atom1 %s atom2 %s eq %g \n",i,
             build_intra->cbond_typ_pow[i].atm1,
             build_intra->cbond_typ_pow[i].atm2,
             mdbond->eq_pow[i]);
   }/*endfor*/
mal_verify(1);PRINTF("Enter an integer ");SCANF("%d",&iii);

PRINTF("------\n");
PRINTF("con bonds typ: %d\n ",mdbond->ntyp_con);
PRINTF("------\n");
  for(i=1;i<=mdbond->ntyp_con;i++){
     PRINTF("%d atom1 %s atom2 %s eq %g \n",i,
             build_intra->cbond_typ_con[i].atm1,
             build_intra->cbond_typ_con[i].atm2,
             mdbond->eq_con[i]);
  }/*endfor*/
mal_verify(1);PRINTF("Enter an integer ");SCANF("%d",&iii);

/*==========================================================================*/
/* ii) Bends */

PRINTF("------\n");
PRINTF("pow bends typ: %d\n ",mdbend->ntyp_pow);
PRINTF("------\n");
   for(i=1;i<=mdbend->ntyp_pow;i++){
     PRINTF("%d atom1 %s atom2 %s atom3 %s eq %g \n",i,
             build_intra->cbend_typ_pow[i].atm1,
             build_intra->cbend_typ_pow[i].atm2,
             build_intra->cbend_typ_pow[i].atm3,
             mdbend->eq_pow[i]);
   }/*endfor*/
mal_verify(1);PRINTF("Enter an integer ");SCANF("%d",&iii);

PRINTF("------\n");
PRINTF("con bends typ: %d\n ",mdbend->ntyp_con);
PRINTF("------\n");
   for(i=1;i<=mdbend->ntyp_con;i++){
     PRINTF("%d atom1 %s atom2 %s atom3 %s eq %g \n",i,
             build_intra->cbend_typ_con[i].atm1,
             build_intra->cbend_typ_con[i].atm2,
             build_intra->cbend_typ_con[i].atm3,
             mdbend->eq_con[i]);
   }/*endfor*/
mal_verify(1);PRINTF("Enter an integer ");SCANF("%d",&iii);

/*==========================================================================*/
/* ii) Bend_bnds */

PRINTF("------\n");
PRINTF("bend_bnd typ: %d\n",mdbend_bnd->ntyp);
PRINTF("------\n");
   for(i=1;i<=mdbend_bnd->ntyp;i++){
     PRINTF("%d atom1 %s atom2 %s atom3 %s eq_bo %g eq_be %g \n",i,
             build_intra->cbend_bnd_typ[i].atm1,
             build_intra->cbend_bnd_typ[i].atm2,
             build_intra->cbend_bnd_typ[i].atm3,
             mdbend_bnd->eq_bond[i],
             mdbend_bnd->eq_bend[i]);
   }/*endfor*/
mal_verify(1);PRINTF("Enter an integer ");SCANF("%d",&iii);

/*==========================================================================*/
/* ii) Tors */

PRINTF("------\n");
PRINTF("pow torsions typ: %d\n",mdtors->ntyp_pow);
PRINTF("------\n");
   for(i=1;i<=mdtors->ntyp_pow;i++){
     PRINTF("%d atom1 %s atom2 %s atom3 %s atom4 %s eq %g \n",i,
             build_intra->ctors_typ_pow[i].atm1,
             build_intra->ctors_typ_pow[i].atm2,
             build_intra->ctors_typ_pow[i].atm3,
             build_intra->ctors_typ_pow[i].atm4,
             mdtors->eq_con[i]);
   }/*endfor*/
mal_verify(1);PRINTF("Enter an integer ");SCANF("%d",&iii);

PRINTF("------\n");
PRINTF("con torsions: %d\n",mdtors->ncon);
PRINTF("------\n");
   for(i=1;i<=mdtors->ntyp_con;i++){
     PRINTF("%d atom1 %s atom2 %s atom3 %s atom4 %s eq %g \n",i,
             build_intra->ctors_typ_con[i].atm1,
             build_intra->ctors_typ_con[i].atm2,
             build_intra->ctors_typ_con[i].atm3,
             build_intra->ctors_typ_con[i].atm4,
             mdtors->eq_con[i]);
   }/*endfor*/
mal_verify(1);PRINTF("Enter an integer ");SCANF("%d",&iii);

/*==========================================================================*/
/* ii) Onfos */

PRINTF("------\n");
PRINTF("onfos typs: %d\n",mdonfo->ntyp);
PRINTF("------\n");
   for(i=1;i<=mdonfo->ntyp;i++){
     PRINTF("%d atom1 %s atom2 %s feps %g s6 %g \n",i,
             build_intra->confo_typ[i].atm1,
             build_intra->confo_typ[i].atm2,
             mdonfo->feps[i],mdonfo->s6[i]);
   }/*endfor*/
mal_verify(1);PRINTF("Enter an integer ");SCANF("%d",&iii);

/*==========================================================================*/
/* i) Grp con 21 */

PRINTF("------\n");mal_verify(1);
PRINTF("Grp con 21 typ: %d \n",mdgrp_bond_con->ntyp_21);
PRINTF("------\n");
   for(i=1;i<=mdgrp_bond_con->ntyp_21;i++){
     PRINTF("%d atom1 %s atom2 %s eqs %g \n",i,
             build_intra->cgrp_typ_21[i].atm1,
             build_intra->cgrp_typ_21[i].atm2,
             mdgrp_bond_con->eq_21[1][i]);
   }/*endfor*/
mal_verify(1);PRINTF("Enter an integer ");SCANF("%d",&iii);

/*==========================================================================*/
/* i) Grp con 23 */

PRINTF("------\n");mal_verify(1);
PRINTF("Grp con 23 typ: %d \n",mdgrp_bond_con->ntyp_23);
PRINTF("------\n");
   for(i=1;i<=mdgrp_bond_con->ntyp_23;i++){
     PRINTF("%d atom1 %s atom2 %s atom3 %s eqs %g %g\n",i,
             build_intra->cgrp_typ_23[i].atm1,
             build_intra->cgrp_typ_23[i].atm2,
             build_intra->cgrp_typ_23[i].atm3,
             mdgrp_bond_con->eq_23[1][i],
             mdgrp_bond_con->eq_23[2][i]);
   }/*endfor*/
mal_verify(1);PRINTF("Enter an integer ");SCANF("%d",&iii);

/*==========================================================================*/
/* i) Grp con 33 */

PRINTF("------\n");mal_verify(1);
PRINTF("Grp con 33 typ: %d \n",mdgrp_bond_con->ntyp_33);
PRINTF("------\n");
   for(i=1;i<=mdgrp_bond_con->ntyp_33;i++){
     PRINTF("%d atom1 %s atom2 %s atom3 %s eqs %g %g %g\n",i,
             build_intra->cgrp_typ_33[i].atm1,
             build_intra->cgrp_typ_33[i].atm2,
             build_intra->cgrp_typ_33[i].atm3,
             mdgrp_bond_con->eq_33[1][i],
             mdgrp_bond_con->eq_33[2][i],
             mdgrp_bond_con->eq_33[3][i]);
   }/*endfor*/
mal_verify(1);PRINTF("Enter an integer ");SCANF("%d",&iii);

/*==========================================================================*/
/* i) Grp con 43 */

PRINTF("------\n");mal_verify(1);
PRINTF("Grp con 43 typ: %d \n",mdgrp_bond_con->ntyp_43);
PRINTF("------\n");
   for(i=1;i<=mdgrp_bond_con->ntyp_43;i++){
     PRINTF("%d atom1 %s atom2 %s atom3 %s atom4 %s\n",i,
             build_intra->cgrp_typ_43[i].atm1,
             build_intra->cgrp_typ_43[i].atm2,
             build_intra->cgrp_typ_43[i].atm3,
             build_intra->cgrp_typ_43[i].atm4);
     PRINTF("%d eqs %g %g %g \n",i,
             mdgrp_bond_con->eq_43[1][i],
             mdgrp_bond_con->eq_43[2][i],
             mdgrp_bond_con->eq_43[3][i]);
   }/*endfor*/
mal_verify(1);PRINTF("Enter an integer ");SCANF("%d",&iii);

/*==========================================================================*/
/* i) Grp con 46 */

PRINTF("------\n");mal_verify(1);
PRINTF("Grp con 46 typ: %d \n",mdgrp_bond_con->ntyp_46);
PRINTF("------\n");
   for(i=1;i<=mdgrp_bond_con->ntyp_46;i++){
     PRINTF("%d atom1 %s atom2 %s atom3 %s atom4 %s\n",i,
             build_intra->cgrp_typ_46[i].atm1,
             build_intra->cgrp_typ_46[i].atm2,
             build_intra->cgrp_typ_46[i].atm3,
             build_intra->cgrp_typ_46[i].atm4);
     PRINTF("%d eqs %g %g %g %g %g %g \n",i,
             mdgrp_bond_con->eq_46[1][i],
             mdgrp_bond_con->eq_46[2][i],
             mdgrp_bond_con->eq_46[3][i],
             mdgrp_bond_con->eq_46[4][i],
             mdgrp_bond_con->eq_46[5][i],
             mdgrp_bond_con->eq_46[6][i]);
   }/*endfor*/
mal_verify(1);PRINTF("Enter an integer ");SCANF("%d",&iii);

/*-----------------------------------------------------------------------*/
}/*end routine */
/*=======================================================================*/

#endif

