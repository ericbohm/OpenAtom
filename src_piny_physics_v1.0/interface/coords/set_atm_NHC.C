/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: set_atm_NHC.c                                */
/*                                                                          */
/* This subprogram sets up the atm/vol NHC for a MD on a                    */ 
/* LD-classical potential energy surface (LD-PES)                           */
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
#include "../class_defs/typedefs_par.h"
#include "../proto_defs/proto_coords_entry.h"
#include "../proto_defs/proto_friend_lib_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_atm_NHC(GENENSOPTS *genensopts,
    GENSTATEPOINT *genstatepoint, 
    GENSIMOPTS *gensimopts,
    MDCLATOMS_INFO *mdclatoms_info,
    MDCLATOMS_PIMD *mdclatoms_pimd,
    MDGHOST_ATOMS *mdghost_atoms,
    MDATOM_MAPS *mdatom_maps,
    MDCONSTRNT *mdconstrnt,
    MDTHERM_INFO *mdtherm_info_bead,
    MDTHERM_INFO *mdtherm_info_class, 
    MDBARO *mdbaro,
    MDPAR_RAHMAN *mdpar_rahman,
    CLASS_PARSE *class_parse,
    double *tot_memory,int pimd_on,
    int hmat_cons_typ)

  /*==========================================================================*/
  /*               Begin subprogram:                                          */
{/*begin routine*/
  /*==========================================================================*/
  /*               Local variables                                            */

  int ind_glob_nhc;          /* Num: Index of global NHC     */
  int jmol_typ;              /* Num: Molecule type for index */
  int jatm;                  /* Num: Atm for index           */
  int jmol;                  /* Num: Mol for index           */
  int i,j;                   /* Num: For indicies            */
  int iend;                  /* Num:                         */
  double now_memory;         /* Num: Memory allocated now    */ 
  double pre_nhc;            /* Num: NHC mass pre-factor     */
  int num_nhc_old;
  int iii,istart;
  int istart_res1,ires;
  int nghost_freeze;
  double nfree_cell; 
  double cst1,cst2,cst3,cst4;
  double text;

  /*    Local pointers   */

  int npt_i                = genensopts->npt_i;
  int npt_f                = genensopts->npt_f;

  double c2_lnv            = 0.0;
  double mass_lnv          = 0.0;
  double c1_hm             = 0.0;
  double mass_hm           = 0.0;

  int pimd_freez_typ       = mdclatoms_pimd->pimd_freez_typ;
  int *ighost_flag         = mdghost_atoms->ighost_flag;  
  int *freeze_flag         = mdconstrnt->freeze_flag; 
  int nfreeze              = mdconstrnt->nfreeze;

  int nmol_typ             = mdatom_maps->nmol_typ;
  int *jatm_jmol_typ_strt  = mdatom_maps->jatm_jmol_typ_strt;
  int *jres_jmol_typ_strt  = mdatom_maps->jres_jmol_typ_strt;
  int *natm_jres_jmol_typ  = mdatom_maps->natm_jres_jmol_typ;
  int *nfree_jres_jmol_typ = mdatom_maps->nfree_jres_jmol_typ;
  int *natm_1mol_jmol_typ  = mdatom_maps->natm_1mol_jmol_typ;
  int *nmol_jmol_typ       = mdatom_maps->nmol_jmol_typ;
  int *nres_1mol_jmol_typ  = mdatom_maps->nres_1mol_jmol_typ;
  int *icons_jmol_typ      = mdconstrnt->icons_jmol_typ;

  int *nfree_1mol_jmol_typ = mdatom_maps->nfree_1mol_jmol_typ;
  int *iatm_mol_typ        = mdatom_maps->iatm_mol_typ;

  int    *imol_nhc_opt     = class_parse->imol_nhc_opt;
  double *tau_nhc_mol      = class_parse->tau_nhc_mol;
  double *text_nhc_mol     = class_parse->text_nhc_mol;
  double tau_nhc_def       = class_parse->tau_nhc_def;
  double tau_vol           = class_parse->tau_vol;
  double tau_vol_nhc       = class_parse->tau_vol_nhc;
  double state_t_ext       = genstatepoint->t_ext;

  int class_num_nhc        = 0;
  int bead_num_nhc         = 0;
  int class_len_nhc        = mdtherm_info_class->len_nhc;
  int bead_len_nhc         = mdtherm_info_bead->len_nhc;

  int nfree                = mdclatoms_info->nfree;
  int pi_beads_true        = mdclatoms_info->pi_beads; 
  int pi_beads             = 1; 
  int natm_tot             = mdclatoms_info->natm_tot;
  int nghost_tot           = mdghost_atoms->nghost_tot;
  double gamma_adb         = mdclatoms_pimd->gamma_adb;


  int anneal_opt           = gensimopts->anneal_opt;
  double ann_start_temp    = gensimopts->ann_start_temp;

  double **class_mass_nhc;  /* assign after malloc */
  double **class_gkt;
  int    *class_inhc_x;
  int    *class_inhc_y;
  int    *class_inhc_z;
  double **bead_mass_nhc;
  double **bead_gkt;
  double *class_text_nhc;
  double *bead_text_nhc;
  int    *bead_inhc_x;
  int    *bead_inhc_y;
  int    *bead_inhc_z;
  double *gkt_vol;
  double *mass_vol_nhc;

  /*==========================================================================*/
  /* 0) Output to screen */

  PRINT_LINE_STAR;
  PRINTF("Set up atom NHC's\n");
  PRINT_LINE_DASH; PRINTF("\n");

  /*==========================================================================*/
  /* I) Calculate the number of atm NHC's by looping over the molecule typs   */

  class_num_nhc= 0;
  ind_glob_nhc   = 0;
  for(jmol_typ=1;jmol_typ<=nmol_typ;jmol_typ++){
    istart = jatm_jmol_typ_strt[jmol_typ];
    iend   = natm_1mol_jmol_typ[jmol_typ] +jatm_jmol_typ_strt[jmol_typ]-1;
    nghost_freeze = 0;
    for(jatm=istart;jatm<=iend;jatm++){
      if((ighost_flag[jatm] > 0) || (freeze_flag[jatm] > 0)){nghost_freeze++;}
    }/*endfor*/
    /*--------------------------------------------------------------------------*/
    /*   1) Couple all molecules of this type to global NHC                     */
    if(imol_nhc_opt[jmol_typ]==1){
      if(ind_glob_nhc == 0){
        class_num_nhc+= 1;
        ind_glob_nhc = class_num_nhc;
      }/*endif*/
    }/*endif*/
    /*--------------------------------------------------------------------------*/
    /*   2) Couple all molecules of this type to the same NHC                   */
    if(imol_nhc_opt[jmol_typ]==2){
      class_num_nhc+= 1;
    }/*endif*/
    /*--------------------------------------------------------------------------*/
    /*   3) Couple each molecule of this type to its own NHC                    */
    if(imol_nhc_opt[jmol_typ]== 3){
      class_num_nhc+= nmol_jmol_typ[jmol_typ];
    }/*endif*/
    /*--------------------------------------------------------------------------*/
    /*   4) Couple each residue of this molecule type to its own NHC            */
    if(imol_nhc_opt[jmol_typ]== 4){
      class_num_nhc+= nmol_jmol_typ[jmol_typ]*nres_1mol_jmol_typ[jmol_typ];
      if(icons_jmol_typ[jmol_typ] > 1){
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("Sorry dude, the residue thermostatting option\n");
        PRINTF("while most excellent is restricted to molecules\n");
        PRINTF("with no internal constraints between residues.\n"); 
        PRINTF("molecule type number %d\n",jmol_typ);
        PRINTF("has some constraints between residues\n"); 
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        FFLUSH(stdout);
        EXIT(1);
      }/*endif*/ 
    }/*endif*/
    /*--------------------------------------------------------------------------*/
    /*    5) Couple each atom in this molecule type to its own NHC              */
    if(imol_nhc_opt[jmol_typ]==5){
      class_num_nhc+= (natm_1mol_jmol_typ[jmol_typ] - nghost_freeze)
        *nmol_jmol_typ[jmol_typ];
      if(icons_jmol_typ[jmol_typ]!=0){
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("Sorry dude, the atomic thermostatting option\n");
        PRINTF("while most excellent is restricted to molecules\n");
        PRINTF("with no internal constraints.\n"); 
        PRINTF("molecule type number %d\n",jmol_typ);
        PRINTF("has some constraints\n"); 
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        FFLUSH(stdout);
        EXIT(1);
      }/*endif*/ 
    }/*endif*/
    /*--------------------------------------------------------------------------*/
    /*    6) Couple each atom in this molecule to its own NHC                   */
    if(imol_nhc_opt[jmol_typ]==6){
      class_num_nhc+= 3*(natm_1mol_jmol_typ[jmol_typ] - nghost_freeze)
        *nmol_jmol_typ[jmol_typ];
      if(icons_jmol_typ[jmol_typ]!=0){
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("sorry dude, the massive thermostatting option\n");
        PRINTF("while most excellent is restricted to molecules\n");
        PRINTF("with no internal constraints.\n"); 
        PRINTF("molecule type number %d\n",jmol_typ);
        PRINTF("has some constraints\n"); 
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        FFLUSH(stdout);
        EXIT(1);
      }/*endif*/
    }/*endif*/

  }/*endfor*/
  num_nhc_old = class_num_nhc;

  /*==========================================================================*/
  /* II) Output the information                                               */

  PRINTF("You have implemented %d nose-hoover chains ",class_num_nhc);
  PRINTF("each of length %d\n",class_len_nhc);

  if(pi_beads>1){
    if(pimd_freez_typ==1){
      bead_num_nhc = 3*(natm_tot - nghost_tot);
    }else{
      bead_num_nhc = 3*(natm_tot - nghost_tot - nfreeze);
    }/*endif*/
    PRINTF("You have implemented %d nose-hoover chains ",bead_num_nhc);
    PRINTF("each of length %d on each bead\n",bead_len_nhc);
  }/*endif*/

  /*==========================================================================*/
  /* III) Allocate the NHC memory and assign local pointers                   */

  now_memory  = ( (class_len_nhc)*(class_num_nhc)*
      (sizeof(double)*5 + sizeof(int)*0 )
      +(natm_tot)*
      (sizeof(double)*0 + sizeof(int)*3 )
      )*1.e-06;
  *tot_memory += now_memory;

  PRINTF("nhc allocation: %g mbytes; total memory %g mbytes\n",
      now_memory,*tot_memory);


#ifdef WHERE_DO_I_COME_FROM
  therm_class->x_nhc  = cmall_mat(1,class_len_nhc,1,class_num_nhc+1,"set_atm_NHC");
  therm_class->v_nhc  = cmall_mat(1,class_len_nhc,1,class_num_nhc+1,"set_atm_NHC");
  therm_class->f_nhc  = cmall_mat(1,class_len_nhc,1,class_num_nhc+1,"set_atm_NHC");
#endif

  mdtherm_info_class->mass_nhc = 
    cmall_mat(1,class_len_nhc,1,class_num_nhc+1,"set_atm_NHC");
  mdtherm_info_class->gkt      = 
    cmall_mat(1,class_len_nhc,1,class_num_nhc+1,"set_atm_NHC");

  mdtherm_info_class->inhc_x = (int *)cmalloc(natm_tot*sizeof(int),"set_atm_NHC")-1;
  mdtherm_info_class->inhc_y = (int *)cmalloc(natm_tot*sizeof(int),"set_atm_NHC")-1;
  mdtherm_info_class->inhc_z = (int *)cmalloc(natm_tot*sizeof(int),"set_atm_NHC")-1;
  mdtherm_info_class->text_nhc = (double *)
    cmalloc((class_num_nhc+1)*sizeof(double),"set_atm_NHC")-1;
  mdtherm_info_class->wdti    = (double *)
    cmalloc(((size_t)25*sizeof(double)),"set_atm_NHC")-1;
  mdtherm_info_class->wdti2    = (double *)
    cmalloc(((size_t)25*sizeof(double)),"set_atm_NHC")-1;
  mdtherm_info_class->wdti4    = (double *)
    cmalloc(((size_t)25*sizeof(double)),"set_atm_NHC")-1;
  mdtherm_info_class->wdti8    = (double *)
    cmalloc(((size_t)25*sizeof(double)),"set_atm_NHC")-1;
  mdtherm_info_class->wdti16   = (double *)
    cmalloc(((size_t)25*sizeof(double)),"set_atm_NHC")-1;

  class_mass_nhc  = mdtherm_info_class->mass_nhc;
  class_gkt       = mdtherm_info_class->gkt;
  class_inhc_x    = mdtherm_info_class->inhc_x;
  class_inhc_y    = mdtherm_info_class->inhc_y;
  class_inhc_z    = mdtherm_info_class->inhc_z;
  class_text_nhc  = mdtherm_info_class->text_nhc;

  if(pimd_on == 1){

#ifdef WHERE_DO_I_COME_FROM
    for(i=1;i<=pi_beads;i++){
      therm_bead[i].x_nhc = cmall_mat(1,bead_len_nhc,1,bead_num_nhc+1,"set_atm_NHC");
      therm_bead[i].v_nhc = cmall_mat(1,bead_len_nhc,1,bead_num_nhc+1,"set_atm_NHC");
      therm_bead[i].f_nhc = cmall_mat(1,bead_len_nhc,1,bead_num_nhc+1,"set_atm_NHC");
    }/*endfor*/
#endif

    mdtherm_info_bead->mass_nhc  = 
      cmall_mat(1,bead_len_nhc,1,bead_num_nhc+1,"set_atm_NHC");
    mdtherm_info_bead->gkt       = 
      cmall_mat(1,bead_len_nhc,1,bead_num_nhc+1,"set_atm_NHC");
    mdtherm_info_bead->inhc_x    = 
      (int *)cmalloc(natm_tot*sizeof(int),"set_atm_NHC")-1;
    mdtherm_info_bead->inhc_y    = 
      (int *)cmalloc(natm_tot*sizeof(int),"set_atm_NHC")-1;
    mdtherm_info_bead->inhc_z    = 
      (int *)cmalloc(natm_tot*sizeof(int),"set_atm_NHC")-1;
    mdtherm_info_bead->text_nhc  = 
      (double *)cmalloc((bead_num_nhc+1)*sizeof(double),"set_atm_NHC")-1;
    mdtherm_info_bead->wdti  = 
      (double *)cmalloc((size_t)9*sizeof(double),"set_atm_NHC")-1;
    mdtherm_info_bead->wdti2  = 
      (double *)cmalloc((size_t)9*sizeof(double),"set_atm_NHC")-1;
    mdtherm_info_bead->wdti4  = 
      (double *)cmalloc((size_t)9*sizeof(double),"set_atm_NHC")-1;
    mdtherm_info_bead->wdti8  = 
      (double *)cmalloc((size_t)9*sizeof(double),"set_atm_NHC")-1;
    mdtherm_info_bead->wdti16 = 
      (double *)cmalloc((size_t)9*sizeof(double),"set_atm_NHC")-1;

    bead_mass_nhc  = mdtherm_info_bead->mass_nhc;
    bead_gkt       = mdtherm_info_bead->gkt;
    bead_inhc_x    = mdtherm_info_bead->inhc_x;
    bead_inhc_y    = mdtherm_info_bead->inhc_y;
    bead_inhc_z    = mdtherm_info_bead->inhc_z;
    bead_text_nhc  = mdtherm_info_bead->text_nhc;

  }/*endif:pimd calculation*/

  mdbaro->x_vol_nhc    = 
    (double *)cmalloc(class_len_nhc*sizeof(double),"set_atm_NHC")-1;
  mdbaro->v_vol_nhc    = 
    (double *)cmalloc(class_len_nhc*sizeof(double),"set_atm_NHC")-1;
  mdbaro->f_vol_nhc    = 
    (double *)cmalloc(class_len_nhc*sizeof(double),"set_atm_NHC")-1;
  mdbaro->gkt_vol      = 
    (double *)cmalloc(class_len_nhc*sizeof(double),"set_atm_NHC")-1;
  mdbaro->mass_vol_nhc = 
    (double *)cmalloc(class_len_nhc*sizeof(double),"set_atm_NHC")-1;

  gkt_vol      = mdbaro->gkt_vol;
  mass_vol_nhc = mdbaro->mass_vol_nhc;

  /*==========================================================================*/
  /* IV) Convert the nhc time scales to au                                    */

  for(jmol_typ=1;jmol_typ<=nmol_typ;jmol_typ++){
    tau_nhc_mol[jmol_typ] /= TIME_CONV;
  }/*endfor*/ 
  tau_nhc_def /= TIME_CONV;
  tau_vol     /= TIME_CONV;
  tau_vol_nhc /= TIME_CONV;

  /*==========================================================================*/
  /* V)Set up the atm nhc's                                                  */

  for(jatm=1;jatm<=natm_tot;jatm++){
    class_inhc_x[jatm] = 0;
    class_inhc_y[jatm] = 0;
    class_inhc_z[jatm] = 0;
  }/*endfor*/

  class_num_nhc     = 0;
  ind_glob_nhc  = 0;
  for(jmol_typ=1;jmol_typ<=nmol_typ;jmol_typ++){
    istart = jatm_jmol_typ_strt[jmol_typ];
    iend   = natm_1mol_jmol_typ[jmol_typ]*nmol_jmol_typ[jmol_typ]
      +jatm_jmol_typ_strt[jmol_typ]-1;
    /*--------------------------------------------------------------------------*/
    /*   1) couple all molecules of this type to global nhc                     */
    if(imol_nhc_opt[jmol_typ]==1){
      if(ind_glob_nhc == 0){
        class_num_nhc+= 1;
        ind_glob_nhc = class_num_nhc;
        class_gkt[1][ind_glob_nhc] = 0.0;
        class_text_nhc[(class_num_nhc)]  = 
          (anneal_opt == 1 ? ann_start_temp : state_t_ext);
        class_mass_nhc[1][ind_glob_nhc] = (tau_nhc_def)*(tau_nhc_def)*state_t_ext;
      }/*endif*/
      class_gkt[1][ind_glob_nhc] +=  (double)((nfree_1mol_jmol_typ[jmol_typ]
            *nmol_jmol_typ[jmol_typ]));
      for(jatm=istart;jatm<=iend;jatm++){
        class_inhc_x[jatm] = ind_glob_nhc;
        class_inhc_y[jatm] = ind_glob_nhc;
        class_inhc_z[jatm] = ind_glob_nhc;
      }/*endfor*/
    }/*endif*/
    /*--------------------------------------------------------------------------*/
    /*   2) couple all molecules of this type to the same nhc                  */
    if(imol_nhc_opt[jmol_typ]==2){
      class_num_nhc+= 1;
      class_gkt[1][(class_num_nhc)] =  (double)((nfree_1mol_jmol_typ[jmol_typ]
            *nmol_jmol_typ[jmol_typ]));
      class_text_nhc[(class_num_nhc)]    = 
        (anneal_opt == 1 ? ann_start_temp : text_nhc_mol[jmol_typ]);
      class_mass_nhc[1][(class_num_nhc)] = tau_nhc_mol[jmol_typ]*
        tau_nhc_mol[jmol_typ]*text_nhc_mol[jmol_typ];
      for(jatm=istart;jatm<=iend;jatm++){
        class_inhc_x[jatm] = class_num_nhc;
        class_inhc_y[jatm] = class_num_nhc;
        class_inhc_z[jatm] = class_num_nhc;
      }/*endfor*/
    }/*endif*/
    /*--------------------------------------------------------------------------*/
    /*   3) couple each molecule of this type to its own nhc                   */
    if(imol_nhc_opt[jmol_typ]== 3){
      for(jmol=1;jmol<=nmol_jmol_typ[jmol_typ];jmol++){
        class_num_nhc+= 1;
        class_gkt[1][(class_num_nhc)] = 
          (double)(nfree_1mol_jmol_typ[jmol_typ]);
        class_text_nhc[(class_num_nhc)]     = 
          (anneal_opt == 1 ? ann_start_temp : text_nhc_mol[jmol_typ]);
        class_mass_nhc[1][(class_num_nhc)]  = tau_nhc_mol[jmol_typ]*
          tau_nhc_mol[jmol_typ]*text_nhc_mol[jmol_typ];
        iend = istart + natm_1mol_jmol_typ[jmol_typ]-1;
        for(jatm=istart;jatm<=iend;jatm++){      
          class_inhc_x[jatm] = class_num_nhc;
          class_inhc_y[jatm] = class_num_nhc;
          class_inhc_z[jatm] = class_num_nhc;
        }/*endfor*/
        istart += natm_1mol_jmol_typ[jmol_typ];
      }/*endfor*/
    }/*endif*/
    /*--------------------------------------------------------------------------*/
    /*   4) couple each residue of this molecule type to its own nhc            */
    if(imol_nhc_opt[jmol_typ]== 4){
      istart_res1   = jres_jmol_typ_strt[jmol_typ]-1;
      for(jmol=1;jmol<=nmol_jmol_typ[jmol_typ];jmol++){
        for(ires=1;ires<=nres_1mol_jmol_typ[jmol_typ];ires++){
          iend = istart+natm_jres_jmol_typ[ires+istart_res1]-1;
          class_num_nhc+= 1;
          class_gkt[1][(class_num_nhc)] = 
            nfree_jres_jmol_typ[ires+istart_res1];
          class_text_nhc[(class_num_nhc)]     = 
            (anneal_opt == 1 ? ann_start_temp : text_nhc_mol[jmol_typ]);
          class_mass_nhc[1][(class_num_nhc)]  = tau_nhc_mol[jmol_typ]*
            tau_nhc_mol[jmol_typ]*text_nhc_mol[jmol_typ];
          for(jatm=istart;jatm<=iend;jatm++){
            class_inhc_x[jatm] = class_num_nhc;
            class_inhc_y[jatm] = class_num_nhc;
            class_inhc_z[jatm] = class_num_nhc;
          }/*endfor: jatm*/
          istart = iend+1;
        }/*endfor: ires*/
      }/*endfor: jmol*/
    }/*endif: res_therm*/
    /*--------------------------------------------------------------------------*/
    /*    5) couple each atom in each molecule to its own nhc                  */
    if(imol_nhc_opt[jmol_typ]==5){
      for(jatm=istart;jatm<=iend;jatm++){
        if((ighost_flag[jatm] == 0) && (freeze_flag[jatm]== 0)) {
          class_num_nhc += 1;
          class_gkt[1][(class_num_nhc)]       = 3.0;
          class_text_nhc[(class_num_nhc)]     = 
            (anneal_opt == 1 ? ann_start_temp : text_nhc_mol[jmol_typ]);
          class_mass_nhc[1][(class_num_nhc)]  = tau_nhc_mol[jmol_typ]*
            tau_nhc_mol[jmol_typ]*text_nhc_mol[jmol_typ];
          class_inhc_x[jatm] = class_num_nhc;
          class_inhc_y[jatm] = class_num_nhc;
          class_inhc_z[jatm] = class_num_nhc;
        }/*endif*/
      }/*endfor*/
    }/*endif*/
    /*--------------------------------------------------------------------------*/
    /*    6) couple each coord of atom in this molecule to its own nhc         */
    if(imol_nhc_opt[jmol_typ]==6){
      for(jatm=istart;jatm<=iend;jatm++){      
        if((ighost_flag[jatm] == 0) && (freeze_flag[jatm]==0)) {
          class_num_nhc+= 1;
          class_gkt[1][(class_num_nhc)]       = 1.0;
          class_text_nhc[(class_num_nhc)]     = 
            (anneal_opt == 1 ? ann_start_temp : text_nhc_mol[jmol_typ]);
          class_mass_nhc[1][(class_num_nhc)]  = tau_nhc_mol[jmol_typ]*
            tau_nhc_mol[jmol_typ]*text_nhc_mol[jmol_typ];
          class_inhc_x[jatm] = class_num_nhc;

          class_num_nhc+= 1;
          class_gkt[1][(class_num_nhc)]       = 1.0;
          class_text_nhc[(class_num_nhc)]     = 
            (anneal_opt == 1 ? ann_start_temp : text_nhc_mol[jmol_typ]);
          class_mass_nhc[1][(class_num_nhc)]  = tau_nhc_mol[jmol_typ]*
            tau_nhc_mol[jmol_typ]*text_nhc_mol[jmol_typ];
          class_inhc_y[jatm] = class_num_nhc;

          class_num_nhc+= 1;
          class_gkt[1][(class_num_nhc)]       = 1.0;
          class_text_nhc[(class_num_nhc)]     = 
            (anneal_opt == 1 ? ann_start_temp : text_nhc_mol[jmol_typ]);
          class_mass_nhc[1][(class_num_nhc)]  = tau_nhc_mol[jmol_typ]*
            tau_nhc_mol[jmol_typ]*text_nhc_mol[jmol_typ];
          class_inhc_z[jatm] = class_num_nhc;
        }/*endif*/
      }/*endfor*/          
    }/*endif*/

  }/*endfor*/

  /*==========================================================================*/
  /* VI) Hook up the ghosts, frozen atoms and non-thermostatted atoms         */
  /*     to the null thermostat                                               */

  for(jmol_typ=1;jmol_typ<=nmol_typ;jmol_typ++){
    if(imol_nhc_opt[jmol_typ]==0){
      istart = jatm_jmol_typ_strt[jmol_typ];
      iend   = natm_1mol_jmol_typ[jmol_typ]*nmol_jmol_typ[jmol_typ]
        +jatm_jmol_typ_strt[jmol_typ]-1;
      for(jatm=istart;jatm<=iend;jatm++){      
        class_inhc_x[jatm] = (class_num_nhc+1);
        class_inhc_y[jatm] = (class_num_nhc+1);
        class_inhc_z[jatm] = (class_num_nhc+1);
      }/*endfor*/
    }/*endif*/
  }/*endfor*/

  for(jatm=1;jatm <= natm_tot;jatm++){
    if((ighost_flag[jatm] > 0) || (freeze_flag[jatm] > 0)){
      class_inhc_x[jatm] = class_num_nhc+1;
      class_inhc_y[jatm] = class_num_nhc+1;
      class_inhc_z[jatm] = class_num_nhc+1;
    }/*endif*/
  }/*endfor*/

  /*==========================================================================*/
  /* VII)Expand out the nhc to len_nhc and tidy up                            */

  if(mdtherm_info_class->therm_typ == 1){
    for(j=2;j<=class_len_nhc;j++){
      for(i=1;i<=class_num_nhc;i++){
        class_gkt[j][i]      = class_text_nhc[i]/BOLTZ;
        class_mass_nhc[j][i] = class_mass_nhc[1][i]/BOLTZ;
      }/*endfor*/
    }/*endfor*/

    /*==========================================================================*/
    /* Get GGMT thermostat masses*/

  } else {
    for(i=1;i<=class_num_nhc;i++){ 
      cst1 = class_gkt[1][i]+2.0;
      cst2 = cst1*(class_gkt[1][i]+4.0);
      text = class_text_nhc[i]/BOLTZ;
      class_mass_nhc[2][i] = class_mass_nhc[1][i]*class_gkt[1][i]
        *text*text*(1.0+cst2/(cst1*cst1));

      if(class_len_nhc == 3){
        cst3 = cst2*(class_gkt[1][i]+6.0);
        cst4 = cst3*(class_gkt[1][i]+8.0);
        class_mass_nhc[3][i] = class_mass_nhc[1][i]*class_gkt[1][i]
          *text*text*text*text
          *(1.0+cst3/(cst1*cst2)+cst4/(cst2*cst2));
      }/*endif*/
    }/*endfor*/
  }/*endelse*/
  /*==========================================================================*/
  /* First thermostat masses */

  for(i=1;i<=class_num_nhc;i++){
    class_mass_nhc[1][i] = class_mass_nhc[1][i]*class_gkt[1][i]/BOLTZ;

    class_gkt[1][i]      = class_gkt[1][i]*class_text_nhc[i]/BOLTZ;
  }/*endfor*/ 

  /*==========================================================================*/
  /* VIII) Check */

  if(class_num_nhc==0 && pi_beads==1){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("dude, if the particles are not coupled to nhcs\n");
    PRINTF("you are just doing an nve calculation \n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }/*endif*/

  if(class_num_nhc!=num_nhc_old){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("internal error in thermostat set up routine \n");
    PRINTF("%d thermos found %d expected \n",class_num_nhc,num_nhc_old);
    PRINTF("*contact your technical support team*\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }/*endif*/

  for(i=1;i<=class_num_nhc;i++){
    if((class_mass_nhc[1][i]<=0.0) || (class_gkt[1][i]<=0)){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("interal error in thermostat set up routine \n");
      PRINTF("mass of thermostat and \n");
      PRINTF("# of degrees of freedom must be >0\n");
      PRINTF("%f %f\n",class_mass_nhc[1][i],class_gkt[1][i]);
      PRINTF("at thermostat number %d\n",i);
      PRINTF("*contact your technical support team*\n");
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }/*endif*/
  }/*endfor*/

  for(jatm=2;jatm<=natm_tot;jatm++){      
    if((class_inhc_x[jatm]<0)               ||
        (class_inhc_x[jatm]>class_num_nhc+1) ||
        (class_inhc_y[jatm]<0)               ||
        (class_inhc_y[jatm]>class_num_nhc+1) ||
        (class_inhc_z[jatm]<0)               || 
        (class_inhc_z[jatm]>class_num_nhc+1) ){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("interal error in atm nhc set up routine \n");
      PRINTF("unassigned pointer to atom number, %d\n",jatm);
      PRINTF("*contact your technical support team*\n");
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }/*endif*/
  }/*endfor*/

  /*==========================================================================*/
  /* IX)Set up the pressure nhc's                                             */

  /*--------------------------------------------------------------------------*/
  /*  1) isotropic */
  if(npt_i==1){
    c2_lnv          = 1.0+3.0/((double)(nfree));
    mass_lnv        = (state_t_ext/BOLTZ)*(tau_vol)*(tau_vol)*
      (double)(nfree+3);
    gkt_vol[1]      = (state_t_ext/BOLTZ);
    mass_vol_nhc[1] = gkt_vol[1]*(tau_vol_nhc)*(tau_vol_nhc);
    for(i=2;i<=class_len_nhc;i++){
      gkt_vol[i]      = state_t_ext/BOLTZ;
      mass_vol_nhc[i] = gkt_vol[i]*(tau_vol_nhc)*(tau_vol_nhc);
    }/*endfor*/
  }/*endif*/
  /*--------------------------------------------------------------------------*/
  /*  2) flexible */
  if(npt_f==1){
    c2_lnv          = 1.0+3.0/((double)(nfree));
    c1_hm           = 1.0/((double)(nfree));
    mass_lnv        = (state_t_ext/BOLTZ)*(tau_vol)*(tau_vol)*
      (double)(nfree+3);
    mass_hm         = (state_t_ext/BOLTZ)*(tau_vol)*(tau_vol)*
      (double)(nfree+3)/3.0;
    if(hmat_cons_typ==0){nfree_cell=6.0;}
    if(hmat_cons_typ==1){nfree_cell=3.0;}
    if(hmat_cons_typ==2){nfree_cell=4.0;}
    gkt_vol[1]      = nfree_cell*(state_t_ext/BOLTZ);
    mass_vol_nhc[1] = gkt_vol[1]*(tau_vol_nhc)*(tau_vol_nhc);
    for(i=2;i<=class_len_nhc;i++){
      gkt_vol[i]      = state_t_ext/BOLTZ;
      mass_vol_nhc[i] = gkt_vol[i]*(tau_vol_nhc)*(tau_vol_nhc);
    }/*endfor*/
  }/*endif*/

  /*==========================================================================*/
  /* X) set up the bead thermostats                                           */

  if(pi_beads>1){

    num_nhc_old = bead_num_nhc;
    for(jatm=1;jatm<=natm_tot;jatm++){
      bead_inhc_x[jatm]  = 0;
      bead_inhc_y[jatm]  = 0;
      bead_inhc_z[jatm]  = 0;
    }/*endfor*/

    bead_num_nhc  = 0;
    ind_glob_nhc  = 0;
    pre_nhc       = (BOLTZ*BOLTZ)/(gamma_adb*gamma_adb*((double)(pi_beads)));

    /*--------------------------------------------------------------------------*/
    /*    1) Couple each coord of the atom bead to its own NHC         */
    for(jatm=1;jatm<=natm_tot;jatm++){      
      jmol_typ = iatm_mol_typ[jatm];
      if(ighost_flag[jatm] == 0) {
        if(pimd_freez_typ==1 || freeze_flag[jatm]==0){
          bead_num_nhc += 1;
          bead_gkt[1][(bead_num_nhc)]       = 1.0;
          bead_mass_nhc[1][(bead_num_nhc)]  = 
            pre_nhc/(text_nhc_mol[jmol_typ]*text_nhc_mol[jmol_typ]);
          bead_text_nhc[(bead_num_nhc)]     = 
            (anneal_opt == 1 ? ann_start_temp : text_nhc_mol[jmol_typ]);
          bead_inhc_x[jatm]                 = bead_num_nhc;

          bead_num_nhc += 1;
          bead_gkt[1][(bead_num_nhc)]       = 1.0;
          bead_mass_nhc[1][(bead_num_nhc)]  = 
            pre_nhc/(text_nhc_mol[jmol_typ]*text_nhc_mol[jmol_typ]);
          bead_text_nhc[(bead_num_nhc)]     = 
            (anneal_opt == 1 ? ann_start_temp : text_nhc_mol[jmol_typ]);
          bead_inhc_y[jatm]                 = bead_num_nhc;

          bead_num_nhc += 1;
          bead_gkt[1][(bead_num_nhc)]       = 1.0;
          bead_mass_nhc[1][(bead_num_nhc)]  = 
            pre_nhc/(text_nhc_mol[jmol_typ]*text_nhc_mol[jmol_typ]);
          bead_text_nhc[(bead_num_nhc)]  = 
            (anneal_opt == 1 ? ann_start_temp : text_nhc_mol[jmol_typ]);
          bead_inhc_z[jatm] = bead_num_nhc;
        }/*endif*/
      }/*endif*/
    }/*endfor:jatm*/          

    /*--------------------------------------------------------------------------*/
    /*    2) Ghost atoms/frozen atoms connected to null thermostat             */
    for(jatm=1;jatm <= natm_tot;jatm++){
      if(ighost_flag[jatm] > 0 ){
        bead_inhc_x[jatm]  = bead_num_nhc+1;
        bead_inhc_y[jatm]  = bead_num_nhc+1;
        bead_inhc_z[jatm]  = bead_num_nhc+1;
      }/*endif*/
      if(pimd_freez_typ==2 && freeze_flag[jatm]>0){
        bead_inhc_x[jatm]  = bead_num_nhc+1;
        bead_inhc_y[jatm]  = bead_num_nhc+1;
        bead_inhc_z[jatm]  = bead_num_nhc+1;
      }/*endif*/
    }/*endfor*/

    /*--------------------------------------------------------------------------*/
    /*   3)Expand out the NHC to len_nhc and tidy up                            */

    for(j=2;j<=bead_len_nhc;j++){
      for(i=1;i<=bead_num_nhc;i++){
        bead_gkt[j][i]      = bead_text_nhc[i]/BOLTZ;
        bead_mass_nhc[j][i] = bead_mass_nhc[1][i]*bead_text_nhc[i]/BOLTZ;
      }/*endfor*/
    }/*endfor*/

    for(i=1;i<=bead_num_nhc;i++){
      bead_mass_nhc[1][i] = bead_mass_nhc[1][i]*bead_text_nhc[i]/BOLTZ;
      bead_gkt[1][i]      = bead_text_nhc[i]/BOLTZ;
    }/*endfor*/     

    /*--------------------------------------------------------------------------*/
    /*   4) Checks                                                             */

    if(bead_num_nhc!=num_nhc_old){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("Internal error in bead thermostat set up routine \n");
      PRINTF("%d thermos found %d expected \n",bead_num_nhc,num_nhc_old);
      PRINTF("*Contact your technical support team*\n");
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }/*endif*/
    for(i=1;i<=bead_num_nhc;i++){
      if((bead_mass_nhc[1][i]<=0.0) || (bead_gkt[1][i]<=0)){
        PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("Internal error in bead thermostat set up routine \n");
        PRINTF("mass of thermostat and \n");
        PRINTF("# of degrees of freedom must be >0\n");
        PRINTF("%g %g\n",bead_mass_nhc[1][i],bead_gkt[1][i]);
        PRINTF("%g %g\n",bead_mass_nhc[1][i],bead_gkt[1][i]);
        PRINTF("at thermostat number %d\n",i);
        PRINTF("*Contact your technical support team*\n");
        PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        FFLUSH(stdout);
        EXIT(1);
      }/*endif*/
    }/*endfor*/
    for(jatm=2;jatm<=natm_tot;jatm++){      
      if((bead_inhc_x[jatm]<0)              || 
          (bead_inhc_x[jatm]>bead_num_nhc+1) || 
          (bead_inhc_y[jatm]<0)              || 
          (bead_inhc_y[jatm]>bead_num_nhc+1) || 
          (bead_inhc_z[jatm]<0)              || 
          (bead_inhc_z[jatm]>bead_num_nhc+1) ){
        PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("Interal error in bead NHC set up routine \n");
        PRINTF("unassigned pointer to atom number, %d\n",jatm);
        PRINTF("*Contact your technical support team*\n");
        PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        FFLUSH(stdout);
        EXIT(1);
      }/*endif*/
    }/*endfor*/

  }/*endif:more than one bead*/       

  /*==========================================================================*/
  /* XI) Tuck things back in the structures                                   */

  mdtherm_info_class->num_nhc = class_num_nhc;
  mdtherm_info_bead->num_nhc  = bead_num_nhc;

  class_parse->tau_nhc_def  = tau_nhc_def;
  class_parse->tau_vol      = tau_vol;     
  class_parse->tau_vol_nhc  = tau_vol_nhc;

  mdbaro->c2_lnv              = c2_lnv;
  mdbaro->mass_lnv            = mass_lnv;
  mdpar_rahman->c1_hm         = c1_hm;
  mdpar_rahman->mass_hm       = mass_hm;

  /*==========================================================================*/
  /* XII) Output to screen                                                    */

  PRINTF("\n");PRINT_LINE_DASH;
  PRINTF("Completed set up of atom NHC's\n");
  PRINT_LINE_STAR;PRINTF("\n");

  /*--------------------------------------------------------------------------*/
}/* end routine */
/*===========================================================================*/

