/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: set_mol_parms.c                              */
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

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_def_base_params(FILENAME_PARSE *filename_parse,
    DICT_WORD def_base_dict[], int num_def_base_dict)
{
  /*==========================================================================*/
  /*-----------------------------------------------------------------------*/ 
  /* 1) /inter_file{} */
  strcpy(filename_parse->def_inter_name,def_base_dict[1].keyarg);
  /*-----------------------------------------------------------------------*/ 
  /* 2) /vps_file{} */
  strcpy(filename_parse->def_vps_name,def_base_dict[2].keyarg);
  /*-----------------------------------------------------------------------*/ 
  /* 3) /bond_file{} */
  strcpy(filename_parse->def_intra_name[1],def_base_dict[3].keyarg);
  /*-----------------------------------------------------------------------*/ 
  /* 4) /bend_file{} */
  strcpy(filename_parse->def_intra_name[2],def_base_dict[4].keyarg);
  /*-----------------------------------------------------------------------*/ 
  /* 5) /tors_file{} */
  strcpy(filename_parse->def_intra_name[3],def_base_dict[5].keyarg);
  /*-----------------------------------------------------------------------*/ 
  /* 6) /tors_file{} */
  strcpy(filename_parse->def_intra_name[4],def_base_dict[6].keyarg);
  /*-----------------------------------------------------------------------*/ 
  /* 7) /surf_file{} */
  strcpy(filename_parse->def_surf_name,def_base_dict[7].keyarg);
  /*-----------------------------------------------------------------------*/ 
  /* 8) /lj_explicit_file{} */
  strcpy(filename_parse->def_lj_name,def_base_dict[8].keyarg);
  /*-----------------------------------------------------------------------*/ 
  /* 9) /grimme_ile{} */
  strcpy(filename_parse->def_grimme_name,def_base_dict[9].keyarg);
  /*=======================================================================*/
}/*end routine*/
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_user_base_params(FILENAME_PARSE *filename_parse,
    DICT_WORD user_base_dict[], int num_user_base_dict)
{
  /*=======================================================================*/
  /*---------------------------------------------------------------------*/ 
  /* 1) /inter_file{} */
  strcpy(filename_parse->user_inter_name,user_base_dict[1].keyarg);
  /*-----------------------------------------------------------------------*/ 
  /* 2) /vps_file{} */
  strcpy(filename_parse->user_vps_name,user_base_dict[2].keyarg);
  /*-----------------------------------------------------------------------*/ 
  /* 3) /bond_file{} */
  strcpy(filename_parse->user_intra_name[1],user_base_dict[3].keyarg);
  /*-----------------------------------------------------------------------*/ 
  /* 4) /bend_file{} */
  strcpy(filename_parse->user_intra_name[2],user_base_dict[4].keyarg);
  /*-----------------------------------------------------------------------*/ 
  /* 5) /tors_file{} */
  strcpy(filename_parse->user_intra_name[3],user_base_dict[5].keyarg);
  /*-----------------------------------------------------------------------*/ 
  /* 6) /onfo_file{} */
  strcpy(filename_parse->user_intra_name[4],user_base_dict[6].keyarg);
  /*-----------------------------------------------------------------------*/ 
  /* 7) /surf_file{} */
  strcpy(filename_parse->user_surf_name,user_base_dict[7].keyarg);
  /*-----------------------------------------------------------------------*/ 
  /* 8) /lj_explicit_file{} */
  strcpy(filename_parse->user_lj_name,user_base_dict[8].keyarg);
  /*-----------------------------------------------------------------------*/ 
  /* 9) /grimme_file{} */
  strcpy(filename_parse->user_grimme_name,user_base_dict[9].keyarg);
  /*=======================================================================*/
}/*end routine*/
/*==========================================================================*/
