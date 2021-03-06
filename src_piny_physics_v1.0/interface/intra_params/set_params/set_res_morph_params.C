/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  set_res_name_params:set up the resname params                           */
/*==========================================================================*/

#include "standard_include.h"
#include "../class_defs/allclass_mdintegrate.h"
#include "../class_defs/allclass_mdatoms.h"
#include "../class_defs/allclass_mdinter.h"
#include "../class_defs/allclass_mdintra.h"
#include "../class_defs/allclass_cp.h"
#include "../class_defs/typedefs_par.h"
#include "../proto_defs/proto_intra_params_local.h"
#include "../proto_defs/proto_handle_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_res_morph_params(DICT_WORD res_name_dict[],int num_res_name_dict,
    char fun_key[],char file_name[],NAME res_typ[],
    int natm_jres_jmol_typ[],int jmol_typ,
    int jres,int jres_off,
    int ires_typ_jres_jmol_typ[])
  /*==========================================================================*/
{
  int i,index,index2;
  /*=======================================================================*/
  /* I) Check for missing key words*/
  for(i=1;i<=1;i++){
    if(res_name_dict[i].iuset==0 && res_name_dict[i].key_type==1){
      keyword_miss(res_name_dict,file_name,fun_key,i);}
  }    /*endfor*/
  /*======================================================================*/
  /* II) Set the keywords */
  /*-------------------------------------------------------------------------*/
  /* 1) \residue_name{} */

  index = 1;  
  index2 = ires_typ_jres_jmol_typ[(jres+jres_off)];
  if(strcasecmp(res_typ[index2],res_name_dict[1].keyarg)!=0){
    keyarg_barf(res_name_dict,file_name,fun_key,index);}

  /*---------------------------------------------------------------------*/
}  /*end routine*/
/*==========================================================================*/
