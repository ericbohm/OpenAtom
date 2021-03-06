/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  set_bend_params:set up intramolecular interaction                       */
/*                    key word dictionary                                   */
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

void set_bend_bnd_params(DICT_WORD intra_dict[],int num_intra_dict,
    char *fun_key,char *file_name,int jmol_typ,
    MDCLATOMS_INFO *clatoms_info,MDATOM_MAPS *atommaps,
    MDBEND_BND *bend_bnd,MDBEND *bend,
    NULL_INTER_PARSE *null_inter_parse,
    BUILD_INTRA *build_intra, int iresidue,int ires_off)

  /*==========================================================================*/
{
  int num,index,ifound,igo;
  int itype1,itype2,itype3;
  int iatm_ind1,iatm_ind2,iatm_ind3;
  int imask1,imask2,imask3;
  int i,itype;

  /*========================================================================*/
  /* I) Check for missing key words*/
  for(i=1;i<=3;i++){
    if(intra_dict[i].iuset==0 && intra_dict[i].key_type==1){
      keyword_miss(intra_dict,file_name,fun_key,i);
    }
  } /*endfor*/

  /*=======================================================================*/
  /* II) Fill the dictionary with words */
  /*-----------------------------------------------------------------------*/
  /*  1) \atom1{}    */

  index = 1;
  sscanf(intra_dict[1].keyarg,"%d",&num);
  iatm_ind1 = num;
  if(iatm_ind1>build_intra->natmind_1res_now||iatm_ind1<0){
    keyarg_barf(intra_dict,file_name,fun_key,index);}
  imask1 = build_intra->mask_atm[iatm_ind1];  
  if(imask1>0)iatm_ind1 = build_intra->index_atm[iatm_ind1];

  /*------------------------------------------------------------------------*/
  /*  2) \atom2{}    */
  index = 2;
  sscanf(intra_dict[2].keyarg,"%d",&num);
  iatm_ind2 = num;
  if(iatm_ind2>build_intra->natmind_1res_now||iatm_ind2<0){
    keyarg_barf(intra_dict,file_name,fun_key,index);
  }
  imask2 = build_intra->mask_atm[iatm_ind2];  
  if(imask2>0)iatm_ind2 = build_intra->index_atm[iatm_ind2];

  /*-----------------------------------------------------------------------*/
  /*  3) \atom3{}    */
  index = 3;
  sscanf(intra_dict[3].keyarg,"%d",&num);
  iatm_ind3 = num;
  if(iatm_ind3>build_intra->natmind_1res_now||iatm_ind3<0){
    keyarg_barf(intra_dict,file_name,fun_key,index);}
  imask3 = build_intra->mask_atm[iatm_ind3];  
  if(imask3>0)iatm_ind3 = build_intra->index_atm[iatm_ind3];
  /*-----------------------------------------------------------------------*/
  /*  5) \modifier{} */
  index  = 5;
  ifound = 0;
  if(strcasecmp(intra_dict[5].keyarg,"on")==0) {ifound = 1;}
  if(strcasecmp(intra_dict[5].keyarg,"con")==0){ifound = 2;}
  if(strcasecmp(intra_dict[5].keyarg,"off")==0){ifound = 3;}
  if(ifound==0){
    keyarg_barf(intra_dict,file_name,fun_key,index);
  }
  if(ifound == 2){ 
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Constrained bends are toast\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }/*endif*/
  /*-----------------------------------------------------------------------*/
  /*  6) bend type    */
  igo = imask1*imask2*imask3;
  if(igo==1){
    itype1 = atommaps->iatm_atm_typ[(clatoms_info->natm_tot+iatm_ind1)];
    itype2 = atommaps->iatm_atm_typ[(clatoms_info->natm_tot+iatm_ind2)];
    itype3 = atommaps->iatm_atm_typ[(clatoms_info->natm_tot+iatm_ind3)];
    strcpy(build_intra->cbend_bnd_typ_now->atm1,atommaps->atm_typ[itype1]);
    strcpy(build_intra->cbend_bnd_typ_now->atm2,atommaps->atm_typ[itype2]);
    strcpy(build_intra->cbend_bnd_typ_now->atm3,atommaps->atm_typ[itype3]);
    strcpy(build_intra->cbend_bnd_typ_now->label,intra_dict[6].keyarg);  
  }/*endif*/

  /*=======================================================================*/
  /* III) Spread bends_bnd */
  if(ifound==1&&igo==1){
    /* A) Add more space */
    if(bend_bnd->num+1 > build_intra->nbend_bnd_max){
      build_intra->nbend_bnd_max += NMEM_MIN;
      bend_bnd->j1 =(int*)crealloc(&(bend_bnd->j1)[1], 
          build_intra->nbend_bnd_max*sizeof(int),"set_bend_bnd_params")-1;
      bend_bnd->j2 =(int*)crealloc(&(bend_bnd->j2)[1],
          build_intra->nbend_bnd_max*sizeof(int),"set_bend_bnd_params")-1;
      bend_bnd->j3 =(int*)crealloc(&(bend_bnd->j3)[1],
          build_intra->nbend_bnd_max*sizeof(int),"set_bend_bnd_params")-1;
      bend_bnd->jtyp=(int*)crealloc(&(bend_bnd->jtyp)[1],
          build_intra->nbend_bnd_max*sizeof(int),"set_bend_bnd_params")-1;
    }/*endif*/
    /*--------------------------------------------------------------------*/
    /* C) Check type */
    itype = bend_bnd->ntyp+1;
    for(i=1;i<=bend_bnd->ntyp;i++){

      if((strcasecmp(build_intra->cbend_bnd_typ[i].atm1,
              build_intra->cbend_bnd_typ_now->atm1)==0)
          &&(strcasecmp(build_intra->cbend_bnd_typ[i].atm2,
              build_intra->cbend_bnd_typ_now->atm2)==0)
          &&(strcasecmp(build_intra->cbend_bnd_typ[i].atm3,
              build_intra->cbend_bnd_typ_now->atm3)==0)
          &&(strcasecmp(build_intra->cbend_bnd_typ[i].label,
              build_intra->cbend_bnd_typ_now->label)==0)) {itype=i;}

      if((strcasecmp(build_intra->cbend_bnd_typ[i].atm1,
              build_intra->cbend_bnd_typ_now->atm3)==0)
          &&(strcasecmp(build_intra->cbend_bnd_typ[i].atm2,
              build_intra->cbend_bnd_typ_now->atm2)==0)
          &&(strcasecmp(build_intra->cbend_bnd_typ[i].atm3,
              build_intra->cbend_bnd_typ_now->atm1)==0)
          &&(strcasecmp(build_intra->cbend_bnd_typ[i].label,
              build_intra->cbend_bnd_typ_now->label)==0)) {itype=i;}
    }      /*endfor*/

    /*---------------------------------------------------------------------*/
    /* D) Add space */
    if(itype>build_intra->nbend_bnd_typ_max){
      build_intra->nbend_bnd_typ_max += NMEM_MIN;
      build_intra->cbend_bnd_typ      = 
        (CBEND *) crealloc(&build_intra->cbend_bnd_typ[1],
            build_intra->nbend_bnd_typ_max*sizeof(CBEND),"set_bend_bnd_params")-1;
    }/*endif*/

    /*----------------------------------------------------------------------*/
    /* E) Add a type */
    if(itype==bend_bnd->ntyp+1){
      bend_bnd->ntyp++;
      strcpy(build_intra->cbend_bnd_typ[itype].atm1,
          build_intra->cbend_bnd_typ_now->atm1);
      strcpy(build_intra->cbend_bnd_typ[itype].atm2,
          build_intra->cbend_bnd_typ_now->atm2);
      strcpy(build_intra->cbend_bnd_typ[itype].atm3,
          build_intra->cbend_bnd_typ_now->atm3);
      strcpy(build_intra->cbend_bnd_typ[itype].label,
          build_intra->cbend_bnd_typ_now->label);
    }/*endif*/

    /*---------------------------------------------------------------------*/
    /* B) Spread */

    bend_bnd->num += 1;
    bend_bnd->j1[bend_bnd->num] = iatm_ind1 + clatoms_info->natm_tot;
    bend_bnd->j2[bend_bnd->num] = iatm_ind2 + clatoms_info->natm_tot;
    bend_bnd->j3[bend_bnd->num] = iatm_ind3 + clatoms_info->natm_tot;
    bend_bnd->jtyp[bend_bnd->num] = itype;
  }/*endif*/
  /*=======================================================================*/
  /* V) Null uri-bradleys into bends */
  if(ifound==3&&igo==1){
    /*---------------------------------------------------------------------*/
    /* A) Add more space */
    if(null_inter_parse->nbend_nul+1 > build_intra->nbend_nul_max){

      build_intra->nbend_nul_max += NMEM_MIN;
      null_inter_parse->jbend1_nul     = 
        (int *) crealloc(&(null_inter_parse->jbend1_nul)[1],
            build_intra->nbend_nul_max*sizeof(int),"set_bend_bnd_params")-1;
      null_inter_parse->jbend2_nul     = 
        (int *) crealloc(&(null_inter_parse->jbend2_nul)[1],
            build_intra->nbend_nul_max*sizeof(int),"set_bend_bnd_params")-1;
      null_inter_parse->jbend3_nul     = 
        (int *) crealloc(&(null_inter_parse->jbend3_nul)[1],
            build_intra->nbend_nul_max*sizeof(int),"set_bend_bnd_params")-1;
    }     /*endif*/
    /*---------------------------------------------------------------------*/
    /* assign */
    null_inter_parse->nbend_nul += 1;
    (null_inter_parse->jbend1_nul)[null_inter_parse->nbend_nul] = 
      iatm_ind1 + clatoms_info->natm_tot;
    (null_inter_parse->jbend2_nul)[null_inter_parse->nbend_nul] = 
      iatm_ind2 + clatoms_info->natm_tot;
    (null_inter_parse->jbend3_nul)[null_inter_parse->nbend_nul] = 
      iatm_ind3 + clatoms_info->natm_tot;
  }/*endif*/
  /*----------------------------------------------------------------------*/
}/*end routine*/
/*========================================================================*/
