/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: interface_hand                               */
/*                                                                          */
/* Subprograms that search the data bases                                   */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../class_defs/typedefs_par.h"
#include "../proto_defs/proto_search_entry.h"
#include "../proto_defs/proto_search_local.h"
#include "../proto_defs/proto_handle_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* search_base_vps: Searches vps data bases                                 */
/*==========================================================================*/

void search_base_vps(char filename[],CVPS *cvps_typ,
    DICT_WORD fun_dict[],int num_fun_dict,
    DICT_WORD *vps_dict_tmp[],
    DICT_WORD vps_dict[],int num_vps_dict_ret,int *ifound)

  /*==========================================================================*/
{/*begin routine*/
  /*=======================================================================*/
  /*             Local variable declarations                               */

  int nline,nkey,i,num;
  int ifirst,nfun_key;
  int num_vps_dict = num_vps_dict_ret;
  NAME fun_key;
  DICT_WORD word;
  FILE *fp;

  /*========================================================================*/

  fp       = cfopen((const char *)filename,"r");
  *ifound  = 0;
  ifirst   = 0;
  nline    = 0;
  nkey     = 0;
  nfun_key = 0;
  while(get_fun_key(fp,fun_key,&nline,&nfun_key,filename)){
    get_fun_key_index(fun_key,num_fun_dict,fun_dict,nline,
        nfun_key,filename,&num);
    if(num==6){
      set_potvps_dict(vps_dict_tmp,&num_vps_dict,ifirst);
      while(get_word(fp,&word,&nline,&nkey,nfun_key,filename)){
        put_word_dict(&word,*vps_dict_tmp,num_vps_dict,
            fun_key,nline,nkey,nfun_key,filename);
      }/*endwhile*/
      if(strcasecmp((*vps_dict_tmp)[1].keyarg,cvps_typ->atm1)==0){
        *ifound = 1;
        dict_save(*vps_dict_tmp,vps_dict,num_vps_dict);
      }
    }else{
      close_fun_key_cnt(fp,fun_key,&nline,nfun_key,filename);      
    }/*endif*/
  }/*end while*/
  fclose(fp);
  if(*ifound==1){
    for(i=1;i<=num_vps_dict;i++){
      if((vps_dict[i].iuset==0)&&(vps_dict[i].key_type==1)){
        keyword_miss(vps_dict,filename,fun_key,i);}
    }/*endfor*/
  }/*endif*/

  /*==========================================================================*/
} /*end routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* search_base_vps: Searches vps data bases                                 */
/*==========================================================================*/
void search_base_grimme(char filename[],CVPS *cgrimme_typ,
    DICT_WORD fun_dict[],int num_fun_dict,
    DICT_WORD *grimme_dict_tmp[],
    DICT_WORD grimme_dict[],int num_grimme_dict_ret,int *ifound)

  /*==========================================================================*/
{/*begin routine*/
  /*=======================================================================*/
  /*             Local variable declarations                               */

  int nline,nkey,i,num;
  int ifirst,nfun_key;
  int num_grimme_dict = num_grimme_dict_ret;
  NAME fun_key;
  DICT_WORD word;
  FILE *fp;

  /*========================================================================*/

  fp       = cfopen((const char *)filename,"r");
  *ifound  = 0;
  ifirst   = 0;
  nline    = 0;
  nkey     = 0;
  nfun_key = 0;
  while(get_fun_key(fp,fun_key,&nline,&nfun_key,filename)){
    get_fun_key_index(fun_key,num_fun_dict,fun_dict,nline,
        nfun_key,filename,&num);
    if(num==10){
      set_grimme_dict(grimme_dict_tmp,&num_grimme_dict,ifirst);
      while(get_word(fp,&word,&nline,&nkey,nfun_key,filename)){
        put_word_dict(&word,*grimme_dict_tmp,num_grimme_dict,
            fun_key,nline,nkey,nfun_key,filename);
      }/*endwhile*/
      if(strcasecmp((*grimme_dict_tmp)[1].keyarg,cgrimme_typ->atm1)==0){
        *ifound = 1;
        dict_save(*grimme_dict_tmp,grimme_dict,num_grimme_dict);
      }
    }else{
      close_fun_key_cnt(fp,fun_key,&nline,nfun_key,filename);      
    }/*endif*/
  }/*end while*/
  fclose(fp);
  if(*ifound==1){
    for(i=1;i<=num_grimme_dict;i++){
      if((grimme_dict[i].iuset==0)&&(grimme_dict[i].key_type==1)){
        keyword_miss(grimme_dict,filename,fun_key,i);}
    }/*endfor*/
  }/*endif*/

/*==========================================================================*/
  } /*end routine*/
/*==========================================================================*/



