/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: fetch_residue.c                              */
/*                                                                          */
/* This subprogram reads in molecular and residue parameter files           */
/* and sets molecular and intramolecular data sets                          */
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
#include "../proto_defs/proto_handle_entry.h"
#include "../proto_defs/proto_friend_lib_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void check_parmfile(char filename[],int *num_fun_dict,
                    DICT_WORD *fun_dict[],char fun_key[],
                    int nres_expt, int *nres_bond,
                    int mol_or_res)

/*==========================================================================*/
{/*begin Routine */
/*==========================================================================*/
/*     Local Variables  */

  int nmol_name_now,natom_now,nbond_now,nbend_now,ntors_now;
  int nonfo_now,nbend_bnd_now,nres_def_now,nres_name_now,nres_bond_now;
  int nres_morph_now,natom_morph_now;
  int ifirst,nline=1,num,nfun_key=0,iii;
  FILE *fp;

/*==========================================================================*/
/* I) Count the stuff */

  nmol_name_now = 0;
  natom_now     = 0;
  nbond_now     = 0;
  nbend_now     = 0;
  ntors_now     = 0;
  nonfo_now     = 0;
  nbend_bnd_now = 0;
  nres_def_now  = 0;
  nres_name_now = 0;
  nres_morph_now = 0;
  natom_morph_now = 0;
  nres_bond_now = 0;

  fp = cfopen((const char*) filename,"r");
  ifirst=0;
  set_intra_fun_dict(fun_dict,num_fun_dict,ifirst);
  while(get_fun_key_cnt(fp,fun_key,&nline,&nfun_key,filename)){
    get_fun_key_index(fun_key,*num_fun_dict,*fun_dict,nline,nfun_key,
                      filename,&num);
    if(!strcasecmp(fun_key,"molecule_name_def"))     nmol_name_now++;
    if(!strcasecmp(fun_key,"atom_def"))         natom_now++;
    if(!strcasecmp(fun_key,"atom_create")) {    natom_now++;natom_morph_now++;}
    if(!strcasecmp(fun_key,"atom_destroy")){    natom_now++;natom_morph_now++;}
    if(!strcasecmp(fun_key,"atom_morph"))  {    natom_now++;natom_morph_now++;}
    if(!strcasecmp(fun_key,"bond_def"))         nbond_now++;
    if(!strcasecmp(fun_key,"bend_def"))         nbend_bnd_now++;
    if(!strcasecmp(fun_key,"torsion_def"))      ntors_now++;
    if(!strcasecmp(fun_key,"onefour_def"))      nonfo_now++;
    if(!strcasecmp(fun_key,"bend_bnd_def"))     nbend_bnd_now++;
    if(!strcasecmp(fun_key,"residue_def"))      nres_def_now++;
    if(!strcasecmp(fun_key,"residue_name_def")) nres_name_now++;
    if(!strcasecmp(fun_key,"residue_morph")){   nres_name_now++;
                                                nres_morph_now++;}
    if(!strcasecmp(fun_key,"residue_bond_def")) nres_bond_now++;
  }/*endwhile*/
  fclose(fp);

  (*nres_bond) = nres_bond_now;

/*==========================================================================*/
/* II) Error check  */


  if(nmol_name_now !=0 && mol_or_res == 2){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Functional keywords %s not expected \n","molecule_name_def");
    PRINTF("in file %s\n",filename);
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }

  if(nmol_name_now != 1 && mol_or_res == 1){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");    
    PRINTF("Molecule name not specified in file %s\n",filename);
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");    
    FFLUSH(stdout);
    EXIT(1);
  }

  if(nres_name_now !=0 && mol_or_res == 1){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Functional keywords %s not expected \n","res_name_def");
    PRINTF("in file %s\n",filename);
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }

  if(nres_name_now != 1 && mol_or_res == 2){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");    
    PRINTF("Residue name not specified in file %s\n",filename);
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");    
    FFLUSH(stdout);
    EXIT(1);
  }

  if(nres_def_now > 0 || nres_bond_now > 0){
    if(natom_now > 0 || nbond_now > 0 || nbend_now > 0 || ntors_now > 0 ||
       nonfo_now > 0 || nbend_bnd_now > 0){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("You have specified atoms,bonds,bends,torsions,onefour \n");
      PRINTF("or bond_bnd with residue definitions or bonding\n");
      PRINTF("in file %s\n",filename);
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }
  }

  if(natom_now == 0 && nres_def_now == 0){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("No atoms or residues specified in file %s\n",filename);
    PRINTF("in file %s\n",filename);
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }
  
  if(nres_def_now == 0 && nres_bond_now != 0){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Residue bond specified without any residues\n");
    PRINTF("in file %s\n",filename);
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");    
    FFLUSH(stdout);
    EXIT(1);
  }

  if(nres_expt != nres_def_now && mol_or_res ==1){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");    
    PRINTF("Number of residue's found %d not equal to number expected %d\n",
           nres_def_now,nres_expt);
    PRINTF("in file %s\n",filename);
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");    
    FFLUSH(stdout);
    EXIT(1);
  }

  if(nres_name_now == 1 && nres_def_now > 0 && nres_bond_now > 0){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");    
    PRINTF("Nested residue definitions not allowed in file %s\n",filename);
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");    
    FFLUSH(stdout);
    EXIT(1);
  }

  if(nres_morph_now !=1 && natom_morph_now>0){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");    
    PRINTF("Atom modifications not permitted in  %s\n",filename);
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");    
    FFLUSH(stdout);
    EXIT(1);
  }


/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void fetch_molname(char filename[],DICT_INTRA *dict_intra,MDATOM_MAPS *atommaps,
               char fun_key[],int jmol_typ,int nres_expt)

/*==========================================================================*/
{ /* begin routine */
/*==========================================================================*/

  int num,nline=1,nfun_key=0,nkey;
  int ifirst,ifound=0;
  FILE *fp;

/*==========================================================================*/

  ifirst=0;
  set_intra_fun_dict(&(dict_intra->fun_dict),&(dict_intra->num_fun_dict),
                     ifirst);
  fp = cfopen((const char*) filename,"r");
  while(get_fun_key(fp,fun_key,&nline,&nfun_key,filename)){
    get_fun_key_index(fun_key,dict_intra->num_fun_dict,
                      dict_intra->fun_dict,nline,nfun_key,
                      filename,&num);

    if(num == 1){
      nkey = 0;
      set_mol_name_dict(&(dict_intra->mol_name_dict),
                        &(dict_intra->num_mol_name_dict),ifirst);
      while(get_word(fp,dict_intra->word,&nline,&nkey,nfun_key,filename)){
        put_word_dict(dict_intra->word,dict_intra->mol_name_dict,
                      dict_intra->num_mol_name_dict,
                  fun_key,nline,nkey,nfun_key,filename);
      } /* end while */
      set_mol_name_params(dict_intra->mol_name_dict,
                          dict_intra->num_mol_name_dict,
                          fun_key,filename,atommaps->mol_typ,
                          atommaps->natm_1mol_jmol_typ,jmol_typ,nres_expt);
      ifound = 1;
      break; /* kick out of while loop */
    } else {
      close_fun_key_cnt(fp,fun_key,&nline,nfun_key,filename);
    } /* end if*/
  } /*end while */

  if(ifound != 1){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Required functional keyword %s not found\n",
          dict_intra->fun_dict[1].keyword);
    PRINTF("in the %dth molecule type parmeter file %s\n",jmol_typ,filename);
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }

  fclose(fp);
/*--------------------------------------------------------------------------*/
} /*end routine*/
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void fetch_residue_defs(char filename[],DICT_INTRA *dict_intra,
                        MDATOM_MAPS *atommaps,FILENAME_PARSE *filename_parse,
                        char fun_key[],int nres_now,int jmol_typ,
                        BUILD_INTRA *build_intra)

/*==========================================================================*/
{/*begin routine */
/*==========================================================================*/

  int nline=1,nfun_key=0,num;
  int i,ifirst,ires_off,nkey;
  FILE *fp;

/*==========================================================================*/
/* I) Get the residues                                                      */

  for(i=1;i<=nres_now;i++){build_intra->ires_ind_chk[i]=0;}
  ires_off = atommaps->jres_jmol_typ_strt[jmol_typ] - 1;
  ifirst=0;
  set_intra_fun_dict(&(dict_intra->fun_dict),&(dict_intra->num_fun_dict),
                     ifirst);
  fp = cfopen((const char*) filename,"r");
  while(get_fun_key(fp,fun_key,&nline,&nfun_key,filename)){
    get_fun_key_index(fun_key,dict_intra->num_fun_dict,
                      dict_intra->fun_dict,nline,nfun_key,
                      filename,&num);
    if(num == 12){
      nkey = 0;
      set_res_def_dict(&(dict_intra->res_def_dict),
                        &(dict_intra->num_res_def_dict),ifirst);
      while(get_word(fp,dict_intra->word,&nline,&nkey,nfun_key,filename)){
        put_word_dict(dict_intra->word,dict_intra->res_def_dict,
                      dict_intra->num_res_def_dict,
                      fun_key,nline,nkey,nfun_key,filename);
      } /* end while */
      set_res_def_params(filename,fun_key,dict_intra->res_def_dict,
                   dict_intra->num_res_def_dict,
                   atommaps,filename_parse->res_param_name,
                   nres_now,build_intra,ires_off,
                         filename_parse);
    } else {
      close_fun_key_cnt(fp,fun_key,&nline,nfun_key,filename);
    } /* end if*/
  } /*end while */

/*==========================================================================*/
/* II) Check the residue indicies                                           */

  for(i=1;i<=nres_now;i++){
    if(build_intra->ires_ind_chk[i]!=1){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("Residue number %d specified %d times in file %s \n",
            i,build_intra->ires_ind_chk[i],filename);
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }/*endif*/
  }/*endfor*/
  fclose(fp);
/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void fetch_residue_def0(MDATOM_MAPS *atommaps,BUILD_INTRA *build_intra, 
                        int jmol_typ)

/*==========================================================================*/
 {/*begin routine*/

int ires_off;

/*==========================================================================*/

   atommaps->nres_typ++;
   if(atommaps->nres_typ > build_intra->nres_typ_max){
      build_intra->nres_typ_max+=NMEM_MIN;
      atommaps->res_typ = (NAME *)
         crealloc(&(atommaps->res_typ[1]),
                    (build_intra->nres_typ_max)*sizeof(NAME),"fetch_residue")-1;
   }/*endif*/
   strcpy(atommaps->res_typ[atommaps->nres_typ],atommaps->mol_typ[jmol_typ]);
   ires_off = atommaps->jres_jmol_typ_strt[jmol_typ]-1;
   atommaps->ires_typ_jres_jmol_typ[1+ires_off] = atommaps->nres_typ;
/*==========================================================================*/
}/*end routine*/
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void fetch_residue_bonds(char filename[],DICT_INTRA *dict_intra,
                         char fun_key[],RESBOND_PRM *resbond_prm,
                         MDATOM_MAPS *atommaps,int nres_now,
                         int nres_bond_now,int jmol_typ,int pi_beads)

/*==========================================================================*/
{/*begin routine */
/*==========================================================================*/

  int ifirst,nline=1,nkey,nfun_key=0,num,ires_bond,ires_off;
  FILE *fp;

/*==========================================================================*/
/* I) Get the residue bonds                                                */

  ires_bond = 0;
  ires_off  = atommaps->jres_jmol_typ_strt[jmol_typ]-1;
  ifirst=0;
  set_intra_fun_dict(&(dict_intra->fun_dict),&(dict_intra->num_fun_dict),
                     ifirst);
  fp = cfopen((const char*) filename,"r");
  while(get_fun_key(fp,fun_key,&nline,&nfun_key,filename)){
    get_fun_key_index(fun_key,dict_intra->num_fun_dict,
                      dict_intra->fun_dict,nline,nfun_key,
                      filename,&num);
    if(num == 13){
      nkey = 0;
      set_res_bond_dict(&(dict_intra->res_bond_dict),
                        &(dict_intra->num_res_bond_dict),ifirst);
      while(get_word(fp,dict_intra->word,&nline,&nkey,nfun_key,filename)){
        put_word_dict(dict_intra->word,dict_intra->res_bond_dict,
                      dict_intra->num_res_bond_dict,
                      fun_key,nline,nkey,nfun_key,filename);
      } /* end while */
      ires_bond++;
      set_res_bond_params(filename,fun_key,ires_bond,dict_intra->res_bond_dict,
                   dict_intra->num_res_bond_dict,resbond_prm,nres_now,
                   atommaps,ires_off,pi_beads);
    } else {
      close_fun_key_cnt(fp,fun_key,&nline,nfun_key,filename);
    } /* end if*/
  } /*end while */

/*==========================================================================*/
/* II) Error check */

  if(ires_bond != nres_bond_now){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("Internal error fetching the residue bonds in file %s\n",
            filename);
      PRINTF("Contact technical support\n");
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
  }/*endif*/
  fclose(fp);
/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void fetch_residue_name(char filename[],DICT_INTRA *dict_intra,
                   MDATOM_MAPS *atommaps,char fun_key[],BUILD_INTRA *build_intra,
                   int jmol_typ,int jres, int jres_off, int iparm,int mol_only)

/*==========================================================================*/
{ /* begin routine */
/*==========================================================================*/

  int num,nline=1,nfun_key=0,nkey;
  int ifirst,ifound=0;
  FILE *fp;

/*==========================================================================*/

  ifound = 0;

  ifirst=0;
  set_intra_fun_dict(&(dict_intra->fun_dict),&(dict_intra->num_fun_dict),
                     ifirst);
  fp = cfopen((const char*) filename,"r");
  while(get_fun_key(fp,fun_key,&nline,&nfun_key,filename)){
    get_fun_key_index(fun_key,dict_intra->num_fun_dict,
                      dict_intra->fun_dict,nline,nfun_key,
                      filename,&num);
    nkey = 0;
    if(num==1||num==14||num==7){
         set_res_name_dict(&(dict_intra->res_name_dict),
                        &(dict_intra->num_res_name_dict),ifirst);
         while(get_word(fp,dict_intra->word,&nline,&nkey,nfun_key,filename)){
               put_word_dict(dict_intra->word,dict_intra->res_name_dict,
                             dict_intra->num_res_name_dict,
                             fun_key,nline,nkey,nfun_key,filename);
         } /*end while*/
         if(num==1 && iparm==1 && mol_only==1){
            set_res_name_params(dict_intra->res_name_dict,
                                dict_intra->num_res_name_dict,
                                fun_key,filename,atommaps->res_typ,
                                atommaps->natm_jres_jmol_typ,jmol_typ,jres,
                                jres_off,atommaps->ires_typ_jres_jmol_typ);
            ifound++;
         }/*endif*/
         if(num==14 && iparm==1){
            set_res_name_params(dict_intra->res_name_dict,
                                dict_intra->num_res_name_dict,
                                fun_key,filename,atommaps->res_typ,
                                atommaps->natm_jres_jmol_typ,jmol_typ,jres,
                                jres_off,atommaps->ires_typ_jres_jmol_typ);
           ifound++;
         }/*endif*/
         if(num==7 && iparm==0){
            set_res_morph_params(dict_intra->res_name_dict,
                                 dict_intra->num_res_name_dict,
                                 fun_key,filename,atommaps->res_typ,
                                 atommaps->natm_jres_jmol_typ,jmol_typ,jres,
                                 jres_off,atommaps->ires_typ_jres_jmol_typ);
            ifound++;
         }/*endif*/
    }else{
      close_fun_key_cnt(fp,fun_key,&nline,nfun_key,filename);
    }/*endif*/
  } /*end while */

  if(ifound != 1){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    if(iparm==1){
      PRINTF("Required functional keyword ~%s not found or repeated \n",
              dict_intra->fun_dict[14].keyword);
      PRINTF("in the residue parameter file %s\n",filename);
    }else{
      PRINTF("Required functional keyword ~%s not found or repeated \n",
              dict_intra->fun_dict[7].keyword);
      PRINTF("in the residue morph file %s\n",filename);
    }
    PRINTF("of the %d residue of the %d molecule\n",jres,jmol_typ);
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
}
  fclose(fp);
/*--------------------------------------------------------------------------*/
} /*end routine*/
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void fetch_residue_atm_prms(char filename[],DICT_INTRA *dict_intra,
                     CPATOM_MAPS *cpatom_maps,MDATOM_MAPS *atommaps,
                     MDCLATOMS_INFO *clatoms_info,MDGHOST_ATOMS *ghost_atoms,
                     MDCONSTRNT *mdconstrnt,char fun_key[],
                     BUILD_INTRA *build_intra,int jmol_typ,int jres,
                     int jres_off,int iparm)

/*==========================================================================*/
    { /* begin routine */
/*==========================================================================*/

  int num,nline=1,nfun_key=0,nkey,ifirst,iii;
  FILE *fp;

/*                                                                          */
/*==========================================================================*/
/* II) Loop over the functional keywords                                    */

  ifirst=0;
  set_intra_fun_dict(&(dict_intra->fun_dict),&(dict_intra->num_fun_dict),
                     ifirst);
  fp = cfopen((const char*) filename,"r");
  while(get_fun_key(fp,fun_key,&nline,&nfun_key,filename)){
    get_fun_key_index(fun_key,dict_intra->num_fun_dict,
                      dict_intra->fun_dict,
                      nline,nfun_key,filename,&num);
    if(num==2||num==9||num==10){
      set_atm_dict(&(dict_intra->atm_dict),&(dict_intra->num_atm_dict),ifirst);
      while(get_word(fp,dict_intra->word,&nline,&nkey,nfun_key,filename)){
              put_word_dict(dict_intra->word,dict_intra->atm_dict,
                            dict_intra->num_atm_dict,
                            fun_key,nline,nkey,nfun_key,filename);
      }/*endwhile*/
      if(num==2||num==9){
        set_atm_params(dict_intra->atm_dict,dict_intra->num_atm_dict,
                       fun_key,filename,jmol_typ,jres,jres_off,
                       cpatom_maps,atommaps,build_intra,clatoms_info,
                       ghost_atoms,mdconstrnt);
      }/*endif*/
      if(num==10){
        set_atm_morph(dict_intra->atm_dict,dict_intra->num_atm_dict,
                      fun_key,filename,jmol_typ,cpatom_maps,atommaps,
                      build_intra,clatoms_info,ghost_atoms);
      }
    }else{
      close_fun_key_cnt(fp,fun_key,&nline,nfun_key,filename);
    }/*endif*/

  }/*endwhile*/
  fclose(fp);
/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/






/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void fetch_residue_atm_masks(char filename[],DICT_INTRA *dict_intra,
                      MDATOM_MAPS *atommaps,char fun_key[],
                      BUILD_INTRA *build_intra,int jmol_typ,int jres, 
                      int jres_off,int iparm)

/*==========================================================================*/
{ /* begin routine */
/*==========================================================================*/

  int num,nline=1,nfun_key=0,nkey,ifirst;
  int iii;
  FILE *fp;

/*==========================================================================*/
/* II) Loop over the functional keywords                                    */

  ifirst=0;
  set_intra_fun_dict(&(dict_intra->fun_dict),&(dict_intra->num_fun_dict),
                     ifirst);
  fp = cfopen((const char*) filename,"r");
  while(get_fun_key(fp,fun_key,&nline,&nfun_key,filename)){
    get_fun_key_index(fun_key,dict_intra->num_fun_dict,dict_intra->fun_dict,
                      nline,nfun_key,filename,&num);
    if(num==10&&iparm==1){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("Atom morphs not permitted \n");
      PRINTF("in the residue parameter file %s\n",filename);
      PRINTF("of the %d residue of the %d molecule\n",jres,jmol_typ);
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }
    if(num==2||num==8||num==9){
      set_atm_dict(&(dict_intra->atm_dict),&(dict_intra->num_atm_dict),ifirst);
      while(get_word(fp,dict_intra->word,&nline,&nkey,nfun_key,filename)){
              put_word_dict(dict_intra->word,dict_intra->atm_dict,
                            dict_intra->num_atm_dict,
                            fun_key,nline,nkey,nfun_key,filename);
      }/*endwhile*/
      if(iparm==1){
        set_atm_mask(dict_intra->atm_dict,
                     dict_intra->num_atm_dict,fun_key,filename,
                     build_intra,num,jres,jmol_typ);
      }else{
        set_atm_mask_rb(dict_intra->atm_dict,
                        dict_intra->num_atm_dict,fun_key,filename,
                        build_intra,num,jres,jmol_typ);
      }/*endif*/
    }else{
      close_fun_key_cnt(fp,fun_key,&nline,nfun_key,filename);
    }/*endif*/

  }/*endwhile*/
  fclose(fp);
/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void fetch_residue_connectivity(char filename[],DICT_INTRA *dict_intra,
             MDATOM_MAPS *atommaps, char fun_key[],
             MDCLATOMS_INFO *clatoms_info, BUILD_INTRA *build_intra,
             MDINTRA *mdintra,NULL_INTER_PARSE *null_inter_parse,
             int jmol_typ,int jres, int jres_off,int iparm,
             int mol_hydrog_con_opt)

/*==========================================================================*/
  { /* begin routine */
/*==========================================================================*/

#include "../class_defs/allclass_strip_mdintra.h"   

  int num,nline=1,nfun_key=0,nkey,ifirst,iii;
  int ires_off;
  FILE *fp;

/*==========================================================================*/
/* II) Loop over the functional keywords                                    */

  ires_off = atommaps->jres_jmol_typ_strt[jmol_typ]-1;
  ifirst=0;
  set_intra_fun_dict(&(dict_intra->fun_dict),&(dict_intra->num_fun_dict),
                     ifirst);
  fp = cfopen((const char*) filename,"r");
  while(get_fun_key(fp,fun_key,&nline,&nfun_key,filename)){
    get_fun_key_index(fun_key,dict_intra->num_fun_dict,dict_intra->fun_dict,
                      nline,nfun_key,filename,&num);
    if(num==3||num==4||num==5||num==6||num==11||num==15){
       set_intra_dict(&(dict_intra->intra_dict),
                      &(dict_intra->num_intra_dict),ifirst);
       while(get_word(fp,dict_intra->word,&nline,&nkey,nfun_key,filename)){
               put_word_dict(dict_intra->word,dict_intra->intra_dict,
                            dict_intra->num_intra_dict,
                             fun_key,nline,nkey,nfun_key,filename);
       }/*endwhile*/

       if(num==3){
         set_bond_params(dict_intra->intra_dict,dict_intra->num_intra_dict,
                         fun_key,filename,jmol_typ,mdconstrnt,mdghost_atoms,
                         clatoms_info,atommaps,mdbond,null_inter_parse,
                         build_intra,jres,ires_off,mol_hydrog_con_opt);
       }/*endif*/

       if(num==5){
         set_tors_params(dict_intra->intra_dict,dict_intra->num_intra_dict,
                          fun_key,filename,jmol_typ,clatoms_info,atommaps,
                          mdtors,null_inter_parse,build_intra,jres,ires_off);
       }/*endif*/
   
       if(num==6){
         set_onfo_params(dict_intra->intra_dict,dict_intra->num_intra_dict,
                          fun_key,filename,jmol_typ,clatoms_info,atommaps,
                          mdonfo,null_inter_parse,build_intra,jres,ires_off);
       }/*endif*/

       if(num==11 || num==4){
         set_bend_bnd_params(dict_intra->intra_dict,dict_intra->num_intra_dict,
                             fun_key,filename,jmol_typ,clatoms_info,atommaps,
                             mdbend_bnd,mdbend,null_inter_parse,
                             build_intra,jres,ires_off);
       }/*endif*/

       if(num==15){
         set_grp_bond_params(dict_intra->intra_dict,dict_intra->num_intra_dict,
                             fun_key,filename,jmol_typ,clatoms_info,atommaps,
                             mdghost_atoms,mdgrp_bond_con,mdgrp_bond_watts,
                             mdconstrnt,null_inter_parse,
                             build_intra,jres,ires_off);
       }/*endif*/

    }else{
      close_fun_key_cnt(fp,fun_key,&nline,nfun_key,filename);
    }/*endif*/

  }/*endwhile*/
  fclose(fp);
/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/













