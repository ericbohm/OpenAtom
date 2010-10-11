//========================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//========================================================================
//                                                                          
//                         PI_MD:                                           
//             The future of simulation technology                          
//             ------------------------------------                         
//                     Module: control_inter_params                         
//                                                                          
// This reads in and sets up the electron-atom interaction pseudopotential  
//                                                                          
//========================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//========================================================================

#include "standard_include.h"
#include "../class_defs/allclass_gen.h"
#include "../class_defs/allclass_cp.h"
#include "../class_defs/allclass_mdintegrate.h"
#include "../class_defs/allclass_mdatoms.h"
#include "../class_defs/allclass_mdinter.h"
#include "../class_defs/allclass_mdintra.h"
#include "../class_defs/typedefs_par.h"
#include "../proto_defs/proto_inter_params_entry.h"
#include "../proto_defs/proto_inter_params_local.h"
#include "../proto_defs/proto_intra_params_local.h"
#include "../proto_defs/proto_search_entry.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_handle_entry.h"



//========================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//========================================================================

void control_inter_params(MDINTERACT *mdinteract,SPLINE_PARSE *spline_parse,
                          FILENAME_PARSE *filename_parse,
                          double alp_ewd,int ncharge,int lj_coul_explicit_opt,
                          int natm_tot,NAME atm_typ[],int iatm_typ[],
                          int iperd,int ishift_pot,double *tot_memory,
                          int int_res_ter)

//====================================================================

{ // begin routine
//====================================================================
//          Local variable declarations                                
  

  DATA_BASE_ENTRIES *inter_base;         // Lst: Database parameters  
  CLJ *clj_typ;
  double *eps,*sig;                      // Lst: Lennard-Jones params  
  double *awill,*bwill,*c6m,*c8m,*c10m;  // Lst: Williams params       
  double *cwill ,*rm_swit, *c9m;         // Lst: Aziz-chen params      
  double *temp_cutoff,*temp_cutoff_res,*temp_cutti;
  int *inter_label;

  int *ifound,*isearch,*igood;           // Lst: found,search goodness flags
  int ifound_explicit_lj;                // Num: Found explicit LJ flag 
  int i,j,iii;                           // Num: Counters               
  int ninter;                            // Num: Number of interactions 
  double now_mem;                        // Num: Memory allocated here  
  char typ[6];
  int nbase,nbase2,ibase_want;
  CATM_LAB *cinter,*cinter_base;
  char *fun_key;
  DICT_WORD *fun_dict;

  DICT_WORD *explicit_lj_dict,*explicit_lj_dict_tmp;
  int num_explicit_lj_dict;
  int nsearch,natm_srch;
  int num_fun_dict;
  int ifirst,ityp;
  int nsplin_mall_tot;
  int ninter_mall;
  int *ifound_lj;

  int ninter_unique;                     /* Num: number of interactions with 
                                             unqiue paramters */
  int ninter_unique_mall;

  int natm_typ = mdinteract->natm_typ;
 
//====================================================================
// 0) Write to screen                                                   

    ninter = natm_typ*(natm_typ + 1)/2;
    PRINTF("\n");
    PRINT_LINE_STAR;
    PRINTF("Searching the data bases (both user defined and default)\n");
    PRINTF("for the %d intermolecular interaction sets\n",ninter);
    PRINT_LINE_DASH;PRINTF("\n");

//====================================================================
//  I) Malloc the memory                                                

    clj_typ = (CLJ *) cmalloc(1*sizeof(CLJ),"control_inter_params");

  ninter = natm_typ*(natm_typ + 1)/2;
  ninter_mall = ninter;

  //  if((ninter_mall!=0)&&((ninter_mall %2)==0)){ninter_mall +=1;} // Maybe later

// Mallocs for explicit (grind it out) interactions (these should NOT be freed at the bottom)

  if(lj_coul_explicit_opt == 1){
    mdinteract->lj_sigma    = (double *) cmalloc(natm_typ*sizeof(double),"control_inter_params")-1;
    mdinteract->lj_epsilon  = (double *) cmalloc(natm_typ*sizeof(double),"control_inter_params")-1;
    mdinteract->cut_atm_typ = (double *) cmalloc(natm_typ*sizeof(double),"control_inter_params")-1;
    mdinteract->cut_atm_typ_res = (double *) cmalloc(natm_typ*sizeof(double),"control_inter_params")-1;
    ifound_lj = (int *) cmalloc(natm_typ*sizeof(int),"control_inter_params")-1;
  } else {
// Mallocs for spline interactions
    inter_label = (int *) cmalloc(ninter_mall*sizeof(int),"control_inter_params")-1;
    eps        = (double *) cmalloc(ninter_mall*sizeof(double),"control_inter_params")-1;
    sig        = (double *) cmalloc(ninter_mall*sizeof(double),"control_inter_params")-1;
    awill      = (double *) cmalloc(ninter_mall*sizeof(double),"control_inter_params")-1;
    bwill      = (double *) cmalloc(ninter_mall*sizeof(double),"control_inter_params")-1;
    cwill      = (double *) cmalloc(ninter_mall*sizeof(double),"control_inter_params")-1;
    rm_swit    = (double *) cmalloc(ninter_mall*sizeof(double),"control_inter_params")-1;
    c6m        = (double *) cmalloc(ninter_mall*sizeof(double),"control_inter_params")-1;
    c8m        = (double *) cmalloc(ninter_mall*sizeof(double),"control_inter_params")-1;
    c9m        = (double *) cmalloc(ninter_mall*sizeof(double),"control_inter_params")-1;
    c10m       = (double *) cmalloc(ninter_mall*sizeof(double),"control_inter_params")-1;
    cinter     = (CATM_LAB *)cmalloc(ninter*sizeof(CATM_LAB),"control_inter_params")-1;  
    ifound     = (int *)cmalloc(ninter*sizeof(int),"control_inter_params")-1;
    isearch    = (int *)cmalloc(ninter*sizeof(int),"control_inter_params")-1;
    igood      = (int *)cmalloc(ninter*sizeof(int),"control_inter_params")-1;
  
    temp_cutoff     = (double *) cmalloc(ninter_mall*sizeof(double),"control_inter_params")-1;
    temp_cutoff_res = (double *) cmalloc(ninter_mall*sizeof(double),"control_inter_params")-1;
    temp_cutti      = (double *) cmalloc(ninter_mall*sizeof(double),"control_inter_params")-1;

    mdinteract->cutoff     = (double *) cmalloc(ninter_mall*sizeof(double),"control_inter_params")-1;
    mdinteract->cutoff_res = (double *) cmalloc(ninter_mall*sizeof(double),"control_inter_params")-1;
    mdinteract->cutti      = (double *) cmalloc(ninter_mall*sizeof(double),"control_inter_params")-1;

    mdinteract->cutskin    = (double *) cmalloc(ninter_mall*sizeof(double),"control_inter_params")-1;
    mdinteract->cutskin_res= (double *) cmalloc(ninter_mall*sizeof(double),"control_inter_params")-1;
    mdinteract->cutskin_root   = (double *) cmalloc(ninter_mall*sizeof(double),"control_inter_params")-1;
    mdinteract->cutskin_root_res= (double *) cmalloc(ninter_mall*sizeof(double),"control_inter_params")-1;
  
    mdinteract->inter_map_index = (int *) cmalloc(ninter_mall*sizeof(int),"control_inter_params")-1;

  }// endif explicit LJ

  fun_key    = (char *)cmalloc(MAXWORD*sizeof(char),"control_inter_params");  

  ifirst =1;
  set_potfun_dict(&fun_dict,&num_fun_dict,ifirst);

  if(lj_coul_explicit_opt == 1){
    set_pot_explicit_lj_dict(&explicit_lj_dict,&num_explicit_lj_dict,ifirst);
    set_pot_explicit_lj_dict(&explicit_lj_dict_tmp,&num_explicit_lj_dict,ifirst);

    for(i=1;i<=natm_typ;i++){

//------------------------------------------------------------------------
//     A) First search the user defined data base                         

      ifound_explicit_lj = 0;

      strcpy(clj_typ->atm1,atm_typ[i]);

      if(strcasecmp(filename_parse->user_lj_name,"")!=0) {
         search_base_explicit_lj(filename_parse->user_lj_name,
                                 clj_typ,fun_dict,num_fun_dict,
                                 &explicit_lj_dict_tmp,explicit_lj_dict,
                                 num_explicit_lj_dict,&ifound_explicit_lj);
         if(ifound_explicit_lj==1){
            set_explicit_lj_params(explicit_lj_dict,
                          filename_parse->user_lj_name,fun_key,
                          &(mdinteract->lj_epsilon[i]),
                	  &(mdinteract->lj_sigma[i]),
                          &(mdinteract->cut_atm_typ[i]),
                          &(mdinteract->cut_atm_typ_res[i]));
         }//endif
     }//endif

//--------------------------------------------------------------------------
//     B) If you haven't found it search the default data base              

      if(ifound_explicit_lj == 0) {

         search_base_explicit_lj(filename_parse->def_lj_name,
                                 clj_typ,fun_dict,num_fun_dict,
                                 &explicit_lj_dict_tmp,explicit_lj_dict,
                                 num_explicit_lj_dict,&ifound_explicit_lj);

       if(ifound_explicit_lj==1){

	 set_explicit_lj_params(explicit_lj_dict,
                               filename_parse->def_lj_name,fun_key,
                               &(mdinteract->lj_epsilon[i]),
                	       &(mdinteract->lj_sigma[i]),
                               &(mdinteract->cut_atm_typ[i]),
                               &(mdinteract->cut_atm_typ_res[i]));

       }//endif
     }//endif
      ifound_lj[i] = (ifound_explicit_lj == 1 ? 0 : i);
     
    }// endfor atom types

      int ierror = 0;
    for(i=1; i <= natm_typ; i++){
      if(ifound_lj[i] != 0){
	PRINTF("The Lennard-Jones parameters were not found for atom %s in the database \n",atm_typ[i]);
        ierror = 1;
      }
    }//end for
    if(ierror){EXIT(1);}

// ------------------------------------------------------------
// Calculate Long-Range Correction for Explicit LJ 

    get_clong_lj_explicit(natm_tot,natm_typ,iatm_typ,
                          mdinteract->lj_epsilon,
                          mdinteract->lj_sigma,
                          mdinteract->cut_atm_typ,
                          mdinteract->cut_atm_typ_res,
                         &(mdinteract->clong),
                         &(mdinteract->clong_res),
                          mdinteract->iswit_vdw,
                          mdinteract->rheal_res);
// ------------------------------------------------------------

  } else {

  ityp = 0;
  for(i=1;i <= natm_typ; i++) {
    for(j=i;j <= natm_typ; j++) {
      ityp++;
      strcpy(cinter[ityp].atm1,atm_typ[i]);
      strcpy(cinter[ityp].atm2,atm_typ[j]);
    }//endfor
  }//endfor
  for(i=1;i<=ninter;i++){ifound[i]=0;}
  for(i=1;i<=ninter;i++){igood[i]=6;}

//====================================================================
/*  III) Search the user defined data base                              */

  natm_srch = 2;

  if(strlen(filename_parse->user_inter_name) != 0){
    nsearch = 1;
    ibase_want = 1;
    count_data_base(filename_parse->user_inter_name,fun_dict,num_fun_dict,
                    &nbase,ibase_want);
    if(nbase>0){
      nbase2 = 2*nbase;
      inter_base  = (DATA_BASE_ENTRIES *)
                       cmalloc(nbase2*sizeof(DATA_BASE_ENTRIES),"control_inter_params")-1;
      cinter_base = (CATM_LAB *)cmalloc(nbase2*sizeof(CATM_LAB),"control_inter_params")-1;
      read_data_base(filename_parse->user_inter_name,fun_dict,num_fun_dict,
                     inter_base,cinter_base,ibase_want,nbase);
      search_base(nbase,nbase2,cinter_base,ninter,cinter,igood,ifound,
                  isearch,nsearch,natm_srch,filename_parse->user_inter_name);
      assign_base_inter(inter_base,nbase,ifound,ninter,
                        sig,eps,awill,bwill,cwill,rm_swit,c6m,c8m,c9m,c10m,
                        inter_label,mdinteract->cutoff,mdinteract->cutoff_res,
                        mdinteract->cutti,
                        isearch,nsearch,cinter,cinter_base);
      cfree(&inter_base[1],"control_inter_params");
      cfree(&cinter_base[1],"control_inter_params");
    }//endif
  }//endif
//====================================================================
/*  IV) Search the default defined data base                            */

  if(strlen(filename_parse->def_inter_name) != 0){
    nsearch = 2;
    ibase_want = 1;
    count_data_base(filename_parse->def_inter_name,fun_dict,num_fun_dict,
                    &nbase,ibase_want);
    if(nbase>0){
      nbase2 = 2*nbase;
      inter_base = (DATA_BASE_ENTRIES *)
                    cmalloc(nbase*sizeof(DATA_BASE_ENTRIES),"control_inter_params")-1;
      cinter_base = (CATM_LAB *)cmalloc(nbase2*sizeof(CATM_LAB),"control_inter_params")-1;
      read_data_base(filename_parse->def_inter_name,fun_dict,num_fun_dict,
                     inter_base,cinter_base,ibase_want,nbase);
      search_base(nbase,nbase2,cinter_base,ninter,cinter,igood,ifound,
                  isearch,nsearch,natm_srch,filename_parse->def_inter_name);
      assign_base_inter(inter_base,nbase,ifound,ninter,
                        sig,eps,awill,bwill,cwill,rm_swit,c6m,c8m,c9m,c10m,
                        inter_label,mdinteract->cutoff,mdinteract->cutoff_res,
                        mdinteract->cutti,
                        isearch,nsearch,cinter,cinter_base);
      cfree(&inter_base[1],"control_inter_params");
      cfree(&cinter_base[1],"control_inter_params");
    }//endif
  }//endif

//====================================================================
// V) Check list for missing entries                                   

  strcpy(typ,"inter");
  atmlst_not_found(ninter,cinter,ifound,natm_srch,typ);

//====================================================================
//Find unique values of epsilon,sigma,rcut for LJ and null interactions
// This involves rearranging all the intermolecular interactions.      


 for(i=1; i<= ninter; i++){
   temp_cutoff[i]     = mdinteract->cutoff[i];
   temp_cutoff_res[i] = mdinteract->cutoff_res[i];
   temp_cutti[i]      = mdinteract->cutti[i];
 }//endfor

#define BYPASS_OFF
#ifdef BYPASS

 ninter_unique = ninter;
 for(i=1;i<=ninter;i++){
   mdinteract->inter_map_index[i] = i;
 }//endfor

#else
 sort_inter_params(eps,sig,
                   awill,bwill,cwill,
                   rm_swit,
                   c6m,c8m,c9m,c10m,
                   temp_cutoff,temp_cutoff_res,temp_cutti,
                   inter_label,mdinteract->inter_map_index,
                   &ninter_unique,ninter);
#endif


 ninter_unique_mall = ninter_unique;
 //  if((ninter_unique_mall!=0)&&((ninter_unique_mall %2)==0)){ninter_unique_mall++;} // Maybe later

  mdinteract->inter_map_mat  = cmall_int_mat(1,natm_typ,1,natm_typ,"control_inter_params.C");
  mdinteract->inter_map_mat0 = cmall_int_mat(1,natm_typ,1,natm_typ,"control_inter_params.C");

  ityp = 0;
  for(i=1;i <= natm_typ; i++) {
    for(j=i;j <= natm_typ; j++) {
      ityp++;
      mdinteract->inter_map_mat[i][j] = mdinteract->inter_map_index[ityp];
      mdinteract->inter_map_mat[j][i] = mdinteract->inter_map_index[ityp];
      mdinteract->inter_map_mat0[i][j] = ityp;
      mdinteract->inter_map_mat0[j][i] = ityp;
    }//endfor
  }//endfor

//====================================================================
// VI) Allocate spline arrays                                          

  mdinteract->nter_typ = ninter;
  mdinteract->nsplin_tot = mdinteract->nsplin*ninter;
  mdinteract->nsplin_tot_unique = ninter_unique*mdinteract->nsplin;
  mdinteract->ninter_unique = ninter_unique;
  now_mem = ((double)(mdinteract->nsplin*ninter_unique*(sizeof(double)*4
                    +sizeof(int)*0)+
                    + sizeof(int)*0)  +
                  ninter_unique*(sizeof(double)*4 + sizeof(int)*0 )
           )*1.e-06;

  *tot_memory += now_mem;
  
   PRINTF("There are %d unique interactions \n",ninter_unique);
   PRINTF("Intermolecular allocation: %g Mbytes; Total memory: %g Mbytes\n",
           now_mem,*tot_memory);

  nsplin_mall_tot = mdinteract->nsplin_tot_unique;

  //  if((nsplin_mall_tot!=0)&&((nsplin_mall_tot % 2)==0)){nsplin_mall_tot += 1;} // Maybe later

  mdinteract->cv0    = (double *) cmalloc(nsplin_mall_tot*sizeof(double),"control_inter_params")-1;
  mdinteract->cdv0   = (double *) cmalloc(nsplin_mall_tot*sizeof(double),"control_inter_params")-1;
  mdinteract->cv0_c  = (double *) cmalloc(nsplin_mall_tot*sizeof(double),"control_inter_params")-1;
  mdinteract->cdv0_c = (double *) cmalloc(nsplin_mall_tot*sizeof(double),"control_inter_params")-1;

  mdinteract->vcut_coul  = (double *) cmalloc(ninter_unique_mall*sizeof(double),"control_inter_params")-1;
  mdinteract->rmin_spl   = (double *) cmalloc(ninter_unique_mall*sizeof(double),"control_inter_params")-1;
  mdinteract->dr_spl     = (double *) cmalloc(ninter_unique_mall*sizeof(double),"control_inter_params")-1;
  mdinteract->dri_spl    = (double *) cmalloc(ninter_unique_mall*sizeof(double),"control_inter_params")-1;

//=====================================================================
/* VII) Set up the splines for the real-space intermolecular potential   */
/*     energy and forces                                                 */

/* Spline the intermolecular interaction */

 set_inter_splin(sig,eps,awill,bwill,cwill,rm_swit,c6m,c8m,c9m,c10m,
                  alp_ewd,ishift_pot,inter_label,
                  spline_parse,mdinteract,
                  temp_cutti,temp_cutoff,
                  ninter_unique,ncharge,iperd);

//=====================================================================
/*  VIII) Get the long range correction                                   */

  mdinteract->clong = 0.0;
  mdinteract->clong_res = 0.0;      
  if(iperd == 3) {

    get_clong(natm_tot,natm_typ,ninter,iatm_typ,c6m,
              mdinteract->inter_map_index,
             &(mdinteract->clong),
             &(mdinteract->clong_res),temp_cutoff,temp_cutoff_res,
             mdinteract->iswit_vdw,mdinteract->rheal_res);

  } /* endif */
  
   PRINTF("Dispersion long range parameter %.15g \n",mdinteract->clong);
   if(int_res_ter==1){
    PRINTF("Dispersion long range parameter(RESPA) %g \n",
         mdinteract->clong_res);
   }//endif

  } // endif explicit LJ opt 

/*=======================================================================*/

  if(lj_coul_explicit_opt == 0){
    cfree(&inter_label[1],"control_inter_params");
    cfree(&eps[1],"control_inter_params");
    cfree(&sig[1],"control_inter_params");
    cfree(&awill[1],"control_inter_params");
    cfree(&bwill[1],"control_inter_params");
    cfree(&cwill[1],"control_inter_params");
    cfree(&rm_swit[1],"control_inter_params");
    cfree(&c6m[1],"control_inter_params");
    cfree(&c8m[1],"control_inter_params");
    cfree(&c9m[1],"control_inter_params");
    cfree(&c10m[1],"control_inter_params");
    cfree(&temp_cutoff[1],"control_inter_params");
    cfree(&temp_cutoff_res[1],"control_inter_params");
    cfree(&temp_cutti[1],"control_inter_params");

    cfree(&cinter[1],"control_inter_params");
    cfree(&ifound[1],"control_inter_params");
    cfree(&isearch[1],"control_inter_params");
    cfree(&igood[1],"control_inter_params");
  }else{
    cfree(clj_typ,"control_inter_params");
    cfree(&ifound_lj[1],"control_inter_params");
  }// endif

  cfree(&fun_dict[1],"control_inter_params");
  cfree(fun_key,"control_inter_params");

//=====================================================================
/* X) Write to screen                                                    */

  PRINTF("\n");
  PRINT_LINE_DASH;
  PRINTF("Completed the data bases searches\n");
  PRINT_LINE_STAR;
  PRINTF("\n");

//--------------------------------------------------------------------
} /* end routine */
//========================================================================



//========================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//========================================================================

void inter_coef(DICT_WORD *dict,char filename[],char fun_key[],
               DATA_BASE_ENTRIES *inter_base,CATM_LAB *cinter_base,int ibase)

//========================================================================
/*               Begin subprogram:                                          */
      {/*begin routine*/
//========================================================================
/*               Local variable declarations:                               */

  int index,iii;
  double bohr2,bohr6,bohr8,bohr9,bohr10;
  double cutti,cutoff,cutoff_res;
  double sig,eps;
  double c6m,c8m,c9m,c10m,awill,bwill,cwill,rm_swit;

//========================================================================
/* I) Fill atom types and label part of the data base     */

  strcpy(cinter_base[ibase].atm1,dict[1].keyarg);
  strcpy(cinter_base[ibase].atm2,dict[2].keyarg);
  strcpy(cinter_base[ibase].label,"");

//=====================================================================
/* II) Set up */

  inter_base[ibase].inter_label = -1;
  bohr2  = BOHR*BOHR;
  bohr6  = bohr2*bohr2*bohr2;
  bohr8  = bohr6*bohr2;
  bohr9  = bohr8*BOHR;
  bohr10 = bohr8*bohr2;

//=====================================================================
/* III) Convert and assign cutoffs */
  
  sscanf(dict[4].keyarg,"%lf",&cutti);
  sscanf(dict[5].keyarg,"%lf",&cutoff);
  sscanf(dict[6].keyarg,"%lf",&cutoff_res);
  if(cutti<0){
    index=4;
    keyarg_barf(dict,filename,fun_key,index);
  }//endif
  if(cutoff<0){
    index=5;
    keyarg_barf(dict,filename,fun_key,index);
  }//endif
  if(cutoff_res<0){
    index=6;
    keyarg_barf(dict,filename,fun_key,index);
  }//endif
  inter_base[ibase].cutti      = cutti/BOHR;
  inter_base[ibase].cutoff     = cutoff/BOHR;
  inter_base[ibase].cutoff_res = cutoff_res/BOHR;

//=====================================================================
/* IV) Convert and assign Lennard-Jones                                  */
  
  if(strcasecmp(dict[3].keyarg,"lennard-Jones") == 0) {
    inter_base[ibase].inter_label = 2;
    sscanf(dict[7].keyarg,"%lg",&sig);
    sscanf(dict[8].keyarg,"%lg",&eps);
    if(sig<0){
      index=7;
      keyarg_barf(dict,filename,fun_key,index);
    }//endif
    if(eps<0){
      index=8;
      keyarg_barf(dict,filename,fun_key,index);
    }//endif
    eps /= BOLTZ;
    sig /= BOHR;
    inter_base[ibase].eps        = eps;
    inter_base[ibase].sig        = sig;
    inter_base[ibase].c6m        = (4.0*eps)*(sig*sig*sig*sig*sig*sig);
  }//endif

//=====================================================================
/*  3) Convert and assign Williams                                        */
  
  if(strcasecmp(dict[3].keyarg,"williams") == 0) {
    inter_base[ibase].inter_label = 3;
    sscanf(dict[9].keyarg,"%lg",&c6m);
    sscanf(dict[10].keyarg,"%lg",&c8m);
    sscanf(dict[11].keyarg,"%lg",&c10m);
    sscanf(dict[12].keyarg,"%lg",&awill);
    sscanf(dict[13].keyarg,"%lg",&bwill);
    if(awill<0){
      index=12;
      keyarg_barf(dict,filename,fun_key,index);
    }//endif
    if(bwill<0){
      index=13;
      keyarg_barf(dict,filename,fun_key,index);
    }//endif
    inter_base[ibase].awill      = awill/BOLTZ;
    inter_base[ibase].bwill      = bwill*BOHR;
    inter_base[ibase].c6m        = c6m/(BOLTZ*bohr6);
    inter_base[ibase].c8m        = c8m/(BOLTZ*bohr8);
    inter_base[ibase].c10m       = c10m/(BOLTZ*bohr10);
  }//endif

//=====================================================================
/*  5) Convert and assign Null                                           */
  
  if(strcasecmp(dict[3].keyarg,"null") == 0) {
    inter_base[ibase].inter_label = 4;
    inter_base[ibase].c6m = 0.0;
  }//endif

//=====================================================================
/*  4) Convert and assign Aziz-Chen                                      */
  
  if(strcasecmp(dict[3].keyarg,"aziz-chen") == 0) {
    inter_base[ibase].inter_label = 6;
    sscanf(dict[9].keyarg,"%lg",&c6m);
    sscanf(dict[10].keyarg,"%lg",&c8m);
    sscanf(dict[16].keyarg,"%lg",&c9m);
    sscanf(dict[11].keyarg,"%lg",&c10m);
    sscanf(dict[12].keyarg,"%lg",&awill);
    sscanf(dict[13].keyarg,"%lg",&bwill);
    sscanf(dict[14].keyarg,"%lg",&cwill);
    sscanf(dict[15].keyarg,"%lg",&rm_swit);
    if(awill<0){
      index=12;
      keyarg_barf(dict,filename,fun_key,index);
    }//endif
    if(bwill<0){
      index=13;
      keyarg_barf(dict,filename,fun_key,index);
    }//endif
    if(cwill<0){
      index=14;
      keyarg_barf(dict,filename,fun_key,index);
    }//endif
    if(rm_swit<0){
      index=15;
      keyarg_barf(dict,filename,fun_key,index);
    }//endif
    inter_base[ibase].awill      = awill/BOLTZ;
    inter_base[ibase].bwill      = bwill*BOHR;
    inter_base[ibase].cwill      = cwill*bohr2;
    inter_base[ibase].rm_swit    = rm_swit/BOHR; 
    inter_base[ibase].c6m        = c6m/(BOLTZ*bohr6);
    inter_base[ibase].c8m        = c8m/(BOLTZ*bohr8);
    inter_base[ibase].c9m        = c9m/(BOLTZ*bohr9);
    inter_base[ibase].c10m       = c10m/(BOLTZ*bohr10);
  }//endif

//=====================================================================
/*  6) Check Potential type                                              */

  if(inter_base[ibase].inter_label==-1){
    index=3;
    keyarg_barf(dict,filename,fun_key,index);
  }//endif


//------------------------------------------------------------------------
   }/*end routine*/
//========================================================================



//========================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//========================================================================
/* set_inter_splin:This subroutine fits the necessary terms to a spline     */
//------------------------------------------------------------------------

void assign_base_inter(DATA_BASE_ENTRIES *inter_base,int nbase,int *ifound,
                      int ninter,
                       double *sig, double *eps,double *awill,double *bwill,
                       double *cwill,double *rm_swit,double *c6m,
                       double *c8m,double *c9m,double *c10m,
                       int *inter_label,double *cutoff,
                       double *cutoff_res,
                       double *cutti,
                       int *isearch,int nsearch,CATM_LAB *cinter,
                       CATM_LAB *cinter_base)

//=====================================================================
/*            Begin subprogram:                                          */
  {/*begin routine*/
//=====================================================================
/*             Local variable declarations                                */
  int ibase,i,iii;
//=====================================================================

  for(i=1;i<=ninter;i++){
    if(ifound[i] > 0 && isearch[i]==nsearch){
      ibase = ifound[i];

      eps[i]         = 0.0;
      sig[i]         = 0.0;
      awill[i]       = 0.0;
      bwill[i]       = 0.0;
      cwill[i]       = 0.0;
      rm_swit[i]     = 0.0;
      c8m[i]         = 0.0;  
      c9m[i]         = 0.0;  
      c10m[i]        = 0.0;  

      inter_label[i] = inter_base[ibase].inter_label;
      cutti[i]       = inter_base[ibase].cutti;
      cutoff[i]      = inter_base[ibase].cutoff;     
      cutoff_res[i]  = inter_base[ibase].cutoff_res;

      switch(inter_label[i]){
        case 2:  eps[i]     = inter_base[ibase].eps;
                 sig[i]     = inter_base[ibase].sig;
                 c6m[i]     = inter_base[ibase].c6m;  
               break;
        case 3:  awill[i]   = inter_base[ibase].awill;
                 bwill[i]   = inter_base[ibase].bwill;
                 c6m[i]     = inter_base[ibase].c6m;
                 c8m[i]     = inter_base[ibase].c8m;
                 c10m[i]    = inter_base[ibase].c10m; 
               break;
        case 4:  c6m[i]     =  0.0;  
               break;
        case 6:  awill[i]   = inter_base[ibase].awill;
                 bwill[i]   = inter_base[ibase].bwill;
                 cwill[i]   = inter_base[ibase].cwill;
                 rm_swit[i] = inter_base[ibase].rm_swit;
                 c6m[i]     = inter_base[ibase].c6m;
                 c8m[i]     = inter_base[ibase].c8m;
                 c9m[i]     = inter_base[ibase].c9m;
                 c10m[i]    = inter_base[ibase].c10m; 
               break;
      }/*endswitch*/
    }//endif
  }//endfor

//------------------------------------------------------------------------
   }/*end routine*/
//========================================================================
//========================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//========================================================================
/* set_inter_splin:This subroutine fits the necessary terms to a spline     */
//------------------------------------------------------------------------

void set_inter_splin(double sig[],double eps[],double a[],double b[],
                     double c[],double rm[],
                     double c6[],double c8[],double c9[],double c10[],
                     double alp_ewd,int ishift,int inter_label[],
                     SPLINE_PARSE *spline_parse,
                     MDINTERACT *mdinteract,
                     double cutti[],double cutoff[],
                     int ninter,int ncharge,int iperd)

//=====================================================================
/*            Begin subprogram:                                          */
  {/*begin routine*/
//=====================================================================
/*             Local variable declarations                                */

  double sigt,epst;                      /* Num: LJ potential params    */
  double at,bt,c6t,c8t,c9t,c10t;             /* Num: Williams pot params    */
  double ct,rmt;                         /* Num: Aziz-Chen pot params   */
  int ioff;                              /* Num: Offset for spline      */
  int i,iii;                             /* Num: For loop counter       */
  int itype;                             /* Num: Potential type label   */
  int iiperd;                            /* Num: Periodicity            */
  int ishift_now;
  double rmax_c,rmax;                    /* Num: Potential cutoff radii */
  double qijp;                           /* Num: product of charges     */
  double dr_tmp,rmin;
  double dri_tmp;

  double *rmin_spl  = mdinteract->rmin_spl;
  double *vcut_coul = mdinteract->vcut_coul;

//====================================================================
/*  I) Spline each interaction                                          */

  ioff = 0;
  qijp = 0.0;
  rmax_c = 0.0;

  for(i=1;i <= ninter;i++) {

    itype       = inter_label[i];
    rmin_spl[i] = cutti[i];
    rmax        = cutoff[i];
    rmax_c      = MAX(rmax,rmax_c);
    iiperd      = iperd;

    spline_vdv(mdinteract->rmin_spl[i],rmax,&(mdinteract->cv0)[ioff],
             &(mdinteract->cdv0)[ioff],mdinteract->nsplin,&(mdinteract->dr_spl)[i],
             &(mdinteract->dri_spl)[i],sig[i],eps[i],a[i],b[i],c[i],rm[i],
             c6[i],c8[i],c9[i],c10[i],alp_ewd,qijp,iiperd,itype,ishift,
             mdinteract->rheal_res,mdinteract->dielectric_opt,
             mdinteract->dielectric_rheal,mdinteract->dielectric_cut,
             mdinteract->dielectric_eps);

    ioff += mdinteract->nsplin;
  } /* endfor:spline each interaction */

  mdinteract->cutoff_max = rmax_c;

//====================================================================
/*  II) Zero the coulomb shift                                          */

  if(ishift == 0) {
    for(i=1;i <= ninter;i++) {
      vcut_coul[i] = 0.0;
    }/* endfor */
  }/* endif: shift off */
  
//====================================================================
/*  III) Spline the ewald coulomb stuff                                 */

  if (ncharge > 0 ){
   if( iperd >0 || mdinteract->dielectric_opt==1 ){
    itype = 1;
    qijp = 1.0;
    iiperd = iperd;
    sigt = 1.0;
    epst = 0.0;
    at = 0.0;
    bt = 0.0;
    ct = 0.0;
    rmt = 0.0;
    c6t = 0.0;
    c8t = 0.0;
    c9t = 0.0;
    c10t = 0.0;
    if(iperd >0){
     if((alp_ewd*rmax_c)<3.4){
       PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       PRINTF("Ewald convergence parameter too small !\n");
       PRINTF("%g*rmax requested, at least 3.4*rmax required\n", 
             (alp_ewd*rmax_c));
       PRINTF("If this criteria is not met, energy will drift.\n");
       PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       FFLUSH(stdout);
      }//endif
      /* EXIT(1);*/
    }/*endif for iperd>0*/    
    ioff = 0;
    ishift_now = ishift;
    for(i=1;i <= ninter;i++) {
       rmin = cutti[i];
       rmax = cutoff[i];
       iiperd = iperd;
       spline_vdv(rmin,rmax,&(mdinteract->cv0_c)[ioff],&(mdinteract->cdv0_c)[ioff],
             mdinteract->nsplin,&dr_tmp,&dri_tmp,sigt,epst,at,bt,ct,rmt,c6t,
             c8t,c9t,c10t,alp_ewd,qijp,iiperd,itype,ishift_now,
             mdinteract->rheal_res,mdinteract->dielectric_opt,
             mdinteract->dielectric_rheal,mdinteract->dielectric_cut,
             mdinteract->dielectric_eps);
       ioff += mdinteract->nsplin;
    }//endfor
   } /* endif:spline special coulomb */
  }/*endif for ncharge*/

//====================================================================
/*  V) Shift straight coulomb                                           */

 if(ncharge > 0 && iperd == 0 &&  ishift == 1) {
     for(i=1;i <= ninter;i++) {
       vcut_coul[i] = 1.0/cutoff[i];
     }/* endfor */
 }/* endif */

//------------------------------------------------------------------
  }/* end routine */
//========================================================================


//========================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//========================================================================
/* sort_inter_params:This subroutine sorts intermolecular parameters        */
/*  and finds the unique eps,sigma,and rcut for lj/null parameters          */
//------------------------------------------------------------------------

void sort_inter_params(double *eps,double *sig,
                       double *awill,double *bwill,double *cwill,
                       double *rm_swit,
                       double *c6m,double *c8m,double *c9m,double *c10m,
                       double *cutoff,double *cutoff_res,double *cutti,
                       int *inter_label,int *inter_map_index,
                       int *pninter_unique,int ninter)

//=====================================================================
/*            Begin subprogram:                                          */
{/*begin routine*/
//=====================================================================
/*             Local variable declarations                                */

  double *epst,*sigt;  
  double *awillt,*bwillt,*cwillt,*rm_switt;
  double *c6mt,*c8mt,*c9mt,*c10mt;
  double *cutofft,*cutoff_rest,*cuttit;
  int *ilabelt;

  int nunique,ioff,icount,index;
  int nlj,nnull,n3,n5,n6,ntot,nother;
  int i,j,k,iii;
  int jlow,jup,iup;

//=====================================================================
//-----------------------------------------------
/*I) Malloc local inter parameters and assign them */
//-----------------------------------------------

  epst     = (double *) cmalloc(ninter*sizeof(double),"control_inter_params")-1; 
  sigt     = (double *) cmalloc(ninter*sizeof(double),"control_inter_params")-1; 
  awillt   = (double *) cmalloc(ninter*sizeof(double),"control_inter_params")-1; 
  bwillt   = (double *) cmalloc(ninter*sizeof(double),"control_inter_params")-1; 
  cwillt   = (double *) cmalloc(ninter*sizeof(double),"control_inter_params")-1; 
  rm_switt = (double *) cmalloc(ninter*sizeof(double),"control_inter_params")-1; 

  c6mt     = (double *) cmalloc(ninter*sizeof(double),"control_inter_params")-1; 
  c8mt     = (double *) cmalloc(ninter*sizeof(double),"control_inter_params")-1; 
  c9mt     = (double *) cmalloc(ninter*sizeof(double),"control_inter_params")-1; 
  c10mt    = (double *) cmalloc(ninter*sizeof(double),"control_inter_params")-1; 

  cutofft     = (double *) cmalloc(ninter*sizeof(double),"control_inter_params")-1; 
  cutoff_rest = (double *) cmalloc(ninter*sizeof(double),"control_inter_params")-1; 
  cuttit      = (double *) cmalloc(ninter*sizeof(double),"control_inter_params")-1; 

  ilabelt = (int    *) cmalloc(ninter*sizeof(int),"control_inter_params")-1; 

  for(i=1; i<= ninter; i++){
    epst[i]     = eps[i];
    sigt[i]     = sig[i];
    awillt[i]   = awill[i];
    bwillt[i]   = bwill[i];
    cwillt[i]   = cwill[i];
    rm_switt[i] = rm_swit[i];
    c6mt[i]     = c6m[i];
    c8mt[i]     = c8m[i];
    c9mt[i]     = c9m[i];
    c10mt[i]    = c10m[i];

    cutofft[i]     = cutoff[i];
    cutoff_rest[i] = cutoff_res[i];
    cuttit[i]      = cutti[i];
    ilabelt[i]     = inter_label[i];
  }

  for(i=1; i<= ninter; i++){
   inter_map_index[i] = i;
  }

//=====================================================================
//---------------------------------------------------------------
/*II) Count up the number of different intermolecular interactions */ 
//---------------------------------------------------------------

  nlj   = 0; 
  n3    = 0;
  nnull = 0; 
  n5    = 0;
  n6    = 0;
 for(i=1; i<= ninter; i++){
  switch(inter_label[i]){
   case 2: nlj++; break; 
   case 3: n3++;  break;
   case 4: nnull++;  break;
   case 5: n5++;  break;
   case 6: n6++;  break;
  }/*endswitch*/
 }//endfor

   ntot = nlj + n3 + nnull + n5 + n6;

  if( ntot != ninter){
   PRINTF("@@@@@@@@@@@@@@@@@@@@@ ERROR @@@@@@@@@@@@@@@@@@@@@@ \n");
   PRINTF("TOTAL NUMBER OF INTERACTIONS NOT EQUAL TO NINTER\n");
   PRINTF("CONTACT TECHNICAL SUPPORT \n");
   PRINTF("@@@@@@@@@@@@@@@@@@@@@ ERROR @@@@@@@@@@@@@@@@@@@@@@ \n");
   EXIT(1);
  }//endif

//=====================================================================
//-----------------------------------------------
/* III) SORT  intermolecular  interactions by type */
//-----------------------------------------------

//-------------------------------
/* 1. inter_label=6 Aziz-Chen last */
//-------------------------------

 for(i = ninter; i > 1; i--){
  if(ilabelt[i] != 6){
   for(j=i-1; j >=  1; j--){
     if(ilabelt[j] == 6){
       switchij(&epst[i],&epst[j]);
       switchij(&sigt[i],&sigt[j]);
       switchij(&awillt[i],&awillt[j]);
       switchij(&bwillt[i],&bwillt[j]);
       switchij(&cwillt[i],&cwillt[j]);
       switchij(&rm_switt[i],&rm_switt[j]);
       switchij(&c6mt[i],&c6mt[j]);
       switchij(&c8mt[i],&c8mt[j]);
       switchij(&c9mt[i],&c9mt[j]);
       switchij(&c10mt[i],&c10mt[j]);
       switchij(&cutofft[i],&cutofft[j]);
       switchij(&cutoff_rest[i],&cutoff_rest[j]);
       switchij(&cuttit[i],&cuttit[j]);
       iswitchij(&ilabelt[i],&ilabelt[j]);
       break;
     }//endif
   }//endfor
  }//endif
 }//endfor

//-------------------------------------------
/* 2. inter_label=5 Williams-LJ second to last */
//-------------------------------------------

   ioff = n6;
   iup  = ninter - ioff; 
 for(i = iup; i > 1; i--){
  if(ilabelt[i] != 5){
   for(j=i-1; j >=  1; j--){
     if(ilabelt[j] == 5){
       switchij(&epst[i],&epst[j]);
       switchij(&sigt[i],&sigt[j]);
       switchij(&awillt[i],&awillt[j]);
       switchij(&bwillt[i],&bwillt[j]);
       switchij(&cwillt[i],&cwillt[j]);
       switchij(&rm_switt[i],&rm_switt[j]);
       switchij(&c6mt[i],&c6mt[j]);
       switchij(&c8mt[i],&c8mt[j]);
       switchij(&c9mt[i],&c9mt[j]);
       switchij(&c10mt[i],&c10mt[j]);
       switchij(&cutofft[i],&cutofft[j]);
       switchij(&cutoff_rest[i],&cutoff_rest[j]);
       switchij(&cuttit[i],&cuttit[j]);
       iswitchij(&ilabelt[i],&ilabelt[j]);
       break;
     }//endif
   }//endfor
  }//endif
 }//endfor

//----------------------------------------
/* 3. inter_label= 3 Williams third to last */
//----------------------------------------

   ioff = n6 + n5;
   iup  = ninter - ioff; 
 for(i = iup; i > 1; i--){
  if(ilabelt[i] != 3){
   for(j=i-1; j >=  1; j--){
     if(ilabelt[j] == 3){
       switchij(&epst[i],&epst[j]);
       switchij(&sigt[i],&sigt[j]);
       switchij(&awillt[i],&awillt[j]);
       switchij(&bwillt[i],&bwillt[j]);
       switchij(&cwillt[i],&cwillt[j]);
       switchij(&rm_switt[i],&rm_switt[j]);
       switchij(&c6mt[i],&c6mt[j]);
       switchij(&c8mt[i],&c8mt[j]);
       switchij(&c9mt[i],&c9mt[j]);
       switchij(&c10mt[i],&c10mt[j]);
       switchij(&cutofft[i],&cutofft[j]);
       switchij(&cutoff_rest[i],&cutoff_rest[j]);
       switchij(&cuttit[i],&cuttit[j]);
       iswitchij(&ilabelt[i],&ilabelt[j]);
       break;
     }//endif
   }//endfor
  }//endif
 }//endfor

//------------------------------------------
/* 2. PLACE inter_label=4 null fourth to last */
//------------------------------------------

   ioff = n6 + n5 + n3;
   iup  = ninter - ioff; 
 for(i = iup; i > 1; i--){
  if(ilabelt[i] != 4){
   for(j=i-1; j >=  1; j--){
     if(ilabelt[j] == 4){
       switchij(&epst[i],&epst[j]);
       switchij(&sigt[i],&sigt[j]);
       switchij(&awillt[i],&awillt[j]);
       switchij(&bwillt[i],&bwillt[j]);
       switchij(&cwillt[i],&cwillt[j]);
       switchij(&rm_switt[i],&rm_switt[j]);
       switchij(&c6mt[i],&c6mt[j]);
       switchij(&c8mt[i],&c8mt[j]);
       switchij(&c9mt[i],&c9mt[j]);
       switchij(&c10mt[i],&c10mt[j]);
       switchij(&cutofft[i],&cutofft[j]);
       switchij(&cutoff_rest[i],&cutoff_rest[j]);
       switchij(&cuttit[i],&cuttit[j]);
       iswitchij(&ilabelt[i],&ilabelt[j]);
       break;
     }//endif
   }//endfor
  }//endif
 }//endfor

 
//=====================================================================
//------------------------------------------
/*IV) SORT LJ and NULL Parameters             */
/*    same eps sig rcut in order              */
//------------------------------------------

 ioff   = n3+n5+n6; 
 i = 1;
 while(i< ninter-ioff){
  k = i;
  for(j=i+1; j <=  ninter-ioff; j++){
    /* Bubble down matching interactions*/
    if( (epst[i]        == epst[j]) && 
        (sigt[i]        == sigt[j]) &&
        (cutofft[i]     == cutofft[j]) &&
        (cutti[i]       == cutti[j]) && 
        (cutoff_rest[i] == cutoff_rest[j])){
      k++;
      switchij(&epst[k],&epst[j]);
      switchij(&sigt[k],&sigt[j]);
      switchij(&awillt[k],&awillt[j]);
      switchij(&bwillt[k],&bwillt[j]);
      switchij(&cwillt[k],&cwillt[j]);
      switchij(&rm_switt[k],&rm_switt[j]);
      switchij(&c6mt[k],&c6mt[j]);
      switchij(&c8mt[k],&c8mt[j]);
      switchij(&c9mt[k],&c9mt[j]);
      switchij(&c10mt[k],&c10mt[j]);
      switchij(&cutofft[k],&cutofft[j]);
      switchij(&cutoff_rest[k],&cutoff_rest[j]);
      switchij(&cuttit[k],&cuttit[j]);
      iswitchij(&ilabelt[k],&ilabelt[j]);
    }//endif
  }//endfor
  i = k+1;      /* The next unique guy is at k+1 */
 }/*end while */

//=====================================================================
//------------------------------------------
/*V) Determine number of unique LJ parameters */
/*    and reassign in new order               */
//------------------------------------------


/*Assign unique LJ/Null parameters */

 nother = n3+n5+n6;
 i = 1;
 nunique = 0;
 while(i<=ninter-ioff){
   nunique++;
   j = i;
   while( ((epst[i]        == epst[j])        &&
           (sigt[i]        == sigt[j])        &&
           (cutofft[i]     == cutofft[j])     &&
           (cutoff_rest[i] == cutoff_rest[j]) &&
           (cuttit[i]      == cuttit[j])     )){
      epst[nunique]        = epst[j];
      sigt[nunique]        = sigt[j];
      awillt[nunique]      = awillt[j];
      bwillt[nunique]      = bwillt[j];
      cwillt[nunique]      = cwillt[j];
      rm_switt[nunique]    = rm_switt[j];
      c6mt[nunique]        = c6mt[j];
      c8mt[nunique]        = c8mt[j];
      c9mt[nunique]        = c9mt[j];
      c10mt[nunique]       = c10mt[j];

      cutofft[nunique]     = cutofft[j];
      cutoff_rest[nunique] = cutoff_rest[j];
      cuttit[nunique]      = cuttit[j];
      ilabelt[nunique]     = ilabelt[j];
      j++;
      if(j==ninter-ioff+1){break;}
    }/*while*/
    i = j;
  }//endfor


/*Reassign all other inter parameters*/
     icount = 0;

 for(i=ninter-nother+1; i<= ninter; i++){
   icount++; 
   index = nunique+icount;
     epst[index]      = epst[i];
     sigt[index]      = sigt[i];
   awillt[index]      = awillt[i];
   bwillt[index]      = bwillt[i];
   cwillt[index]      = cwillt[i];
   rm_switt[index]    = rm_switt[i];
     c6mt[index]      = c6mt[i];
     c8mt[index]      = c8mt[i];
     c9mt[index]      = c9mt[i];
    c10mt[index]      = c10mt[i];
   cutofft[index]     = cutofft[i];
   cutoff_rest[index] = cutoff_rest[i];
   cuttit[index]      = cuttit[i];
   ilabelt[index]     = ilabelt[i];

 }//endfor


   nunique += nother;
  *pninter_unique = nunique;

  if(nunique>ninter){
   PRINTF("@@@@@@@@@@@@@@@@@@@@@ ERROR @@@@@@@@@@@@@@@@@@@@@@ \n");
   PRINTF("Total number of unique interactions > then ninter!!\n");
   PRINTF("Contact technical support : %d vs %d\n",ninter,nunique);
   PRINTF("@@@@@@@@@@@@@@@@@@@@@ ERROR @@@@@@@@@@@@@@@@@@@@@@ \n");
   EXIT(1);
  }

  if(nunique==0 && ninter!=0){
   PRINTF("@@@@@@@@@@@@@@@@@@@@@ ERROR @@@@@@@@@@@@@@@@@@@@@@ \n");
   PRINTF("Total number of unique interactions  = 0!!\n");
   PRINTF("Contact technical support");
   PRINTF("@@@@@@@@@@@@@@@@@@@@@ ERROR @@@@@@@@@@@@@@@@@@@@@@ \n");
   EXIT(1);
  }


//=====================================================================
//--------------------------------------------
/*VI) Assign lj_ind base on unique eps sig rcut */
//--------------------------------------------

  for(i=1; i<= ninter; i++){
   inter_map_index[i] = 0;
  }

  for(i=1; i<= ninter; i++){
   if(inter_label[i] == 2 || 
      inter_label[i] == 4){

    for(j=1; j<= (nunique-nother); j++){
     if( (eps[i] == epst[j]) && 
         (sig[i] == sigt[j]) && 
         (cutoff[i] == cutofft[j]) &&
         (cutti[i] == cuttit[j]) &&
         (cutoff_res[i] == cutoff_rest[j]) ){
       inter_map_index[i] = j; 
     }//endif 
    }//endfor

   }//endif
  }//endfor

//---------------------------------------------------------------
/* Assign  Williams  compare awill bwill c6m c8m c10m and cutoffs  */
//---------------------------------------------------------------


   jlow = nunique - nother + 1;
   jup  = jlow + n3;

  for(i=1; i<= ninter; i++){
   if(inter_label[i] == 3 ){
    for(j=jlow; j<= jup; j++){ 
     if((awill[i] == awillt[j]) && 
        (bwill[i] == bwillt[j]) &&
        (c6m[i]==c6mt[j]) && 
        (c8m[i]==c8mt[j]) && 
        (c10m[i] == c10mt[j]) &&
        (cutoff[i] == cutofft[j]) && 
        (cutti[i] == cuttit[j]) &&
        (cutoff_res[i] == cutoff_res[j]) ){
       inter_map_index[i] = j; 
     }//endif 
    }//endfor
   }//endif
  }//endfor


//-------------------------------------------------------------------------
// Assign  Williams-LJ compare awill bwill c6m c8m c10m eps sig and cutoffs 
//-------------------------------------------------------------------------
   jlow = nunique - nother + n3 + 1; 
   jup  = jlow + n5;

  for(i=1; i<= ninter; i++){
   if(inter_label[i] == 5 ){
    for(j=jlow; j<= jup; j++){
     if( (eps[i] == epst[j]) && 
         (sig[i] == sigt[j]) && 
         (awill[i] == awillt[j]) && 
         (bwill[i] == bwillt[j]) &&
         (c6m[i] == c6mt[j]) && 
         (c8m[i] == c8mt[j]) && 
         (c10m[i] == c10mt[j]) &&
         (cutoff[i] == cutofft[j]) && 
         (cutti[i] == cuttit[j]) &&
         (cutoff_res[i] == cutoff_rest[j])  ){
       inter_map_index[i] = j; 
     }//endif 
    }//endfor
   }//endif
  }//endfor

//------------------------------------------------------------------------
// Assign Aziz-Chen compare awill bwill cwill rm_swit c6m c8m c9m c10m  
// and cutoffs  
//------------------------------------------------------------------------

   jlow = nunique - nother + n3 + n5 + 1; 
   jup  = jlow + n6;


  for(i=1; i<= ninter; i++){
   if(inter_label[i] == 6 ){

    for(j=jlow; j<= jup; j++){
     if((awill[i] == awillt[j]) &&
        (bwill[i] == bwillt[j]) &&
        (cwill[i] == cwillt[j]) && 
        (rm_swit[i] == rm_switt[j]) && 
        (c6m[i] == c6mt[j]) && 
        (c8m[i] == c8mt[j]) && 
        (c9m[i] == c9mt[j]) && 
        (c10m[i] == c10mt[j]) &&
        (cutoff[i] == cutofft[j]) && 
        (cutti[i] == cuttit[j]) &&
        (cutoff_res[i] == cutoff_rest[j])  ){
          inter_map_index[i] = j; 
      }//endif 
    }//endfor

   }//endif
  }//endfor

//---------------------------------------------
//  Reassign the unique atom types               
//---------------------------------------------

   for(i=1; i<= nunique; i++){
     eps[i]        = epst[i];
     sig[i]        = sigt[i];
     awill[i]      = awillt[i];
     bwill[i]      = bwillt[i];
     cwill[i]      = cwillt[i];
     rm_swit[i]    = rm_switt[i];
     c6m[i]        = c6mt[i];
     c8m[i]        = c8mt[i];
     c9m[i]        = c9mt[i];
     c10m[i]       = c10mt[i];
     cutoff[i]     = cutofft[i];
     cutti[i]      = cuttit[i];
     cutoff_res[i] = cutoff_rest[i];
    inter_label[i] = ilabelt[i];
   }//endfor

//---------------------------------------------
// free locally assigned memory                 
//---------------------------------------------

  cfree(&epst[1],"sort_inter_params");
  cfree(&sigt[1],"sort_inter_params");
  cfree(&awillt[1],"sort_inter_params");
  cfree(&bwillt[1],"sort_inter_params");
  cfree(&cwillt[1],"sort_inter_params");
  cfree(&rm_switt[1],"sort_inter_params");

  cfree(&c6mt[1],"sort_inter_params");
  cfree(&c8mt[1],"sort_inter_params");
  cfree(&c9mt[1],"sort_inter_params");
  cfree(&c10mt[1],"sort_inter_params");

  cfree(&cutofft[1],"sort_inter_params");
  cfree(&cutoff_rest[1],"sort_inter_params");
  cfree(&cuttit[1],"sort_inter_params");

  cfree(&ilabelt[1],"sort_inter_params");

//========================================================================
}// end routine
//========================================================================



//========================================================================
  void switchij(double *valuei,double *valuej)
//========================================================================
    {// begin routine
//========================================================================
 double temp;

  temp     = *valuei;
  *valuei  = *valuej;
  *valuej  = temp;

//========================================================================
  }// end routine
//========================================================================



//========================================================================
  void iswitchij(int *valuei,int *valuej)
//========================================================================
  {// begin routine
//========================================================================
 int temp;

  temp     = *valuei;
  *valuei  = *valuej;
  *valuej  = temp;
 
//========================================================================
  }// end routine
//========================================================================



//========================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//========================================================================

void set_explicit_lj_params(DICT_WORD explicit_lj_dict[],
                            char *filename,char *fun_key,
                            double *epsilon,double *sigma,
                            double *rcut,double *rcut_res)

//========================================================================
//               Begin subprogram:                                          
 { // begin routine

  double temp;
  int index;

//------------------------------------------------------------------------
//  0) Set up                                                            

         strcpy(fun_key,"explicit_lj_parm");

//------------------------------------------------------------------------
//  1) Assign epsilon                                                    

	 sscanf(explicit_lj_dict[2].keyarg,"%lg",&temp);
	 *epsilon = (double) temp;
         *epsilon /= BOLTZ;
         if(*epsilon < 0) {
            index = 2;
            keyarg_barf(explicit_lj_dict,filename,fun_key,index);
         }

//------------------------------------------------------------------------
//  2) Assign sigma                                                     

	 sscanf(explicit_lj_dict[3].keyarg,"%lg",&temp);
	 *sigma = (double) temp;
         *sigma /= BOHR;
         if(*sigma < 0) {
            index = 3;
            keyarg_barf(explicit_lj_dict,filename,fun_key,index);
         }

//------------------------------------------------------------------------
//  3) Assign rcut                                                      

	 sscanf(explicit_lj_dict[4].keyarg,"%lg",&temp);
	 *rcut = (double) temp;
         *rcut /= BOHR;
         if(*rcut < 0) {
            index = 4;
            keyarg_barf(explicit_lj_dict,filename,fun_key,index);
         }
//------------------------------------------------------------------------
//  5) Assign rcut_res                                                    

	 sscanf(explicit_lj_dict[5].keyarg,"%lg",&temp);
	 *rcut_res = (double) temp;
         *rcut_res /= BOHR;
         if(*rcut_res < 0) {
            index = 5;
            keyarg_barf(explicit_lj_dict,filename,fun_key,index);
         }

//========================================================================
 }/*end routine*/
//========================================================================
