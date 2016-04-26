/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: control_vps_params                           */
/*                                                                          */
/* This reads in and sets up the electron-atom interaction pseudopotential  */ 
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
#include "../proto_defs/proto_vps_params_entry.h"
#include "../proto_defs/proto_search_entry.h"
#include "../proto_defs/proto_intra_params_local.h"
#include "../proto_defs/proto_vps_params_local.h"
#include "../proto_defs/proto_handle_entry.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "proto_math.h"

#define DEBUG_DKNY_OFF

#define JUERG_FACTOR_ON
#ifdef  JUERG_FACTOR_ON
#define JUERG_FACTOR 0.72
#else
#define JUERG_FACTOR 1.0
#endif

typedef struct vps_file{
  char name[MAXWORD];
}VPS_FILE;

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void control_vps_params(CPPSEUDO *cppseudo,GENCELL *gencell,
    FILENAME_PARSE *filename_parse,
    SPLINE_PARSE *spline_parse,int *iatm_atm_typ,
    int natm_typ,NAME *atm_typ,
    double *tot_memory,int natm_tot,int natm_ab_init,
    int cp_ptens_calc,int cp_dual_grid_opt,double ecut_cp)

  /*==========================================================================*/
  /*               Begin subprogram:                                          */
{/*begin routine*/
  /*==========================================================================*/
  /*               Local variable declarations                                */
  int i;                       /* Num:  For loop counters             */
  int ifound;                  /* Num:  Data base match flag          */
  int ishift,ishift2,ishift3;  /* Num:  Angular momentum shifts       */
  int ngh_now;                 /* Num:  Number of ngh points read in  */
  int ngh_max;                 /* Num:  Max number of ngh points      */

  VPS_FILE *vps_file;          /* Fle:  Pseudopotential file          */
  double now_mem;              /* Num:  Current memory usage          */

  /* Dictionary memory */
  char *filename;              /* Char: temp file name                */
  CVPS *cvps_typ;
  char *fun_key;
  DICT_WORD *word;
  DICT_WORD *fun_dict;
  int num_fun_dict;
  DICT_WORD *vps_dict,*vps_dict_tmp;
  int num_vps_dict,ifirst,iii;
  int natm_typ_mall,natm_mall,nsplin_mall,norm_mall,nlist_mall;
  double dummy0,dummy1,dummy2,dummy3;
  int nmall_gh,natm_typ_gh = 0;   /* starting malloc value for gauss-hermite */
  /* Volume and grid variables */
  double alpha_conv_dual = cppseudo->alpha_conv_dual;
  double vol_cp   = gencell->vol_cp;

  int nsplin_g,n_rad_max,n_ang_max1;
  int ngh_tot,nsplin_g_tot,nlist,norm_size;
  int natm_nonloc_now;
  double *q_typ = cppseudo->q_typ;

  /*==========================================================================*/
  /* 0) Output                                                                */

  PRINTF("\n");
  PRINT_LINE_STAR
    PRINTF("Searching the data bases both user defined and default\n");
  PRINTF("for the electron-atom pseudopotentials\n");
  PRINT_LINE_DASH;PRINTF("\n");

  /*==========================================================================*/
  /* I) Toast bad spline points : (note conversion back to Ry)              */

  if( (cppseudo->nsplin_g< 4000) && (2.0*ecut_cp > 60.0) ){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Now, dude, lets clean up those files and use a \n");
    PRINTF("reasonable number of psuedo spline points at large \n");
    PRINTF("cutoffs. Its clear, %d points at %g Ry, won't do!\n",
        cppseudo->nsplin_g,2.0*ecut_cp);
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }/*endif*/

  /*==========================================================================*/
  /* II) Convert alpha_conv_dual and Allocate some temporary character arrays  */

  alpha_conv_dual        /= (pow(vol_cp,1.0/3.0));
  cppseudo->alpha_conv_dual = alpha_conv_dual;

  filename_parse->vps_name = (NAME *) 
    cmalloc(natm_typ*sizeof(NAME),"control_vps_params")-1;
  vps_file  = (VPS_FILE *) 
    cmalloc(natm_typ*sizeof(VPS_FILE),"control_vps_params")-1;
  fun_key   = (char *)cmalloc(MAXWORD*sizeof(char),"control_vps_params");  
  filename  = (char *)cmalloc(MAXWORD*sizeof(char),"control_vps_params");  
  word      = (DICT_WORD *)cmalloc(sizeof(DICT_WORD),"control_vps_params")-1;
  cvps_typ  = (CVPS *)cmalloc(sizeof(CVPS),"control_vps_params");  
  ifirst    = 1;
  set_potfun_dict(&fun_dict,&num_fun_dict,ifirst);
  set_potvps_dict(&vps_dict,&num_vps_dict,ifirst);      
  set_potvps_dict(&vps_dict_tmp,&num_vps_dict,ifirst);      


  /*==========================================================================*/
  /* III) Malloc up the vps stuff of isze natm_typ                            */ 

  natm_typ_mall      = natm_typ;

  if( (natm_typ_mall % 2)==0){natm_typ_mall++;}

  cppseudo->n_ang      = (int *) 
    cmalloc(natm_typ_mall*sizeof(int),"control_vps_params")-1;
  cppseudo->loc_opt    = (int *)
    cmalloc(natm_typ_mall*sizeof(int),"control_vps_params")-1;
  cppseudo->ivps_label = (int *)
    cmalloc(natm_typ_mall*sizeof(int),"control_vps_params")-1;
  cppseudo->rcut_nl    = (double *)
    cmalloc(natm_typ_mall*sizeof(double),"control_vps_params")-1;
  cppseudo->q_pseud    = (double *)
    cmalloc(natm_typ_mall*sizeof(double),"control_vps_params")-1;

  cppseudo->nrad_0 = (int *)
    cmalloc(natm_typ_mall*sizeof(double),"control_vps_params")-1;
  cppseudo->nrad_1 = (int *)
    cmalloc(natm_typ_mall*sizeof(double),"control_vps_params")-1;
  cppseudo->nrad_2 = (int *)
    cmalloc(natm_typ_mall*sizeof(double),"control_vps_params")-1;
  cppseudo->nrad_3 = (int *)
    cmalloc(natm_typ_mall*sizeof(double),"control_vps_params")-1;

  now_mem = (natm_typ_mall*(sizeof(double)*2 + sizeof(int))*7)*1.e-06;
  *tot_memory += now_mem;

  /*==========================================================================*/
  /* IV) Initial malloc and zeroing of the nonlocal list                      */ 
  /*     List will be reallocated when the exact number of nonlocal atoms     */
  /*     is determined                                                        */

  cppseudo->iatm_nonloc = (int *) 
    cmalloc(natm_ab_init*sizeof(int),"control_vps_params")-1;
  for(int i=1;i<=natm_ab_init;i++){
    cppseudo->iatm_nonloc[i] = 0;
  }
  natm_nonloc_now = 0;

  /*==========================================================================*/
  /* V) Loop over all unique atom types and get the vps stuff               */ 

  cppseudo->n_rad_max    = 0;
  cppseudo->n_ang_max    = 0;
  cppseudo->n_ang_max_kb = 0;
  cppseudo->n_ang_max_gh = 0;
  ngh_max = 0;
  for(i=1;i<=natm_typ;i++) {
    /*--------------------------------------------------------------------------*/
    /*     A) First search the user defined data base                           */

    ifound = 0;
    strcpy(cvps_typ->atm1,atm_typ[i]);
    if(strcasecmp(filename_parse->user_vps_name,"")!=0) {
      search_base_vps(filename_parse->user_vps_name,
          cvps_typ,fun_dict,num_fun_dict,
          &vps_dict_tmp,vps_dict,num_vps_dict,&ifound);
      if(ifound==1){
        set_vps_params(vps_dict,
            filename_parse->user_vps_name,fun_key,
            &(cppseudo->ivps_label[i]),filename,
            &(cppseudo->loc_opt[i]),&(cppseudo->n_ang[i]),
            &(cppseudo->rcut_nl[i]),&ngh_now,
            &(cppseudo->nrad_0[i]),
            &(cppseudo->nrad_1[i]),&(cppseudo->nrad_2[i]),
            &(cppseudo->nrad_3[i]));
      }/*endif*/
    }/*endif*/
    /*--------------------------------------------------------------------------*/
    /*     B) If you haven't found it search the default data base              */

    if(ifound == 0) {
      search_base_vps(filename_parse->def_vps_name,
          cvps_typ,fun_dict,num_fun_dict,
          &vps_dict_tmp,vps_dict,num_vps_dict,&ifound);
      if(ifound==1){
        set_vps_params(vps_dict,
            filename_parse->def_vps_name,fun_key,
            &(cppseudo->ivps_label[i]),filename,
            &(cppseudo->loc_opt[i]),&(cppseudo->n_ang[i]),
            &(cppseudo->rcut_nl[i]),&ngh_now,
            &(cppseudo->nrad_0[i]),
            &(cppseudo->nrad_1[i]),&(cppseudo->nrad_2[i]),
            &(cppseudo->nrad_3[i]));
      }/*endif*/
    }/*endif*/
    /*--------------------------------------------------------------------------*/
    /*     C) Make sure you have now found this puppy, if not exit              */

    if(ifound == 0) {
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("Electron pseudopotential interaction with\n"); 
      PRINTF("%s\n",atm_typ[i]);
      PRINTF("not found in default interaction data base\n");
      PRINTF("pi_md.vps\n");
      if(strlen(filename_parse->user_vps_name) > 0)  {
        PRINTF("or in user defined pseudopot data base\n");
        PRINTF("%s\n",filename_parse->user_vps_name);
        /*endif*/}
      PRINTF("\n");
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }/*endif*/

    /*--------------------------------------------------------------------------*/
    /*     D) Find maximum angular momentum component                           */
    cppseudo->n_ang_max = (cppseudo->n_ang_max > cppseudo->n_ang[i] ? 
        cppseudo->n_ang_max : cppseudo->n_ang[i]);

    if(cppseudo->ivps_label[i] != 2){
      cppseudo->n_ang_max_kb = (cppseudo->n_ang_max_kb > cppseudo->n_ang[i] ? 
          cppseudo->n_ang_max_kb : cppseudo->n_ang[i]);
    }else{
      cppseudo->n_ang_max_gh = (cppseudo->n_ang_max_gh > cppseudo->n_ang[i] ? 
          cppseudo->n_ang_max_gh : cppseudo->n_ang[i]);
      natm_typ_gh++;
    }/*endif*/

    cppseudo->n_rad_max = (cppseudo->n_rad_max > cppseudo->nrad_0[i] ? 
        cppseudo->n_rad_max : cppseudo->nrad_0[i]);
    cppseudo->n_rad_max = (cppseudo->n_rad_max > cppseudo->nrad_1[i] ?
        cppseudo->n_rad_max : cppseudo->nrad_1[i]);
    cppseudo->n_rad_max = (cppseudo->n_rad_max > cppseudo->nrad_2[i] ?
        cppseudo->n_rad_max : cppseudo->nrad_2[i]);
    cppseudo->n_rad_max = (cppseudo->n_rad_max > cppseudo->nrad_3[i] ?
        cppseudo->n_rad_max : cppseudo->nrad_3[i]);
    ngh_max             = MAX(ngh_max,ngh_now);

    strcpy(vps_file[i].name,filename);
    strcpy(filename_parse->vps_name[i],filename);
  }/*endfor natm_typ*/

  /*==========================================================================*/
  /*  VI) Allocate more memory for pseudopotentials                          */

  /*--------------------------------------------------------------------------*/
  /*    0) Store malloc sizes in class  */

  nsplin_g     = cppseudo->nsplin_g;
  n_rad_max    = cppseudo->n_rad_max;
  n_ang_max1   = cppseudo->n_ang_max+1;

  norm_size    = n_ang_max1*natm_typ*n_rad_max*n_rad_max;
  nsplin_g_tot = n_ang_max1*n_rad_max*nsplin_g*natm_typ;
  nlist        = n_ang_max1*natm_tot;
  ngh_tot      = ngh_max*natm_typ_gh*n_ang_max1;


  cppseudo->n_ang_max1   = n_ang_max1;
  cppseudo->natm_typ     = natm_typ;
  cppseudo->natm_tot     = natm_tot;
  cppseudo->nsplin_g_tot = nsplin_g_tot;
  cppseudo->norm_size    = norm_size;
  cppseudo->nlist        = nlist;
  cppseudo->natm_typ_gh  = natm_typ_gh;
  cppseudo->ngh          = ngh_max;

  nsplin_mall  = nsplin_g_tot;
  norm_mall    = norm_size;
  natm_mall    = natm_tot;
  nlist_mall   = nlist;
  if((nsplin_mall % 2)==0){nsplin_mall++;}
  if((norm_mall   % 2)==0){norm_mall++;}
  if((natm_mall   % 2)==0){natm_mall++;}
  if((nlist_mall  % 2)==0){nlist_mall++;}

  /*--------------------------------------------------------------------------*/
  /* i) Malloc Pseudo spline and other stuff */

  cppseudo->vps0    = 
    (double *) cmalloc(nsplin_mall*sizeof(double),"control_vps_params")-1;
  cppseudo->vps1    =
    (double *) cmalloc(nsplin_mall*sizeof(double),"control_vps_params")-1;
  cppseudo->vps2    =
    (double *) cmalloc(nsplin_mall*sizeof(double),"control_vps_params")-1;
  cppseudo->vps3    =
    (double *) cmalloc(nsplin_mall*sizeof(double),"control_vps_params")-1;

  now_mem   = ( nsplin_mall*4 *sizeof(double))*1.e-06;
  *tot_memory += now_mem;

  if(cp_ptens_calc == 1){
    cppseudo->dvps0    =
      (double *) cmalloc(nsplin_mall*sizeof(double),"control_vps_params")-1;
    cppseudo->dvps1    =
      (double *) cmalloc(nsplin_mall*sizeof(double),"control_vps_params")-1;
    cppseudo->dvps2    =
      (double *) cmalloc(nsplin_mall*sizeof(double),"control_vps_params")-1;
    cppseudo->dvps3    =
      (double *) cmalloc(nsplin_mall*sizeof(double),"control_vps_params")-1;
    now_mem   = ( nsplin_mall*4 *sizeof(double))*1.e-06;
    *tot_memory += now_mem;
  }/* endif */
  cppseudo->vpsnorm =
    (double *) cmalloc(norm_mall*sizeof(double),"control_vps_params")-1;
  cppseudo->gzvps   =
    (double *) cmalloc(natm_typ_mall*sizeof(double),"control_vps_params")-1;
  cppseudo->gzvps0  = 
    (double *) cmalloc(natm_typ_mall*cppseudo->n_rad_max*
        sizeof(double),"control_vps_params")-1;
  cppseudo->nrad_max_l =
    (int *)cmalloc(natm_typ_mall*5*sizeof(int),"control_vps_params")-1;

  now_mem = ((norm_mall+natm_typ_mall
        +natm_typ_mall*cppseudo->n_rad_max)*sizeof(double)
      +(natm_typ_mall*5)*sizeof(int))*1.e-06;

  *tot_memory += now_mem;

  /*--------------------------------------------------------------------------*/
  /* iii) Pseudo list : MAJOR HACKET JOB */

  natm_mall = natm_tot;
  if((natm_mall % 2)==0){natm_mall++;}
  nlist_mall = (cppseudo->n_ang_max1)*natm_tot;
  if((nlist_mall % 2)==0){nlist_mall++;}

  cppseudo->np_nl        = (int *)
    cmalloc((cppseudo->n_ang_max1)*sizeof(int),"control_vps_params")-1;
  cppseudo->np_nl_gh     = (int *)
    cmalloc((cppseudo->n_ang_max1)*sizeof(int),"control_vps_params")-1;
  cppseudo->ip_nl        = (int *)
    cmalloc(nlist_mall*sizeof(int),"control_vps_params")-1;
  cppseudo->ip_nl_gh     = (int *)
    cmalloc(nlist_mall*sizeof(int),"control_vps_params")-1;
  cppseudo->ip_nl_rev    = (int *)
    cmalloc(nlist_mall*sizeof(int),"control_vps_params")-1;
  cppseudo->ip_nl_rev_gh = (int *)
    cmalloc(nlist_mall*sizeof(int),"control_vps_params")-1;

  cppseudo->map_nl = (int *) 
    cmalloc(natm_mall*sizeof(int),"control_vps_params")-1;
  cppseudo->ip_loc_cp_box = (int *)
    cmalloc(natm_mall *sizeof(int),"control_vps_params")-1;

  cppseudo->np_nl_rad_str  = 
    cmall_int_mat(1,cppseudo->n_ang_max1,
        1,cppseudo->n_rad_max,"control_vps_params.h");
  cppseudo->np_nl_rad_end  = 
    cmall_int_mat(1,cppseudo->n_ang_max1,
        1,cppseudo->n_rad_max,"control_vps_params.h");

  cppseudo->rgh = (double *)
    cmalloc(ngh_max*sizeof(double),"control_vps_params")-1;

  nmall_gh =  ngh_max*(cppseudo->natm_typ_gh)*(cppseudo->n_ang_max1);

  cppseudo->wgh = (double *)
    cmalloc((nmall_gh)*sizeof(double),"control_vps_params")-1;

  for(i=1; i<= nmall_gh; i++){
    cppseudo->wgh[i] = 0.0;
  }/*endfor*/ 

  now_mem = ((natm_mall*3 + 2*ngh_max)*sizeof(double)
      + nmall_gh*sizeof(double)
      +((cppseudo->n_ang_max1)+nlist_mall+natm_mall*2
        + 2*(cppseudo->n_ang_max1)*cppseudo->n_rad_max )
      *sizeof(int))*1.e-06;

  *tot_memory += now_mem;

  /*--------------------------------------------------------------------------*/
  /* iv) Output */

  PRINTF("Pseudopotential allocation: %g Mbytes; Total memory: %g Mbytes\n",
      now_mem,*tot_memory);

  /*==========================================================================*/
  /*  VII) Spline up the stuff                                                 */


  for(i=1;i<=natm_typ;i++){
    cppseudo->nrad_max_l[1]= MAX(cppseudo->nrad_max_l[1],cppseudo->nrad_0[i]);
    cppseudo->nrad_max_l[2]= MAX(cppseudo->nrad_max_l[2],cppseudo->nrad_1[i]);
    cppseudo->nrad_max_l[3]= MAX(cppseudo->nrad_max_l[3],cppseudo->nrad_2[i]);
    cppseudo->nrad_max_l[4]= MAX(cppseudo->nrad_max_l[4],cppseudo->nrad_3[i]);
  }/*endfor*/

  cppseudo->dg_spl = ((cppseudo->gmax_spl)-(cppseudo->gmin_spl))
    /((double)(cppseudo->nsplin_g));

  for(i=1;i<=natm_typ*(cppseudo->n_rad_max);i++){
    (cppseudo->gzvps0)[i] = 0;
  }/*endfor*/

  for(i=1;i<=natm_typ;i++) {
    ishift  = (i-1)*(cppseudo->n_ang_max1)*(cppseudo->nsplin_g)
      *(cppseudo->n_rad_max);
    ishift2 = (i-1)*(cppseudo->n_ang_max1)*(cppseudo->n_rad_max)
      *(cppseudo->n_rad_max);
    ishift3 = (i-1)*(cppseudo->n_rad_max);
    strcpy(filename,vps_file[i].name);
    if(cp_ptens_calc == 1){
      make_vps_splin(filename,cppseudo->loc_opt[i],cppseudo->n_ang[i],
          cppseudo->ivps_label[i],
          cppseudo->nsplin_g,cppseudo->dg_spl,
          cppseudo->gmin_spl,
          cppseudo->gmax_spl,cppseudo->gmin_true,
          &(cppseudo->vps0)[ishift],&(cppseudo->vps1)[ishift],
          &(cppseudo->vps2)[ishift],&(cppseudo->vps3)[ishift],
          &(cppseudo->dvps0)[ishift],&(cppseudo->dvps1)[ishift],
          &(cppseudo->dvps2)[ishift],&(cppseudo->dvps3)[ishift],
          &(cppseudo->gzvps)[i],&(cppseudo->gzvps0)[ishift3],
          &(cppseudo->q_pseud[i]),
          &(cppseudo->vpsnorm)[ishift2],
          (cppseudo->nrad_0[i]),
          (cppseudo->nrad_1[i]),(cppseudo->nrad_2[i]),
          (cppseudo->nrad_3[i]),
          cp_ptens_calc,
          cp_dual_grid_opt,alpha_conv_dual,cppseudo->n_rad_max,
          cppseudo->n_ang_max_gh,
          &(cppseudo->ngh),cppseudo->rgh,cppseudo->wgh);
    } else {
      make_vps_splin(filename,cppseudo->loc_opt[i],cppseudo->n_ang[i],
          cppseudo->ivps_label[i],
          cppseudo->nsplin_g,cppseudo->dg_spl,
          cppseudo->gmin_spl,
          cppseudo->gmax_spl,cppseudo->gmin_true,
          &(cppseudo->vps0)[ishift],&(cppseudo->vps1)[ishift],
          &(cppseudo->vps2)[ishift],&(cppseudo->vps3)[ishift],
          &dummy0,&dummy1,
          &dummy2,&dummy3,
          &(cppseudo->gzvps)[i],&(cppseudo->gzvps0)[ishift3],
          &(cppseudo->q_pseud[i]),
          &(cppseudo->vpsnorm)[ishift2],
          (cppseudo->nrad_0[i]),
          (cppseudo->nrad_1[i]),(cppseudo->nrad_2[i]),
          (cppseudo->nrad_3[i]),
          cp_ptens_calc,
          cp_dual_grid_opt,alpha_conv_dual,cppseudo->n_rad_max,
          cppseudo->n_ang_max_gh,
          &(cppseudo->ngh),cppseudo->rgh,cppseudo->wgh);
    }/* endif */
  }/*endfor*/


  PRINTF("Electron-atom interactions succesfully assigned\n");
  if( (cppseudo->n_rad_max>1) || (cppseudo->natm_typ_gh > 0)){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("One radial channel and no Gauss Hermite puppies.");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    EXIT(1);
  }//endif

  /*==========================================================================*/
  /*  IX) Free                                                 */

  cfree(&vps_file[1],"control_vps_params");
  cfree(fun_key,"control_vps_params");
  cfree(filename,"control_vps_params");
  cfree(&word[1],"control_vps_params");
  cfree(&fun_dict[1],"control_vps_params");
  cfree(&vps_dict[1],"control_vps_params");
  cfree(&vps_dict_tmp[1],"control_vps_params");
  cfree(cvps_typ,"control_vps_params");

  /*==========================================================================*/
  /* X) Create the non-local list */

  create_non_local_list(cppseudo,natm_tot,iatm_atm_typ,natm_typ);

  double *q_pseud=cppseudo->q_pseud;
  /*==========================================================================*/
  /* Test the q_pseud */

  for(i=1;i<=natm_typ;i++){
    if(fabs(q_typ[i]-q_pseud[i])>1.e-05){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("Charges assigned to the atom types must\n");
      PRINTF("match the pseudo potential charge : %d : %g %g!\n",
          i,q_typ[i],q_pseud[i]);
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }//endif
  }//endfor

  /*==========================================================================*/
  /*  XII) Output                                                             */

  PRINTF("\n");PRINT_LINE_DASH
    PRINTF("Completed the pseudopotential data bases searches\n");
  PRINT_LINE_STAR;PRINTF("\n");

  /*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/






/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_vps_params(DICT_WORD vps_dict[],char *filename,char *fun_key,
    int *ivps_label,char *vps_file,
    int *loc_opt,int *n_ang,double *rcut_nl,int *pngh,
    int *nrad_0, int *nrad_1, int *nrad_2, int *nrad_3)

  /*==========================================================================*/
  /*               Begin subprogram:                                          */
{/*begin routine*/
  int index,ierr,ngh;
  double tmp;
  /*--------------------------------------------------------------------------*/
  /*  0) Set up                                                               */
  strcpy(fun_key,"vps_parm");
  /*--------------------------------------------------------------------------*/
  /*  1) Assign vps file                                                      */
  strcpy(vps_file,vps_dict[3].keyarg);
  /*--------------------------------------------------------------------------*/
  /*  2) Assign radial channel info                                           */

  sscanf(vps_dict[7].keyarg,"%lg",&tmp);
  (*nrad_0)  = (int )(tmp);
  sscanf(vps_dict[8].keyarg,"%lg",&tmp);
  (*nrad_1)  = (int )(tmp);
  sscanf(vps_dict[9].keyarg,"%lg",&tmp);
  (*nrad_2)  = (int )(tmp);
  sscanf(vps_dict[10].keyarg,"%lg",&tmp);
  (*nrad_3)  = (int )(tmp);

  if((*nrad_0)<0){index=7;
    keyarg_barf(vps_dict,filename,fun_key,index);}
  if((*nrad_1)<0){index=7;
    keyarg_barf(vps_dict,filename,fun_key,index);}
  if((*nrad_2)<0){index=7;
    keyarg_barf(vps_dict,filename,fun_key,index);}
  if((*nrad_3)<0){index=7;
    keyarg_barf(vps_dict,filename,fun_key,index);}

  /*--------------------------------------------------------------------------*/
  /*  4) Assign vps type                                                      */

  (*ivps_label) = -1;
  if(strcasecmp(vps_dict[2].keyarg,"loc") == 0)            {(*ivps_label) = 0;}
  if(strcasecmp(vps_dict[2].keyarg,"local") == 0)          {(*ivps_label) = 0;}
  if(strcasecmp(vps_dict[2].keyarg,"kb")  == 0)            {(*ivps_label) = 1;}
  if(strcasecmp(vps_dict[2].keyarg,"gauss_hermite")  == 0) {(*ivps_label) = 2;}
  if(strcasecmp(vps_dict[2].keyarg,"vdb") == 0)            {(*ivps_label) = 3;}
  if(strcasecmp(vps_dict[2].keyarg,"null")== 0)            {(*ivps_label) = 4;}
  if(strcasecmp(vps_dict[2].keyarg,"goedecker")== 0)       {(*ivps_label) = 5;}
  if((*ivps_label)==-1){index=2;
    keyarg_barf(vps_dict,filename,fun_key,index);}

  /*--------------------------------------------------------------------------*/
  /*  2) Assign non-local info                                                */

  sscanf(vps_dict[4].keyarg,"%lg",&tmp);
  (*n_ang)  = (int )(tmp);
  if((*n_ang)<0){index=4;
    keyarg_barf(vps_dict,filename,fun_key,index);}
  if((*n_ang)>4){index=4;
    keyarg_barf(vps_dict,filename,fun_key,index);}

  sscanf(vps_dict[6].keyarg,"%lg",&tmp);
  (*rcut_nl) = tmp;
  if((*rcut_nl)<0){index=6;
    keyarg_barf(vps_dict,filename,fun_key,index);}

  /*--------------------------------------------------------------------------*/
  /*  4) Local option                                                        */

  sscanf(vps_dict[5].keyarg,"%lg",&tmp);
  (*loc_opt) = (int )(tmp);

  if((*ivps_label)!=5){
    if(((*loc_opt)<0)||((*loc_opt)>(*n_ang))){index=5;
      keyarg_barf(vps_dict,filename,fun_key,index);}
  }/*endif*/

  if(*ivps_label == 5){
    if( (*loc_opt)!=(*n_ang+1) ){
      if( ((*loc_opt)==0) && ((*nrad_0)==0) ){
        index=1;
      }else{
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("Dude, the local option must be set to \n"); 
        PRINTF("the number of available angular momentum channels+1\n");
        PRINTF("%d vs %d \n",(*loc_opt),(*n_ang+1));
        PRINTF("for a Goedecker-type pseudopotential\n");
        PRINTF("unless the lmax=0 and n_rad_0=0,\n");
        PRINTF("that is, unless Stephan has made a local pseudo\n");
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        FFLUSH(stdout);
        EXIT(1);
      }/*endif*/
    }/*endif*/
    ierr = 0;
    if( ((*loc_opt)==1)&&((*nrad_0)==0) ){ierr++;}
    if( ((*loc_opt)==2)&&((*nrad_1)==0) ){ierr++;}
    if( ((*loc_opt)==3)&&((*nrad_2)==0) ){ierr++;}
    if( ((*loc_opt)==4)&&((*nrad_3)==0) ){ierr++;} 
    if(ierr>0){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("You have set up the last angular channel with \n");
      PRINTF("no radial channels for a Goedecker-type pseudopotential\n");
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }/*endif*/
    if( (*loc_opt)!=0 ){*n_ang = (*n_ang + 1);}
  }/*endif*/

  /*--------------------------------------------------------------------------*/
  /*  5) Gauss hermite pseudopotentials                                      */

  index = 11;
  sscanf(vps_dict[11].keyarg,"%lg",&tmp);
  ngh = (int) tmp;
  if(*ivps_label == 2 && ngh > 180) {
    keyarg_barf(vps_dict,filename,fun_key,index);
  }
  if(*ivps_label == 2 && ngh == 0){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("You have set a gauss-hermite nonlocal atom\n");
    PRINTF("but have zero integration points: ngh =  %d.\n",ngh);
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }/*endif*/

  *pngh = ngh;
  /*--------------------------------------------------------------------------*/
  /* Set Radial channels if not Goedecker */

  if(*ivps_label != 5){
    *nrad_0=1;
    *nrad_1=1;
    *nrad_2=1;
    *nrad_3=1;
  }/*endif*/

  /*-----------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/







/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void make_vps_splin(char *vps_file,int loc_opt,int n_ang,
    int ivps_label,int nsplin_g,double dg_spl,
    double gmin_spl,double gmax_spl,double gmin_true,
    double *vps0,double *vps1,double *vps2,double *vps3,
    double *dvps0,double *dvps1,double *dvps2,double *dvps3,
    double *gzvps,double *gzvps0,double *q_pseud,
    double *vpsnorm, int nrad_0,int nrad_1,int nrad_2,
    int nrad_3,int cp_ptens_calc,int cp_dual_grid_opt,
    double alpha_conv_dual,
    int n_rad_max,int n_ang_max_gh,
    int *pngh,double *rgh,double *wgh)

  /*==========================================================================*/
  /*               Begin subprogram:                                          */
{/*begin routine*/
  /*==========================================================================*/
  /*               Local variable declarations                                */

  double rmax;                   /* Num: Maximum spline radius       */
  double z_1,z_2;                /* Num: Charges for long range piece */
  double alpha_1,alpha_2;        /* Num: Ewald alpha's for long range */ 
  double dr;                     /* Num: delta r for spline           */
  double amat;                   /* Num: Useful temporary for spline  */
  double v_now;                  /* Num: current value of pp          */
  double rphi_now;               /* Num: current value of pseudo wf   */
  double zpol;                   /* Num: Polarization charge          */
  double gamma;                  /* Num: Useful temporary for long range */
  double dummy;                  /* Num: Really dum for pressure tensor*/
  int i;                         /* Num: Generic counter             */
  int ilong;                     /* Num: Flag for long range piece   */
  int n_ang_now;                 /* Num: Number of angular momentum  
                                    components                  */
  int nchan_tot;                 /* Num: Total number of radial channels*/
  int n_rad0_now,n_rad1_now;     /* Num: Number of radial channels for */
  int n_rad2_now,n_rad3_now;     /*   for angular momentum channel     */
  int iii;
  int iang;                      /* Num: Angular momentum counter    */
  int ishift_now;                /* Num: Angular momentum shift      */
  int iang_now;                  /* Num: Angular momentum counter    */
  int irad,jrad;                 /* Num: Radial channel counter      */
  int ioff_rad;                  /* Num: Radial channel shift        */
  int *nrad_l_tmp;               /* Temp array  # radial channels for
                                    given angular momentum channel */
  double *g;                     /* Lst: temporary spline array      */
  double *v_rphi,*v_loc,*r;      /* Lst: temporary spline arrays     */
  int ir;                        /* Num: counter for above arrays    */
  int nr;                        /* Num: Length of above arrays      */ 
  FILE *fp_vps_file;             /* File pointer to pseudo pot file  */

  int ierr;
  double tmp;
  int n_rad_max_sq = n_rad_max*n_rad_max;
  int ngh = *pngh;

  /* Temp gauss-hermite spline arrays */
  static int natm_typ_gh = 0;
  double *wgh_tmp;
  double *c0,*c1,*c2,*c3;
  double *dvl;

  /*==========================================================================*/
  /*  I) Allocate g-space vector for spline                                   */

  g = (double *) cmalloc((nsplin_g)*sizeof(double),"make_vps_splin")-1;

  /*==========================================================================*/
  /*  II) Open the electron-atom pseudopotential file for reading             */

  fp_vps_file = cfopen((const char *) vps_file,"r");

  /*==========================================================================*/
  /*  III) Set up g vectors                                                   */

  for(i=1;i <= nsplin_g;i++){ 
    g[i] = dg_spl*((double) (i-1)) + gmin_spl;
  }/*endfor*/

  /*==========================================================================*/
  /*  IV) Set up a local pseudo potential                                     */

  if(ivps_label <= 4) {

    if(fscanf(fp_vps_file,"%d %lf %d\n",&nr,&rmax,&n_ang_now) != 3) 
    {vps_read_error(vps_file);}
    if(fscanf(fp_vps_file,"%lf %lf %lf %lf\n",
          &z_1,&alpha_1,&z_2,&alpha_2) != 4) 
    {vps_read_error(vps_file);}
    if(fscanf(fp_vps_file,"%lf %lf\n",&zpol,&gamma) != 2)      
    {vps_read_error(vps_file);}


  }/*endif ivps_label */
  /*--------------------------------------*/
  /* ivps_label == 5  GOEDECKER Potential */

  if(ivps_label == 5){
    if(fscanf(fp_vps_file,"%d %lf %d\n",&nr,&rmax,&n_ang_now) != 3) 
    {vps_read_error(vps_file);}

    if(fscanf(fp_vps_file,"%d %d %d %d\n",&n_rad0_now,
          &n_rad1_now,
          &n_rad2_now, 
          &n_rad3_now) != 4) 
    {vps_read_error(vps_file);}
    if(fscanf(fp_vps_file,"%lf %lf %lf %lf\n",
          &z_1,&alpha_1,&z_2,&alpha_2) != 4) 
    {vps_read_error(vps_file);}
    if(fscanf(fp_vps_file,"%lf %lf\n",&zpol,&gamma) != 2) 
    {vps_read_error(vps_file);}

    /*-------------------------------------*/
    /* Read in Radial Channel Matrices     */
    /*-------------------------------------*/

    /*-----------  S CHANNEL --------------*/
    ioff_rad = 0;
    for(irad=1; irad <= n_rad0_now; irad++){
      for(jrad=1; jrad <= n_rad0_now; jrad++){
        if(fscanf(fp_vps_file,"%lf ",&(vpsnorm[jrad+ioff_rad])))
        {vps_read_error(vps_file);}
      }/*endfor*/
      ioff_rad += n_rad_max;
    }/*endfor*/

    /*-----------  P CHANNEL --------------*/
    ioff_rad = n_rad_max*n_rad_max;
    for(irad=1; irad <= n_rad1_now; irad++){
      for(jrad=1; jrad <= n_rad1_now; jrad++){
        if(fscanf(fp_vps_file,"%lf ",&(vpsnorm[jrad+ioff_rad])))
        {vps_read_error(vps_file);}
      }/*endfor*/
      ioff_rad += n_rad_max;
    }/*endfor*/

    /*-----------  D CHANNEL --------------*/
    ioff_rad = 2*n_rad_max*n_rad_max;
    for(irad=1; irad <= n_rad2_now; irad++){
      for(jrad=1; jrad <= n_rad2_now; jrad++){
        if(fscanf(fp_vps_file,"%lf ",&(vpsnorm[jrad+ioff_rad])))
        {vps_read_error(vps_file);}
      }/*endfor*/
      ioff_rad += n_rad_max;
    }/*endfor*/

    /*-----------  F CHANNEL --------------*/
    ioff_rad = 3*n_rad_max*n_rad_max;
    for(irad=1; irad <= n_rad3_now; irad++){
      for(jrad=1; jrad <= n_rad3_now; jrad++){
        if(fscanf(fp_vps_file,"%lf ",&(vpsnorm[jrad+ioff_rad])))
        {vps_read_error(vps_file);}
      }/*endfor*/
      ioff_rad += n_rad_max;
    }/*endfor*/



  }/*endif ivps_label goedecker*/

  v_rphi = (double *) cmalloc(nr*sizeof(double),"make_vps_splin")-1;
  v_loc  = (double *) cmalloc(nr*sizeof(double),"make_vps_splin")-1;
  r      = (double *) cmalloc(nr*sizeof(double),"make_vps_splin")-1;

  *q_pseud = z_1 + z_2;

  /*--------------------------------------------------------------------------*/
  /*   A) Allocate real space arrays                                          */

  dr = rmax/((double ) nr);
  for(i=1;i <= nr;i++){
    r[i] = ((double ) (i-1))*dr;
  }/*endfor*/

  if( (n_ang_now < n_ang) && (ivps_label!=5) ) {
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Dude, like given the number of angular momentum\n"); 
    PRINTF("components you've specified in the\n");
    PRINTF("pseudopotential file %s\n",vps_file);
    PRINTF("proceeding with the simulation would be a\n");
    PRINTF("pointless exercise leading to most bogus results\n");
    PRINTF("%d vs %d\n",n_ang_now,n_ang);
    PRINTF("\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }/*endif*/

  if( (n_ang_now != n_ang-1) && (ivps_label==5)) {
    if( (n_ang_now==0)&&(loc_opt==0) && (n_rad0_now==0)){
      iii = 1;
    }else{
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("Dude, the number of angular momentum\n"); 
      PRINTF("components you've specified in the\n");
      PRINTF("pseudopotential file %s\n",vps_file);
      PRINTF("does not match does the value in the controller,\n");
      PRINTF("%d vs %d\n,",n_ang_now,n_ang-1);
      PRINTF("for a Goedecker-type pseudopotential. Stephan would\n");
      PRINTF("not approve and neither do we!\n");
      PRINTF("\n");
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }/*endif*/
  }/*endif*/

  if(ivps_label==5) {   
    ierr = 0;
    if( (nrad_0 != n_rad0_now) && (n_ang-1>=0)){ierr++;}
    if( (nrad_1 != n_rad1_now) && (n_ang-1>=1)){ierr++;}
    if( (nrad_2 != n_rad2_now) && (n_ang-1>=2)){ierr++;}
    if( (nrad_3 != n_rad3_now) && (n_ang-1>=3)){ierr++;}
    if(ierr>0){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("Dude, the number of radial channel \n"); 
      PRINTF("components you've specified in the\n");
      PRINTF("pseudopotential file %s\n",vps_file);
      PRINTF("does not match does the value in the controller,\n");
      PRINTF("for a Goedecker-type pseudopotential. Stephan would\n");
      PRINTF("not approve and neither do we!\n");
      if( (nrad_0 != n_rad0_now) && (n_ang-1>0)){
        PRINTF("l=0 : %d vs %d\n",nrad_0,n_rad0_now);
      }/*endif*/
      if( (nrad_1 != n_rad1_now) && (n_ang-1>1)){
        PRINTF("l=1 : %d vs %d\n",nrad_1,n_rad1_now);
      }/*endif*/
      if( (nrad_2 != n_rad2_now) && (n_ang-1>2)){
        PRINTF("l=2 : %d vs %d\n",nrad_2,n_rad2_now);
      }/*endif*/
      if( (nrad_3 != n_rad3_now) && (n_ang-1>3)){
        PRINTF("l=3 : %d vs %d\n",nrad_3,n_rad3_now);
      }/*endif*/
      PRINTF("\n");
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }/*endif*/
  }/*endif*/


  /*--------------------------------------------------------------------------*/
  /*   B) Find v_loc                                                         */

  if(ivps_label < 5){ /* NOT GOEDECKER TYPE */


    for(iang=1;iang <= (loc_opt + 1);iang++) {
      for(ir=1;ir <= nr;ir++) {          
        if(fscanf(fp_vps_file,"%lf %lf\n",&v_now,&rphi_now) != 2) 
        {vps_read_error(vps_file);}
        v_loc[ir] = v_now;
      }/* endfor */ 
    } /* endfor */

  }else{ /*GOEDECKER TYPE */

    nchan_tot = n_rad0_now+n_rad1_now+n_rad2_now+n_rad3_now;
    for(iang=1; iang<=(nchan_tot+1); iang++){
      for(ir=1;ir <= nr;ir++) {          
        if(fscanf(fp_vps_file,"%lf %lf\n",&v_now,&rphi_now) != 2) 
        {vps_read_error(vps_file);}
        v_loc[ir] = v_now;
      }/* endfor */ 
    }/*endfor*/


  }/*endif ivps_label*/


  for(ir=1;ir <= nr;ir++) {v_rphi[ir] = v_loc[ir]*r[ir];}

  ilong      = 1;
  iang_now   = 0;
  ishift_now = (loc_opt)*n_rad_max*nsplin_g;
  if(cp_ptens_calc == 1){
    slow_bess_vps(v_rphi,nr,dr,r,&(vps0)[ishift_now],
        &(dvps0)[ishift_now],
        nsplin_g,g,gmin_true,
        z_1,alpha_1,z_2,alpha_2,zpol,gamma,ilong,iang_now,
        gzvps,&gzvps0[1],cp_ptens_calc,cp_dual_grid_opt,
        alpha_conv_dual);
  } else {
    slow_bess_vps(v_rphi,nr,dr,r,&(vps0)[ishift_now],
        &dummy,
        nsplin_g,g,gmin_true,
        z_1,alpha_1,z_2,alpha_2,zpol,gamma,ilong,iang_now,
        gzvps,&gzvps0[1],cp_ptens_calc,cp_dual_grid_opt,
        alpha_conv_dual);
  }/* endif */
  spline_fit(&(vps0)[ishift_now],&(vps1)[ishift_now],
      &(vps2)[ishift_now],&(vps3)[ishift_now],g,nsplin_g);
  if(cp_ptens_calc == 1){
    spline_fit(&(dvps0)[ishift_now],&(dvps1)[ishift_now],
        &(dvps2)[ishift_now],&(dvps3)[ishift_now],g,nsplin_g);
  }/* endif */

  /*==========================================================================*/
  /*  V) Set up a null local pseudo potential                                 */

  if(ivps_label == 4) {

    loc_opt  = 0;
    n_ang    = 0;
    *q_pseud = 0.0;
    for(i=1;i <= nsplin_g; i++) {
      vps0[i] = 0.0;
      vps1[i] = 0.0;
      vps2[i] = 0.0;
      vps3[i] = 0.0;
    }/*endfor*/
    if(cp_ptens_calc == 1){
      for(i=1;i <= nsplin_g; i++) {
        dvps0[i] = 0.0;
        dvps1[i] = 0.0;
        dvps2[i] = 0.0;
        dvps3[i] = 0.0;
      } /* endfor */
    }/* endif */

  }/* endif:null potential */

  /*==========================================================================*/
  /* VI) KB nonlocal potential                                                */

  if(ivps_label == 1) {
    /*--------------------------------------------------------------------------*/
    /*    A) rewind                                                      */

    rewind(fp_vps_file);      
    if(fscanf(fp_vps_file,"%d %lf %d\n",&nr,&rmax,&n_ang_now) != 3) 
    {vps_read_error(vps_file);}
    if(fscanf(fp_vps_file,"%lf %lf %lf %lf\n",
          &z_1,&alpha_1,&z_2,&alpha_2)!= 4) 
    {vps_read_error(vps_file);}
    if(fscanf(fp_vps_file,"%lf %lf\n",&zpol,&gamma) != 2) 
    {vps_read_error(vps_file);}
    /*--------------------------------------------------------------------------*/
    /*    B) Spline projection operator                                         */

    vpsnorm[(loc_opt*n_rad_max_sq + 1)] = 0.0;
    for(iang = 1; iang <= n_ang + 1; iang++) {
      amat = 0.0;
      for(ir=1; ir <= nr; ir++) {
        if(fscanf(fp_vps_file,"%lf %lf\n",&v_now,&rphi_now) != 2) 
        {vps_read_error(vps_file); }
        v_rphi[ir] = (v_now-v_loc[ir])*rphi_now;
        amat += rphi_now*v_rphi[ir]*dr;
      } /* endfor */
      if(iang != loc_opt+1) {
        vpsnorm[((iang-1)*n_rad_max_sq+1)] = (1.0/amat);
        ishift_now = (iang-1)*n_rad_max*nsplin_g;
        ilong = 0;
        iang_now = iang-1;
        if(cp_ptens_calc == 1){
          slow_bess_vps(v_rphi,nr,dr,r,&(vps0[ishift_now]),
              &(dvps0[ishift_now]),
              nsplin_g,g,gmin_true,
              z_1,alpha_1,z_2,alpha_2,zpol,gamma,ilong,iang_now,
              gzvps,&gzvps0[1],cp_ptens_calc,cp_dual_grid_opt,
              alpha_conv_dual);
        } else {
          slow_bess_vps(v_rphi,nr,dr,r,&(vps0[ishift_now]),
              &dummy,
              nsplin_g,g,gmin_true,
              z_1,alpha_1,z_2,alpha_2,zpol,gamma,ilong,iang_now,
              gzvps,&gzvps0[1],cp_ptens_calc,cp_dual_grid_opt,
              alpha_conv_dual);
        }/* endif */
        spline_fit(&(vps0[ishift_now]),&(vps1[ishift_now]),
            &(vps2[ishift_now]),&(vps3[ishift_now]),g,nsplin_g);

        if(cp_ptens_calc == 1){
          spline_fit(&(dvps0[ishift_now]),&(dvps1[ishift_now]),
              &(dvps2[ishift_now]),&(dvps3[ishift_now]),g,nsplin_g);
        }/* endif */

      } /* endif:loc_opt */
    } /* endfor:channels */
  } /* endif:KB nonlocal */

  /*==========================================================================*/
  /* VII) Gauss-Hermite nonlocal pseudopotentials                             */

  if(ivps_label == 2){
    natm_typ_gh++;
    wgh_tmp = (double *) cmalloc(ngh*sizeof(double),"make_vps_splin")-1;

    weight_node_gauss_hermite(ngh,rgh,wgh_tmp);

    /*Only use Gauss-Hermite integration points with a weight > 1.0e-10 */   
    /* CHECK WHETHER THIS LIMIT STATEMENT IS NEEDED */

    limit_gauss_hermite(&(ngh),rgh,wgh_tmp);  

    *pngh = ngh;


    /*--------------------------------------------------------------------------*/
    /* Read in the pseudo potential                                             */

    /*    A) rewind                                                             */

    rewind(fp_vps_file);      
    if(fscanf(fp_vps_file,"%d %lf %d\n",&nr,&rmax,&n_ang_now) != 3) 
    {vps_read_error(vps_file);}
    if(fscanf(fp_vps_file,"%lf %lf %lf %lf\n",
          &z_1,&alpha_1,&z_2,&alpha_2)!= 4) 
    {vps_read_error(vps_file);}
    if(fscanf(fp_vps_file,"%lf %lf\n",&zpol,&gamma) != 2) 
    {vps_read_error(vps_file);}
    /*--------------------------------------------------------------------------*/
    /* This only works when vloc is the highest angular momentum channel        */

    for(iang = 1; iang < (loc_opt + 1); iang++) {
      for(ir=1; ir <= nr; ir++) {
        if(fscanf(fp_vps_file,"%lf %lf\n",&v_now,&rphi_now) != 2) 
        {vps_read_error(vps_file); }
        v_rphi[ir] = (v_now-v_loc[ir]); 
      } /* endfor */


      /*--------------------------------------------------------------------------*/
      /* i) make r                                                                */
      c0  = (double *) cmalloc(nr*sizeof(double),"make_vps_splin")-1;  
      c1  = (double *) cmalloc(nr*sizeof(double),"make_vps_splin")-1;  
      c2  = (double *) cmalloc(nr*sizeof(double),"make_vps_splin")-1;  
      c3  = (double *) cmalloc(nr*sizeof(double),"make_vps_splin")-1;  
      r   = (double *) cmalloc(nr*sizeof(double),"make_vps_splin")-1;  
      dvl = (double *) cmalloc(ngh*sizeof(double),"make_vps_splin")-1;

      dr = rmax/((double) nr);

      for(i=1; i<= nr; i++){
        r[i]  = (double)(i-1)*dr;
        c0[i] = v_rphi[i];
      }/*endfor*/

      /*--------------------------------------------------------------------------*/
      /* ii) Fit to spline                                                         */

      spline_fit(c0,c1,c2,c3,r,nr);

      /*--------------------------------------------------------------------------*/
      /* ii) Fetch dvl                                                           */

      get_dvl(ngh,rgh,dvl,c0,c1,c2,c3,dr,nr);

      /*--------------------------------------------------------------------------*/
      /* iii) Make the general weight                                             */

      make_weight_gen(wgh_tmp,wgh,rgh,dvl,ngh,
          iang,natm_typ_gh,n_ang_max_gh);
    }/* endfor:channels */

    cfree(&(wgh_tmp[1]),"make_vps_splin");
    cfree(&(c0[1]),"make_vps_splin");
    cfree(&(c1[1]),"make_vps_splin");
    cfree(&(c2[1]),"make_vps_splin");
    cfree(&(c3[1]),"make_vps_splin");
    cfree(&(r[1]),"make_vps_splin");
    cfree(&(dvl[1]),"make_vps_splin");

  }/*endif GAUSS-HERMITE*/

  /*==========================================================================*/
  /* VIIII) GOEDECKER nonlocal potential                                         */

  if(ivps_label == 5) {
    /*--------------------------------------------------------------------------*/
    /*    A) rewind                                                      */
    rewind(fp_vps_file);      
    if(fscanf(fp_vps_file,"%d %lf %d\n",&nr,&rmax,&n_ang_now) != 3) 
    {vps_read_error(vps_file);}
    if(fscanf(fp_vps_file,"%d %d %d %d\n",&n_rad0_now,
          &n_rad1_now,
          &n_rad2_now,
          &n_rad3_now) != 4) 
    {vps_read_error(vps_file);}
    if(fscanf(fp_vps_file,"%lf %lf %lf %lf\n",
          &z_1,&alpha_1,&z_2,&alpha_2) != 4) 
    {vps_read_error(vps_file);}
    if(fscanf(fp_vps_file,"%lf %lf\n",&zpol,&gamma) != 2) 
    {vps_read_error(vps_file);}
    /*-----------  Skip the matrix  --------------*/
    for(irad=1; irad <= (n_rad0_now*n_rad0_now); irad++){
      if(fscanf(fp_vps_file,"%lf ",&tmp)){}
    }/*endfor*/
    for(irad=1; irad <= (n_rad1_now*n_rad1_now); irad++){
      if(fscanf(fp_vps_file,"%lf ",&tmp)){}
    }/*endfor*/
    for(irad=1; irad <= (n_rad2_now*n_rad2_now); irad++){
      if(fscanf(fp_vps_file,"%lf ",&tmp)){}
    }/*endfor*/
    for(irad=1; irad <= (n_rad3_now*n_rad3_now); irad++){
      if(fscanf(fp_vps_file,"%lf ",&tmp)){}
    }/*endfor*/

    /*--------------------------------------------------------------------------*/
    /*    B) Spline projection operator                                         */

    nrad_l_tmp     = (int *) cmalloc(4*sizeof(int),"make_vps_splin");
    nrad_l_tmp[0] = n_rad0_now;
    nrad_l_tmp[1] = n_rad1_now;
    nrad_l_tmp[2] = n_rad2_now;
    nrad_l_tmp[3] = n_rad3_now;

    /* loc_opt=n_ang so loop stops at n_ang*/
    for(iang = 1; iang <= n_ang ; iang++) {
      for(irad = 1; irad <= nrad_l_tmp[(iang-1)]; irad++){ 
        for(ir=1; ir <= nr; ir++) {
          if(fscanf(fp_vps_file,"%lf %lf\n",&rphi_now,&tmp) != 2) 
          {vps_read_error(vps_file); }
          v_rphi[ir] = rphi_now;
        } /* endfor */
        ishift_now = (iang-1)*n_rad_max*nsplin_g + (irad-1)*nsplin_g;
        ilong = 0;
        iang_now = iang-1;
        if(cp_ptens_calc == 1){
          slow_bess_vps(v_rphi,nr,dr,r,&(vps0[ishift_now]),
              &(dvps0[ishift_now]),
              nsplin_g,g,gmin_true,
              z_1,alpha_1,z_2,alpha_2,zpol,gamma,ilong,iang_now,
              gzvps,&gzvps0[irad],cp_ptens_calc,cp_dual_grid_opt,
              alpha_conv_dual);
        } else {
          slow_bess_vps(v_rphi,nr,dr,r,&(vps0[ishift_now]),
              &dummy,
              nsplin_g,g,gmin_true,
              z_1,alpha_1,z_2,alpha_2,zpol,gamma,ilong,iang_now,
              gzvps,&gzvps0[irad],cp_ptens_calc,cp_dual_grid_opt,
              alpha_conv_dual);
        }/* endif */
        spline_fit(&(vps0[ishift_now]),&(vps1[ishift_now]),
            &(vps2[ishift_now]),&(vps3[ishift_now]),g,nsplin_g);

        if(cp_ptens_calc == 1){
          spline_fit(&(dvps0[ishift_now]),&(dvps1[ishift_now]),
              &(dvps2[ishift_now]),&(dvps3[ishift_now]),g,nsplin_g);
        }/* endif */

      }/*endfor: radial channels */
    }/* endfor:angular channels */

    cfree(nrad_l_tmp,"make_vps_splin");

  } /* endif:GOEDECKER nonlocal */

  /*==========================================================================*/
  /*   VII) Free memory                                                       */

  if(ivps_label != 4) {
    cfree(&(v_rphi[1]),"make_vps_splin");
    cfree(&(v_loc[1]),"make_vps_splin");
    if(ivps_label != 2) cfree(&(r[1]),"make_vps_splin");
  }/*endif*/
  cfree(&(g[1]),"make_vps_splin");


  /*==========================================================================*/
  /*   VIII) Close the file and done                                          */

  fclose(fp_vps_file);

  /*-----------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/







/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void slow_bess_vps(double v_rphi[],int nr,double dr,double r[],
    double fv_rphi[],double fdv_rphi[],int nsplin_g,double g[],
    double gmin_true,
    double z_1,double alpha_1,double z_2,double alpha_2,
    double zpol,double gamma,int ilong,int iang,
    double *gzvps,double *gzvps0,int cp_ptens_calc,
    int cp_dual_grid_opt, double alpha_conv_dual)

  /*==========================================================================*/
  /*               Begin subprogram:                                          */
{/*begin routine*/
  /*==========================================================================*/
  /*               Local variable declarations                                */

  double c[8],angamma[8],angamma2[8];  /* Lst: Useful temporary arrays */
  double phi;                    /* Num: Angle used in 
                                    polarization corrections*/
  double g2,rj0,rj1,rj2,rj3,gzero,arg,fpi,pi,tpi,fpidr;
  /* Num: Useful constants and
     temporaries             */
  double r2dj0,r2dj1,r2dj2,r2dj3,j4,j3,j2,j1,j0;
  int ig,ir,i;              /* Num: Counters                */
  int iii;
  double falpha_12,falpha_22;

  double gm1,gamm_plu,gamm_min,pre_c;
  double ztot,falpha_conv_dual;


  /*==========================================================================*/
  /* II) Get some constants                                                    */

  pi  = M_PI;
  tpi = 2.0*pi;
  fpi = 4.0*pi;
  fpidr = fpi*dr;
  for(ig=1;ig <= nsplin_g;ig++) {
    fv_rphi[ig] = 0.0;
  }
  if(cp_ptens_calc ==1 ){
    for(ig=1;ig <= nsplin_g;ig++) {
      fdv_rphi[ig] = 0.0;
    }/* endfor */
  }/* endif */

  /*==========================================================================*/
  /*==========================================================================*/

  switch(iang) {

    /*==========================================================================*/
    /* I) L=0 Term */

    case 0:

      /*--------------------------------------------------------------------------*/
      /*    i) g = 0 : fdv = 0                                                    */

      gzero = 0.0;
      for(ir=2;ir <= nr;ir++) {gzero += fpidr*r[ir]*v_rphi[ir];}
      if(ilong == 1) {*gzvps  = gzero;}
      if(ilong != 1) {*gzvps0 = gzero;}

      /*--------------------------------------------------------------------------*/
      /*    ii) g ne 0                                                            */

      for(ig=1;ig <= nsplin_g; ig++) {
        for(ir=2;ir <= nr; ir++) {
          arg = r[ir]*g[ig];
          rj0 = (sin(arg)/arg)*r[ir];
          fv_rphi[ig]  += fpidr*rj0*v_rphi[ir];
        } /* endfor */
        if(cp_ptens_calc ==1 ){
          for(ir=2;ir <= nr; ir++) {
            arg = r[ir]*g[ig];
            j1 = (sin(arg)/arg - cos(arg))/arg;
            r2dj0 = -(j1)*r[ir]*r[ir];
            fdv_rphi[ig] += fpidr*r2dj0*v_rphi[ir]/g[ig];
          } /* endfor */  
        }/* endif */
      } /* endfor */


      break;

      /*==========================================================================*/
      /* II) L=1 Term */

    case 1:

      /*--------------------------------------------------------------------------*/
      /*    i) g ne 0                                                             */

      for(ig=1;ig <= nsplin_g; ig++) {
        for(ir=2;ir <= nr; ir++) {
          arg = r[ir]*g[ig];
          rj1 = ((sin(arg)/arg - cos(arg))/arg)*r[ir];
          fv_rphi[ig] += fpidr*rj1*v_rphi[ir];
        } /* endfor */
        if(cp_ptens_calc ==1 ){
          for(ir=2;ir <= nr; ir++) {
            arg = r[ir]*g[ig];
            j0 = 1.0*(sin(arg)/arg);
            j2 = ((3.0/(arg*arg)-1.0)*sin(arg)-3.0*cos(arg)/arg)/arg;
            r2dj1 = (j0-2.0*j2 )*r[ir]*r[ir]/3.0;
            fdv_rphi[ig] += fpidr*r2dj1*v_rphi[ir]/g[ig];
          } /* endfor */
        }/* endif */
      } /* endfor */

      break;

      /*==========================================================================*/
      /* III) L=2 Term */

    case 2:

      /*--------------------------------------------------------------------------*/
      /*    i) g ne 0                                                             */

      for(ig=1;ig <= nsplin_g; ig++) {
        for(ir=2;ir <= nr; ir++) {
          arg = r[ir]*g[ig];
          rj2 = (((3.0/(arg*arg)-1.0)*sin(arg)-3.0*cos(arg)/arg)/arg)*r[ir];
          fv_rphi[ig] += fpidr*rj2*v_rphi[ir];
        } /* endfor */
        if(cp_ptens_calc ==1 ){
          for(ir=2;ir <= nr; ir++) {
            arg = r[ir]*g[ig];
            j1  = (sin(arg)/arg - cos(arg))/arg;
            j3  = ((15.0/(arg*arg) - 6.0)*sin(arg)/arg + 
                (1.0 - 15.0/(arg*arg))*cos(arg))/arg;
            r2dj2 = ( 2.0*j1 - 3.0*j3 )*r[ir]*r[ir]/5.0;
            fdv_rphi[ig] += fpidr*r2dj2*v_rphi[ir]/g[ig];
          } /* endfor */
        }/* endif */
      } /* endfor */

      break;

      /*==========================================================================*/
      /* IV) L=3 Term */

    case 3:

      /*--------------------------------------------------------------------------*/
      /*    i) g ne 0                                                             */

      for(ig=1;ig <= nsplin_g; ig++) {
        for(ir=2;ir <= nr; ir++) {
          arg = r[ir]*g[ig];
          rj3 = (((15.0/(arg*arg) - 6.0)*sin(arg)/arg + 
                (1.0 - 15.0/(arg*arg))*cos(arg))/arg)*r[ir];
          fv_rphi[ig] += fpidr*rj3*v_rphi[ir];
        } /* endfor */
        if(cp_ptens_calc ==1 ){
          for(ir=2;ir <= nr; ir++) {
            arg = r[ir]*g[ig];
            j3 = ((15.0/(arg*arg) - 6.0)*sin(arg)/arg + 
                (1.0 - 15.0/(arg*arg))*cos(arg))/arg;
            j2 = ((3.0/(arg*arg)-1.0)*sin(arg)-3.0*cos(arg)/arg)/arg;
            j4 = (7.0/arg)*j3 - j2;
            r2dj3 = ( 3.0*j2-4.0*j4 )*r[ir]*r[ir]/7.0;
            fdv_rphi[ig] += fpidr*r2dj3*v_rphi[ir]/g[ig];
          } /* endfor */
        }/* endif */
      } /* endfor */

      break;

  } /* end switch */

  /*==========================================================================*/
  /*==========================================================================*/


  /*==========================================================================*/
  /* V) Add in long range part if necessary                                  */

  if(ilong == 1) {
    /*--------------------------------------------------------------------------*/
    /*   A) Coulomb                                                             */
    if((z_1 != 0.0) || (z_2 != 0.0)) {
      ztot = 0.0;
      if(cp_dual_grid_opt==2){ztot = z_1 + z_2;}
      for(ig=1;ig <= nsplin_g;ig++) {
        g2 = g[ig]*g[ig];
        fv_rphi[ig] -= 
          (z_1*fpi*exp((-0.25*g2/(alpha_1*alpha_1)))/g2
           +  z_2*fpi*exp((-0.25*g2/(alpha_2*alpha_2)))/g2
           -  ztot*fpi*(
             exp((-0.25*g2/(alpha_conv_dual*alpha_conv_dual))))/g2
          );
      }/* endfor */
      falpha_12        = 4.0*alpha_1*alpha_1;
      falpha_22        = 4.0*alpha_2*alpha_2; 
      falpha_conv_dual = 4.0*alpha_conv_dual*alpha_conv_dual; 
      (*gzvps) += fpi*(z_1/falpha_12 + z_2/falpha_22
          -ztot/falpha_conv_dual);
      if(cp_ptens_calc ==1 ){
        for(ig=1;ig <= nsplin_g;ig++) {
          g2 = g[ig]*g[ig];
          fdv_rphi[ig] += 
            ( (2.0*z_1*fpi*exp((-g2/falpha_12))/g2)
              *(1.0/(falpha_12)+1.0/g2)
              +  (2.0*z_2*fpi*exp((-g2/falpha_22))/g2)
              *(1.0/(falpha_22)+1.0/g2)
            );
        } /* endfor */
      }/* endif */
    } /* endif */
    /*--------------------------------------------------------------------------*/
    /*   B) Polarization corrections: */
    if(zpol != 0) {
      for(i=1;i<=7;i++) {
        angamma[i]  = ((double) (i-1))*gamma;
        angamma2[i] = angamma[i]*angamma[i];
      } /* endfor */
      c[1] =  1.0;  c[2] = -6.0;  c[3] = 15.0;
      c[4] =-20.0;  c[5] = 15.0;  c[6] = -6.0;
      c[7] =  1.0;
      for(i=2;i<=7;i++) {
        (*gzvps)  -= 0.50*zpol*fpi*c[i]*angamma[i]*log(angamma[i]);
      } /* endfor */
      if(cp_ptens_calc ==1 ){
        for(i=1;i<=7;i++) {
          for(ig=1;ig<=nsplin_g;ig++) {
            g2            = g[ig]*g[ig];  
            gm1           = 1.0/g[ig];
            gamm_plu      = angamma2[i]+g2;
            gamm_min      = angamma2[i]-g2;
            phi           = atan2(g[ig],angamma[i]);
            pre_c         = 0.50*zpol*tpi*c[i];
            fv_rphi[ig]  -= pre_c*(gm1*gamm_min*phi+angamma[i]*log(gamm_plu));
            fdv_rphi[ig] -= pre_c*(gm1*gamm_min*(-phi*gm1+angamma[i]/gamm_plu)
                -2.0*phi+2.0*g[ig]*angamma[i]/gamm_plu)*gm1;
          } /* endfor ig */
        } /* endfor i */
      }else{
        for(i=1;i<=7;i++) {
          for(ig=1;ig<=nsplin_g;ig++) {
            g2            = g[ig]*g[ig];  
            gm1           = 1.0/g[ig];
            gamm_plu      = angamma2[i]+g2;
            gamm_min      = angamma2[i]-g2;
            phi           = atan2(g[ig],angamma[i]);
            pre_c         = 0.50*zpol*tpi*c[i];
            fv_rphi[ig]  -= pre_c*(gm1*gamm_min*phi+angamma[i]*log(gamm_plu));
          } /* endfor ig */
        } /* endfor i */
      }/*endif*/
    } /* endif zpol */

  } /* endif ilong */

  /*-----------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/




/*==========================================================================*/

void vps_read_error(char *vps_file)

  /*==========================================================================*/
{/*Begin routine*/
  /*-----------------------------------------------------------------------*/

  PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
  PRINTF("Error reading input from electron-atom\n"); 
  PRINTF("pseudopotential file %s\n",vps_file);
  PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
  FFLUSH(stdout);
  EXIT(1);

  /*-----------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/



/*==========================================================================*/
/*==========================================================================*/

void get_dvl(int ngh,double *rgh,double *dvl,double *c0,
    double *c1, double *c2, double *c3,double dr,int nr)

  /*==========================================================================*/
{/*Begin routine*/

  int igh,iii;
  double r,h,h0;  
  /*-----------------------------------------------------------------------*/

  for(igh=1;igh<=ngh;igh++){
    r = rgh[igh]*JUERG_FACTOR;
    iii = (int)(r/dr) + 1;
    iii = MIN(iii,nr);
    iii = MAX(iii,1);
    h0  = (double)(iii-1)*dr;
    h = r-h0;
    dvl[igh] = ((c3[iii]*h+c2[iii])*h+c1[iii])*h+c0[iii];
  }/* endfor */

  /*-----------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/



/*==========================================================================*/
/*==========================================================================*/

void make_weight_gen(double *wgh_tmp,double *wgh,double *rgh, double *dvl,
    int ngh,int iang,int natm_typ_gh,
    int n_ang_max_gh)

  /*==========================================================================*/
{/*Begin routine*/

  int igh,iii,ioff;
  double juerg_wght,arg;
  double r;
  double fpisq;

  /*-----------------------------------------------------------------------*/

  fpisq = (4.0*M_PI)*(4.0*M_PI);

  ioff = (natm_typ_gh-1)*n_ang_max_gh*ngh + (iang-1)*ngh;

  for(igh=1;igh<=ngh;igh++){
    r   = rgh[igh]*JUERG_FACTOR;
    arg = rgh[igh]*rgh[igh];
    juerg_wght = JUERG_FACTOR*exp(arg);

    /* FIX UP FOR DUAL GRID OPTION */

#ifdef  JUERG_FACTOR_ON
    wgh[igh+ioff] = (wgh_tmp[igh]*fpisq*dvl[igh]*r*r*juerg_wght);  

#define CFPI_SQ_OFF
#ifdef  CFPI_SQ
    if(iang == 1){
      wgh[igh+ioff] = fpisq;
    }else{
      wgh[igh+ioff] = fpisq;
    }
#endif

#else
    wgh[igh+ioff] = (wgh_tmp[igh]*fpisq*dvl[igh]*r*r);
#endif


  }/* endfor */

  /*-----------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void make_cp_atom_list(CPATOM_MAPS *cpatom_maps,CPPSEUDO *cppseudo,
    int natm_tot, int natm_typ, double *q,int *iatm_atm_typ)
  //========================================================================
{//Begin subprogram: 
  //========================================================================
  // Create the list of cp-atoms

  int *cp_atm_flag = cpatom_maps->cp_atm_flag;
  int nab_initio   = 0;
  for(int i=1;i<=natm_tot;i++){
    if(cp_atm_flag[i]==1){nab_initio++;}
  }//endfor

  int *cp_atm_lst  = (int *)cmalloc(nab_initio*sizeof(int),"make_cp_atom_list")-1;
  int j = 1;
  for(int i=1;i<=natm_tot;i++){
    if(cp_atm_flag[i]==1){cp_atm_lst[j]=i;j++;}
  }//endfor

  cpatom_maps->cp_atm_lst = cp_atm_lst;
  cpatom_maps->nab_initio = nab_initio;

  if(nab_initio!=natm_tot){
    PRINTF("All Atoms must be Abinitio\n");
    EXIT(1);
  }//endif

  //------------------------------------------------------------------------
  // Set up the charge of each atom type : For now, all the same!!

  double *q_typ = (double *)cmalloc(natm_typ*sizeof(double),"cp_atm_list")-1;

  for(int i=1;i<=natm_tot;i++){q_typ[iatm_atm_typ[i]] = q[i];}

  for(int i=1;i<=natm_tot;i++){
    if(q[i]!=q_typ[iatm_atm_typ[i]]){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("Charges of all atoms with the same type label must\n");
      PRINTF("match for now: q[%d]=%g qtyp[%d]=%g\n",
          i,q[i],iatm_atm_typ[i],q_typ[iatm_atm_typ[i]]);
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }//endif
  }//endfor

  cppseudo->q_typ    = q_typ;
  cppseudo->natm_typ = natm_typ;

  //------------------------------------------------------------------------
  // Create external list

  int *natm_eext = (int *)cmalloc((natm_typ+1)*sizeof(int),"control_vps_params")-1;
  for(int i=1;i<=natm_typ;i++){natm_eext[i]=0;}
  for(int i=1;i<=natm_tot;i++){natm_eext[iatm_atm_typ[i]]+=1;}

  int natm_eext_max = 0;
  for(int i=1;i<=natm_typ;i++){natm_eext_max = MAX(natm_eext_max,natm_eext[i]);}

  int **map_eext = cmall_int_mat(1,natm_typ,1,natm_eext_max,"ctrl_vps_params");
  for(int i=1;i<=natm_typ;i++){natm_eext[i]=0;}
  for(int i=1;i<=natm_tot;i++){
    natm_eext[iatm_atm_typ[i]]+=1;
    int ic = natm_eext[iatm_atm_typ[i]];
    map_eext[iatm_atm_typ[i]][ic] = i;
  }//endfor
  natm_eext[(natm_typ+1)]=natm_tot;

  cppseudo->natm_eext     = natm_eext;
  cppseudo->natm_eext_max = natm_eext_max;
  cppseudo->map_eext      = map_eext;

  //------------------------------------------------------------------------
}//end routine 
//==========================================================================




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void create_non_local_list(CPPSEUDO *cppseudo,int natm_tot,
    int *iatm_atm_typ, int natm_typ)
  /*==========================================================================*/
{/*begin routine*/
  /*==========================================================================*/

  int ityp,lang,i,ic,jc,lang1;
  int n_ang_max1;

  int n_ang_max_all= cppseudo->n_ang_max;
  int n_ang_max    = cppseudo->n_ang_max_kb;
  int *ivps_label  = cppseudo->ivps_label;
  int *n_ang       = cppseudo->n_ang;
  int *loc_opt     = cppseudo->loc_opt;
  double *vpsnorm  = cppseudo->vpsnorm;
  int ees_on       = cppseudo->nonlocal.ees_on;
  int ees_eext_on  = cppseudo->nonlocal.ees_eext_on;
  int n_interp     = cppseudo->nonlocal.n_interp;

  /* malloced, setup and stored below*/
  int natm_nl;           /* number of total non-local atoms */
  int *natm_typ_lang;    /* # of non-local atm typs in lth channel */
  int **iatm_typ_lang;   /* atom type lookup: [ityp][lang]         */

  int *map_nl;           /* map[i]=j jth atm is ith non-local atom */
  int *natm_lang;        /* # atms of this type: iatm_typ_lang     */
  int *iatm_str_lang;    /* where atm is in list: iatm_typ_lang    */


  /*=========================================================================*/
  /* I) Some redudant checking */

  for(ityp=1;ityp<=natm_typ;ityp++){
    if(ivps_label[ityp]==1 && n_ang[ityp]==0 && loc_opt[ityp]==0){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("KB pseudo with max-l-channel=0 and loc_opt l=0?\n");
      PRINTF("This is a local pseudpotential. Please redefine it.\n");
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }/*endif*/
  }/*endfor*/

  /*=========================================================================*/
  /* II) Find the KB non-local atom types sorted by channel                  */

  n_ang_max1    = n_ang_max+1;
  natm_typ_lang=(int *)cmalloc(n_ang_max1*sizeof(int),
      "create_nonlocal_list")-1;
  iatm_typ_lang=cmall_int_mat(1,natm_typ,1,n_ang_max1,"create_nonlocal_list");

  for(lang=0,lang1=1;lang<=n_ang_max;lang++,lang1++){
    natm_typ_lang[lang1]=0;
    for(ityp=1;ityp<=natm_typ;ityp++){
      if(ivps_label[ityp]==1 && loc_opt[ityp]!=lang && n_ang[ityp]>=lang){
        natm_typ_lang[lang1]+=1;
        iatm_typ_lang[natm_typ_lang[lang1]][lang1]=ityp;
      }/*endif */
    }/*endfor : types*/
  }/*endfor : lang*/

  for(lang=0,lang1=1;lang<=n_ang_max;lang++,lang1++){
    PRINTF("   In the l=%d channel: there are %d active atom types\n",
        lang,natm_typ_lang[lang1]);
    if( (natm_typ_lang[lang1]>0) && (lang!=0) ){
      if(ees_on==1 && lang>=2){
        PRINTF("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");
        PRINTF("l>=2 non-locality experimental at present\n");
        PRINTF("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");
      }//endif
      if(ees_on==0){
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("l!=0 non-locality only supported with ees\n");
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        EXIT(1);
      }//endif
    }//endif
  }/*endfor*/

  /*=========================================================================*/
  /* III) Find the atoms */

  natm_nl = 0;
  for(i=1;i<=natm_tot;i++){
    ityp = iatm_atm_typ[i];
    if(ivps_label[ityp]==1){natm_nl++;}
  }/*endfor*/

  map_nl       = (int *)cmalloc(natm_nl*sizeof(int),"create_nonlocal_list")-1;
  natm_lang    = (int *)cmalloc(natm_typ*sizeof(int),"create_nonlocal_list")-1;
  iatm_str_lang= (int *)cmalloc(natm_typ*sizeof(int),"create_nonlocal_list")-1;

  ic = 0;
  for(ityp=1;ityp<=natm_typ;ityp++){
    natm_lang[ityp]     = 0;
    iatm_str_lang[ityp] = ic+1;
    if(ivps_label[ityp]==1){
      for(i=1;i<=natm_tot;i++){
        if(iatm_atm_typ[i]==ityp){
          ic++; natm_lang[ityp]++;
          map_nl[ic] = i;
        }/*endif*/
      }/*endfor*/
    }/*endif*/
  }/*endfor*/

  jc = 0; 
  for(ityp=1;ityp<=natm_typ;ityp++){jc+=natm_lang[ityp];}

  if( (ic!=natm_nl) || (jc!=natm_nl) ){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Internal error setting up non-local lists\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }/*endif*/

  PRINTF("   There are %d nonlocal atoms\n",natm_nl);

  /*===================================================================*/
  // Fill the norm for l=0 with 1 radial channel

  double *vnorm_0  = (double  *)cmalloc(natm_nl*sizeof(double),"nonloc_list")-1;
  for(i=1;i<=natm_nl;i++){vnorm_0[i]=0.0;}

  lang = 0; lang1 = 1;
  int natm_typ_now = natm_typ_lang[lang1];
  if(natm_typ_now>0){

    for(ityp=1;ityp<=natm_typ_now;ityp++){  // l=0 atm types
      int iatm_typ = iatm_typ_lang[ityp][lang1];// true atom type index
      int natm_now = natm_lang[iatm_typ];       // # atms of this type 
      int iatm_str = iatm_str_lang[iatm_typ];   // where atms strt in list
      int ind_now  = (iatm_typ-1)*(n_ang_max_all+1)+lang+1;
      for(int jatm=1;jatm<=natm_now;jatm++){
        int iatm = jatm+iatm_str-1;
        vnorm_0[iatm] = vpsnorm[ind_now];
      }//endfor
    }//endfor
  }//endif
  /*===================================================================*/
  // put it back
  cppseudo->nonlocal.vnorm_0       = vnorm_0; 
  cppseudo->nonlocal.natm_tot      = natm_tot;
  cppseudo->nonlocal.natm          = natm_nl;
  cppseudo->nonlocal.natm_typ      = natm_typ;
  cppseudo->nonlocal.nl_max        = n_ang_max;
  cppseudo->nonlocal.nl_max1       = (n_ang_max+1);

  cppseudo->nonlocal.natm_typ_lang = natm_typ_lang;
  cppseudo->nonlocal.iatm_typ_lang = iatm_typ_lang;
  cppseudo->nonlocal.map_nl        = map_nl;
  cppseudo->nonlocal.natm_lang     = natm_lang;
  cppseudo->nonlocal.iatm_str_lang = iatm_str_lang;
  if(ees_on==1 || ees_eext_on==1){
    nlEesSetIter(cppseudo);
  }

  /*-----------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void nlEesSetIter(CPPSEUDO *cppseudo){
  //==========================================================================

  int nl_max         = cppseudo->nonlocal.nl_max;
  int *natm_typ_lang = cppseudo->nonlocal.natm_typ_lang;
  int **iatm_typ_lang = cppseudo->nonlocal.iatm_typ_lang;
  int *natm_lang     = cppseudo->nonlocal.natm_lang;

  //==========================================================================
  // Count the number of iterations

  int nl_iter   = 0;
  int ntot_zmat = 0;  
  int nmax_zmat = 0;  
  for(int lang=0,lang1=1;lang<=nl_max;lang++,lang1++){// l-channels      
    int natm_typ = natm_typ_lang[lang1];// num atom types in this l-channel 
    if(natm_typ>0){
      for(int mang=-lang;mang<=lang;mang++){// m-channels in this l-channel 
        for(int ityp=1;ityp<=natm_typ;ityp++){// atm types in this channel  
          int iatm_typ = iatm_typ_lang[ityp][lang1]; // true atom type      
          int natm     = natm_lang[iatm_typ];        // # atms of this type 
          nl_iter     += 1;
          ntot_zmat   += natm;
          nmax_zmat   = MAX(nmax_zmat,natm);
        }//endfor
      }//endfor
    }//endfor
  }//endfor

  //==========================================================================
  // Malloc

  int *lang_v    = (int *)cmalloc((nl_iter+1)*sizeof(int),"nlEesSet")-1;
  int *mang_v    = (int *)cmalloc((nl_iter+1)*sizeof(int),"nlEesSet")-1;
  int *ityp_v    = (int *)cmalloc((nl_iter+1)*sizeof(int),"nlEesSet")-1;
  int *n_zmat    = (int *)cmalloc((nl_iter+1)*sizeof(int),"nlEesSet")-1;
  int *ioff_zmat = (int *)cmalloc((nl_iter+1)*sizeof(int),"nlEesSet")-1;

  //==========================================================================
  // Store the {l,m,ityp} of each iteration 

  int i=0;
  for(int lang=0,lang1=1;lang<=nl_max;lang++,lang1++){// l-channels      
    int natm_typ = natm_typ_lang[lang1];// num atom types in this l-channel 
    if(natm_typ>0){
      for(int mang=-lang;mang<=lang;mang++){// m-channels in this l-channel 
        for(int ityp=1;ityp<=natm_typ;ityp++){// atm types in this channel  
          int iatm_typ = iatm_typ_lang[ityp][lang1]; // true atom type      
          int natm     = natm_lang[iatm_typ];        // # atms of this type 
          i++;
          lang_v[i]    = lang;
          mang_v[i]    = mang;
          ityp_v[i]    = ityp;
          n_zmat[i]    = natm;
        }//endfor
      }//endfor
    }//endif
  }//endfor

  ioff_zmat[1] = 0;
  for(i=2;i<=nl_iter;i++){
    ioff_zmat[i] = ioff_zmat[(i-1)] + n_zmat[i];
  }//endfor

  //==========================================================================
  // Pack the data

  cppseudo->nonlocal.nl_iter   = nl_iter;
  cppseudo->nonlocal.ntot_zmat = ntot_zmat;
  cppseudo->nonlocal.nmax_zmat = nmax_zmat;
  cppseudo->nonlocal.lang_v    = lang_v;
  cppseudo->nonlocal.mang_v    = mang_v;
  cppseudo->nonlocal.ityp_v    = ityp_v;
  cppseudo->nonlocal.n_zmat    = n_zmat;
  cppseudo->nonlocal.ioff_zmat = ioff_zmat;

  //--------------------------------------------------------------------------
}//end routine
//==========================================================================

//==========================================================================
// spherical harmonic constants : Easiest to do here. They never change
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void set_ylm_cons(CPYLM_CONS *ylm_cons)
  //==========================================================================
{// begin routine*/
  //==========================================================================

  double tpi           = 2.0*M_PI;
  double fpi           = 4.0*M_PI;

  ylm_cons->rt_fpi     = 1.0/sqrt(fpi);
  ylm_cons->rt_thrfpi  = sqrt(3.0/fpi);
  ylm_cons->rt_threpi  = sqrt(1.50/fpi);
  ylm_cons->hrt_fivfpi = 0.50*sqrt(5.0/fpi);
  ylm_cons->rt_fiftepi = sqrt(7.50/fpi);
  ylm_cons->hrt_sevfpi = 0.50*sqrt(7.0/fpi);
  ylm_cons->hrt_toepi  = 0.50*sqrt(10.50/fpi)/sqrt(2.0);
  ylm_cons->hrt_ohffpi = 0.50*sqrt(105.0/fpi)/sqrt(2.0);
  ylm_cons->hrt_tfepi  = 0.50*sqrt(17.50/fpi)/sqrt(2.0);

  //--------------------------------------------------------------------------
}//end routine
//==========================================================================



PSSCRATCH::PSSCRATCH(PSNONLOCAL *_psnonlocal, CPATOM_MAPS *cpatom_maps)
{

  psnonlocal=_psnonlocal;

  int natm=psnonlocal->natm;
  int natm_tot=psnonlocal->natm_tot;
  int natm_nl = natm_tot;  // too big, but we don't really care.
  int ityp=0;

  /*===================================================================*/
  // Scratch nicely initialized for pupping

  x          = (double  *)cmalloc(natm_nl*sizeof(double),"nonloc_list")-1;
  y          = (double  *)cmalloc(natm_nl*sizeof(double),"nonloc_list")-1;
  z          = (double  *)cmalloc(natm_nl*sizeof(double),"nonloc_list")-1;

  // reused by hartree so a bit bigger
  index_atm  = (int *)cmalloc(natm_tot*sizeof(int),"nonloc_list")-1;
  vtemp      = (double *)cmalloc(natm_tot*sizeof(double),"nonloc_list")-1;
  ei_inc     = (complex *)cmalloc(natm_tot*sizeof(complex),"nonloc_list")-1;
  ti_inc     = (complex *)cmalloc(natm_tot*sizeof(complex),"nonloc_list")-1;
  for(int i=1;i<=natm_nl;i++){x[i]=0.0;}
  for(int i=1;i<=natm_nl;i++){y[i]=0.0;}
  for(int i=1;i<=natm_nl;i++){z[i]=0.0;}
  for(int i=1;i<=natm_tot;i++){ei_inc[i]=complex(0.0,0.0);}
  for(int i=1;i<=natm_tot;i++){ti_inc[i]=complex(0.0,0.0);}
  int n_interp=psnonlocal->n_interp;

  if(psnonlocal->ees_on==1 || psnonlocal->ees_eext_on==1){
    int natm = natm_tot;
    aj      = (double *)cmalloc(n_interp*sizeof(double),"psnl_pup")-1;
    rn      = (double *)cmalloc(n_interp*sizeof(double),"psnl_pup")-1;
    rn1     = (double *)cmalloc(n_interp*sizeof(double),"psnl_pup")-1;
    index_a = (int *)cmalloc(n_interp*sizeof(int),"psnl_pup")-1;
    index_b = (int *)cmalloc(n_interp*sizeof(int),"psnl_pup")-1;
    igrid_at= (int *)cmalloc(n_interp*sizeof(int),"psnl_pup")-1;
    igrid_bt= (int *)cmalloc(n_interp*sizeof(int),"psnl_pup")-1;
    frac_a  = (double *)cmalloc(natm*sizeof(double),"psnl_pup");
    frac_b  = (double *)cmalloc(natm*sizeof(double),"psnl_pup");
    frac_c  = (double *)cmalloc(natm*sizeof(double),"psnl_pup");
    iatemp  = (int *)cmalloc(natm*sizeof(int),"psnl_pup");
    ibtemp  = (int *)cmalloc(natm*sizeof(int),"psnl_pup");
    ictemp  = (int *)cmalloc(natm*sizeof(int),"psnl_pup");
    igrid_a = cmall_int_mat(1,n_interp,0,natm,"psnl_pup");
    igrid_b = cmall_int_mat(1,n_interp,0,natm,"psnl_pup");
    igrid_c = cmall_int_mat(1,n_interp,0,natm,"psnl_pup");
    mn_a    = cmall_mat(1,n_interp,0,natm,"psnl_pup");
    mn_b    = cmall_mat(1,n_interp,0,natm,"psnl_pup");
    mn_c    = cmall_mat(1,n_interp,0,natm,"psnl_pup");
    ua      = cmall_mat(1,n_interp,0,natm,"psnl_pup");
    ub      = cmall_mat(1,n_interp,0,natm,"psnl_pup");
    uc      = cmall_mat(1,n_interp,0,natm,"psnl_pup");
    dmn_a   = cmall_mat(1,n_interp,0,natm,"psnl_pup");
    dmn_b   = cmall_mat(1,n_interp,0,natm,"psnl_pup");
    dmn_c   = cmall_mat(1,n_interp,0,natm,"psnl_pup");
  }//endif
//----------------------------------------------------------------------------
  }
//============================================================================


//============================================================================
//****************************************************************************
//============================================================================
void control_grimme_params(FILENAME_PARSE *filename_parse,
            int *iatm_atm_typ, int natm_typ,NAME *atm_typ,
 	    double *tot_memory,int natm_tot,double *c6Grimme,double *r0Grimme)
//============================================================================
  {// begin routine
//============================================================================
  char *filename;              /* Char: temp file name                */
  CVPS *cgrimme_typ;
  char *fun_key;
  DICT_WORD *word;
  DICT_WORD *fun_dict;
  int num_fun_dict;
  DICT_WORD *grimme_dict,*grimme_dict_tmp;
  int num_grimme_dict;
  double *r0, *c6;
//============================================================================
// Tell the user your looking up the Grimme parameters
  PRINTF("\n");
  PRINT_LINE_STAR
  PRINTF("Searching the data bases both user defined and default\n");
  PRINTF("for the Grimme parameters\n");
  PRINT_LINE_DASH;PRINTF("\n");
//============================================================================
// Malloc the memory
  fun_key     = (char *)cmalloc(MAXWORD*sizeof(char),"control_grimme_params");  
  filename    = (char *)cmalloc(MAXWORD*sizeof(char),"control_grimme_params");  
  word        = (DICT_WORD *)cmalloc(sizeof(DICT_WORD),"control_grimme_params")-1;
  cgrimme_typ = (CVPS *)cmalloc(sizeof(CVPS),"control_vps_params");  
//============================================================================
// set up the dictionary

  int ifirst    = 1;
  set_potfun_dict(&fun_dict,&num_fun_dict,ifirst);
  set_grimme_dict(&grimme_dict,&num_grimme_dict,ifirst);      
  set_grimme_dict(&grimme_dict_tmp,&num_grimme_dict,ifirst);      

  c6 = (double *) cmalloc(natm_typ*sizeof(double),"control_vps_params")-1;
  r0 = (double *) cmalloc(natm_typ*sizeof(double),"control_vps_params")-1;

//============================================================================
//  Loop over all atom types and find the grimme parameters
  for(int i=1;i<=natm_typ;i++) {
//--------------------------------------------------
    int ifound = 0;
    strcpy(cgrimme_typ->atm1,atm_typ[i]);
   //--------------------------------------------------------------
   // A) Search the user define file if there is one
    if(strcasecmp(filename_parse->user_grimme_name,"")!=0) {
      search_base_grimme(filename_parse->user_grimme_name,
          cgrimme_typ,fun_dict,num_fun_dict,
          &grimme_dict_tmp,grimme_dict,num_grimme_dict,&ifound);
      if(ifound==1){
        set_grimme_params(grimme_dict,filename_parse->user_grimme_name,fun_key,
                          &c6[i],&r0[i]);
      }/*endif*/
    }/*endif*/
    /*--------------------------------------------------------------------------*/
    /* B) If you haven't found it search the default data base              */
    if(ifound == 0) {
      search_base_grimme(filename_parse->def_grimme_name,
          cgrimme_typ,fun_dict,num_fun_dict,
          &grimme_dict_tmp,grimme_dict,num_grimme_dict,&ifound);
      if(ifound==1){
        set_grimme_params(grimme_dict,filename_parse->def_grimme_name,fun_key,
                          &c6[i],&r0[i]);
      }/*endif*/
    }/*endif*/
    /*--------------------------------------------------------------------------*/
    /*     C) Make sure you have now found this puppy, if not exit              */
    if(ifound == 0){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("Grimme atom type \n"); 
      PRINTF("%s\n",atm_typ[i]);
      PRINTF("not found in default data base\n");
      PRINTF("file: %s\n",filename_parse->def_grimme_name);
      if(strlen(filename_parse->user_grimme_name) > 0)  {
        PRINTF("or in user defined pseudopot data base\n");
        PRINTF("file: %s\n",filename_parse->user_grimme_name);
        /*endif*/}
      PRINTF("\n");
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }/*endif*/
//--------------------------------------------------
  }//endfor : atom types
//============================================================================
// Make the grimme parameters a list of length natm_tot
  for(int i=1;i<=natm_tot;i++){
    c6Grimme[i] =  c6[iatm_atm_typ[i]];
    r0Grimme[i] =  r0[iatm_atm_typ[i]];
  }//endfor
//============================================================================
// clean up memory
  PRINTF("\n");PRINT_LINE_DASH
  PRINTF("Completed the Grimme data bases searches\n");
  PRINT_LINE_STAR;PRINTF("\n");
//============================================================================
// clean up memory
  cfree(fun_key, "control_grimme_params");
  cfree(filename,"control_grimme_params");
  cfree(&word[1],"control_grimme_params");
  cfree(&fun_dict[1],"control_grimme_params");
  cfree(&grimme_dict[1],"control_grimme_params");
  cfree(&grimme_dict_tmp[1],"control_grimme_params");
  cfree(cgrimme_typ,"control_grimme_params");
  cfree(&c6[1],"control_grimme_params");
  cfree(&r0[1],"control_grimme_params");
//----------------------------------------------------------------------------
  }//end routine
//============================================================================


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void set_grimme_params(DICT_WORD grimme_dict[],char *filename,char *fun_key,
                       double *c6, double *r0)
/*==========================================================================*/
/*               Begin subprogram:                                          */
  {/*begin routine*/
  /*--------------------------------------------------------------------------*/
  int index;
  double tmp;
  double Hart_per_J_Mol_m1 = 1.0/2625500.0;
  double Bohr_per_nm   = 10.0/BOHR;
  /*--------------------------------------------------------------------------*/
  /*  0) Set up                                                               */
  strcpy(fun_key,"grimme_parm");
  /*--------------------------------------------------------------------------*/
  /*  2) c6                                                */
  index=2;
  sscanf(grimme_dict[index].keyarg,"%lg",&tmp);
  if(tmp<0.0){keyarg_barf(grimme_dict,filename,fun_key,index);}
  c6[0] = (tmp*Hart_per_J_Mol_m1)*pow(Bohr_per_nm,6.0);
  /*--------------------------------------------------------------------------*/
  /*  3) r0                                                */
  index=3;
  sscanf(grimme_dict[index].keyarg,"%lg",&tmp);
  if(tmp<=0.0){keyarg_barf(grimme_dict,filename,fun_key,index);}
  r0[0] = tmp/BOHR;
 /*-----------------------------------------------------------------------*/
   }/*end routine*/
/*==========================================================================*/
