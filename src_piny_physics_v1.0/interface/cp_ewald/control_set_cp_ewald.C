/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: set_ewald                                    */
/*                                                                          */
/* This reads in and sets up the k-space                                    */
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
#include "../proto_defs/proto_cp_ewald_entry.h"
#include "../proto_defs/proto_cp_ewald_local.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_friend_lib_entry.h"

/*==========================================================================*/
#define DEBUG_CLUS_CORR_OFF
#define CHECK_CLUS_OFF

#define MAX_INT 12.0

#define ORIG_OFF
#define PME
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Control ewald/cp g-space initialization */
/*==========================================================================*/

void control_set_cp_ewald(GENSIMOPTS *simopts,GENCELL *cell,
                          CPCOEFFS_INFO *cpcoeffs_info, CPOPTS *cpopts,
                          GENEWALD *ewald, CPEWALD *cpewald,CP_PARSE *cp_parse,
                          double *gmin_true,double *gmin_spl,double *gmax_spl,
                          int kmax_ewd,int kmax_res,
                          double *tot_memory,int int_res_ter,
                          MDPART_MESH *part_mesh,MDECOR *ecor, 
                          int cp_lsda,int cp_min_diis,int cp_dual_grid_opt_on) 

/*=======================================================================*/
/*            Begin subprogram:                                          */
   {/*begin routine*/
/*=======================================================================*/
/*            Local variable declarations:                               */

   double ecut_now;                    /* Num: Energy cutoff                 */
   double ecut_dens_cp_box_now;        /* Num: Energy cutoff                 */
   double ecut_res,ecut_tmp;
   double ecut_lg;                     /* Num: Energy cutoff for dens        */
   double ecut_sm;                     /* Num: Energy cutoff for wf          */
   double *hmati_ewd;                  /* Num: Inverse ewald h matrix        */
   double *hmati_ewd_cp;               /* Num: Inverse cp ewald h matrix     */
   double deth,deth_cp,side;           /* Num: Volumes and sizes             */
   double now_mem;                     /* Num: Memory used here              */
   int idum1=1,idum2=1,idum3=1;
   int iii;                            /* Num: Debug tool                    */
   int nmall;
   int cp_on;                          /* Num: CP on flag                    */
   int nk1,nk2,nk3;                    /* Num: Number of k vec on small grid */
   int nkf1_dens_cp_box,
       nkf2_dens_cp_box,
       nkf3_dens_cp_box;               /* Num: Number of k vec on large  grid*/
   int nkf1,nkf2,nkf3;                 /* Num: Number of k vec on mega  grid */

   int nktot,nktot_res,nktot_sm;       /* Num: Spherically cutoff k vecs     */
   int nktot_dens_cp_box;
   int ncoef,ncoef_dens_cp_box,ncoef_l;

   int ngrid_a_res,ngrid_b_res,ngrid_c_res;
   int ngrid_a, ngrid_b, ngrid_c;

   int *kmaxv;                         /* Lst: K-vector ranges               */
   int *kmax_cp_tmp,*kmaxv_res,cp_on_tmp;
   int *kmax_cp;
   int *kmaxv_dens_cp_box;
   int *kmax_cp_dens_cp_box;

   int pme_b_opt;
   double *bfact_r, *bfact_i;

   double gmin_spl_tmp,gmin_true_tmp;  /* Num : Min/Max g-vectors            */
   double gmax_spl_tmp;

   /*--------------------------------------*/
   /* Local pointers                       */

   int pme_on          = part_mesh->pme_on;
   int pme_res_on      = part_mesh->pme_res_on;
   int n_interp        = part_mesh->n_interp;
   int n_interp_res    = part_mesh->n_interp_res;
   int iperd           = cell->iperd;

   double *hmat_ewd    = cell->hmat_ewd;
   double *hmat_ewd_cp = cell->hmat_ewd_cp;

   int box_rat      =   cpewald->box_rat;

/*=======================================================================*/
/* 0) Output to screen                                                   */

     PRINTF("\n");PRINT_LINE_STAR;
     PRINTF("Setting up reciprocal space\n");
     PRINT_LINE_DASH;PRINTF("\n");

/*=======================================================================*/
/* I) Set cp switch and initialize respa kvectors                        */

   cp_on = simopts->cp_min 
          +simopts->cp_wave_min
          +simopts->cp
          +simopts->cp_wave
          +simopts->debug_cp
          +simopts->cp_pimd
          +simopts->debug_cp_pimd+simopts->cp_wave_pimd
          +simopts->cp_wave_min_pimd;

   if(int_res_ter==0){ewald->nktot_res=0;ecor->nktot_res=0;}

/*=======================================================================*/
/* II) Allocate simple memory                                            */

   hmati_ewd      = (double *) cmalloc((size_t)9*sizeof(double),"control_set_cp_ewald")-1;
   hmati_ewd_cp   = (double *) cmalloc((size_t)9*sizeof(double),"control_set_cp_ewald")-1;
   kmaxv          =    (int *) cmalloc((size_t)3*sizeof(int),"control_set_cp_ewald")-1;
   kmaxv_res      =    (int *) cmalloc((size_t)3*sizeof(int),"control_set_cp_ewald")-1;
   kmax_cp_tmp    =    (int *) cmalloc((size_t)3*sizeof(int),"control_set_cp_ewald")-1;

   if(cp_on == 1){
     cpewald->kmax_cp = (int *) cmalloc((size_t)3*sizeof(int),"control_set_cp_ewald")-1;
     kmax_cp          = cpewald->kmax_cp;
     cpewald->kmax_cp_dens_cp_box = (int *) cmalloc((size_t)3*sizeof(int),"control_set_cp_ewald")-1;
     kmax_cp_dens_cp_box          = cpewald->kmax_cp_dens_cp_box;
     if(cp_dual_grid_opt_on >= 1){
       kmaxv_dens_cp_box = (int *) cmalloc((size_t)3*sizeof(int),"control_set_cp_ewald")-1;
     }/*endif*/
   }/* endif cp_on */


/*==========================================================================*/
/* III) Get inverse cell matrix and convert ewald_alpha                     */

   gethinv(hmat_ewd,hmati_ewd,&deth,iperd);    
   gethinv(hmat_ewd_cp,hmati_ewd_cp,&deth_cp,iperd);    

   side  = pow(deth,(1.0/3.0));  
   (ewald->alp_ewd) /= side;
  
/*==========================================================================*/
/* IV) Calculate cutoff, count number k vectors, Malloc and Fill            */

/*----------------------------------------------------------------------*/
/* A) Calculate cutoff, count number k vectors, malloc and fill        */
/*    Without the dual box this is standard grid for cp/ewald          */
/*    With the dual box this is the small box calculation              */

   calc_cutoff(kmax_ewd,&ecut_now,&(cp_parse->cp_ecut),cp_on,
                kmax_cp,kmaxv,hmati_ewd_cp,deth_cp);  

   countkvec3d(&(ewald->nktot),ecut_now,kmaxv,hmati_ewd_cp); 

   nktot                   = ewald->nktot;
   cpcoeffs_info->ecut     = ecut_now;
   cpcoeffs_info->ncoef_l  = nktot+1;
   cpewald->nktot_lg       = nktot;
   ncoef_l                 = nktot+1;
   ecor->ecut              = 4.0*ecut_now;
   ewald->ecut             = 4.0*ecut_now;
   ewald->nkc_max          = kmaxv[3];

/*----------------------------------------------------------------------*/
/* A.1) For dualing : Calculate cutoff and count kvectors for the large */
/*      box and save the small box.                                     */
if(cp_on == 1){
   switch(cp_dual_grid_opt_on){
    case 0:  
      ecut_dens_cp_box_now   = ecut_now;
      kmax_cp_dens_cp_box[1] = kmax_cp[1];
      kmax_cp_dens_cp_box[2] = kmax_cp[2];
      kmax_cp_dens_cp_box[3] = kmax_cp[3];
    break;
    case 1:  
      ecut_dens_cp_box_now               = ecut_now;
      nktot_dens_cp_box                  = nktot;
      cpewald->nktot_dens_cp_box         = nktot;
      cpcoeffs_info->ncoef_l_dens_cp_box = nktot + 1;
      ncoef_dens_cp_box                  = nktot + 1;
      cpcoeffs_info->ecut_dens_cp_box    = ecut_now;

      kmaxv_dens_cp_box[1]   = kmaxv[1];
      kmaxv_dens_cp_box[2]   = kmaxv[2];
      kmaxv_dens_cp_box[3]   = kmaxv[3];

      kmax_cp_dens_cp_box[1] = kmax_cp[1];
      kmax_cp_dens_cp_box[2] = kmax_cp[2];
      kmax_cp_dens_cp_box[3] = kmax_cp[3];

      nk1 = 4*(kmax_cp[1]+1);
      nk2 = 4*(kmax_cp[2]+1);
      nk3 = 4*(kmax_cp[3]+1);
/*!!!!  DERIVE THIS EXPRESSION !!!! */
/* How mamy g-vectors do you need to get the density on the large grid */
/* correctly if you fft to g-space with a spherical cut-off and fft back */
/* Glenn's analysis indicates thou shalt not truncate g space  at all   */
/*    for this option , great since this makes the pme a more important */
/*    method  */

      kmaxv[1] = 2*(box_rat*(kmax_cp[1]+1)-1);
      kmaxv[2] = 2*(box_rat*(kmax_cp[2]+1)-1);
      kmaxv[3] = 2*(box_rat*(kmax_cp[3]+1)-1);

      countkvec3d(&(ewald->nktot),ecut_now,kmaxv,hmati_ewd);

      nktot                   = ewald->nktot;
      cpcoeffs_info->ecut     = ecut_now;
      cpcoeffs_info->ncoef_l  = nktot+1;
      ncoef_l                 = nktot+1;
      ecor->ecut              = 4.0*ecut_now;
      ewald->ecut             = 4.0*ecut_now;
      ewald->nkc_max          = kmaxv[3];
    break;
#ifdef  ORIG
    case 2:
      ecut_dens_cp_box_now               = ecut_now;
      nktot_dens_cp_box                  = nktot;
      cpewald->nktot_dens_cp_box         = nktot;
      cpcoeffs_info->ncoef_l_dens_cp_box = nktot + 1;
      ncoef_dens_cp_box                  = nktot + 1;
      cpcoeffs_info->ecut_dens_cp_box    = ecut_now;

      kmaxv_dens_cp_box[1]   = kmaxv[1];
      kmaxv_dens_cp_box[2]   = kmaxv[2];
      kmaxv_dens_cp_box[3]   = kmaxv[3];

      kmax_cp_dens_cp_box[1] = kmax_cp[1];
      kmax_cp_dens_cp_box[2] = kmax_cp[2];
      kmax_cp_dens_cp_box[3] = kmax_cp[3];

      nk1 = 4*(kmax_cp[1]+1);
      nk2 = 4*(kmax_cp[2]+1);
      nk3 = 4*(kmax_cp[3]+1);
      kmaxv[1] = 2*(box_rat*(kmax_cp[1]+1)-1);
      kmaxv[2] = 2*(box_rat*(kmax_cp[2]+1)-1);
      kmaxv[3] = 2*(box_rat*(kmax_cp[3]+1)-1);

      countkvec3d(&(ewald->nktot),ecut_now,kmaxv,hmati_ewd);

      nktot                   = ewald->nktot;
      cpcoeffs_info->ecut     = ecut_now;
      cpcoeffs_info->ncoef_l  = nktot+1;
      ncoef_l                 = nktot+1;
      ecor->ecut              = 4.0*ecut_now;
      ewald->ecut             = 4.0*ecut_now;
      ewald->nkc_max          = kmaxv[3];
    break;
#endif
#ifdef  PME
    case 2:
      ecut_dens_cp_box_now               = ecut_now;
      nktot_dens_cp_box                  = nktot;
      cpewald->nktot_dens_cp_box         = nktot;
      cpcoeffs_info->ncoef_l_dens_cp_box = nktot + 1;
      ncoef_dens_cp_box                  = nktot + 1;
      cpcoeffs_info->ecut_dens_cp_box    = ecut_now;

      kmaxv_dens_cp_box[1]   = kmaxv[1];
      kmaxv_dens_cp_box[2]   = kmaxv[2];
      kmaxv_dens_cp_box[3]   = kmaxv[3];

      kmax_cp_dens_cp_box[1] = kmax_cp[1];
      kmax_cp_dens_cp_box[2] = kmax_cp[2];
      kmax_cp_dens_cp_box[3] = kmax_cp[3];

      calc_cutoff(kmax_ewd,&ecut_now,&(cp_parse->cp_ecut_dual_grid),cp_on,
                  kmax_cp,kmaxv,hmati_ewd,deth);  
      countkvec3d(&(ewald->nktot),ecut_now,kmaxv,hmati_ewd); 

      nktot                   = ewald->nktot;
      cpcoeffs_info->ecut     = ecut_now;
      cpcoeffs_info->ncoef_l  = nktot+1;
      ncoef_l                 = nktot+1;
      ecor->ecut              = 4.0*ecut_now;
      ewald->ecut             = 4.0*ecut_now;
      ewald->nkc_max          = kmaxv[3];
    break;
#endif
   }/*end switch*/
}/*endif cp_on*/
 
/*----------------------------------------------------------------------*/
/* B) Malloc                                                            */

   nmall =  nktot+1;   if((nmall % 2)==0){nmall++;}
   now_mem    = (nmall*(sizeof(double)*0 + sizeof(int)*5))*1.e-06;
   if(int_res_ter!=0){
      now_mem  += (nmall*(sizeof(double)*0 + sizeof(int)*6 ))*1.e-06;
   }/*endif*/
   (*tot_memory) += now_mem;

   ewald->kastr = (int *) cmalloc(nmall*sizeof(int),"control_set_cp_ewald")-1;
   ewald->kbstr = (int *) cmalloc(nmall*sizeof(int),"control_set_cp_ewald")-1;
   ewald->kcstr = (int *) cmalloc(nmall*sizeof(int),"control_set_cp_ewald")-1;
   ewald->ibrk1 = (int *) cmalloc(nmall*sizeof(int),"control_set_cp_ewald")-1;
   ewald->ibrk2 = (int *) cmalloc(nmall*sizeof(int),"control_set_cp_ewald")-1;
   if(int_res_ter != 0) {
      ewald->kastr_res      = (int *) cmalloc(nmall*sizeof(int),"control_set_cp_ewald")-1;
      ewald->kbstr_res      = (int *) cmalloc(nmall*sizeof(int),"control_set_cp_ewald")-1;
      ewald->kcstr_res      = (int *) cmalloc(nmall*sizeof(int),"control_set_cp_ewald")-1;
      ewald->ibrk1_res      = (int *) cmalloc(nmall*sizeof(int),"control_set_cp_ewald")-1;
      ewald->ibrk2_res      = (int *) cmalloc(nmall*sizeof(int),"control_set_cp_ewald")-1;
      ewald->ibrk3          = (int *) cmalloc(nmall*sizeof(int),"control_set_cp_ewald")-1;
   }/*endif*/

   PRINTF("Ewald allocation: %g Mbytes; Total memory %g Mbytes\n",
           now_mem,*tot_memory);

/*------------------------------------------------------------------------*/

   if( cp_dual_grid_opt_on >= 1 && cp_on == 1){
     nmall      =  nktot_dens_cp_box+1;   if((nmall % 2)==0){nmall++;}
     now_mem    = (nmall*(sizeof(double)*0 + sizeof(int)*5))*1.e-06;
     (*tot_memory) += now_mem;

     cpewald->kastr_dens_cp_box    = (int *) cmalloc(nmall*sizeof(int),"control_set_cp_ewald")-1;
     cpewald->kbstr_dens_cp_box    = (int *) cmalloc(nmall*sizeof(int),"control_set_cp_ewald")-1;
     cpewald->kcstr_dens_cp_box    = (int *) cmalloc(nmall*sizeof(int),"control_set_cp_ewald")-1;
     cpewald->ibrk1_dens_cp_box    = (int *) cmalloc(nmall*sizeof(int),"control_set_cp_ewald")-1;
     cpewald->ibrk2_dens_cp_box    = (int *) cmalloc(nmall*sizeof(int),"control_set_cp_ewald")-1;

     PRINTF("cp grid allocation: %g Mbytes; Total memory %g Mbytes\n",
            now_mem,*tot_memory);

   }/*endif cp_dual_grid_opt_on*/

/*------------------------------------------------------------------------*/
/* C) Fill                                                                */
   setkvec3d(nktot,ecut_now,kmaxv,hmati_ewd,
             ewald->kastr,ewald->kbstr,ewald->kcstr,
             ewald->ibrk1,ewald->ibrk2,cp_on,
             gmin_spl,gmin_true,gmax_spl);
/*------------------------------------------------------------------------*/
/* C) Fill DENS_CP_BOX                                                    */

   if( cp_dual_grid_opt_on >= 1 && cp_on == 1){
     setkvec3d(nktot_dens_cp_box,ecut_dens_cp_box_now,
               kmaxv_dens_cp_box,hmati_ewd_cp,
               cpewald->kastr_dens_cp_box,
               cpewald->kbstr_dens_cp_box,cpewald->kcstr_dens_cp_box,
               cpewald->ibrk1_dens_cp_box,cpewald->ibrk2_dens_cp_box,cp_on,
               &gmin_spl_tmp,&gmin_true_tmp,&gmax_spl_tmp);
   }/*endif cp_dual_grid_opt_on*/

   if( cp_dual_grid_opt_on == 2 && cp_on == 1){
     *gmin_spl  = gmin_spl_tmp;
     *gmin_true = gmin_true_tmp;
     *gmax_spl  = gmax_spl_tmp;
   }/*endif*/

/*=======================================================================*/
/* V) Setup PME                                                          */

   if(pme_on==1){
/*-----------------------------------------------------------------------*/
/*  A) Get grid size */

      part_mesh->nktot_pme = nktot;
      set_pme_grid(ecut_now,deth,hmati_ewd,kmaxv,
                  &(part_mesh->ngrid_a),&(part_mesh->ngrid_b),
                  &(part_mesh->ngrid_c),n_interp,
                  part_mesh->kmax_pme);
      part_mesh->ecut     = 4.0*ecut_now;
      ngrid_a = part_mesh->ngrid_a;
      ngrid_b = part_mesh->ngrid_b;
      ngrid_c = part_mesh->ngrid_c;
      if((ngrid_a<n_interp) || (ngrid_b<n_interp) || (ngrid_c<n_interp)){
        PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("The PME n_interp parameter > number of grid points \n");
        PRINTF("This is not allowed\n");      
        PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
        FFLUSH(stdout);
        EXIT(1);      
      }/*endif*/
/*-----------------------------------------------------------------------*/
/*  B) Malloc */
      now_mem = ( (nktot)*(sizeof(double))
                +3*(n_interp) *(sizeof(double))  )*1.e-06;
      if(pme_res_on == 1) {
        now_mem += ((nktot) *(sizeof(double)) )*1.e-06;
      }/*endif*/
     *tot_memory += now_mem;

   /*------------------------*/
   /*  Malloc the bweights   */
     nmall                     =  nktot;
     part_mesh->bweight_tot    = (double *) cmalloc(nmall*sizeof(double),"control_set_cp_ewald")-1;
     bfact_r =  part_mesh->bweight_tot;
     bfact_i =  part_mesh->bweight_tot;
     if(pme_res_on == 1) {
      part_mesh->bweight_tot_res = (double *) cmalloc(nmall*sizeof(double),"control_set_cp_ewald")-1;
     }/*endif*/

   /*------------------------*/
   /*  Malloc some scratch  */
      nmall                          = n_interp;

      part_mesh->aj    = (double *)cmalloc(nmall*sizeof(double),"control_set_cp_ewald")-1;
      part_mesh->rn    = (double *)cmalloc(nmall*sizeof(double),"control_set_cp_ewald")-1;
      part_mesh->rn1   = (double *)cmalloc(nmall*sizeof(double),"control_set_cp_ewald")-1;

   /*-----------------------*/
   /*  Output               */

      PRINTF("PME allocation: %g Mbytes; Total memory %g Mbytes\n",
                                           now_mem,*tot_memory);

/*-----------------------------------------------------------------------*/
/*  C) Create PME Bweight */
      pme_b_opt = 1;
      set_pme_wght(nktot,ewald->kastr,ewald->kbstr,ewald->kcstr,
                   ngrid_a,ngrid_b,ngrid_c,
                   idum1,idum2,idum3,
                   pme_b_opt,bfact_r,bfact_i,
                   part_mesh->bweight_tot,n_interp,
                   part_mesh->aj,part_mesh->rn,part_mesh->rn1);
   }/*endif*/

/*=======================================================================*/
/* VI) Set up RESPA k vectors and RESPA PME                              */

   if(int_res_ter != 0) {
/*-----------------------------------------------------------------------*/
/*  A) Repsa K vectors */
      ecut_tmp = 0.0;cp_on_tmp=0;

      calc_cutoff(kmax_res,&ecut_res,&ecut_tmp,cp_on_tmp,
                  kmax_cp_tmp,kmaxv_res,hmati_ewd,deth);

      ecor->ecut_res   = 4.0*ecut_res;
      ewald->ecut_res  = 4.0*ecut_res;
      countkvec3d(&(ewald->nktot_res),ecut_res,kmaxv_res,hmati_ewd);
      nktot_res         = ewald->nktot_res;
      ecor->nktot_res   = nktot_res;
      setkvec3d(nktot_res,ecut_res,kmaxv_res,hmati_ewd,
                ewald->kastr_res,ewald->kbstr_res,ewald->kcstr_res,
                ewald->ibrk1_res,ewald->ibrk2_res,cp_on_tmp,
                &gmin_spl_tmp,&gmin_true_tmp,&gmax_spl_tmp);
      setkvec3d_res(kmax_res,hmati_ewd,
                    ewald->kastr,ewald->kbstr,ewald->kcstr,ewald->ibrk3,
                    nktot,nktot_res);
      if(pme_res_on==1){
/*-----------------------------------------------------------------------*/
/*  B) Set the PME GRID */
       part_mesh->nktot_pme_res = nktot_res;
       set_pme_grid(ecut_res,deth,hmati_ewd,kmaxv_res,
                  &(part_mesh->ngrid_a_res),&(part_mesh->ngrid_b_res),
                  &(part_mesh->ngrid_c_res),n_interp_res,
                  part_mesh->kmax_pme_res);
       part_mesh->ecut_res     = 4.0*ecut_res;
       ngrid_a_res = part_mesh->ngrid_a_res;
       ngrid_b_res = part_mesh->ngrid_b_res;
       ngrid_c_res = part_mesh->ngrid_c_res;
       if((ngrid_a_res<n_interp_res) || (ngrid_b_res<n_interp_res) || 
          (ngrid_c_res<n_interp_res)){
         PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
         PRINTF("The RESPA PME n_interp parameter > number of grid points \n");
         PRINTF("This is not allowed\n");      
         PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
        FFLUSH(stdout);
        EXIT(1);      
       }/*endif*/
/*-----------------------------------------------------------------------*/
/*  C) Create PME Bweight */
       pme_b_opt = 1;
       set_pme_wght(nktot_res,ewald->kastr_res,ewald->kbstr_res,
                  ewald->kcstr_res,ngrid_a_res,ngrid_b_res,ngrid_c_res,
                  idum1,idum2,idum3,
                  pme_b_opt,bfact_r,bfact_i,
                  part_mesh->bweight_tot_res,n_interp_res,
                  part_mesh->aj,part_mesh->rn,part_mesh->rn1);
      }/*endif*/
    }/*endif*/

/*=======================================================================*/
/* VII) Setup CP coefs and k vectors                                     */

   if(cp_on == 1) {

/*--------------------------------------------------------------------*/
/*  A)  Count the k-vectors                                           */
      ecut_sm = ecut_dens_cp_box_now;
      countkvec3d_sm(&(cpewald->nktot_sm),ecut_sm,
                     kmax_cp_dens_cp_box,hmati_ewd_cp);
      nktot_sm = cpewald->nktot_sm;
      cpcoeffs_info->ncoef   = nktot_sm+1;
      ncoef                  = nktot_sm+1;
/*--------------------------------------------------------------------*/
/*  B)  Malloc                                                       */
      nmall =  (nktot_sm+1);if((nmall % 2)==0){nmall++;}
      now_mem = (nmall*(sizeof(double)*0 + sizeof(int)*5 ))*1.e-06;
      *tot_memory += now_mem;
      cpewald->kastr_sm = (int *) cmalloc(nmall*sizeof(int),"control_set_cp_ewald")-1;
      cpewald->kbstr_sm = (int *) cmalloc(nmall*sizeof(int),"control_set_cp_ewald")-1;
      cpewald->kcstr_sm = (int *) cmalloc(nmall*sizeof(int),"control_set_cp_ewald")-1;
      cpewald->ibrk1_sm = (int *) cmalloc(nmall*sizeof(int),"control_set_cp_ewald")-1;
      cpewald->ibrk2_sm = (int *) cmalloc(nmall*sizeof(int),"control_set_cp_ewald")-1;

      nmall =  ncoef; if((nmall % 2)==0){nmall++;}
      cpcoeffs_info->cmass = (double *)cmalloc(nmall*sizeof(double),"control_set_cp_ewald")-1;

        PRINTF("CP allocation: %g Mbytes; Total memory %g Mbytes\n",
                 now_mem,*tot_memory);

/*--------------------------------------------------------------------*/
/*  C)  Fill and check                                                */

      setkvec3d_sm(nktot_sm,ecut_sm,kmax_cp_dens_cp_box,hmati_ewd_cp,
                   cpewald->kastr_sm,cpewald->kbstr_sm,cpewald->kcstr_sm,
                   cpewald->ibrk1_sm,cpewald->ibrk2_sm,
                   &(cpewald->gw_gmin),&(cpewald->gw_gmax));

    if(cp_dual_grid_opt_on == 0 && cp_on == 1){
      check_kvec(ewald->nktot,ewald->kastr,ewald->kbstr,ewald->kcstr,nktot_sm,
                  cpewald->kastr_sm,cpewald->kbstr_sm,cpewald->kcstr_sm);
    }

      if(cp_dual_grid_opt_on >= 1 && cp_on == 1){
        check_kvec(cpewald->nktot_dens_cp_box,cpewald->kastr_dens_cp_box,
                   cpewald->kbstr_dens_cp_box,
                   cpewald->kcstr_dens_cp_box,nktot_sm,
                   cpewald->kastr_sm,cpewald->kbstr_sm,cpewald->kcstr_sm);
      }/*endif cp_dual_grid_opt_on */   
/*--------------------------------------------------------------------*/
/*  D) Set up the cp masses                                           */

      set_cpmass(ncoef,cpewald->kastr_sm,
                 cpewald->kbstr_sm,cpewald->kcstr_sm,
                 cpcoeffs_info->cmass,hmati_ewd_cp,
                 &(cp_parse->cp_mass_tau_def),cp_parse->cp_mass_cut_def,
                 &(cpcoeffs_info->icmass_unif));

      cpopts->cp_hess_cut  = cp_parse->cp_mass_cut_def;
      cpopts->cp_hess_tau = cp_parse->cp_mass_tau_def;

   }/*endif:cpon*/

/*=======================================================================*/
/* VIII) Output time has arrived                                         */

 /*-----------------------------------------------------------------------*/
 /* A) CP output */
     if(cp_on == 1) {
       switch(cp_dual_grid_opt_on){
         case 0:
           nkf1 = 4*(kmax_cp[1]+1);
           nkf2 = 4*(kmax_cp[2]+1);
           nkf3 = 4*(kmax_cp[3]+1);
         break;
         case 1:
           nkf1 = 4*box_rat*(kmax_cp_dens_cp_box[1] + 1);
           nkf2 = 4*box_rat*(kmax_cp_dens_cp_box[2] + 1);
           nkf3 = 4*box_rat*(kmax_cp_dens_cp_box[3] + 1);

           nkf1_dens_cp_box = 4*(kmax_cp_dens_cp_box[1] + 1);
           nkf2_dens_cp_box = 4*(kmax_cp_dens_cp_box[2] + 1);
           nkf3_dens_cp_box = 4*(kmax_cp_dens_cp_box[3] + 1);
          break;
#ifdef ORIG
         case 2:
           nkf1 = 4*box_rat*(kmax_cp_dens_cp_box[1] + 1);
           nkf2 = 4*box_rat*(kmax_cp_dens_cp_box[2] + 1);
           nkf3 = 4*box_rat*(kmax_cp_dens_cp_box[3] + 1);

           nkf1_dens_cp_box = 4*(kmax_cp_dens_cp_box[1] + 1);
           nkf2_dens_cp_box = 4*(kmax_cp_dens_cp_box[2] + 1);
           nkf3_dens_cp_box = 4*(kmax_cp_dens_cp_box[3] + 1);
         break;
#endif
#ifdef  PME
         case 2:
           nkf1 = 4*(kmax_cp[1] + 1);
           nkf2 = 4*(kmax_cp[2] + 1);
           nkf3 = 4*(kmax_cp[3] + 1);

           nkf1_dens_cp_box = 4*(kmax_cp_dens_cp_box[1] + 1);
           nkf2_dens_cp_box = 4*(kmax_cp_dens_cp_box[2] + 1);
           nkf3_dens_cp_box = 4*(kmax_cp_dens_cp_box[3] + 1);
         break;
#endif
       }/*end switch*/

       nk1  = 2*kmax_cp_dens_cp_box[1];
       nk2  = 2*kmax_cp_dens_cp_box[2];
       nk3  = 2*kmax_cp_dens_cp_box[3];

       ecut_lg = 4.0*ecut_now;
       PRINTF("Your large cp-fft grid is  %d by %d by %d\n",nkf1,nkf2,nkf3);
       PRINTF("There are  %d total k-vectors ",ncoef_l); 
       PRINTF("upon spherical truncation. \n");
       PRINTF("The large energy cutoff is 4*Ecut= %f Ryd\n",2.0*ecut_lg);
       PRINTF("\n");

       if(cp_dual_grid_opt_on >= 1){
         PRINTF("Your large density grid for the cp-fft box is %d by %d by %d",
                 nkf1_dens_cp_box,nkf2_dens_cp_box,nkf3_dens_cp_box);
         PRINTF("\nThere are  %d total k-vectors ",ncoef_dens_cp_box); 
         PRINTF("upon spherical truncation. \n");
         PRINTF("The large energy cutoff is 4*Ecut= %f Ryd\n",8.0*ecut_sm);
         PRINTF("\n");
       }/*endif cp_dual_grid_opt_on */

       PRINTF("Your small cp-fft grid is  %d by %d by %d\n",nk1,nk2,nk3);
       PRINTF("There are %d total k-vectors ",ncoef); 
       PRINTF("upon spherical truncation. \n");
       PRINTF("The small energy cutoff is Ecut= %f Ryd\n",2.0*ecut_sm);

     }else {
 /*-----------------------------------------------------------------------*/
 /* B) Non-CP output */
       PRINTF("You are using %d k-vectors in your ewald sum\n",nktot);
       PRINTF("Your reciprocal space shape: (-%d,%d) by (-%d,%d) by (-%d,%d)",
              kmaxv[1],kmaxv[1],kmaxv[2],kmaxv[2],kmaxv[3],kmaxv[3]);
       PRINTF("\n");
     } /* endif:cp_on */

 /*-----------------------------------------------------------------------*/
 /* C) PME output  */
     if(pme_on==1){
       PRINTF("Your particle mesh grid shape: %d by %d by %d\n",
                                                   ngrid_a,ngrid_b,ngrid_c);
       if((pme_res_on==1)&&(int_res_ter != 0)&&(kmax_res > 0)){
         PRINTF("Your respa particle mesh grid shape: %d by %d by %d\n",
                                       ngrid_a_res,ngrid_b_res,ngrid_c_res);
       }/*endif*/
     }/*endif*/

     PRINTF("\n");PRINT_LINE_DASH;
     PRINTF("Completed reciprocal space set up\n");
     PRINT_LINE_STAR;PRINTF("\n");


/*=======================================================================*/
/* IX) Free excess memory                                                */

   cfree(&(hmati_ewd)[1],"control_set_cp_ewald");
   cfree(&(hmati_ewd_cp)[1],"control_set_cp_ewald");
   cfree(&(kmaxv)[1],"control_set_cp_ewald");
   cfree(&(kmaxv_res)[1],"control_set_cp_ewald");
   cfree(&(kmax_cp_tmp)[1],"control_set_cp_ewald");

/*-----------------------------------------------------------------------*/ 
  }/* end routine */
/*=======================================================================*/
