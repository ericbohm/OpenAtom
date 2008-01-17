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
                          int cp_lsda,int cp_min_diis,int cp_dual_grid_opt_on,
                          PSNONLOCAL *psnonlocal, CPPSEUDO * cppseudo) 

/*=======================================================================*/
/*            Begin subprogram:                                          */
   {/*begin routine*/
/*=======================================================================*/
/*            Local variable declarations:                               */

   double ecut_now;                    /* Num: Energy cutoff                 */
   double ecut_lg;                     /* Num: Energy cutoff for dens        */
   double ecut_sm;                     /* Num: Energy cutoff for wf          */
   double *hmati_ewd;                  /* Num: Inverse ewald h matrix        */
   double *hmati_ewd_cp;               /* Num: Inverse cp ewald h matrix     */
   double deth,deth_cp,side;           /* Num: Volumes and sizes             */
   int iii;                            /* Num: Debug tool                    */
   int cp_on;                          /* Num: CP on flag                    */
   int nk1,nk2,nk3;                    /* Num: Number of k vec on small grid */
   int nkf1_dens_cp_box,
       nkf2_dens_cp_box,
       nkf3_dens_cp_box;               /* Num: Number of k vec on large  grid*/
   int nkf1,nkf2,nkf3;                 /* Num: Number of k vec on mega  grid */

   int nktot,nktot_sm;                 /* Num: Spherically cutoff k vecs     */
   int nktot_dens_cp_box;
   int ncoef,ncoef_l;

   int *kmaxv;                         /* Lst: K-vector ranges               */
   int *kmaxv_res,cp_on_tmp;
   int *kmax_cp;
   int *kmax_cp_dens_cp_box;

   /*--------------------------------------*/
   /* Local pointers                       */

   int iperd           = cell->iperd;
   double *hmat_ewd    = cell->hmat_ewd;
   double *hmat_ewd_cp = cell->hmat_ewd_cp;

   int  fft_opt        =  simopts->fftopt;

   int cp_nonlocal_ees_opt = psnonlocal->ees_on;
   int cp_eext_ees_opt     = psnonlocal->ees_eext_on;
   int *nfft;
   int *nfft_dens;

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

   if(int_res_ter==0){
     ewald->nktot_res=0;
     ecor->nktot_res=0;
   }//endif

   if(cp_on==0 || int_res_ter!=0){
     PRINTF("Sorry only CP for now\n"); EXIT(1);
   }//endif

/*=======================================================================*/
/* II) Allocate simple memory                                            */

   hmati_ewd      = (double *)cmalloc((size_t)9*sizeof(double),"control_set_cp_ewald")-1;
   hmati_ewd_cp   = (double *)cmalloc((size_t)9*sizeof(double),"control_set_cp_ewald")-1;

   kmaxv          =    (int *)cmalloc((size_t)3*sizeof(int),"control_set_cp_ewald")-1;
   kmax_cp        =    (int *)cmalloc((size_t)3*sizeof(int),"control_set_cp_ewald")-1;
   kmax_cp_dens_cp_box = (int *)cmalloc((size_t)3*sizeof(int),"control_set_cp_ewald")-1;

   nfft        =    (int *)cmalloc((size_t)3*sizeof(int),"control_set_cp_ewald")-1;
   nfft_dens   = (int *)cmalloc((size_t)3*sizeof(int),"control_set_cp_ewald")-1;

   cpewald->kmax_cp             = kmax_cp;
   cpewald->kmax_cp_dens_cp_box = kmax_cp_dens_cp_box;

   cpewald->nfft             = nfft;
   cpewald->nfft_dens        = nfft_dens;

/*==========================================================================*/
/* III) Get inverse cell matrix and convert ewald_alpha                     */

   gethinv(hmat_ewd,hmati_ewd,&deth,iperd);    
   gethinv(hmat_ewd_cp,hmati_ewd_cp,&deth_cp,iperd);    

   side  = pow(deth,(1.0/3.0));  
   (ewald->alp_ewd) /= side;
  
/*==========================================================================*/
/* IV) Calculate cutoff, count number k vectors, Malloc and Fill : Density grid */

/*----------------------------------------------------------------------*/
/* A) Calculate cutoff, count number k vectors, malloc and fill        */
/*    Without the dual box this is standard grid for cp/ewald          */
/*    With the dual box this is the small box calculation              */

   calc_cutoff(kmax_ewd,&ecut_now,&(cp_parse->cp_ecut),cp_on,kmax_cp,kmaxv,
               hmati_ewd_cp,deth_cp,nfft,fft_opt);  
   countkvec3d(&nktot,ecut_now,kmaxv,hmati_ewd_cp,gmin_spl,gmin_true,gmax_spl);

   ewald->nktot            = ewald->nktot;
   ewald->ecut             = 4.0*ecut_now;
   ewald->nka_max          = kmaxv[1];
   ewald->nkb_max          = kmaxv[2];
   ewald->nkc_max          = kmaxv[3];
   ecor->ecut              = 4.0*ecut_now;

   ncoef_l                 = nktot+1;
   ecut_sm                 = ecut_now;
   cpcoeffs_info->ecut     = ecut_now;
   cpcoeffs_info->ncoef_l  = ncoef_l;
   cpewald->nktot_lg       = nktot;

   kmax_cp_dens_cp_box[1] = kmax_cp[1];
   kmax_cp_dens_cp_box[2] = kmax_cp[2];
   kmax_cp_dens_cp_box[3] = kmax_cp[3];

   nfft_dens[1] = nfft[1];
   nfft_dens[2] = nfft[2];
   nfft_dens[3] = nfft[3];
 
/*=======================================================================*/
/* VII) Setup CP coefs and k vectors  : wavefunction grid                */

/*--------------------------------------------------------------------*/
/*  A)  Count the small k-vectors                                     */

   countkvec3d_sm(&nktot_sm,ecut_sm,kmax_cp,hmati_ewd_cp,
                  &(cpewald->gw_gmin),&(cpewald->gw_gmax));
   ncoef                  = nktot_sm+1;
   cpewald->nktot_sm      = nktot_sm;
   cpcoeffs_info->ncoef   = ncoef;

/*--------------------------------------------------------------------*/
/*  E) Set up the pme-nonlocal                                        */

   if(cp_nonlocal_ees_opt==1){
      init_nonlocal_ees(kmax_cp,ecut_sm,psnonlocal,fft_opt);
   }/*endif*/

   if(cp_eext_ees_opt==1){
      init_eext_ees(kmax_cp,cppseudo,fft_opt);
   }/*endif*/

/*=======================================================================*/
// Hessian options

   cpopts->cp_hess_cut = cp_parse->cp_mass_cut_def;
   cpopts->cp_hess_tau = cp_parse->cp_mass_tau_def;

/*=======================================================================*/
/* VIII) Output time has arrived                                         */

    nkf1 = nfft[1];
    nkf2 = nfft[2]; 
    nkf3 = nfft[3];

    nk1  = 2*kmax_cp_dens_cp_box[1];
    nk2  = 2*kmax_cp_dens_cp_box[2];
    nk3  = 2*kmax_cp_dens_cp_box[3];

    ecut_lg = 4.0*ecut_now;

    PRINTF("Your large cp-fft grid is  %d by %d by %d\n",nkf1,nkf2,nkf3);
    PRINTF("There are  %d total k-vectors ",ncoef_l); 
    PRINTF("upon spherical truncation. \n");
    PRINTF("The large energy cutoff is 4*Ecut= %f Ryd\n",2.0*ecut_lg);
    PRINTF("\n");

    PRINTF("Your small cp-fft grid is  %d by %d by %d\n",nk1,nk2,nk3);
    PRINTF("There are %d total k-vectors ",ncoef); 
    PRINTF("upon spherical truncation. \n");
    PRINTF("The small energy cutoff is Ecut= %f Ryd\n",2.0*ecut_sm);

/*=======================================================================*/

    PRINTF("\n");PRINT_LINE_DASH;
    PRINTF("Completed reciprocal space set up\n");
    PRINT_LINE_STAR;PRINTF("\n");

/*=======================================================================*/
/* IX) Free excess memory                                                */

    cfree(&(hmati_ewd)[1],"control_set_cp_ewald");
    cfree(&(hmati_ewd_cp)[1],"control_set_cp_ewald");
    cfree(&(kmaxv)[1],"control_set_cp_ewald");

/*-----------------------------------------------------------------------*/ 
  }/* end routine */
/*=======================================================================*/

