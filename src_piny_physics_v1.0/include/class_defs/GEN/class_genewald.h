//==========================================================================
//                 Indexing info for CP coeffs and Ewald                    
//             {Variables needed for mem allocation:                        
//                   nktot_sm,nktot,nktot_res}                              
//                                                                          

#ifndef _GENEWALD_
#define _GENEWALD_

#include "../proto_defs/proto_friend_lib_entry.h"

class GENEWALD {

 //----------------
 public:
  int cp_any_on;              // opt : any cp opt on
  int ewald_on;               // opt : ewald sum on
  int int_res_ter;            // opt : intermolecular respa
  int nsplin_g;               // NUm : spline pts for cp stuff
  int nktot;                  // Num: # of PW coeff on large sphere
                              //        cutoff g-space grid =ncoef_l    
  int nktot1;                 // Num: nktot+1
  int nktot_res;              // Num: # of PW coeffon RESPA sphere
                              //       cutoff g-space grid <=ncoef_l    
  int nktot_res1;             // Num: nktot_res+1

  int iperd;
  int nfix_para;              // nkvec_perd_fix
  int nfix_perp;              // nkvec_perd_expnd
  int nka_max;                // Num: Max value of ka, calculated by all 
  int nkb_max;                // Num: Max value of kb, calculated by all 
  int nkc_max;                // Num: Max value of kc, calculated by
			      // all 
  int ncorr_c;
  int ncorr_b;

  double alp_ewd;             // Num: Convergence param of Ewald sum  
  double self_erf;            
  double ecut,ecut_res;
  double ecut_clus;           // Num: Energy cutoff for cluster BCs  
  double alp_clus;            // Num: Ewald alpha for cluster BCs 
  double gs_max;              // Num: Max (gx^2+gy^2)^1/2 calc by all 

  int *ka_corr_c;
  int *kb_corr_c;
  int *kc_corr_c;
  int *ka_corr_b;
  int *kb_corr_b;
  int *kc_corr_b;
  double *kernel_corr_c;
  double *kernel_corr_b;

 //----------------
 //con-destruct:
   GENEWALD(){
    cp_any_on  = 0;
    ewald_on   = 0;
    int_res_ter= 0;
    nsplin_g   = 0;
    nktot      = 0;       
    nktot1     = 0;      
    nktot_res  = 0;   
    nktot_res1 = 0;  
    nfix_para  = 0;
    nfix_perp  = 1;
    nka_max    = 0;
    nkb_max    = 0;     
    nkc_max    = 0;     
    ncorr_c    = 0;
    ncorr_b    = 0;
    iperd      = 3;
   };
  ~GENEWALD(){
    if(ncorr_c>0){
      free(ka_corr_c);
      free(kb_corr_c);
      free(kc_corr_c);
      free(kernel_corr_c);
    }//endif
    if(ncorr_b>0){
      free(ka_corr_b);
      free(kb_corr_b);
      free(kc_corr_b);
      free(kernel_corr_b);
    }//endif
  };

#ifdef PUP_ON
 //----------------
 //pupping
  void pup(PUP::er &p){
    //pupping ints
        p | cp_any_on;
        p | ewald_on;
        p | int_res_ter;
        p | nsplin_g;
        p | nktot;
        p | nktot1;
        p | nktot_res;
        p | nktot_res1;
        p | nfix_para;
        p | nfix_perp;
        p | nka_max;
        p | nkb_max;
        p | nkc_max;
        p | ncorr_b;
        p | ncorr_c;
        p | iperd;
    //pupping dbles
        p | alp_ewd;
        p | self_erf;
        p | ecut;
        p | ecut_res;
        p | ecut_clus;
        p | alp_clus;
        p | gs_max;
    //pupping arrays that start at 0
	if(ncorr_c>0){
          if(p.isUnpacking()) {
            ka_corr_c    = (int *)cmalloc(ncorr_c*sizeof(int),"genewald");
            kb_corr_c    = (int *)cmalloc(ncorr_c*sizeof(int),"genewald");
            kc_corr_c    = (int *)cmalloc(ncorr_c*sizeof(int),"genewald");
            kernel_corr_c = (double *)cmalloc(ncorr_c*sizeof(double),"genewald");
	  }//endif
          PUParray(p,ka_corr_c,ncorr_c);
          PUParray(p,kb_corr_c,ncorr_c);
          PUParray(p,kc_corr_c,ncorr_c);
          PUParray(p,kernel_corr_c,ncorr_c);
	}//endif
	if(ncorr_b>0){
          if(p.isUnpacking()) {
            ka_corr_b    = (int *)cmalloc(ncorr_b*sizeof(int),"genewald");
            kb_corr_b    = (int *)cmalloc(ncorr_b*sizeof(int),"genewald");
            kc_corr_b    = (int *)cmalloc(ncorr_b*sizeof(int),"genewald");
            kernel_corr_b = (double *)cmalloc(ncorr_b*sizeof(double),"genewald");
	  }//endif
          PUParray(p,ka_corr_b,ncorr_b);
          PUParray(p,kb_corr_b,ncorr_b);
          PUParray(p,kc_corr_b,ncorr_b);
          PUParray(p,kernel_corr_b,ncorr_b);
	}//endif
#ifdef _PARALLEL_DEBUG_        
        if (p.isUnpacking())
          state_class_out ();
#endif
  } // end pup
#endif

  void state_class_out(){

     char fileName [255];
     sprintf (fileName, "%d_genewald.state", CkMyPe());
     FILE *fp;  fp = fopen(fileName,"w");

     fprintf(fp,"nsplin_g %d\n",nsplin_g);
     fprintf(fp,"nktot %d\n",nktot);
     fprintf(fp,"nktot1 %d\n",nktot1);
     fprintf(fp,"nktot_res %d\n",nktot_res);
     fprintf(fp,"nktot_res1 %d\n",nktot_res1);
     fprintf(fp,"nkc_max %d\n",nkc_max);
     fprintf(fp,"nkb_max %d\n",nkb_max);
    //pupping dbles
     fprintf(fp,"alp_ewd %g\n",alp_ewd);
     fprintf(fp,"self_erf %g\n",self_erf);
     fprintf(fp,"ecut %g\n",ecut);
     fprintf(fp,"ecut_res %g\n",ecut_res);
     fprintf(fp,"ecut_clus %g\n",ecut_clus);
     fprintf(fp,"alp_clus %g\n",alp_clus);
     fprintf(fp,"gs_max %g\n",gs_max);
   fclose(fp);

  } // end routine

}; // GENEWALD;

#ifdef PUP_ON
PUPmarshall(GENEWALD);
#endif

#endif

//==========================================================================
