//==========================================================================
//               PME for particles under dual gridding: should reuse other code
//             {Variables needed for mem allocation: }                      
//                                                                          

#ifndef _CPDUAL_PME_
#define _CPDUAL_PME_

class CPDUAL_PME{

 //----------------
 public:
 int cp_dual_grid_opt;        //Opt: on = 1 off = 0 
 int n_interp_dual_pme;
 int nkf1_cp_box;
 int nkf2_cp_box;
 int nkf3_cp_box;
 int nfft;

 int *iatemp,*ibtemp,*ictemp;      //Lth: nkf1(2,3)_dens_cp_box   
 int **igrid_a,**igrid_b,**igrid_c;//Lth: n_interp_dual_pme 
                                   //  X nkf1(2,3)_dens_cp_box
 int **igrid_now;                  //Lth: n_interp_dual_pme X nkf1_cp_box 


 double *a_pme,*b_pme,*c_pme;     //arrays to hold scaled grid points
                                  //Lth: nkf1(2,3)_dens_cp_box       
 double *frac_a,*frac_b,*frac_c;  //Lth: nkf1_cp_box
 double *aj,*rn,*rn1;             //Lth: n_interp_pme_dual
 double *bw_r,*bw_i;              //Lth: nfft  weighting factors             
                                  //  for r->g space pme                     
 double **mn_a,**mn_b,**mn_c;     //Lth: n_interp_pme_dual X nkf1(2,3)_cp_box
 double **ua,**ub,**uc;           //Lth: n_interp_pme_dual X nkf1(2,3)_cp_box

 //----------------
 //con-destruct:
   CPDUAL_PME(){
    cp_dual_grid_opt  = 0;
    n_interp_dual_pme = 0;
    nkf1_cp_box       = 0;
    nkf2_cp_box       = 0;
    nkf3_cp_box       = 0;
    nfft              = 0;
   };
  ~CPDUAL_PME(){};

#ifdef PUP_ON
 //----------------
 //pupping
  void pup(PUP::er &p){
    //pupping int
     p | cp_dual_grid_opt;
     p | n_interp_dual_pme;
     p | nkf1_cp_box;
     p | nkf2_cp_box;
     p | nkf3_cp_box;
     p | nfft;
    //pupping int arrays
#ifdef CPDUAL_IMPLEMENTED
     if(cp_dual_grid_opt>0){
      pup1d_int(p,&iatemp,nkf1_cp_box);
      pup1d_int(p,&ibtemp,nkf2_cp_box);
      pup1d_int(p,&ictemp,nkf3_cp_box);
      pup2d_int(p,&igrid_a,n_interp_dual_pme,nkf1_cp_box);
      pup2d_int(p,&igrid_b,n_interp_dual_pme,nkf2_cp_box);
      pup2d_int(p,&igrid_c,n_interp_dual_pme,nkf3_cp_box);
      pup2d_int(p,&igrid_now,n_interp_dual_pme,nkf1_cp_box);
    //pupping dble arrays
      pup1d_dbl(p,&a_pme,nkf1_cp_box);
      pup1d_dbl(p,&b_pme,nkf2_cp_box);
      pup1d_dbl(p,&c_pme,nkf3_cp_box);
      pup1d_dbl(p,&frac_a,nkf1_cp_box);
      pup1d_dbl(p,&frac_b,nkf2_cp_box);
      pup1d_dbl(p,&frac_c,nkf3_cp_box);
      pup1d_dbl(p,&aj,n_interp_dual_pme);
      pup1d_dbl(p,&rn,n_interp_dual_pme);
      pup1d_dbl(p,&rn1,n_interp_dual_pme);
      pup1d_dbl(p,&bw_r,nfft);
      pup1d_dbl(p,&bw_i,nfft);
      pup2d_dbl(p,&uc,n_interp_dual_pme,nkf1_cp_box);
      pup2d_dbl(p,&ua,n_interp_dual_pme,nkf2_cp_box);
      pup2d_dbl(p,&ub,n_interp_dual_pme,nkf3_cp_box);
      pup2d_dbl(p,&mn_a,n_interp_dual_pme,nkf1_cp_box);
      pup2d_dbl(p,&mn_b,n_interp_dual_pme,nkf2_cp_box);
      pup2d_dbl(p,&mn_c,n_interp_dual_pme,nkf3_cp_box);
     }// endif
#endif
#ifdef _PARALLEL_DEBUG_        
    if (p.isUnpacking())
     state_class_out ();
#endif     
  } // end pup
#endif

  void state_class_out(){

     char fileName [255];
     sprintf (fileName, "%d_cpdual_pme.state", CkMyPe());
     FILE *fp;  fp = fopen(fileName,"w");

     fprintf(fp," cp_dual_grid_opt %d\n", cp_dual_grid_opt);
     fprintf(fp," n_interp_dual_pme %d\n", n_interp_dual_pme);
     fprintf(fp," nkf1_cp_box %d\n", nkf1_cp_box);
     fprintf(fp," nkf2_cp_box %d\n", nkf2_cp_box);
     fprintf(fp," nkf3_cp_box %d\n", nkf3_cp_box);
     fprintf(fp," nfft %d\n", nfft);
   fclose(fp);

  }// end routine

}; //CPDUAL_PME

#ifdef PUP_ON
PUPmarshall(CPDUAL_PME);
#endif



#endif

