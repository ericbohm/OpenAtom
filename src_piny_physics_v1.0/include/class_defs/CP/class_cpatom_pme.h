//==========================================================================
//               PME for particles under dual gridding: should reuse other code
//             {Variables needed for mem allocation: }                      
//                                                                          

#ifndef _CPATOM_PME_
#define _CPATOM_PME_

class CPATOM_PME{

 //----------------
 public:
  int pme_on;                      //Opt: PME on 
  int n_interp;                    //Num: Order of interpolation 
  int nlen_pme;                    //Num: Scr lngth
  int nktot_pme;                   //Num: equal to ewald->nktot 
  int nkf1,nkf2,nkf3;              //Num: PME mesh same as large sparse grid 
  int nktot;                       // Num : spherically cutoff
  int ngrid_a,ngrid_b,ngrid_c;     //Num: PME mesh same as large sparse grid   

  int *iatemp,*ibtemp,*ictemp;     //Lst: Lth: nlen_pme                  
  int *nc,*ioff_c;                 //Lst: Lth nkf3                       

  int **igrid_a,**igrid_b,**igrid_c;//Lst: Lth: ninterp*nlen_pme
  int **igrid_now;                 //Lst: Lth: ninterp*nlen_pme


  double *frac_a,*frac_b,*frac_c;  //Lst: Lth:nlen_pme 
  double *aj,*rn,*rn1;             //Lst: Lth: ninterp

  double **ua,**ub,**uc;           //Lst: Lth:ninterp*nlen_pme
  double **mn_a,**mn_b,**mn_c;     //Lst: Lth: ninterp*nlen_pme
  double **dmn_a,**dmn_b,**dmn_c;  //Lst: Lth: ninterp*nlen_pme
  double **qgrid_now;              //Lst: Lth: ninterp*nlen_pme
  double *bweight_tot;             //Lst: Lth: nktot : 
  double *bw_r,*bw_i;              //Lth: nktot  weighting factors
                                   //  for r->g space pme
 //----------------
 //con-destruct:
   CPATOM_PME(){
    pme_on    = 0;                    
    n_interp  = 0;                  
    nlen_pme  = 0;                  
    nktot_pme = 0;                 
    nktot     = 0;  
    nkf1      = 0;
    nkf2      = 0;
    nkf3      = 0;
    nktot     = 0;                     
    ngrid_a   = 0;
    ngrid_b   = 0;
    ngrid_c   = 0;
   };
  ~CPATOM_PME(){};

#ifdef PUP_ON
 //----------------
 //pupping
  void pup(PUP::er &p){
    //pupping ints
      p | pme_on;
      p | n_interp;
      p | nlen_pme;
      p | nktot;
      p | nktot_pme;
      p | nkf1;
      p | nkf2;
      p | nkf3;
      p | ngrid_a;
      p | ngrid_b;
      p | ngrid_c;
    //pupping int arrays
    if(pme_on==1){
#ifdef PME_CP_IMPLEMENTED 
      pup1d_int(p,&iatemp,nlen_pme);
      pup1d_int(p,&ibtemp,nlen_pme);
      pup1d_int(p,&ictemp,nlen_pme);
      pup1d_int(p,&nc,nkf3);
      pup1d_int(p,&ioff_c,nkf3);
      pup2d_int(p,&igrid_a,n_interp,nlen_pme);
      pup2d_int(p,&igrid_b,n_interp,nlen_pme);
      pup2d_int(p,&igrid_c,n_interp,nlen_pme);
      pup2d_int(p,&igrid_now,n_interp,nlen_pme);
    //pupping dble arrays
      pup1d_dbl(p,&frac_a,nlen_pme);
      pup1d_dbl(p,&frac_b,nlen_pme);
      pup1d_dbl(p,&frac_c,nlen_pme);
      pup1d_dbl(p,&aj,n_interp);
      pup1d_dbl(p,&rn,n_interp);
      pup1d_dbl(p,&rn1,n_interp);
      pup1d_dbl(p,&bweight_tot,nktot);
      pup1d_dbl(p,&bw_r,nktot);
      pup1d_dbl(p,&bw_i,nktot);
      pup2d_dbl(p,&ua,n_interp,nlen_pme);
      pup2d_dbl(p,&ub,n_interp,nlen_pme);
      pup2d_dbl(p,&uc,n_interp,nlen_pme);
      pup2d_dbl(p,&mn_a,n_interp,nlen_pme);
      pup2d_dbl(p,&mn_b,n_interp,nlen_pme);
      pup2d_dbl(p,&mn_c,n_interp,nlen_pme);
      pup2d_dbl(p,&dmn_a,n_interp,nlen_pme);
      pup2d_dbl(p,&dmn_b,n_interp,nlen_pme);
      pup2d_dbl(p,&dmn_c,n_interp,nlen_pme);
      pup2d_dbl(p,&qgrid_now,n_interp,nlen_pme);
#endif

    }//endif
#ifdef _PARALLEL_DEBUG_        
    if (p.isUnpacking())
     state_class_out ();
#endif         
  } // end pup
#endif

  void state_class_out(){

    char fileName [255];
    sprintf (fileName, "%d_cpatom_pme.state", CkMyPe());
    FILE *fp;  fp = fopen(fileName,"w");

    fprintf(fp,"pme_on %d\n",pme_on);
    fprintf(fp,"n_interp %d\n",n_interp);
    fprintf(fp,"nlen_pme %d\n",nlen_pme);
    fprintf(fp,"nktot_pme %d\n",nktot_pme);
    fprintf(fp,"nkf1 %d\n",nkf1);
    fprintf(fp,"nkf2 %d\n",nkf2);
    fprintf(fp,"nkf3 %d\n",nkf3);
    fprintf(fp,"ngrid_a %d\n",ngrid_a);
    fprintf(fp,"ngrid_b %d\n",ngrid_b);
    fprintf(fp,"ngrid_c %d\n",ngrid_c);
   fclose(fp);

  }// end routine



}; //CPATOM_PME

#ifdef PUP_ON
PUPmarshall(CPATOM_PME);
#endif

#endif

//==========================================================================
