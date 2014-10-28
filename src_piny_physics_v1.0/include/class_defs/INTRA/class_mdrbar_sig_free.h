//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                                                                          
//                            class_mdrbar_sigma.h                          
//                                                                          
//    Class definition for rbar-sigma free energy calculations              
//                                                                          
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

#ifndef _MDRBAR_SIG_FREE_
#define _MDRBAR_SIG_FREE_

class MDRBAR_SIG_FREE {
  public:
    int nfree;                   /* Num: # of rbar_sigma free energy bonds */
    int nhist_bar,nhist_sig;     /* Num: # of pts in histogram          */
    int iopt;                    /* Num: 1/0 if rbar-sig is on/off      */
    int *j1,*j2;                 /* Num: atm index in ith bond          */
    double fk_bar;               /* Num: rbar  force constant           */ 
    double fk_sigma;             /* Num: sigma  force constant          */ 
    double eq_bar;               /* Num: Equil mean                     */
    double eq_sigma;             /* Num: Equil std                      */
    double del_bar,del_sig;      /* Num: Bin width of histogram         */
    double rmin,rmax;            /* Num: Min/max distance in histogram  */
    double smin,smax;            /* Num: Min/max distance in histogram  */
    double rnfree;               /* Num: # of bonds                     */
    double **hist;               /* Lst:  Free energy histogram    
Lth:  nhist_sig*nhist_bar           */
    double **hist_rn;            /* Lst: Tors free energy histogram    
Lth: nfree*nhist_bar */
    char *file;                  /* Chr: Free energy output file        */

    // Default constructor/destructor

    MDRBAR_SIG_FREE(){
      nfree     = 0;                   
      nhist_bar = 0;
      nhist_sig = 0;     
      iopt      = 0;
    }
    ~MDRBAR_SIG_FREE(){}

    //---------------------------------------------------------------------------
#ifdef PUP_ON
    void pup(PUP::er &p){
      // PUP integers
      p | nfree;
      p | nhist_bar;
      p | nhist_sig;
      p | iopt;
      // PUP doubles
      p | fk_bar;        
      p | fk_sigma;      
      p | eq_bar;        
      p | eq_sigma;      
      p | del_bar;
      p | del_sig;
      p | rmin;
      p | rmax;      
      p | smin;
      p | smax;      
      p | rnfree;         
      // PUP Arrays
      if(nfree>0){
        pup1d_int(p,&j1,nfree);    
        pup1d_int(p,&j2,nfree);    
        pup2d_dbl(p,&hist,nhist_sig,nhist_bar,"mdbarsig");
        pup2d_dbl(p,&hist_rn,nfree,nhist_bar,"mdbarsig");
        pup1d_char(p,&file,MAXWORD);
      }//endif
#ifdef _PARALLEL_DEBUG_        
      if (p.isUnpacking())
        state_class_out ();
#endif 
    }// end pack/unpack
#endif


    void state_class_out(){
      int i;
      char fileName [255];
      sprintf (fileName, "%d_mdrbar_sig_free.state", CkMyPe());
      FILE *fp;  fp = fopen(fileName,"w");


      // integers
      fprintf(fp,"nfree %d\n",nfree);
      fprintf(fp,"nhist_bar %d\n",nhist_bar);
      fprintf(fp,"nhist_sig %d\n",nhist_sig);
      fprintf(fp,"iopt %d\n",iopt);
      // doubles
      fprintf(fp,"fk_bar %g\n",fk_bar);        
      fprintf(fp,"fk_sigma %g\n",fk_sigma);      
      fprintf(fp,"eq_bar %g\n",eq_bar);        
      fprintf(fp,"eq_sigma %g\n",eq_sigma);      
      fprintf(fp,"del_bar,del_sig %g %g\n",del_bar,del_sig);
      fprintf(fp,"rmin,rmax %g %g\n",rmin,rmax);      
      fprintf(fp,"smin,smax %g %g\n",smin,smax);      
      fprintf(fp,"rnfree %g\n",rnfree);         
      // Arrays
      for(i=1;i<=nfree;i++){fprintf(fp,"j1[%d] %d\n",i,j1[i]);}
      for(i=1;i<=nfree;i++){fprintf(fp,"j2[%d] %d\n",i,j2[i]);}
      if(nfree > 0){
        fprintf(fp,"file %s\n",file);
      }// endif
      fclose(fp);
    }// end routine

}; // End class definition

#ifdef PUP_ON
PUPmarshall(MDRBAR_SIG_FREE);
#endif

#endif
//==========================================================================
