/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                            class_mdtors_free.h                             */
/*                                                                          */
/*        Class definition for tors free energy calculations                */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#ifndef _MDTORS_FREE_
#define _MDTORS_FREE_

class MDTORS_FREE {
  public:
    int num;                     /* Num: # of free energy tors          */
    int npow;                    /* Num: Tor power                      */
    int nhist;                   /* Num: # of pts in histogram          */

    int *j1,*j2,*j3,*j4; /* Num: Indices of atms in free E tors;*/


    double fk;                   /* Num: Tors  force constant           */ 
    double *eq;                /* Num: Equil tors angle               */
    double del;                  /* Num: Bin width of histogram         */

    char *file;                  /* Chr: Tors free energy output file   */
    double *hist;                /* Lst: Tors free energy histogram    
Lth:  nhist_tors_free               */
    double **hist_2d;            /* Lst: Tors free energy histogram    
Lth: nhist_tors_free*nhist_tors_free*/

    //============================================================================
    // Default Constructor/Destructor

    MDTORS_FREE(){
      num   = 0;
      npow  = 0;
      nhist = 0;
      j1 = (int *) cmalloc(3*sizeof(int),"class_mdtors_free_constructor")-1;
      j2 = (int *) cmalloc(3*sizeof(int),"class_mdtors_free_constructor")-1;
      j3 = (int *) cmalloc(3*sizeof(int),"class_mdtors_free_constructor")-1;
      j4 = (int *) cmalloc(3*sizeof(int),"class_mdtors_free_constructor")-1;
      eq =(double *)cmalloc(3*sizeof(double),"class_mdtors_free_constructor")-1;
    }
    ~MDTORS_FREE(){
      cfree(&j1[1],"class_mdtors_free_destructor");
      cfree(&j2[1],"class_mdtors_free_destructor");
      cfree(&j3[1],"class_mdtors_free_destructor");
      cfree(&j4[1],"class_mdtors_free_destructor");
      cfree(&eq[1],"class_mdtors_free_destructor");
    }

    //============================================================================
    // Pack/Unpack

#ifdef PUP_ON
    void pup(PUP::er &p){

      // Pupping ints

      p | num;
      p | npow;
      p | nhist;

      // PUP doubles

      p | fk;
      p | del;

      // PUP Arrays

      if(num>0){
        pup1d_int(p,&j1,3);
        pup1d_int(p,&j2,3);
        pup1d_int(p,&j3,3);
        pup1d_int(p,&j4,3);
        pup1d_dbl(p,&eq,3);
        pup1d_char(p,&file,MAXWORD);
        if(num == 1){pup1d_dbl(p,&hist,nhist);}
        if(num == 2){pup2d_dbl(p,&hist_2d,nhist,nhist,"mdtorsfree");}
      }/* endif */
#ifdef _PARALLEL_DEBUG_        
      if (p.isUnpacking())
        state_class_out ();
#endif 
    }/* end Pack/Unpack */
#endif

    //============================================================================
    // Print out state of class 

    void state_class_out(){
      char fileName [255];
      sprintf (fileName, "%d_mdtors_free.state", CkMyPe());
      FILE *fp;  fp = fopen(fileName,"w");    

      fprintf(fp,"mdtors_free:  num %d\n",num);
      fprintf(fp,"mdtors_free:  npow %d\n",npow);
      fprintf(fp,"mdtors_free:  nhist %d\n",nhist);

      fprintf(fp,"mdtors_free fk %.12g\n",fk);
      fprintf(fp,"mdtors_free del %.12g\n",del);

      for(int i=1;i<=3;i++){fprintf(fp,"mdtors_free:  j1[%d] %d\n",i,j1[i]);}
      for(int i=1;i<=3;i++){fprintf(fp,"mdtors_free:  j2[%d] %d\n",i,j2[i]);}
      for(int i=1;i<=3;i++){fprintf(fp,"mdtors_free:  j3[%d] %d\n",i,j3[i]);}
      for(int i=1;i<=3;i++){fprintf(fp,"mdtors_free:  eq[%d] %g\n",i,eq[i]);}
      fprintf(fp,"mdtors_free:  file %s\n",file);

      fclose(fp);
    }/* end member function */

    //---------------------------------------------------------------------------
}; /* end class definition */
//============================================================================

#ifdef PUP_ON
PUPmarshall(MDTORS_FREE);
#endif

#endif
//============================================================================
