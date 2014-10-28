/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                          class_mdbend_free.h                             */
/*                                                                          */
/*        Class definition for bend free energy calculations                */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#ifndef _MDBEND_FREE_
#define _MDBEND_FREE_

class MDBEND_FREE {
  public:
    int num;                     /* Num: # of free energy bends         */
    int j1,j2,j3;                /* Num: Indices of atms in Free E bend */
    int npow;                    /* Num: Bend power                     */
    int nhist;                   /* Num: # of pts in histogram          */

    double fk;                   /* Num: Bend force constant            */ 
    double eq;                   /* Num: Equil bend angle               */
    double del;                  /* Num: Bin width of histogram         */

    char *file;                  /* Chr: Bend free energy out filename  */ 
    double *hist;                /* Lst: Bend free energy histogram
Lth:  nhist                         */
    // Constructor/Destructor

    MDBEND_FREE(){
      num   = 0;                   
      j1    = 0;
      j2    = 0;
      j3    = 0;              
      npow  = 0;                  
      nhist = 0;                 
    }
    ~MDBEND_FREE(){}

#ifdef PUP_ON
    void pup(PUP::er &p){

      // PUP the ints
      p | num;
      p | j1;
      p | j2;
      p | j3;
      p | npow;
      p | nhist;

      // PUP the doubles

      p | fk;
      p | eq;
      p | del;

      // PUP the arrays

      if(num>0){
        pup1d_char(p,&file,MAXWORD);
        pup1d_dbl(p,&hist,nhist);
      }/*endif*/
#ifdef _PARALLEL_DEBUG_        
      if (p.isUnpacking())
        state_class_out ();
#endif 
    }/* end Pack/Unpack */
#endif

    void state_class_out(){
      char fileName [255];
      sprintf (fileName, "%d_mdbend_free.state", CkMyPe());
      FILE *fp;  fp = fopen(fileName,"w");

      // int
      fprintf(fp,"num %d\n",num);
      fprintf(fp,"j1 %d\n",j1);
      fprintf(fp,"j2 %d\n",j2);
      fprintf(fp,"j3 %d\n",j3);
      fprintf(fp,"npow %d\n",npow);
      fprintf(fp,"nhist %d\n",nhist);
      // dbles
      fprintf(fp,"fk %g\n",fk);
      fprintf(fp,"eq %g\n",eq);
      fprintf(fp,"del %g\n",del);

      fprintf(fp,"file %s\n",file);

      fclose(fp);
    }// end routine


}; /* class definition */

#ifdef PUP_ON
PUPmarshall(MDBEND_FREE);
#endif

#endif
