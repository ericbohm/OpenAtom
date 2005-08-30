/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                          class_mdbond_free.h                             */
/*                                                                          */
/*        Class definition for bond free energy calculations                */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#ifndef _MDBOND_FREE_
#define _MDBOND_FREE_

class MDBOND_FREE {
 public:
  int num;                     /* Num: # of free energy bonds         */
  int j1,j2;                   /* Num: Indices of atms in free E bonds*/
  int npow;                    /* Num: Bond power                     */
  int nhist;                   /* Num: # of pts in histogram          */

  double fk;                   /* Num: Bond force constant            */ 
  double eq;                   /* Num: Equil bond length              */
  double del;                  /* Num: Bin width of histogram         */
  double rmin;                 /* Num: Min distance in histogram      */
  double rmax;                 /* Num: Max distance in histogram      */

  char *file;                  /* Chr: Bond free energy outputfile    */ 
  double *hist;                /* Lst: Bond free energy histogram  
                                  Lth:  nhist_bond_free               */
// Constructor/Destructor

   MDBOND_FREE(){
    num   = 0;                  
    j1    = 0;
    j2    = 0;                
    npow  = 0;                 
    nhist = 0;                
   }
  ~MDBOND_FREE(){}

#ifdef PUP_ON
  void pup(PUP::er &p){

    // PUP the ints
    p | num;
    p | j1;
    p | j2;
    p | npow;
    p | nhist;

    // PUP the doubles

    p | fk;
    p | eq;
    p | del;
    p | rmin;
    p | rmax;

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
   sprintf (fileName, "%d_mdbond_free.state", CkMyPe());
   FILE *fp;  fp = fopen(fileName,"w");

    // int
    fprintf(fp,"num %d\n",num);
    fprintf(fp,"j1 %d\n",j1);
    fprintf(fp,"j2 %d\n",j2);
    fprintf(fp,"npow %d\n",npow);
    fprintf(fp,"nhist %d\n",nhist);
    // dbles
    fprintf(fp,"fk %g\n",fk);
    fprintf(fp,"eq %g\n",eq);
    fprintf(fp,"del %g\n",del);
    fprintf(fp,"rmin %g\n",rmin);
    fprintf(fp,"rmax %g\n",rmax);

    fprintf(fp,"file %s\n",file);

   fclose(fp);

  }// end routine

}; /* end class definition */

#ifdef PUP_ON
PUPmarshall(MDBOND_FREE);
#endif

#endif

