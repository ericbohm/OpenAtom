//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                                                                          
//                            class_mdonfo.h                                
//                                                                          
//                Class definition for classical one-four interactions      
//                                                                          
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

#ifndef _MDONFO_
#define _MDONFO_

class MDONFO{
 public:
  int num;                     /* Num: # one-four interaction         */
  int ntyp;                    /* Num: # one-four typs                */ 

  int *j1,*j2;                 /* Lst: Indices of atms in 1-4s;     
                                  Lth: num                            */
  int *jtyp;                   /* Map: index of 1-4 -> 1-4 type;
                                  Lth: num                            */
  double *feps;                /* Lst: List of  one-four epsilons;
                                  Lth: ntyp                           */
  double *s6;                  /* Lst: List of 1-4 sigma^6's;     
                                  Lth: num                            */
  double *sc;                  /* Lst: List of 1-4 scaling factors;
                                  Lth: ntyp                           */
// Default constructor/destructor

   MDONFO(){
     num       = 0;
     ntyp      = 0;
   }
  ~MDONFO(){}

#ifdef PUP_ON
  void pup(PUP::er &p){
  // PUP integers
    p | num;         
    p | ntyp;        
    // PUP Arrays 
    if(num>0){
      pup1d_int(p,&j1,num);    
      pup1d_int(p,&j2,num);
      pup1d_int(p,&jtyp,num);
    }// endif
    if(ntyp>0){
      pup1d_dbl(p,&feps,ntyp);
      pup1d_dbl(p,&s6,ntyp);
      pup1d_dbl(p,&sc,ntyp);
    }// endif
#ifdef _PARALLEL_DEBUG_        
    if (p.isUnpacking())
     state_class_out ();
#endif 
  }/* End pack/unpack */
#endif

//------------------------------------------------------------------------------
// Print out state of class
  void state_class_out(){
     int i;
     char fileName [255];
     sprintf (fileName, "%d_mdonfo.state", CkMyPe());
     FILE *fp;  fp = fopen(fileName,"w");

   // scalars
     fprintf(fp,"mdonfo: num %d\n",num);
     fprintf(fp,"mdonfo: ntyp %d\n",ntyp);
   // Print Integer Arrays 
     for(i=1;i<=num;i++){fprintf(fp,"mdonfo: j1[%d] %d\n",i,j1[i]);}
     for(i=1;i<=num;i++){fprintf(fp,"mdonfo: j2[%d] %d\n",i,j2[i]);}
     for(i=1;i<=num;i++){fprintf(fp,"mdonfo: jtyp[%d] %d\n",i,jtyp[i]);}
   // Print Double Arrays
     for(i=1;i<=ntyp;i++){fprintf(fp,"mdonfo: feps[%d] %.12g\n",i,feps[i]);}
     for(i=1;i<=ntyp;i++){fprintf(fp,"mdonfo: s6[%d] %.12g\n",i,s6[i]);}
     for(i=1;i<=ntyp;i++){fprintf(fp,"mdonfo: sc[%d] %.12g\n",i,sc[i]);}

     fclose(fp);
  }// end member function

//------------------------------------------------------------------------------

}; // End class definition

#ifdef PUP_ON
PUPmarshall(MDONFO);
#endif

#endif
//==========================================================================
