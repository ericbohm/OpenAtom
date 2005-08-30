//==========================================================================
//                  Ensemble options                                        
//             {Variables needed for mem allocation:}                       

#ifndef _GENENSOPTS_
#define _GENENSOPTS_

class GENENSOPTS{
 //----------------
 public:
  int nve,nvt,npt_i,npt_f,nst;// Opt: Stat mech ensembles            

 //----------------
 //con-destruct:
  GENENSOPTS(){
    nve   = 0;
    nvt   = 0;
    npt_i = 0;
    npt_f = 0;
    nst   = 0;
  };
 ~GENENSOPTS(){};

#ifdef PUP_ON
 //----------------
 //pupping
  void pup(PUP::er &p){
    //pupping ints
        p | nve;
        p | nvt;
        p | npt_i;
        p | npt_f;
        p | nst;
#ifdef _PARALLEL_DEBUG_        
    if (p.isUnpacking())
     state_class_out ();
#endif        
  } // end pup
#endif

  void state_class_out(){

     char fileName [255];
     sprintf (fileName, "%d_genensopts.state", CkMyPe());
     FILE *fp;  fp = fopen(fileName,"w");

     fprintf(fp,"nve %d\n",nve);
     fprintf(fp,"nvt %d\n",nvt);
     fprintf(fp,"npt_i %d\n",npt_i);
     fprintf(fp,"npt_f %d\n",npt_f);
     fprintf(fp,"nst %d\n",nst);
   fclose(fp);
  }// end routine



}; // GENENSOPTS;

#ifdef PUP_ON
PUPmarshall(GENENSOPTS);
#endif

#endif

//==========================================================================
