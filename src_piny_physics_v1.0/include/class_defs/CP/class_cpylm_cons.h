//==========================================================================
//                 Spherical Harmonic constants                             
//             {Variables needed for mem allocation:                        
//                                                  }                       
//                                                                          

#ifndef _CPYLM_CONS_
#define _CPYLM_CONS_

class CPYLM_CONS {

 //----------------
 public :
  double rt_fpi,rt_thrfpi,rt_threpi;
  double hrt_fivfpi,rt_fiftepi,hrt_sevfpi;
  double hrt_toepi,hrt_ohffpi,hrt_tfepi;

 //----------------
 //con-destruct:
   CPYLM_CONS(){
      rt_fpi     = 0;
      rt_thrfpi  = 0;
      rt_threpi  = 0;
      hrt_fivfpi = 0;
      rt_fiftepi = 0;
      hrt_sevfpi = 0;
      hrt_toepi  = 0;
      hrt_ohffpi = 0;
      hrt_tfepi  = 0;
   };
  ~CPYLM_CONS(){};

#ifdef PUP_ON
 //----------------
 //pupping
  void pup(PUP::er &p){

    //pupping dbles
      p | rt_fpi;
      p | rt_thrfpi;
      p | rt_threpi;
      p | hrt_fivfpi;
      p | rt_fiftepi;
      p | hrt_sevfpi;
      p | hrt_toepi;
      p | hrt_ohffpi;
      p | hrt_tfepi;
#ifdef _PARALLEL_DEBUG_        
    if (p.isUnpacking())
     state_class_out ();
#endif
  } // end pup
#endif

  void state_class_out(){
     char fileName [255];
     sprintf (fileName, "%d_cpylm_cons.state", CkMyPe());
     FILE *fp;  fp = fopen(fileName,"w");

     fprintf(fp,"rt_fpi %g\n",rt_fpi);
     fprintf(fp,"rt_thrfpi %g\n",rt_thrfpi);
     fprintf(fp,"rt_threpi %g\n",rt_threpi);
     fprintf(fp,"hrt_fivfpi %g\n",hrt_fivfpi);
     fprintf(fp,"rt_fiftepi %g\n",rt_fiftepi);
     fprintf(fp,"hrt_sevfpi %g\n",hrt_sevfpi);
     fprintf(fp,"hrt_toepi %g\n",hrt_toepi);
     fprintf(fp,"hrt_ohffpi %g\n",hrt_ohffpi);
     fprintf(fp,"hrt_tfepi %g\n",hrt_tfepi);
   fclose(fp);
  } // end routine 


}; //CPYLM_CONS

#ifdef PUP_ON
PUPmarshall(CPYLM_CONS);
#endif

#endif


//==========================================================================
