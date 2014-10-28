//==========================================================================
//               Pressure and Kinetic tensors                               

#ifndef _GENPTENS_
#define _GENPTENS_

class GENPTENS {

  //----------------
  public:
    double *tvten;              // Lst: KE  tensor           ; Lth: 9  
    double *pvten,*pvten_tot;   // Lst: PV  tensors          ; Lth: 9  

    //----------------
    //con-destruct:
    GENPTENS(){
      tvten     = (double *)cmalloc(9*sizeof(double),"class_genptens_constructor")-1;
      pvten     = (double *)cmalloc(9*sizeof(double),"class_genptens_constructor")-1;
      pvten_tot = (double *)cmalloc(9*sizeof(double),"class_genptens_constructor")-1;
    };
    ~GENPTENS(){};

#ifdef PUP_ON
    //----------------
    //pupping
    void pup(PUP::er &p){
      //pupping dbles
      pup1d_dbl(p,&tvten,9);
      pup1d_dbl(p,&pvten,9);
      pup1d_dbl(p,&pvten_tot,9);
#ifdef _PARALLEL_DEBUG_        
      if (p.isUnpacking())
        state_class_out ();
#endif       
    } // end pup
#endif

    void state_class_out(){

      char fileName [255];
      sprintf (fileName, "%d_gentens.state", CkMyPe());
      FILE *fp;  fp = fopen(fileName,"w");

      fprintf(fp,"nothing to report\n");
      fclose(fp);
    }// end routine

}; // GENPTENS;

#ifdef PUP_ON
PUPmarshall(GENPTENS);
#endif

#endif

//==========================================================================
