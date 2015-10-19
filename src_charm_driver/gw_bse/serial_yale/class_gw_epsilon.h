//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//
//                  GWBSE simulation options                                   
//
//                class definition for GW_EPSILON
//                                                                          
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
#ifndef _GW_EPSILON_
#define _GW_EPSILON_


class GW_EPSILON{

  //---------------
  public:
    double Ecuteps;              // Num: Epsilon matrix cutoff
    double tol_iter;             // Num: Tolerance of the iterative matrix inversion method


    //----------------
    //con-destruct:
    GW_EPSILON(){
      Ecuteps         = 0;
      tol_iter        = 0;
    };
    ~GW_EPSILON(){};

#ifdef PUP_ON
    //-------------
    //pupping
    void pup(PUP::er &p){
      //pupping dbles
      p | Ecuteps;
      p | tol_iter;
      
#ifdef _PARALLEL_DEBUG_        
      if (p.isUnpacking())
        state_class_out ();
#endif
    } // end pup
#endif

    void state_class_out(){
      char fileName [255];
      sprintf ( fileName, "%d_gw_epsilon.out", CKMYPE());
      FILE *fp; fp = fopen(fileName,"w");
      //dbles
      fprintf(fp,"Ecuteps %lg\n",Ecuteps);
      fprintf(fp,"tol_iter %lg\n",tol_iter);
      fclose(fp);
    }// end routine
    
}; // GW_EPSILON

#ifdef PUP_ON
PUPmarshall(GW_EPSILON);
#endif

#endif
//==========================================================================
