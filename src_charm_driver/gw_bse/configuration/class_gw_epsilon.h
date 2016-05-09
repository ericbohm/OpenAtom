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
    int nspin;                   // Num: number of spins
    int nkpt;                    // Num: number of k points
    int nocc;                    // Num: number of occupied states
    int nunocc;                  // Num: number of unoccupied states
    double*** Eocc;              // Eigenvalues for occupied states
    double*** Eunocc;            // Eigenvalues for unoccupied states
    double EcutFFT;              // Num: Energy cutoff for FFT (Rydberg)
    double Ecuteps;              // Num: Epsilon matrix cutoff (Rydberg)
    double tol_iter;             // Num: Tolerance of the iterative matrix inversion method
    char eigFileName[200];       // No need to pup this, contents stored in Eocc and Eunocc

    //----------------
    //con-destruct:
    GW_EPSILON(){
      EcutFFT  = 0;
      Ecuteps  = 0;
      tol_iter = 0;
      nocc = 0;
      nunocc = 0;
      Eocc = NULL;
      Eunocc = NULL;
    };
    ~GW_EPSILON(){};

#ifdef PUP_ON
    //-------------
    //pupping
    void pup(PUP::er &p){
      //pupping dbles
      p | EcutFFT;
      p | Ecuteps;
      p | tol_iter;
      p | nspin;
      p | nkpt;
      p | nocc;
      p | nunocc;

      // pupping arrays
      if (p.isUnpacking()) {
        Eocc = new double**[nspin];
        Eunocc = new double**[nspin];
        for (int s = 0; s < nspin; s++) {
          Eocc[s] = new double*[nkpt];
          Eunocc[s] = new double*[nkpt];
          for (int k = 0; k < nkpt; k++) {
            Eocc[s][k] = new double[nocc];
            Eunocc[s][k] = new double[nunocc];
          }
        }
      }
      for (int s = 0; s < nspin; s++) {
        for (int k = 0; k < nkpt; k++) {
          PUParray(p, Eocc[s][k], nocc);
          PUParray(p, Eunocc[s][k], nunocc);
        }
      }
      
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
      fprintf(fp,"EcutFFT %lg\n",EcutFFT);
      fprintf(fp,"Ecuteps %lg\n",Ecuteps);
      fprintf(fp,"tol_iter %lg\n",tol_iter);
      fprintf(fp,"nocc %d\n",nocc);
      fprintf(fp,"nunocc %d\n",nunocc);
      for (int s = 0; s < nspin; s++) {
        for (int k = 0; k < nkpt; k++) {
          for (int i = 0; i < nocc; i++) {
            fprintf(fp,"Eocc[%d,%d,%d] = %lg\n", s,k,i, Eocc[s][k][i]);
          }
          for (int i = 0; i < nunocc; i++) {
            fprintf(fp,"Eunocc[%d,%d,%d] = %lg\n", s,k,i, Eunocc[s][k][i]);
          }
        }
      }
      fclose(fp);
    }// end routine
    
}; // GW_EPSILON

#ifdef PUP_ON
PUPmarshall(GW_EPSILON);
#endif

#endif
//==========================================================================
