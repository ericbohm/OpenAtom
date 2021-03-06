#ifndef USRINPUT_H
#define USRINPUT_H

#include <cstdio>
#include <cstdlib>

class USRINPUT{

  public:

    int nstate;          // number of total states
    int nocc;            // number of occupied states
    int nunocc;          // number of unoccupied states
    int nkpt;            // number of k points
    int nspin;           // number of spin channels

    double EcutFFT;      // set FFTsize for both wavefunction and polarizability 
    double EpsCut;       // Epsilon matrix cutoff 
    double shift[3];     // shift vector (for epsilon calculation)

    int    iter_maxiter; // maximum number of iteration for matrix inversion
    double iter_convg;   // convergence of iterative matrix inversion

    // these variables are for reduced order P calculation
    // interpolation method variables
    int nPitp;           // number of spacing for P interpolation 
    // Laplace method variables
    double ptol;          // percentage tolerance 

    //constructor
    USRINPUT(){
      read_usrinput();
#ifdef DEBUG
      USRINPUT_class_out();
#endif
    }

    // member functions
    void read_usrinput();

    // tell user what are the inputs
    void USRINPUT_class_out(){
      FILE *fp; fp = fopen("USRINPUT_class_out","w");
      fprintf(fp,"number of states: %d\n", nstate);
      fprintf(fp,"number of occupied/unoccupied states: %d / %d\n",nocc,nunocc);
      fprintf(fp,"number of k points: %d\n",nkpt);
      fprintf(fp,"number of spin: %d\n",nspin);
      fprintf(fp,"FFT cutoff energy (Rydberg unit): %lg\n",EcutFFT);
      fprintf(fp,"Epsilon matrix cutoff (Hartree unit): %lg\n",EpsCut);
      fprintf(fp,"k point shift for Polarizability calculations: %lg %lg %lg\n",shift[0], shift[1], shift[2]);
      fprintf(fp,"maximum number of iteration for matrix inversion: %d\n",iter_maxiter);
      fprintf(fp,"convergence for iterative matrix inversion: %lg\n",iter_convg);
#ifdef USE_P_INTERPOLATION
      fprintf(fp,"number of spacing for P Interpolation: %d\n",nPitp);
#elif USE_P_LAPLACE
      fprintf(fp,"tolerance for 1/(Ec-Ev) for Laplace method: %lg\n",ptol);
#endif
      fclose(fp);
    }
    
      

};

#endif
