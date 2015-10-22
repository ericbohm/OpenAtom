#ifndef USRINPUT_H
#define USRINPUT_H

class USRINPUT{

  public:

    int nstate;          // number of total states
    int nocc;            // number of occupied states
    int nunocc;          // number of unoccupied states
    int nkpt;            // number of k points
    int nspin;           // number of spin channels

    double wfnFFTsize;   // set FFTsize for both wavefunction and polarizability 
    double EpsCut;       // Epsilon matrix cutoff 
    double shift[3];     // shift vector (for epsilon calculation)

    int    iter_maxiter; // maximum number of iteration for matrix inversion
    double iter_convg;   // convergence of iterative matrix inversion 

};

#endif
