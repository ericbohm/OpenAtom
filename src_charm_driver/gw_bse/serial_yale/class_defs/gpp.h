#ifndef GPPOPTS_H
#define GPPOPTS_H


#include "../include/ckcomplex.h"
#include "sysinfo.h"
#include "../util.h"
// class for generalized plasmon-pole model calculations

// class GPPOPTS reads variables from input file(s) and set the basic information for GPP calculations
class GPPOPTS{

 public:

  // input options
  bool gpp_is_on;      // if gpp is on(true/1) or not(false/0)
  bool qespresso;      // if rho and state data are from Quantum Espresso(true/1) or not(false/0)
  int numEnergy;       // number of points for energy that we want to calculate.
                       // should be read from the input file
  double *Energy;      // values to be evaluated. read in
  double delta;        // broadening factor
  
  int numq;            // number of q vectors
  char rhoFile[1000];  // file name that reads density data
  complex *rhoData;    // density data
  int *ga, *gb, *gc;   // g index for density
  int nr[3];           // number of data points in each direction

  // read rho data from file
  void readRho();
  void fft_R_to_G();

  SYSINFO *sys;        // it is good to have system information

  // constructor
  GPPOPTS(){};
  GPPOPTS(SYSINFO& _sys);

  void readInputFile(char*);
  void readRho(char*);
  
  // say what you have
  void state_class_out(){
    FILE *fp; fp = fopen("GPPOPTS_class.out", "w");
    fprintf(fp,"generalized plasmon-pole model calculation is on/off (1/0): %d\n", gpp_is_on);
    fprintf(fp,"Is the density file from Quantum Espresso?: %d\n",qespresso);
    fprintf(fp,"number of points for energy: %d\n", numEnergy);
    fprintf(fp,"Energy values to be evaluated:\n");
    for (int i=0; i<numEnergy; i++){
        fprintf(fp,"%d   %lg \n", i+1, Energy[i]);
    }
    fprintf(fp,"number of q points: %d\n",numq);
    fprintf(fp,"density file: %s\n", rhoFile);
    fclose(fp);
  }//end routines

};







// class GPPDATA is the actual class that holds the data at each frequency
#include "matrix.h"
#include "vcoulb.h"
class GPPDATA{
 public:

  GPPOPTS *opt;     // GPP input values

  int qIndex;       // this q index
  int ng;           // number of g vectors

  VCOULB *vc;       // coulomb potential
  CMATRIX *S;       // S matrix  (pointer to Epsmat) after the eigendecomposition, S will contain eigenvectors
  
  double *eigval;   // eigenvalues, real numbers
  double *omsq;     // omega square

  CMATRIX **Sw;     // frequency dependent S matrix
  
  // constructor
  GPPDATA(VCOULB* _vc, CMATRIX* epsmat, GPPOPTS gppopts){
    vc = _vc;
    qIndex = _vc->qIndex;
    ng = _vc->ng;
    S = epsmat;      // get Epsmat
    opt = &gppopts;  // get GPPOPTS
    // m'alloc
    eigval = new double [ng];
    omsq = new double [ng];
    Sw = new CMATRIX *[opt->numEnergy];
  }

  // member functions
  void compute(){
    calc_Smtrx();
    eigenDecomp();
    calc_Omsq();
    for(int iw=0; iw<opt->numEnergy; iw++){
#ifdef VERBOSE
      printf(" Calculating S(w) matrix at w index %d\n",iw);
#endif
      Sw[iw] = new CMATRIX( ng, ng );
      calc_Sw( iw );
      delete Sw[iw];
    }
  }

  void calc_Smtrx();
  void eigenDecomp();
  void calc_Omsq();
  void find_rho_gdiff(int [3], int&, bool &);
  void calc_Sw(int);

};



#endif
