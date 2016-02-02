#ifndef __STATES_H__
#define __STATES_H__

#include "states.decl.h"
#include "ckcomplex.h"

class States : public CBase_States {

 public:
  complex *stateCoeff;   // state coefficient in reciprocal space (G space)
  complex *stateCoeffR;  // state in R space (well.. it's not really coefficient in R space)
  int *ga, *gb, *gc;
  int numCoeff;
  int ikpt;        // index for k point
  int ispin;       // index for spin
  int istate;      // index for state

  int ibinary_opt; // binary file option to read state file
  bool doublePack; // if only have gamma point, then true (=1). Otherwise, false(=0)
  char fileName[1000];  // file name for state

  int nfft[3]; // number of fft grid in each direction
  
  /// Constructors ///
  States();
  States(CkMigrateMessage *msg);

  /// Entry Methods ///                                                                         
  void fftGtoR();
  void sendToCache();
  void sendToP();

  /// fftw routines ///
  void fft_G_to_R();

  /// scalar routines ///
  void readState(char *);
  
};

extern /* readonly */ CProxy_States states_proxy;

#endif //__STATES_H__
