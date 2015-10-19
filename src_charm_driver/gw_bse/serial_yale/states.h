#ifndef __STATES_H__
#define __STATES_H__
#include "ckcomplex.h"

class states_occ : public CBase_states_occ {

 public:
  complex *stateCoeff;   // state coefficient in reciprocal space (G space)
  complex *stateCoeffR;  // state in R space (well.. it's not really coefficient in R space)
  int *ga, *gb, *gc;
  int numCoeff;
  int ikpt;        // index for k point
  int ispin;       // index for spin
  int countdebug;

  int ibinary_opt; // binary file option to read state file
  bool doublePack; // if only have gamma point, then true (=1). Otherwise, false(=0)
  char fileName[1000];  // file name for state

  int nfft[3]; // number of fft grid in each direction
  
  /// Constructors ///
  states_occ();
  states_occ(CkMigrateMessage *msg);
  /// Entry Methods ///                                                                         
  void beamoutMyState(int iteration, int qindex);
  void exitfordebugging();

  /// scalar routines ///
  void readState(char *);
  // fftw routines
  void fft_G_to_R();
  void set_fftsize();
  
};


class states_unocc : public CBase_states_unocc {

 public:
  complex *stateCoeff;   // state coefficient in reciprocal space (G space)
  complex *stateCoeffR;  // state in R space (well.. it's not really coefficient in R space)
  int *ga, *gb, *gc;
  int numCoeff;
  int ikpt;        // index for k point
  int ispin;       // index for spin

  int ibinary_opt; // binary file option to read state file
  bool doublePack; // if only have gamma point, then true (=1). Otherwise, false(=0)
  char fileName[1000];  // file name for state

  int nfft[3]; // number of fft grid in each direction
  
  /// Constructors ///                                                                            
  states_unocc();
  states_unocc(CkMigrateMessage *msg);

  /// Entry Methods ///                                                                         
  void beamoutMyState(int iteration, int qindex);

  /// scalar routines ///
  void readState(char *);
  // fftw routines
  void fft_G_to_R();
  void set_fftsize();
};




#endif //__STATES_H__
