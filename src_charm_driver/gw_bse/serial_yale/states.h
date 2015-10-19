#ifndef __STATES_H__
#define __STATES_H__
#include "ckcomplex.h"

// Message sent from psi used to compute f. Some of these are cached on each
// node, and others are streamed in to the PMatrix as needed.
class PsiMessage : public CMessage_PsiMessage {
  public:
    PsiMessage(unsigned s, complex* p) : size(s) {
      std::copy(p, p+size, psi);
    }
    unsigned spin_index, k_index, state_index, size;
    complex* psi;
};


class PsiCache: public CBase_PsiCache {
  public:
    PsiCache();

    void receivePsi(PsiMessage*);
    complex* getPsi(unsigned) const;
  private:
    unsigned psi_count, psi_size, received_psis;
    complex** psis;
};

class states_occ : public CBase_states_occ {

 public:
  complex *stateCoeff;   // state coefficient in reciprocal space (G space)
  complex *stateCoeffR;  // state in R space (well.. it's not really coefficient in R space)
  int *ga, *gb, *gc;
  int numCoeff;
  int ikpt;        // index for k point
  int ispin;       // index for spin
  int istate;      // index for state
  int countdebug;

  int ibinary_opt; // binary file option to read state file
  bool doublePack; // if only have gamma point, then true (=1). Otherwise, false(=0)
  char fileName[1000];  // file name for state

  int nfft[3]; // number of fft grid in each direction
  
  /// Constructors ///
  states_occ();
  states_occ(CkMigrateMessage *msg);
  /// Entry Methods ///                                                                         
  void sendToCache();
  void sendToP();
  void beamoutMyState(int iteration, int qindex);
  void exitfordebugging();

  /// scalar routines ///
  void readState(char *);
  // fftw routines
  void fft_G_to_R();
  
};


class states_unocc : public CBase_states_unocc {

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
  states_unocc();
  states_unocc(CkMigrateMessage *msg);

  /// Entry Methods ///                                                                         
  void sendToCache();
  void sendToP();
  void beamoutMyState(int iteration, int qindex);

  /// scalar routines ///
  void readState(char *);
  // fftw routines
  void fft_G_to_R();
};




#endif //__STATES_H__
