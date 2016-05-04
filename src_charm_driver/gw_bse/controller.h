#ifndef __CONTROLLER_H__
#define __CONTROLLER_H__

#include <cstdlib>
#include "ckcomplex.h"

#include "fft_controller.decl.h"
#include "psi_cache.decl.h"
#include "controller.decl.h"

// Structure keeping track of all timers we report
struct Timers {
  // Setup timers
  double total_setup, chare_creation, fft_states, caches_filled;
  // Phase 1 timers
  double total_phase1, form_p, get_times;
  int fcomp_count, pcomp_count;
  double max_fcomp, avg_fcomp, total_pcomp, avg_pcomp;
  // Phase 2 timers
  double total_phase2, to1D, to2D, fft1, fft2, trans1, trans2;
  // Phase 3 timers
  double total_phase3;
};

class Controller : public CBase_Controller {
  Controller_SDAG_CODE
  public:
    Controller();

  private:
    bool do_output;
    unsigned debug_stages;
    unsigned K, L, M, pipeline_stages;
    unsigned next_K, next_state, total_sent, total_complete;
    Timers timers;
};

// A struct containing the required info for computing a set of f vectors for a
// single unoccupied state. For each f the equation is:
// f[i] = psi_occ[i] * psi_unocc[i].conj() * scaling_factor * umklapp_factor[i]
// f, psi vectors, and umklapp_factor all have 'size' elements
// e_unocc and umklapp_factor are the same for every f
// e_occ has an entry for each f to be computed
// The set of occupied psis are all occupied psis needed for the given unocc psi
struct FComputePacket {
  unsigned int size;
  double e_unocc;
  double* e_occ;
  complex* unocc_psi;
  complex** occ_psis;
  complex* umklapp_factor;
  complex* fs;
};

class PsiCache : public CBase_PsiCache {
  public:
    PsiCache();

    void receivePsi(PsiMessage*);
    void computeFs(PsiMessage*);
    void reportFTime();
    complex* getPsi(unsigned, unsigned, unsigned) const;
    complex* getF(unsigned) const;
  private:
    void kqIndex(unsigned, unsigned&, int*);
    void computeUmklappFactor(int*);

    // Used for CkLoop parameters
    FComputePacket f_packet;

    unsigned K, L, psi_size, received_psis, qindex;
    complex*** psis;
    complex* fs;
    complex* umklapp_factor;

    double total_time;
};

extern /* readonly */ CProxy_Controller controller_proxy;
extern /* readonly */ CProxy_PsiCache psi_cache_proxy;

#endif
