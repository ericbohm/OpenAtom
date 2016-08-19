#ifndef __CONTROLLER_H__
#define __CONTROLLER_H__

#include <cstdlib>
#include "ckcomplex.h"

#include "fft_controller.decl.h"
#include "psi_cache.decl.h"

#include "controller.decl.h"

#define ITERATION 100 //needs to be read from epsilon.in
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
    void prep();
    void createEpsilon();
    void calc_Geps();
    void got_geps(std::vector<int> accept, int epsilon_size);
    void got_vcoulb(std::vector<double> vcoulb_in);
  private:
    bool do_output;
    int msg_received;
    int iter, maxiter;
    int iteration;
    unsigned dimension, rows;
    bool resultInsert;
    double epsCut;
    double alat;
    double vol;
    std::vector<double> vcoulb;
    double shift[3];
    unsigned K, L, M, pipeline_stages;
    unsigned next_K, next_state, total_sent, total_complete;
    unsigned max_sends, next_report_threshold;
    unsigned p_matrix_dimension, num_p_rows;
    int global_inew, global_jnew;
    int max_local_inew;
    int jnew_arr[140];
    int padded_epsilon_size;
    complex tmp_arr[140*140];
    double prev_max;
    std::vector<int> accept_result;
    CLA_Matrix_interface matA, matB, matC, matA2, matB2, matC2, matA3, matB3, matC3;
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
    int getWrote();
  private:
    void kqIndex(unsigned, unsigned&, int*);
    void computeUmklappFactor(int*);

    // Used for CkLoop parameters
    FComputePacket f_packet;

    unsigned K, L, psi_size, received_psis, qindex;
    complex*** psis;
    complex* fs;
    int wrote;
    complex* umklapp_factor;
    
    double total_time;
};

extern /* readonly */ CProxy_Controller controller_proxy;
extern /* readonly */ CProxy_PsiCache psi_cache_proxy;

#endif
