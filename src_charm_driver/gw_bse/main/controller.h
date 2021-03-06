#ifndef __CONTROLLER_H__
#define __CONTROLLER_H__

#include <cstdlib>
#include "ckcomplex.h"

#include "fft_controller.decl.h"
#include "psi_cache.decl.h"
#include "fvector_cache.decl.h"

#include "controller.decl.h"
#define ITERATION 10 //needs to be read from epsilon.in
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
	double total_phase4;
};

class Controller : public CBase_Controller {
  Controller_SDAG_CODE
  public:
    Controller();
    void prep();
    void calc_Geps();
    void got_geps(std::vector<int> accept, int epsilon_size);
    void got_vcoulb(std::vector<double> vcoulb_in);
  private:
    bool do_output;
    int msg_received;
    int iter, maxiter;
    int iteration;
    unsigned index;
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
    int padded_epsilon_size;
    double prev_max;
    std::vector<int> accept_result;
    Timers timers;

    // Epsilon proxies and matrices
    CLA_Matrix_interface matA, matB, matC;
    CLA_Matrix_interface matA2, matB2, matC2;
    CLA_Matrix_interface matA3, matB3, matC3;
    CProxy_EpsMatrix eps_matrix1D_proxy;
    CProxy_EpsMatrix eps_matrix2D_proxy;
    CProxy_EpsMatrix eps_matrix2D_bproxy;
    CProxy_EpsMatrix eps_matrix2D_cproxy;
    CProxy_EpsMatrix eps_matrix2D_bbproxy;
    CProxy_EpsMatrix eps_matrix2D_ccproxy;
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
    complex* getF(unsigned,unsigned) const;
    void setVCoulb(std::vector<double> vcoulb_in);
    std::vector<double> getVCoulb();
    void setRegionData(int start_row, int start_col, int tile_nrows, int tile_ncols);
    void reportInfo();

  private:
    void kqIndex(unsigned, unsigned&, int*);
    void computeUmklappFactor(int*);

    // Used for CkLoop parameters
    FComputePacket f_packet;

    unsigned K, L, psi_size, received_psis, qindex, pipeline_stages, received_chunks;
    // TODO: Flatten arrays?
    complex*** psis;
    complex*** psis_shifted;
    complex* fs;
    std::vector<double> vcoulb;
    complex* umklapp_factor;
    
    double total_time;

    // Used for registering fvector regions
    int min_row, min_col, max_row, max_col;
    CmiNodeLock tile_lock;
};

class FVectorCache : public CBase_FVectorCache {
  FVectorCache_SDAG_CODE
  public:
    FVectorCache();
    void init(int size_xy);
    void putFVec(Phase4Message* msg);
    void computeFTilde(complex* fs_in);
    void putFVec(int n, complex* fs_input);
    complex* getFVec(int n, int l, int start, int size);
    void applyCutoff(int size, int* accept);
    void findIndices();
  private:
    FComputePacket f_packet;
    unsigned L, psi_size, fcount, n_list_size, node_count;
    int ndata, totalSize, data_size_x, data_size_y;
    int eps_chares_x, eps_chares_y, my_chare_count, my_chare_start;
    int eps_start_chare_x, eps_start_chare_y, eps_end_chare_x, eps_end_chare_y;
    int *my_eps_chare_indices_x;
    int *my_eps_chare_indices_y;
    int *data_offset_x;
    int *data_offset_y;
    int count;
    complex sum;
    complex* fs;
    int num_chares, num_chares_x, num_chares_y, chare_factor;
    int *charesX;
    int *charesY;
    int *local_offset;
    int *global_offset;
    int num_rows, num_cols;
    int *offsets;//[num_chares*2];//already calculated
};


extern /* readonly */ CProxy_Controller controller_proxy;
extern /* readonly */ CProxy_PsiCache psi_cache_proxy;
extern /* readonly */ CProxy_FVectorCache fvector_cache_proxy;

#endif
