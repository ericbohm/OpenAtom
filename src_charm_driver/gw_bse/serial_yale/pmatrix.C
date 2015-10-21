#include "standard_include_gwbse.h"
#include "allclass_gwbse.h"
#include "pmatrix.h"
#include "messages.h"
#include "controller.h"
#include "states.h"

extern /* readonly */ CProxy_States states_proxy;
extern /* readonly */ CProxy_PsiCache psi_cache_proxy;

PMatrix::PMatrix() {
  GWBSE* gwbse = GWBSE::get();
  pipeline_stages = gwbse->gw_parallel.pipeline_stages;
  L = gwbse->gw_parallel.L;
  M = gwbse->gw_parallel.M;
  num_rows = gwbse->gw_parallel.rows_per_chare;
  num_cols = gwbse->gw_parallel.n_elems;

  start_row = thisIndex * num_rows;
  start_col = 0;
  done_count = 0;

  data = new complex*[num_rows];
  for (int i = 0; i < num_rows; i++) {
    data[i] = new complex[num_cols];
  }
}

void PMatrix::receivePsi(PsiMessage* msg) {
  GWBSE* gwbse = GWBSE::get();

  // Variables for f = psi(i) * psi(j)
  const unsigned size = msg->size;
  complex* psi_occ;               // Comes from the cache
  complex* psi_unocc = msg->psi;  // Sent directly to us
  complex f[size];                // Formed locally from the two psis

  // Variables for indexing into the eigenvalues arrays
  const unsigned ispin = msg->spin_index;
  const unsigned ikpt = msg->k_index;
  const unsigned m = msg->state_index - L; // Eigenvalues indexed separately for occ and unocc for now

  // Eigenvalues used to scale the entries of f
  double*** e_occ = gwbse->gw_epsilon.Eocc;
  double*** e_unocc = gwbse->gw_epsilon.Eunocc;

  // Loop over all of the cached psis, and compute an f via pointwise
  // multiplication with the received psi. Then compute the outer product of
  // f x f' and accumulate it's contribution in P.
  for (int l = 0; l < L; l++) {
    // Compute f based on each pair of Psis, and the two associated eigenvalues
    psi_occ = psi_cache_proxy.ckLocalBranch()->getPsi(l);
    double scaling_factor = sqrt(4/(e_occ[ispin][ikpt][l] - e_unocc[ispin][ikpt][m]));
    for (int i = 0; i < size; i++) {
      f[i] = psi_occ[i]*psi_unocc[i].conj() * scaling_factor;
    }
    
    // Once we've computed f, compute its contribution to our chunk of P.
    for (int r = 0; r < num_rows; r++) {
      for (int c = 0; c < num_cols; c++) {
        data[r][c] += f[r+start_row]*f[c+start_col].conj();
      }
    }
  }

  // Once all P chares have finished with this psi, the next psi can broadcast.
  if (msg->state_index + pipeline_stages < L + M) {
    contribute(CkCallback(CkReductionTarget(States, sendToP), states_proxy(msg->spin_index, msg->k_index, msg->state_index + pipeline_stages)));
  }
  if (++done_count == M) {
    contribute(CkCallback(CkCallback::ckExit));
  }
}

#include "pmatrix.def.h"
