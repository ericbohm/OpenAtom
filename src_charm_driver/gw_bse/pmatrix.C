#include "pmatrix.h"
#include "psi.h"
#include "gw_bse.h"

PMatrix::PMatrix() {
  num_rows = config.rows_per_chare;
  num_cols = config.n_elems;
  start_row = thisIndex * num_rows;
  start_col = 0;
  done_count = 0;

  data = new double*[num_rows];
  for (int i = 0; i < num_rows; i++) {
    data[i] = new double[num_cols];
  }
}

void PMatrix::receivePsi(PsiMessage* msg) {
  // Compute all f's associated with the received psi, and accumulate their
  // contributions to P.
  CkPrintf("[%i]: Received psi [%i,%i]\n", thisIndex, msg->k_index, msg->state_index);
  const unsigned size = msg->size;
  double* psi1 = msg->psi;
  double* psi2;
  double f[size];

  // Loop over all of the cached psis, and compute an f via pointwise
  // multiplication with the received psi. Then compute the outer product of
  // f x f' and accumulate it's contribution in P.
  for (int l = 0; l < config.L; l++) {
    // Compute f based on each pair of Psis
    // TODO: Figure out the most cache effective way to do this. Should we
    // explicitly compute f or use the psis to compute the addition to P.
    // TODO: Instead of computing f on the P chares should the cache do it?
    psi2 = psicache.ckLocalBranch()->getPsi(l);
    for (int i = 0; i < size; i++) {
      f[i] = psi1[i]*psi2[i];
    }
    
    // Once we've computed f, compute its contribution to our chunk of P.
    for (int r = 0; r < num_rows; r++) {
      for (int c = 0; c < num_cols; c++) {
        data[r][c] += f[r+start_row]*f[c+start_col];
      }
    }
  }

  // Once all P chares have finished with this psi, the next psi can broadcast.
  if (msg->state_index + config.pipeline_stages < config.L + config.M) {
    contribute(CkCallback(CkReductionTarget(Psi, sendToP), psi(msg->k_index, msg->state_index + config.pipeline_stages)));
  }
  if (++done_count == config.M) {
    contribute(CkCallback(CkCallback::ckExit));
  }
}
