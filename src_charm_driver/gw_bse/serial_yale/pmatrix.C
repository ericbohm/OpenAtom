#include "standard_include_gwbse.h"
#include "allclass_gwbse.h"
#include "pmatrix.h"
#include "messages.h"
#include "controller.h"
#include "states.h"

PMatrix::PMatrix() {
  GWBSE* gwbse = GWBSE::get();
  pipeline_stages = gwbse->gw_parallel.pipeline_stages;
  L = gwbse->gw_parallel.L;
  M = gwbse->gw_parallel.M;
  num_rows = gwbse->gw_parallel.rows_per_chare;
  num_cols = gwbse->gw_parallel.n_elems;
  qindex = 1; // let's make set it to 1 for now, but I think it should be something like thisIndex.x (if P is 2D Chare array) 

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
  const unsigned istate = msg->state_index;
  const unsigned m = istate - L; // Eigenvalues indexed separately for occ/unocc

  // index for k+q point
  unsigned ikq;
  int umklapp[3]; // modify wavefunction for psi_occ(ikq) if umklapp process applies
  kqIndex(ikpt, ikq, umklapp);

  // Eigenvalues used to scale the entries of f
  double*** e_occ = gwbse->gw_epsilon.Eocc;
  double*** e_unocc = gwbse->gw_epsilon.Eunocc;

  // Loop over all of the cached psis, and compute an f via pointwise
  // multiplication with the received psi. Then compute the outer product of
  // f x f' and accumulate it's contribution in P.
  for (int l = 0; l < L; l++) {
    // Compute f based on each pair of Psis, and the two associated eigenvalues
    psi_occ = psi_cache_proxy.ckLocalBranch()->getPsi(0, ispin, ikq, l); // k+q index (ikq)
    for (int i = 0; i < size; i++) {
      f[i] = psi_occ[i]*psi_unocc[i].conj();
    }
    
    // Once we've computed f, compute its contribution to our chunk of P.
    double scaling_factor = 4/(e_occ[ispin][ikq][l] - e_unocc[ispin][ikpt][m]);
    for (int r = 0; r < num_rows; r++) {
      for (int c = 0; c < num_cols; c++) {
        data[r][c] += f[r+start_row]*f[c+start_col].conj() * scaling_factor;
      }
    }
  }

  // Tell the controller we've completed work on this psi
  contribute(CkCallback(CkReductionTarget(Controller, psiComplete), controller_proxy));
}

void PMatrix::printRowAndExit(int row) {
  if (row >= start_row && row < start_row + num_rows) {
    FILE* fp;
    char filename[200];
    sprintf(filename, "P_Rspace_row%d.dat", row);
    fp = fopen(filename, "w");
    for (int i = 0; i < num_cols; i++) {
      fprintf(fp, "row %d col %d %lg %lg\n", row, i, data[row-start_row][i].re, data[row-start_row][i].im);
    }
    fclose(fp);
  }
  contribute(CkCallback(CkCallback::ckExit));
}


void PMatrix::kqIndex(unsigned ikpt, unsigned ikq, int (&uklpp)[3]){
  GWBSE* gwbse = GWBSE::get();

  // temporary space to save k/q/k+q vectors
  double *this_k, *this_q;
  double k_plus_q[3], k_plus_q_orig[3];
  
  this_k = gwbse->gwbseopts.kvec[ikpt];
  this_q = gwbse->gwbseopts.qvec[qindex];

  for (int i=0; i<3; i++) {
    // calculate k+q vector 
    k_plus_q[i] = this_k[i] + this_q[i]; // k+q vector
    k_plus_q_orig[i] = k_plus_q[i]; // save it for Umklapp process
    // if not 0 =< k+q [i] <1, adjust k+q so that k+q[i] is in the Brillouine zone 
    if ( k_plus_q[i] >= 1 ) {
      k_plus_q[i] -= 1;
    }
    else if( k_plus_q[i] < 0 ){
      k_plus_q[i] += 1;
    }
  }
    
  // find k+q vector index
  for (int kk=0; kk < gwbse->gwbseopts.nkpt; kk++) {
    double diff = 0;
    //this_k is now a difference between k and k+q
    for (int i=0; i<3; i++) {
      this_k[i] = abs( gwbse->gwbseopts.kvec[kk][i] - k_plus_q[i] );
      diff += this_k[i];
    }
    if (diff == 0) {
      ikq = kk;  //k+q index is found
    }
  }
    // save umklapp scattering information
  for (int i=0; i<3; i++) {
    uklpp[i] = int( k_plus_q_orig[i] - k_plus_q[i] );
  }


}

#include "pmatrix.def.h"
