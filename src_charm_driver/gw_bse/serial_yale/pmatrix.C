#include "standard_include_gwbse.h"
#include "allclass_gwbse.h"
#include "pmatrix.h"
#include "messages.h"
#include "controller.h"
#include "states.h"
#include "my_fftw.h"
#include "fft_routines.h"

PMatrix::PMatrix() {
  GWBSE* gwbse = GWBSE::get();
  L = gwbse->gw_parallel.L;
  num_rows = gwbse->gw_parallel.rows_per_chare;
  num_cols = gwbse->gw_parallel.n_elems;
  nfft = gwbse->gw_parallel.fft_nelems;
  qindex = 0; // Eventually the controller will set this

  start_row = thisIndex * num_rows;
  start_col = 0;
  
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
  /*if (thisIndex == 0) {
    CkPrintf("Incoming kpt: %d, my qindex: %d, k+q: %d\n", ikpt, qindex, ikq);
  }*/
  for (int l = 0; l < L; l++) {
    // Compute f based on each pair of Psis, and the two associated eigenvalues
    psi_occ = psi_cache_proxy.ckLocalBranch()->getPsi(ispin, ikq, l);
    for (int i = 0; i < size; i++) {
      f[i] = psi_occ[i]*psi_unocc[i].conj();
    }
    /*if (thisIndex == 0) {
      CkPrintf("f[%i]: %lg %lg\n", 0, f[0].re, f[0].im);
    }*/
    
    // Once we've computed f, compute its contribution to our chunk of P.
    double scaling_factor = 4/(e_occ[ispin][ikq][l] - e_unocc[ispin][ikpt][m]);
    for (int r = 0; r < num_rows; r++) {
      for (int c = 0; c < num_cols; c++) {
        data[r][c] += f[r+start_row]*f[c+start_col].conj() * scaling_factor;
      }
    }
    /*if (thisIndex == 0) {
      CkPrintf("e_occ: %lg\n", e_occ[ispin][ikq][l]);
      CkPrintf("e_unocc: %lg\n", e_unocc[ispin][ikpt][m]);
      CkPrintf("scaling_factor: %lg\n", scaling_factor);
      CkPrintf("P[0,0]: %lg %lg\n", data[0][0].re, data[0][0].im);
    }*/
  }

  // Tell the controller we've completed work on this psi
  contribute(CkCallback(CkReductionTarget(Controller, psiComplete), controller_proxy));
}

void PMatrix::fftRows() {
  // TODO: Minjung will add the serial code here to fft each row stored in this chare
  // NOTE: The rows are stored in the data array
  // The first row in this chare is start_row, and num_rows is the number of rows in this chare
  // Each row will need to be fft'd.
  int direction = 1; // for rows, for column, direction should be -1
  for (int i=0; i < num_rows; i++){
    setup_fftw_3d(nfft,direction);
    put_into_fftbox(nfft, data[i], in_pointer);
    do_fftw();
    fftbox_to_array(num_cols, out_pointer, data[i], 1);
  }
  // Minjung: Need to call destroy_fftw_stuff here?
}

void PMatrix::doTranspose() {
}

void PMatrix::printRowAndExit(int row) {
  if (row >= start_row && row < start_row + num_rows) {
    FILE* fp;
    char filename[200];
    sprintf(filename, "P_Rspace_q%d_row%d.dat", qindex, row);
    fp = fopen(filename, "w");
    for (int i = 0; i < num_cols; i++) {
      fprintf(fp, "row %d col %d %lg %lg\n", row, i, data[row-start_row][i].re, data[row-start_row][i].im);
    }
    fclose(fp);
  }
  contribute(CkCallback(CkCallback::ckExit));
}


void PMatrix::kqIndex(unsigned ikpt, unsigned& ikq, int* uklapp){
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
    bool match = true;
    this_k = gwbse->gwbseopts.kvec[kk];
    //this_k is now a difference between k and k+q
    for (int i=0; i<3; i++) {
      if (this_k[i] != k_plus_q[i]) {
        match = false;
        break;
      }
    }
    if (match) {
      ikq = kk;
      break;
    }
  }
    // save umklapp scattering information
  for (int i=0; i<3; i++) {
    uklapp[i] = int( k_plus_q_orig[i] - k_plus_q[i] );
  }

}

#include "pmatrix.def.h"
