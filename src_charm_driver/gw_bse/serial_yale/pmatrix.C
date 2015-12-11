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
  qindex = 3; // Eventually the controller will set this

  num_chares = num_cols / num_rows;

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
  complex* f = new complex[size]; // Formed locally from the two psis

  // Variables for indexing into the eigenvalues arrays
  const unsigned ispin = msg->spin_index;
  const unsigned ikpt = msg->k_index;
  const unsigned istate = msg->state_index;
  const unsigned m = istate - L; // Eigenvalues indexed separately for occ/unocc

  // index for k+q point
  unsigned ikq;
  int umklapp[3]; // modify wavefunction for psi_occ(ikq) if umklapp process applies
  kqIndex(ikpt, ikq, umklapp); // TODO: Rather than compute each time, just make a table at startup

  // U-process modification
  bool Uproc = false;
  complex* umklapp_factor = new complex[size];
  // if umklapp is non-zero then it is U-process, not N-process, so Uproc=true
  if (umklapp[0] != 0 || umklapp[1] != 0 || umklapp[2] != 0) {
    Uproc = true;
    getUmklappFactor(umklapp_factor, umklapp);
  }

  // Eigenvalues used to scale the entries of f
  double*** e_occ = gwbse->gw_epsilon.Eocc;
  double*** e_unocc = gwbse->gw_epsilon.Eunocc;

  // Loop over all of the cached psis, and compute an f via pointwise
  // multiplication with the received psi. Then compute the outer product of
  // f x f' and accumulate it's contribution in P.
  for (int l = 0; l < L; l++) {
    // Compute f based on each pair of Psis, taking into account Uproc
    psi_occ = psi_cache_proxy.ckLocalBranch()->getPsi(ispin, ikq, l);
    for (int i = 0; i < size; i++) {
      f[i] = psi_occ[i]*psi_unocc[i].conj();
      // If we are doing U-process, then multiply by our umklapp factor
      if (Uproc) {
        f[i] *= umklapp_factor[i];
      }
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

  delete[] f;
  delete[] umklapp_factor;
}

void PMatrix::fftRows(int direction) {
  // FFT each row stored in this chare
  // TODO: The data for the FFTs shouldn't be stored in global variables
  for (int i=0; i < num_rows; i++){
    setup_fftw_3d(nfft, direction);
    put_into_fftbox(nfft, data[i], in_pointer);
    do_fftw();
    fftbox_to_array(num_cols, out_pointer, data[i], 1);
  }

  // Let the controller know we have completed the fft
  contribute(CkCallback(CkReductionTarget(Controller, fftComplete), controller_proxy));
}

// Print first n rows to file named with prefix
void PMatrix::printRows(int n, const char* prefix) {
  for (int r = start_row; r < n && r < start_row + num_rows; r++) {
    FILE* fp;
    char filename[200];
    sprintf(filename, "row_data/%s_q%d_row%d.dat", prefix, qindex, r);
    fp = fopen(filename, "w");
    for (int i = 0; i < num_cols; i++) {
      fprintf(fp, "row %d col %d %lg %lg\n", r, i, data[r-start_row][i].re, data[r-start_row][i].im);
    }
    fclose(fp);
  }
  contribute(CkCallback(CkReductionTarget(Controller, printingComplete), controller_proxy));
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


void PMatrix::getUmklappFactor(complex* umklapp_factor, int uklpp[3]){

  if (uklpp[0]==0 && uklpp[1]==0 && uklpp[2]==0){
    // do nothing
  }
  else{
    GWBSE *gwbse = GWBSE::get();
    int* nfft;
    nfft = gwbse->gw_parallel.fft_nelems;
    double *a1, *a2, *a3, *b1, *b2, *b3;
    a1 = gwbse->gwbseopts.a1;
    a2 = gwbse->gwbseopts.a2;
    a3 = gwbse->gwbseopts.a3;
    b1 = gwbse->gwbseopts.b1;
    b2 = gwbse->gwbseopts.b2;
    b3 = gwbse->gwbseopts.b3;
    double lattconst = gwbse->gwbseopts.latt;

    double rijk, G0, phase;
    unsigned counter = 0;
    for(int i=0; i<nfft[0]; i++){
      for(int j=0; j<nfft[1]; j++){
        for(int k=0; k<nfft[2]; k++){
          phase = 0;
          for (int l=0; l<3; l++){
            rijk = a1[l]*i/nfft[0] + a2[l]*j/nfft[1] + a3[l]*k/nfft[2];
            G0 = b1[l]*uklpp[0] + b2[l]*uklpp[1] + b3[l]*uklpp[2];
            G0 *= -2*M_PI/lattconst;
            phase += rijk*G0;
          }
          umklapp_factor[counter].re = cos(phase);
          umklapp_factor[counter].im = sin(phase);
          counter += 1;
        }// end k loop
      }// end j loop
    }// end i loop
  }//end if-else statement

}//end function

#include "pmatrix.def.h"
