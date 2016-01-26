/* This routine calculates P matrix in R space first then perform FFT

   P[r,r'] = 4/(Eocc(k+q)-Eunocc(k)) * f[r] * f[r'].conj()

   (TODO: factor 4 has to be changed if nspin != 0.) 

   f[r] = psi_occ(k+q)[r] * psi_unocc(k)[r].conj()
 
           FFT       IFFT        
   P[r,r'] => P[G,r'] => P[G,G']  (here it performs columns first)
   no scaling factor is necessary because we did forward and backward consecutively

*/

#include <cstdlib>
#include <iostream>
#include <fstream>
#include "calc_P.h"
#include "util.h"

#include "include/mylapack.h"

using namespace std;

// calculate Polarizability matrix 
// psiR_k is unoccupied states
// psiR_kq is occupied states

void calc_Pmtrx_Rspace(CMATRIX *P, SYSINFO sys, STATES *psiR_k, STATES *psiR_kq, int uklpp[3], int nfft[3]){

    int ndata = nfft[0]*nfft[1]*nfft[2];
    // check if the numbers match
    if(ndata!=P->nrow){
	Die("number of FFT points do not match with number of row in P matrix");
    }
    if(ndata!=P->ncol){
	Die("number of FFT points do not match with number of column in P matrix");
    }

    // we need this since occupied state wavefunction may change due to U process
    complex *psiR_kq_tmp;
    psiR_kq_tmp = new complex [ndata];
        
    int nv = sys.nocc; // number of valence bands (=occupied states)
    int nc = sys.nunocc; // number of conduction bands (=unoccupied states)

    int vsq; // length square of uklpp vector
    
    // loop over occupied states (valence band)
    for (int iv=0; iv<nv; iv++) {
	 
	// psiR_kq_tmp is exactly the same as psiR_kq[iv] if uklpp is zero.
        for (int i=0; i<ndata; i++) { 
            psiR_kq_tmp[i] = psiR_kq->coeff[iv][i];
        }

	double Evkq = psiR_kq->eig[iv];

	// modify wavefunction if umklapp scattering applies 
        vsq = uklpp[0]*uklpp[0] + uklpp[1]*uklpp[1] + uklpp[2]*uklpp[2];
        if ( vsq != 0 ){
            modify_state_Uproc(psiR_kq_tmp, uklpp, nfft, sys);
        }

	// loop for unoccupied states 
        for (int ic=nv; ic<nv+nc; ic++) {

	    double Eck = psiR_k->eig[ic];

	    // update Pmtrx
            update_Pmtrx( psiR_kq_tmp, Evkq, psiR_k->coeff[ic], Eck, ndata, sys, P);
	}
    }
    delete [] psiR_kq_tmp;
}












// update_Pmtrx
// this function updates P matrix
void update_Pmtrx(complex *psi_vkq, double Evkq, complex *psi_ck, double Eck, int ndata, SYSINFO sys, CMATRIX *P){

    // TO DO: fact depends on nspin. check and modify
#ifdef BGW
    // if P(G,G') is compared with Berkeley GW...
    double fact = (double)(2)/ (double)(sys.nkpt) / sys.vol / (Evkq-Eck);
#else
    double fact = (double)(4)/( Evkq - Eck );
#endif

    complex* f = new complex[ndata];

    // calculate f
    for (int i=0; i<ndata; i++) {
	f[i] = psi_vkq[i] * psi_ck[i].conj();
     }
    

#ifdef USE_LAPACK

    char transformT='C'; // conjugate transpose (Hermitian conjugate)
    char transform='N';  // do nothing
    double Lalpha = fact; // scale A*B by this factor
    double Lbeta = double(1.0);  // scale initial value of C by this factor
    int one = 1;

    //GEMM   C:= alpha*A*B + beta*C
    myGEMM( &transform, &transformT, &ndata, &ndata, &one, &Lalpha, f, &ndata, f, &ndata, &Lbeta, P->m, &ndata);    

#else

    for (int i=0; i<ndata; i++) {
        for (int j=0; j<ndata; j++) {
            P->get(i,j) += fact * f[i] * f[j].conj();
        }
    }
#endif

  delete[] f;
}









// modify wavefunction when U process applies
void modify_state_Uproc(complex *a_v, int uklpp[3], int nFFT[3], SYSINFO sys){

    int ndata = nFFT[0]*nFFT[1]*nFFT[2];
    double rijk, G0, phase;
    complex fact;
    int icount = 0;
    double PI = M_PI;


    for (int i=0; i<nFFT[0]; i++){
        for (int j=0; j<nFFT[1]; j++){
            for (int k=0; k<nFFT[2]; k++){
  	        phase = 0;
	        fact = 0;
	        for (int ii=0; ii<3; ii++){
	            rijk = sys.a1[ii]*i/nFFT[0] + sys.a2[ii]*j/nFFT[1] + sys.a3[ii]*k/nFFT[2];
	            G0 = sys.b1[ii]*uklpp[0] + sys.b2[ii]*uklpp[1] + sys.b3[ii]*uklpp[2];
	            G0 *= -2*PI/sys.alat;
	            phase += rijk*G0;
	        }
	        fact.re = cos(phase);
	        fact.im = sin(phase);
	        a_v[icount] *= fact;

	        icount += 1;
            }
        }
    }

    // Print out error message, but don't die
    if ( ndata != icount ){
        cout << "something went wrong when Umklapp vector is applied." << endl;
    }
 
}  










// Performing FFTW
// P(r,r') -> P(g,g')
void Pmtrx_R_to_G(int iq, CMATRIX* P, int nfft[3], double vol){

    // number of FFT grid
    const int ndata = nfft[0]*nfft[1]*nfft[2];
    

    if (ndata!=P->nrow || ndata!=P->ncol){
	Die("number of data in P matrix does not match with total number of fftgrid");
    }
    
    
    int i, j;
    double real, imag;
    
    // loop over columns 
    setup_fftw_3d(nfft,-1); 

    for (int j=0; j<ndata; j++) {
        //fftin = 0;
        for (int i=0; i<ndata; i++) {
            in_pointer[i][0] = P->get(i,j).re;
            in_pointer[i][1] = P->get(i,j).im;
        }
        // Fourier transform
	do_fftw();
        
        for (int i=0; i<ndata; i++) {        
            P->get(i,j).re = out_pointer[i][0];
	    P->get(i,j).im = out_pointer[i][1];
        }
    }// end loop for

    
    // loop over rows
    setup_fftw_3d(nfft,1); 

    for (int i=0; i<ndata; i++) {
        
        for (int j=0; j<ndata; j++) {
            in_pointer[j][0] = P->get(i,j).re;
            in_pointer[j][1] = P->get(i,j).im;
        }
        // Fourier transform
	do_fftw();
        
        for (int j=0; j<ndata; j++) {        
            P->get(i,j).re = out_pointer[j][0];
	    P->get(i,j).im = out_pointer[j][1];
        }   
    
    }// end loop for rows

    // destroy FFTW stuffs
    //destroy_fftw_stuff();
   
}
