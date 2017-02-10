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

//===================================
//  calc_Pmtrx_Rspace
//===================================
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



//===================================
//  update_Pmtrx
//===================================
// this function updates P matrix
void update_Pmtrx(complex *psi_vkq, double Evkq, complex *psi_ck, double Eck, int ndata, SYSINFO sys, CMATRIX *P){

    // TODO  fact depends on nspin. check and modify
#ifdef BGW
    // if P(G,G') is compared with Berkeley GW...
    double fact = (double)(2)/ (double)(sys.nkpt) / sys.vol / (Evkq-Eck);
#else
    double fact = (double)(4)/( Evkq - Eck );
#endif

    complex f[ndata];

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
}


//===================================
//  modify_state_Uproc
//===================================
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


#ifdef USE_P_INTERPOLATION
//===================================
//  CalcPmtrxRspaceInterpolation
//===================================
#include "class_defs/interpolator.h"
void CalcPmtrxRspaceInterpolation(CMATRIX* P, INTERPOLATOR* G, STATES *psiR_occ, int uklpp[3], int nfft[3], SYSINFO sys){

  int ndata = psiR_occ->ndata;
  int nocc = sys.nocc;
  // TODO  fact depends on nspin. check and modify
#ifdef BGW
  // if P(G,G') is compared with Berkeley GW...
  double fact = (double)(-2)/ (double)(sys.nkpt) / sys.vol ;
#else
  double fact = (double)(-4) ;
#endif
  if (ndata != nfft[0]*nfft[1]*nfft[2]){
    Die("In CalcPmtrxRspaceInterpolation, ndata does not match with nfft");
  }
  
  // modify EoccMin and EoccMax since psiR_occ,k+q has different energy values
  G->EoccMin = psiR_occ->eig[0];
  G->EoccMax = psiR_occ->eig[sys.nocc-1];
  G->make_zlist();
  int nzpt = G->nzpt;


  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  // adding new functionality here
  // allocate colG (dimension: nzpt x nr)
  G->colG = new complex* [nzpt];
  for(int iz=0; iz<nzpt; iz++){
    G->colG[iz] = new complex [G->nr];
  }

  // initialize
  for(int iz=0; iz<nzpt; iz++){
    for(int ir=0; ir<G->nr; ir++){
      G->colG[iz][ir] = 0;
    }
  }
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@




  
  // Okay, now I'm printing out everything
#ifdef DEBUG
  /*
  for(int r1=0; r1<ndata; r1++){
    int r2 = r1;
    G->G_eachband(r1,r2);
  }
  */
#endif





  
  // update P matrix column by column
  // outer loop for column
  
  for(int c=0; c<ndata; c++){

    // 1. calculate G(z) for this column
    for(int z=0; z<nzpt; z++){
      G->G_calculator(c, z);
    }
    
#ifdef DEBUG
    /*if ( c==0 || c==863 || c==1727 ){
      G->print_G( c );
    }*/
#endif

    // Loop over valence bands (occupied states)
    // 2. calculate psiR_occ(:) * G(Ev) * psiR_occ(c)
    for(int occ=0; occ<nocc; occ++){

      // a. interpolate G(Eocc)
      double Eocc = psiR_occ->eig[occ];
      complex colG_Eocc[ndata]; // temporary space to save G at energy Eocc
      G->LinearInterpolation(Eocc, colG_Eocc);
      
      // b. multiply psiR_occ(:)
      // b-1. Modify PsiR_occ when k+q is outside of 1st B.Z.
      complex a_v[ndata]; // temporary space to save psi->coeff at this occ state
      for(int i=0; i<ndata; i++){
	a_v[i] = psiR_occ->coeff[occ][i];
      }
      int vsq = uklpp[0]*uklpp[0] + uklpp[1]*uklpp[1] + uklpp[2]*uklpp[2];
      if ( vsq != 0 ){
        modify_state_Uproc(a_v, uklpp, nfft, sys);
      }
      // TODO: Check if I should call LAPACK here?
      for(int i=0; i<ndata; i++){
#ifdef USE_LAPACK
	// Since P will be transposed when LAPACK is on (and this is always going to be the case),
	// let's save it in row major
      	P->get(c,i) += fact * a_v[i] * colG_Eocc[i] * a_v[c].conj();
#else
	P->get(i,c) += fact * a_v[i] * colG_Eocc[i] * a_v[c].conj();
#endif
      }
    }// P(:,c)_k is done
  }// P_k calculation is done

}
#endif



//**********************************************************************************
// Laplace method starts here
//**********************************************************************************


#ifdef USE_P_LAPLACE
#include "class_defs/laplace.h"
//===================================
//  Laplace and Gauss-Lagerre quadrature
//===================================
void CalcPmtrxLaplace(CMATRIX* P, LAPLACE &L, STATES *psiR_kq, STATES *psiR_k, int uklpp[3], int nfft[3], SYSINFO sys){

  // local variables
  int nwocc = L.nwocc; // number of windows in occupied states
  int nwunocc = L.nwunocc; // number of windows in unoccupied states

  int nocc = sys.nocc;
  int nunocc = sys.nunocc;
  
  int ndata = nfft[0]*nfft[1]*nfft[2];

  int counter;


  // loop for window : occupied states
  for (int iwocc=0; iwocc<nwocc; iwocc++){

    // Occupied state energy in this window
    int this_nocc=0;
    for (int iv=0; iv<nocc; iv++){
      if (L.bin_occ[iwocc] <= psiR_kq->eig[iv]  &&  psiR_kq->eig[iv] < L.bin_occ[iwocc+1])
	this_nocc += 1;
    }
    double Eocc[this_nocc];
    int Eoccidx[this_nocc];
    counter = 0;
    for (int iv=0; iv<nocc; iv++){
      if (L.bin_occ[iwocc] <= psiR_kq->eig[iv]  &&  psiR_kq->eig[iv] < L.bin_occ[iwocc+1]){
	Eocc[counter] = psiR_kq->eig[iv];
	Eoccidx[counter] = iv;
	counter += 1;
      }
    }


    // debugging message
    /*
    printf("\n occupied states window %d , number of occupied states in this window: %d\n",iwocc,this_nocc);
    for (int i=0;i<this_nocc;i++){
      printf("This occupied state: %d  Eocc: %lg , Eoccidx: %d\n",i,Eocc[i],Eoccidx[i]);
    }
    */
    
    

    // loop for window : unoccupied states
    for (int iwunocc=0; iwunocc<nwunocc; iwunocc++){

      // Unoccupied state energy in this window
      int this_nunocc=0;
      for (int ic=0; ic<nunocc; ic++){
	if(L.bin_unocc[iwunocc] <= psiR_k->eig[nocc+ic]  &&  psiR_k->eig[nocc+ic] < L.bin_unocc[iwunocc+1])
	  this_nunocc += 1;
      }
      double Eunocc[this_nunocc];
      int Eunoccidx[this_nunocc];
      counter = 0;
      for (int ic=0; ic<nunocc; ic++){
	if(L.bin_unocc[iwunocc] <= psiR_k->eig[nocc+ic]  &&  psiR_k->eig[nocc+ic] < L.bin_unocc[iwunocc+1]){
	  Eunocc[counter] = psiR_k->eig[nocc+ic];
	  Eunoccidx[counter] = nocc+ic; 
	  counter += 1;
	}
      }

      // debugging message
      /*
      printf("\n unoccupied states window %d , number of unoccupied states in this window: %d\n",iwunocc, this_nunocc);
      for (int i=0;i<this_nunocc;i++){
	printf("This unoccupied state: %d  Eunocc: %lg , Eunoccidx: %d\n",i,Eunocc[i],Eunoccidx[i]);
      }
      */

      // number of nodes and optimized a at this window pair (iwocc,iwunocc)
      int Ng = L.Nnodes[iwocc][iwunocc];
      double opta = L.opta[iwocc][iwunocc];


      // counting number of computations
      L.ncounter += Ng * ( this_nocc + this_nunocc );

      
      //----------------------------------
      // Gaussian quadrature starts here
      //----------------------------------
      
      for (int n=0; n<Ng; n++){

	double Wn = L.w[Ng-1][n];
	double Xn = L.n[Ng-1][n];
	double factor =  (-4.0/opta/sys.nspin) * Wn * exp( -( (L.epsc-L.epsv)/opta -1.0 ) * Xn );
       

	// debugging
	//printf("optimal a: %lg\n",opta);
	//printf("ig: %d\n Xn: %lg Wn: %lg\n",n,Xn,Wn);
	//printf("factor: -4/a Wn*Exp(-(epsc-epsv)/a-1)*Xn: %lg \n",factor);

	// f_occ and f_unocc calculation
	CMATRIX *focc;
	CMATRIX *funocc;
	focc = new CMATRIX(ndata,ndata);
	funocc = new CMATRIX(ndata,ndata);
	
	// Occupied states
	// focc[c] = \sum_v \psi_{v,k+q}(r) x \psi_{v,k+q}(c)* x exp( -(eps^v - Ev,k+q)/a * Xn)
	// we can build a matrix that contains all valence band wavefunctions
	
        CMATRIX *psis;
	psis = new CMATRIX(this_nocc,ndata);
	
	for (int iv=0; iv<this_nocc; iv++){

	  // copy of the occupied state wavefunction as it may change due to U process
	  complex psi_occ[ndata];
	  
	  for (int i=0; i<ndata; i++){
	    psi_occ[i] = psiR_kq->coeff[Eoccidx[iv]][i];
	  }
	  
	  // modify wavefunction if umklapp scattering applies
	  int uklsq = uklpp[0]*uklpp[0] + uklpp[1]*uklpp[1] + uklpp[2]*uklpp[2];
	  if ( uklsq != 0 ){
		modify_state_Uproc(psi_occ, uklpp, nfft, sys);
	  }

	  
	  for (int i=0; i<ndata; i++){
	    psis->get(iv,i) = psi_occ[i] * sqrt( exp(-(L.epsv-Eocc[iv])/opta*Xn) );
	  }
	  
	  /*
	  // construct focc(r,r')
	  // focc(r,r') = focc(r,r') + psi_occ(r) * psi_occ(r').conj() * exp(-(epsv-Ev)/a*Xn)
#ifdef USE_LAPACK
	  char transformT='C'; // conjugate transpose (Hermitian conjugate)
	  char transform='N'; // do nothing
	  double Lalpha = exp( -(L.epsv-Eocc[iv])/opta * Xn ); // scale A*B by this factor
	  double Lbeta = double(1.0); // scale initial value of C by this factor
	  int one = 1;
	  // LAPACK GEMM C = alpha*A*B + beta*C
	  myGEMM( &transform, &transformT, &ndata, &ndata, &one, &Lalpha, psi_occ, &ndata, psi_occ, &ndata, &Lbeta, focc->m, &ndata);
#else
	  for (int row=0; row < ndata ; row++){
	    for (int col=0; col < ndata; col++){
	      focc->get(row,col) += psi_occ[row] * psi_occ[col].conj() * exp( -(L.epsv-Eocc[iv])/opta * Xn );
	    }
	  }
#endif
	  */
	} // end occupied states
	
#ifdef USE_LAPACK
	char transformT = 'C'; // conjugate transpose (Hermitian conjugate)
	char transform = 'N'; // do nothing
	double Lalpha = double(1.0); // scale A*B by this factor
	double Lbeta = double(1.0); // scale initial value of C by this factor
	// LAPACK GEMM: C = alpha*A*B + beta*C
	myGEMM( &transform, &transformT, &ndata, &ndata, &this_nocc, &Lalpha, psis->m, &ndata, psis->m, &ndata, &Lbeta, focc->m, &ndata);
#else
	Die("Without -DUSE_LAPACK flag, polarizability calculation does not work!");
#endif  
	delete psis;
	

	 

	// Unoccupied states
	psis = new CMATRIX(this_nunocc,ndata);
	
	for (int ic=0; ic<this_nunocc; ic++){

	  // copy unoccupied states into the temporary array
	  complex psi_unocc[ndata];
	  for (int i=0; i<ndata; i++){
	    psi_unocc[i] = psiR_k->coeff[Eunoccidx[ic]][i].conj();
	  }

	  for (int i=0; i<ndata; i++){
	    psis->get(ic,i) = psi_unocc[i] * sqrt( exp(-(Eunocc[ic]-L.epsc)/opta * Xn) );
	  }
	  
	  /*
	  // construct funocc(r,r')
	  // funocc(r,r') = funocc(r,r') + psi_unocc(r).conj() * psi_unocc(r') * exp( -(Ec-epsc)/a*Xn) 
#ifdef USE_LAPACK
	  char transformT='C'; // conjugate transpose (Hermitian conjugate)
	  char transform='N'; // do nothing
	  double Lalpha = exp( -(Eunocc[ic]-L.epsc)/opta * Xn ); // scale A*B by this factor
	  double Lbeta = double(1.0); // scale initial value of C by this factor
	  int one = 1;
	  // LAPACK GEMM C = alpha*A*B + beta*C
	  myGEMM( &transform, &transformT, &ndata, &ndata, &one, &Lalpha, psi_unocc, &ndata, psi_unocc, &ndata, &Lbeta, funocc->m, &ndata);
#else
	  for (int row=0; row < ndata ; row++){
	    for (int col=0; col < ndata; col++){
	      funocc->get(row,col) += psi_unocc[row] * psi_unocc[col].conj() * exp( -(Eunocc[ic]-L.epsc)/opta * Xn );
	    }
	  }
#endif	
	  */
	}// end unoccupied state

#ifdef USE_LAPACK
	myGEMM( &transform, &transformT, &ndata, &ndata, &this_nunocc, &Lalpha, psis->m, &ndata, psis->m, &ndata, &Lbeta, funocc->m, &ndata);
#else
	Die("Without -DUSE_LAPACK flag, polarizability calculation does not work!");
#endif  
	delete psis;
	

	
	// Update P
	for (int i=0; i<ndata*ndata; i++){
	  P->m[i] += factor * focc->m[i] * funocc->m[i];
	}
	//for (int r=0; r<ndata; r++){
	//  for (int c=0; c<ndata; c++){
	//    P->get(r,c) += factor * focc->get(r,c) * funocc->get(r,c);
	//  }
	//}
	
	delete focc;
	delete funocc;
	
      }// end Gaussian quadrature loop
      
    }// end iwunocc loop
  }// end iwocc loop
  
}



#endif



//===================================
//  Pmtrx_R_to_G
//===================================
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

