/*
        main program for epsilon calculations

        1. read information from sysinfo.dat
        2. read user input file (epsilon.in)
        3. set fftsize for psi
        4. FFTW for all psi
        5. P matrix in R space is calculated
        6. P matrix in G space is calculated by performing dual FFTs
        7. Epsilon matrix is calculated, and the size of the matrix is about 1/4 of P matrix
        8. Epsilon inverse matrix is calculated by using iterative method

*/

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include "epsilon.h"
#include "constant.h"


int main(){

    // read system and cell information 
    SYSINFO sys = SYSINFO();

    // read user input values 
    USRINPUT usrin = USRINPUT();

    // check consistency between userinput and system input
    check_inputs(usrin, sys);

    // okay, it's not a good way, but I'm just assigning occ/unocc to sys here
    sys.nocc = usrin.nocc;
    sys.nunocc = usrin.nunocc;
 
    // determine the size of FFT box for wavefunction conversion
    // number of fft points in 3D direction
    int nfft[3];
    wfn_fftsize(usrin, sys, nfft);
    
    // define constant integers
    const int nspin = sys.nspin;
    const int nkpt = sys.nkpt;
    const int nstate = sys.nstate;
    const int nocc = usrin.nocc;
    const int nunocc = usrin.nunocc;

    // STATES are defined as a 2 dimensional array
    // that includes spin, and kpt index

    STATES** psi[nspin];  // in momentum space
    STATES** psiR[nspin]; // in real space

    // loop over spin
    for (int is = 0; is < nspin; is++){

	// flag for wavefunction shift
	// this flag is for Polarizability matrix calculations
	// wavefunction is not shifted
	bool shifted(false);

	psi[is] = new STATES *[nkpt];
	psiR[is] = new STATES *[nkpt];
	
	// loop over k points
	for (int ik = 0; ik < nkpt; ik++){

	    psi[is][ik] = new STATES (is, ik, nstate, shifted);
	    psiR[is][ik] = new STATES (is, ik, nstate, shifted);

	    read_states(psi[is][ik]);

	}// end ik loop

	// FFT: psi->psiR
	do_fft_states(psi[is], psiR[is], sys, nfft);

    }// end is loop

    // delete psi
    for (int is=0; is<nspin; is++){
	for (int ik=0; ik<nkpt; ik++){
	    delete[] psi[is][ik];
	}
	delete[] psi[is];
    }

    mymessage("(unshifted) states have been read and FFTed");
    
    /////////////////////////////////////////
    //                                     //
    // P matrix -> EPS INVERSE calculation //
    //                                     //
    /////////////////////////////////////////

    // number of q vectors = number of k points
    int nq = nkpt;
    
    // Polarizability variables
    int ikq; // index for k+q vector
    int uklpp[3]; // umklapp scattering
    CMATRIX **P[nspin]; // Polarizability matrix
    int ndata = nfft[0]*nfft[1]*nfft[2];// number of data for P matrix
    
    // Epsilon matrix variables
    CMATRIX **Epsmat[nspin];
    GSPACE *geps[nkpt];
    bool *accept = new bool [ndata];
    double *vcoulb;

    // additional states for shifted wavefunctions (psi_, psiR_)
    STATES** psi_[nspin];
    STATES** psiR_[nspin];

    // loop over spin
    for (int is=0; is<nspin; is++){

	//m'alloc psi_ and psiR_
	psi_[is] = new STATES* [nkpt];
	psiR_[is] = new STATES* [nkpt];
	//m'alloc P and Epsmat
	P[is] = new CMATRIX* [nq];
	Epsmat[is] = new CMATRIX* [nq];


#ifdef DEBUG
nq = 4;
#endif


        // loop over q 
	for (int iq=0; iq<nq; iq++){
	    // m'alloc P matrix
	    P[is][iq] = new CMATRIX (ndata, ndata);

	    // Polarizability calculation requires shifted wavefunctions
	    // when q=0 & G=0
/*
	    if (iq==0){
		bool shifted(true);
		psi_[is] = new STATES *[nkpt];
		psiR_[is] = new STATES *[nkpt];

		for(int ik=0; ik<nkpt; ik++){
		    psi_[is][ik] = new STATES (is, ik, nocc, shifted);
		    psiR_[is][ik] = new STATES (is, ik, nocc, shifted);
		    read_states(psi_[is][ik]);
		}
		
		// FFT: psi_ --> psiR_
		do_fft_states(psi_[is], psiR_[is], sys, nfft);
		mymessage("(shifted) states has been read and FFTed");

		//delete psi_
		for (int ik=0; ik<nkpt; ik++){
		    delete[] psi_[is][ik];
		}

		// calculte P matrix
		for (int ik=0; ik<nkpt; ik++){

		    // index for k+q
		    get_k_plus_q_index(ik, iq, ikq, sys, uklpp);

		    // calculate P(r,r') matrix
		    calc_Pmtrx_Rspace(P[is][iq], sys, psiR[is][ik], psiR_[is][ikq], uklpp, nfft);
		}//end k loop
	    }//end iq=0

	    if(iq!=0){
*/
		// calculate P matrix
		for (int ik=0; ik<nkpt; ik++){

		    // index for k+q
		    get_k_plus_q_index(ik, iq, ikq, sys, uklpp);

		    // calculate P(r,r') matrix
		    calc_Pmtrx_Rspace(P[is][iq], sys, psiR[is][ik], psiR[is][ikq], uklpp, nfft);
		}
/*
	    }
*/
	    // P(r,r';q) is done

#ifdef USE_LAPACK
	    // we need to transpose P matrix as it is stored in column major since we've called LAPACK routine
	    P[is][iq]->transpose();
#endif


	    // FFT: P(r,r';q) -> P(g,g';q)
	    Pmtrx_R_to_G(iq, P[is][iq], nfft, sys.vol);
	    //------------------------ DONE Polarizability calculations -------------------------//

	    

	    // prep to calculate Epsmat
	    
	    // later on, we are going to do it only once because all q points share the same g space
	    if(is==0){
		geps[iq] = new GSPACE;
		get_geps(geps[iq], usrin, sys, iq, nfft, accept);
	    }

	    Epsmat[is][iq] = new CMATRIX(geps[iq]->ng, geps[iq]->ng);

	    // get shift vectors
	    if(iq==0){
		for(int ii=0; ii<3; ii++){sys.shift[ii] = usrin.shift[ii];}
	    }
	    
	    // calculate epsilon matrix
	    calc_Epsmat(vcoulb, iq, sys, geps[iq], Epsmat[is][iq], P[is][iq], ndata, accept);	    

	    //delete P matrix
	    delete P[is][iq];

	    // calculate epsilon inverse matrix
	    iter_invmtrx(Epsmat[is][iq], usrin, geps[iq]->ng);



	}//end q loop
    }//end spin loop

    delete[] accept;


    TimeStamp();
    return 0;

}//end of main

