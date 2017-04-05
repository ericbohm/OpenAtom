/*
        main program for epsilon calculations
*/

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include "epsilon.h"
#include "constant.h"
#include "sigma.h"


int main(){

  // read system and cell information
  SYSINFO sys;
  sys.construct();

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

  // delete Gspace psi
  for (int is=0; is<nspin; is++){
    for (int ik=0; ik<nkpt; ik++){
      delete[] psi[is][ik];
    }
    delete[] psi[is];
  }
#ifdef VERBOSE
  mymessage("(unshifted) states have been read and FFTed");
#endif

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


  // initial setting if one chooses interpolation or laplace method for P
#ifdef USE_P_INTERPOLATION
  INTERPOLATOR ***G[nspin];
#elif USE_P_LAPLACE
  LAPLACE L = LAPLACE(usrin, sys, psiR);
#endif


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

#ifdef VERBOSE
    printf("  Inside P->EPS loop, spin at %d \n", is);
#endif
    
    //m'alloc psi_ and psiR_
    psi_[is] = new STATES* [nkpt];
    psiR_[is] = new STATES* [nkpt];
    //m'alloc P and Epsmat
    P[is] = new CMATRIX* [nq];
    Epsmat[is] = new CMATRIX* [nq];

#ifdef USE_P_INTERPOLATION
    G[is] = new INTERPOLATOR** [nq];
#endif


#ifdef DEBUG
nq = 1;
#endif

    // loop over q 
    for (int iq=0; iq<nq; iq++){
#ifdef VERBOSE
      printf("    q point at %d\n", iq);
      TimeStamp();
#endif
      // m'alloc P matrix
      P[is][iq] = new CMATRIX (ndata, ndata);

#ifdef USE_P_INTERPOLATION
      G[is][iq] = new INTERPOLATOR* [nkpt];
#endif


      // Polarizability calculation requires shifted wavefunctions

      // when q=0 & G=0

      if (iq==0){

	// reading wavefunctions
	
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
#ifdef VERBOSE
        mymessage("(shifted) states has been read and FFTed");
#endif
        //delete psi_
        for (int ik=0; ik<nkpt; ik++){
          delete[] psi_[is][ik];
        }

	//========================
        // calculte P matrix
	//========================
        for (int ik=0; ik<nkpt; ik++){            

#ifdef VERBOSE
          clock_t startTime = clock();
#endif

          // index for k+q
          get_k_plus_q_index(ik, iq, ikq, sys, uklpp);


          // calculate P(r,r') matrix
#ifdef USE_P_INTERPOLATION

	  G[is][iq][ik] = new INTERPOLATOR(iq, usrin.nPitp, psiR[is][ik], psiR_[is][ikq], sys);
	  
          CalcPmtrxRspaceInterpolation(P[is][iq], G[is][iq][ik], psiR_[is][ikq], uklpp, nfft, sys);

	  // remove G from the memory
	  delete G[is][iq][ik];

#elif USE_P_LAPLACE

	  CalcPmtrxLaplace(P[is][iq], L, psiR_[is][ikq], psiR[is][ik], uklpp, nfft, sys);

#else	  
          calc_Pmtrx_Rspace(P[is][iq], sys, psiR[is][ik], psiR_[is][ikq], uklpp, nfft);
#endif

#ifdef VERBOSE
	  clock_t endTime = clock();
	  clock_t clockTicksTaken = endTime - startTime;
	  double timeInSeconds = clockTicksTaken / (double) CLOCKS_PER_SEC;

	  printf("k index %d Time for calculating partial P matrix = %lg sec\n",ik,timeInSeconds);
#endif
	  

        }//end ik loop

#ifdef USE_P_LAPLACE
	printf("\n number of computation for P calculation with Gauss-Laguerre quadrature: %d\n\n",L.ncounter);
#endif	

      }//end iq=0

      
      if(iq!=0){
	// calculate P matrix
	for (int ik=0; ik<nkpt; ik++){
	  // index for k+q
	  get_k_plus_q_index(ik, iq, ikq, sys, uklpp);
	  // calculate P(r,r') matrix

#ifdef USE_P_INTERPOLATION

	  G[is][iq][ik] = new INTERPOLATOR(iq, usrin.nPitp, psiR[is][ik], psiR[is][ikq], sys);

	  CalcPmtrxRspaceInterpolation(P[is][iq], G[is][iq][ik], psiR[is][ikq], uklpp, nfft, sys);

#elif USE_P_LAPLACE
	  
	  CalcPmtrxLaplace(P[is][iq], L, psiR[is][ikq], psiR[is][ik], uklpp, nfft, sys);
#else
	  calc_Pmtrx_Rspace(P[is][iq], sys, psiR[is][ik], psiR[is][ikq], uklpp, nfft);
#endif

	  

	}//end ik loop
      }
      // P(r,r';q) is done

#ifdef USE_LAPACK
      // we need to transpose P matrix as it is stored in column major since we've called LAPACK routine
      P[is][iq]->transpose();
#endif

      // FFT: P(r,r';q) -> P(g,g';q)      
      Pmtrx_R_to_G(iq, P[is][iq], nfft, sys.vol);
      
//@@@@@@@@@@ print out real values to evaluate interpolation scheme
#ifdef DEBUG
      /*
      printf("We are printing out P matrix in G\n");
      {
      GSPACE g;
      printf("fft grid: %d %d %d\n",nfft[0],nfft[1],nfft[2]);
      g.ng = nfft[0]*nfft[1]*nfft[2];
      printf("number of points: %d\n",g.ng);
      g.ig = new int[g.ng];
      g.jg = new int[g.ng];
      g.kg = new int[g.ng];
      printf("allocation done\n");
      fftidx_to_gidx( g.ig, g.jg, g.kg, nfft );
      printf("Got G index!!!\n");
      */
      
      printf("\n Printing out P matrix in G\n");
      
      char fname[100];
#ifdef USE_P_LAPLACE
      sprintf(fname,"P_q%d_GL",iq);
      
#else
      sprintf(fname,"P_q%d",iq);
#endif
      /*
      int nr = nfft[0]*nfft[1]*nfft[2];
      for(int i=0; i<nr; i++){
	sprintf(fname,"P_row_%d",i);
	P[is][iq]->printRow(i,fname);
      }
      */
      //P[is][iq]->printRow(0,fname);
      //P[is][iq]->printallG(fname,&g);
      
#endif
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      
      
      //------------------------ DONE Polarizability calculations -------------------------//


#ifdef VERIFICATION
      if(iq==1) {
        char fname[100];
        sprintf(fname,"PMatrix_iq%d", iq);
        P[is][iq]->printAllRows(fname);
      }
#endif

#ifdef FULLRUN      
      // prep to calculate Epsmat
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

      if(iq==0){
	complex data = Epsmat[is][iq]->m[0];
	printf("\n Epsilon inverse at q=0 & G=G'=0:  %lg   %lg\n",data.re, data.im);
      }

      /*
      //@@@@@@@@@@@ let's print out values for Epsmat. for iq=1
      if(iq==1){
	for(int irow=0;irow<10;irow++){
	  char fname[100];
	  sprintf(fname,"EpsilonInverseMatrix_row%d",irow);
	  Epsmat[is][iq]->printRowEps(irow,fname, geps[iq]);
	}
      }
      //@@@@@@@@@@@ print out value done

      //@@@@@@@@@@ print out real values to evaluate interpolation scheme
#ifdef DEBUG
      printf("We are printing out epsilon inverse matrix in G space\n");
      {
      GSPACE *g;
      g = geps[iq];
      printf("Got G index!!!\n");
      
      char fname[100];
#ifdef USE_P_INTERPOLATION
      sprintf(fname,"EpsInv_q%d",iq);
#else
      sprintf(fname,"EpsInv_q%d_orig",iq);
#endif
      Epsmat[is][iq]->printallG(fname,g);
      }
#endif
      //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@      
      */
#endif
      
    }//end q loop
  }//end spin loop


  delete[] accept;


  mymessage("epsilon calculation is done");








#ifdef GPP

  /////////////////////////////////////////
  //          GPP calculations           //
  /////////////////////////////////////////

  // Set GPP options
  GPPOPTS gppopts = GPPOPTS(sys);

  // print out GPP options
#ifdef DEBUG
  gppopts.state_class_out();
#endif

  // Set GPP data
  GPPDATA **gppdata[nspin];

  for(int is=0; is<nspin; is++){

    gppdata[is] = new GPPDATA *[nq];
    
    for(int iq=0; iq<nq; iq++){
      // set-up coulomb potential
      VCOULB* vc = new VCOULB(iq, geps[iq], sys);

      // construct gppdata
      gppdata[is][iq] = new GPPDATA(vc, Epsmat[is][iq], gppopts);
      // Calculate S matrix and do eigendecomposition
      gppdata[is][iq]->compute();
      
    }// end q
  }// end spin

#endif // GPP calculation

  //************************************************************
  /* PART 2 
     SIGMA calculation
     No matrix method (see sigma_strategy.pdf) is implemented
     TODO : Maybe later, Sigma calculation will be combined with GPP calculations...?
  *///**********************************************************
  
  /////////////////////////////////////////
  //         Screened Exchange           //
  /////////////////////////////////////////
/*
  // SIGMA options
  SIGMAOPTS sigmaopts;
  // class for sigma dynamic calculations
  SIGMADYN **sigmadyn[nspin];
  
  // loop over spin
  for(int is=0; is<nspin; is++){
    // loop over q
    for(int iq=0; iq<nq; iq++){

      
      
    }
  }
*/



  TimeStamp();
  return 0;

}//end of main

