#include <cstdlib>
#include <iostream>
#include <math.h>
#include "include/ckcomplex.h"
#include "class_defs/laplace.h"
#include "class_defs/states.h"


//=================================================
//   LAPLACE class member function
//=================================================

// initialize max and min eigenvalues 
void LAPLACE::initialize(STATES* psi, int nocc, int nunocc){
  minEocc = psi->eig[0];
  maxEocc = psi->eig[nocc-1];
  minEunocc = psi->eig[nocc];
  maxEunocc = psi->eig[nocc+nunocc-1];
}

// set band gap and band width
void LAPLACE::set_GapBandwidth(STATES* psi, int nocc, int nunocc){

  if( psi->eig[0] < minEocc ){
    minEocc = psi->eig[0];
  }
  if( psi->eig[nocc-1] > maxEocc ){
    maxEocc = psi->eig[nocc-1];
  }
  if( psi->eig[nocc] < minEunocc ){
    minEunocc = psi->eig[nocc];
  }
  if( psi->eig[nocc+nunocc-1] > maxEunocc ){
    maxEunocc = psi->eig[nocc+nunocc-1];
  }

  gap = minEunocc - maxEocc;
  bandwidth = maxEunocc - minEocc;

}

// read Gauss-Laguerre nodes and weights
void LAPLACE::get_nodes_weights(){
  
  // open Gauss-Laguerre node and weight table file
  // I generated this file using matlab code (GLagIntP.m in working directory)
  FILE* fp = fopen("GLagNodeWeight.dat","r");
  int nptmax = 100;
  
  npt = new int [nptmax];
  n = new double *[nptmax];
  w = new double *[nptmax];
  for (int i=0; i<nptmax; i++){

    fscanf(fp,"%d\n",&npt[i]);

    n[i] = new double [i+1];
    w[i] = new double [i+1];
    
    for (int j=0; j<i+1; j++){
      fscanf(fp,"%lf %lf\n",&n[i][j],&w[i][j]);
    }
  }
  fclose(fp);
}


// set windows for given gap and width
// this is just based on 1/(Eunocc-Eocc) function
// I'm using cost estimation function to find windows quickly
// observation was nwocc = 1 is the lowest for any given nwunocc
void LAPLACE::set_windows(){
  
  int maxnw = 10;  // number of maximum windows. I chose to use just 10. This number can be changed
  // open file to read possible combinations for a given number of windows
  int** slice[maxnw];
  // array to save the total number of windowing combinations at each number of windows
  int totNumComb[maxnw];
  for (int nw=0; nw<maxnw; nw++){
    char filename[256];
    sprintf(filename,"windowing/%dw.dat",nw+1);
    FILE* fp = fopen(filename,"r");
    fscanf(fp,"%d\n",&totNumComb[nw]);
    slice[nw] = new int *[totNumComb[nw]];
    for(int i=0; i<totNumComb[nw]; i++){
      slice[nw][i] = new int [nw+1];
      if( nw == 0 ){
        fscanf(fp,"%d\n",&slice[nw][i][0]);
      }
      if( nw == 1 ){
        fscanf(fp,"%d %d\n",&slice[nw][i][0],&slice[nw][i][1]);
      }
      if( nw == 2 ){
        fscanf(fp,"%d %d %d\n",&slice[nw][i][0],&slice[nw][i][1],&slice[nw][i][2]);
      }
      if( nw == 3 ){
        fscanf(fp,"%d %d %d %d\n",&slice[nw][i][0],&slice[nw][i][1],&slice[nw][i][2],&slice[nw][i][3]);
      }
      if( nw == 4 ){
        fscanf(fp,"%d %d %d %d %d\n",&slice[nw][i][0],&slice[nw][i][1],&slice[nw][i][2],&slice[nw][i][3],&slice[nw][i][4]);
      }
      if( nw == 5 ){
        fscanf(fp,"%d %d %d %d %d %d\n",&slice[nw][i][0],&slice[nw][i][1],&slice[nw][i][2],&slice[nw][i][3],&slice[nw][i][4],&slice[nw][i][5]);
      }
      if( nw == 6 ){
        fscanf(fp,"%d %d %d %d %d %d %d\n",&slice[nw][i][0],&slice[nw][i][1],&slice[nw][i][2],&slice[nw][i][3],&slice[nw][i][4],&slice[nw][i][5],&slice[nw][i][6]);
      }
      if( nw == 7 ){
        fscanf(fp,"%d %d %d %d %d %d %d %d\n",&slice[nw][i][0],&slice[nw][i][1],&slice[nw][i][2],&slice[nw][i][3],&slice[nw][i][4],&slice[nw][i][5],&slice[nw][i][6],&slice[nw][i][7]);
      }
      if( nw == 8 ){
        fscanf(fp,"%d %d %d %d %d %d %d %d %d\n",&slice[nw][i][0],&slice[nw][i][1],&slice[nw][i][2],&slice[nw][i][3],&slice[nw][i][4],&slice[nw][i][5],&slice[nw][i][6],&slice[nw][i][7],&slice[nw][i][8]);
      }
      if( nw == 9 ){
        fscanf(fp,"%d %d %d %d %d %d %d %d %d %d\n",&slice[nw][i][0],&slice[nw][i][1],&slice[nw][i][2],&slice[nw][i][3],&slice[nw][i][4],&slice[nw][i][5],&slice[nw][i][6],&slice[nw][i][7],&slice[nw][i][8],&slice[nw][i][9]);
      }
    }
    fclose(fp);
  }


  // local variables
  double Evbw = maxEocc - minEocc;
  double Ecbw = maxEunocc - minEunocc;

  int nEv = usrin->nocc * usrin->nkpt;
  int nEc = usrin->nunocc * usrin->nkpt;

  double Evbin[maxnw+1];
  double Ecbin[maxnw+1];

  // set up bins
  for (int i=0; i<maxnw; i++){
    Evbin[i] = minEocc + i*Evbw/double(maxnw);
    Ecbin[i] = minEunocc + i*Ecbw/double(maxnw);
  }
  Evbin[maxnw] = maxEocc;
  Ecbin[maxnw] = maxEunocc;

  double minNumComp = pow(10,10); // minimum number of computation, and initialization
  int wIdxv; // number of windows that gives minimum number of computations (occupied states)
  int wIdxc; // number of windows that gives minimum number of computations (unoccupied states)

  // array that saves (ipv,ipc) pair (the lowest computation) at (nwv,nwc)
  int minIdxv[maxnw][maxnw]; 
  int minIdxc[maxnw][maxnw];

  /* debugging
  double MINCOST[maxnw][maxnw];
  */
  
  // loop over number of windows for occupied and unoccupied states
  for (int nwv=0; nwv<maxnw; nwv++){
    for (int nwc=0; nwc<maxnw; nwc++){

      // set anchor to make windows
      double anchor_v[nwv+2];
      double anchor_c[nwc+2];
      anchor_v[0] = Evbin[0];
      anchor_v[nwv+1] = Evbin[maxnw];
      anchor_c[0] = Ecbin[0];
      anchor_c[nwc+1] = Ecbin[maxnw];


      // search all possibilities at this nwv and nwc
      
      int ncomb = totNumComb[nwv]*totNumComb[nwc];
      double cost[ncomb]; // save cost
      int icounter = 0; // counter
      double mincost = pow(10,10); // initialize

 
      // for occupied states (valence bands)
      for (int ipv=0; ipv<totNumComb[nwv]; ipv++){
	
	double stepv = Evbw/double(maxnw);
	for (int i=0; i<nwv ; i++){
	  anchor_v[i+1] = anchor_v[i] + ( stepv * slice[nwv][ipv][i] );
	}
      
	// for unoccupied states (conduction bands)
	for(int ipc=0; ipc<totNumComb[nwc]; ipc++){

	  double stepc = Ecbw/double(maxnw);
	  for (int i=0; i<nwc; i++){
	    anchor_c[i+1] = anchor_c[i] + ( stepc * slice[nwc][ipc][i] );
	  }

	  // new loop starts to estimate the cost
	  cost[icounter] = 0.0; // cost for this particular windowing (ipv,ipc) at this (nwv,nwc)
	  for (int i=0; i<nwv+1; i++){
	    for (int j=0; j<nwc+1; j++){
	      double this_vbw = anchor_v[i+1] - anchor_v[i];
	      double this_cbw = anchor_c[j+1] - anchor_c[j];
	      double this_bw = anchor_c[j+1] - anchor_v[i];
	      double this_gap = anchor_c[j] - anchor_v[i+1];
	      cost[icounter] += sqrt(this_bw/this_gap)*(this_cbw/Ecbw*double(nEc) + this_vbw/Evbw*double(nEv));
	    }
	  }
	  
	  if ( cost[icounter] < mincost ){
	    mincost = cost[icounter];
	    // save (ipv,ipc) to find the minimum computation combination
	    minIdxv[nwv][nwc] = ipv;
	    minIdxc[nwv][nwc] = ipc;
	    /* debugging
	    MINCOST[nwv][nwc] = mincost;
            */
	  }
	  icounter += 1;

	}// end ipc loop
      }// end ipv loop


      if ( mincost < minNumComp ){
	minNumComp = mincost;
	wIdxv = nwv;
	wIdxc = nwc;
      }


    }// end nwc loop
  }// end nwv loop

  printf("\n Reporting windowing strategy:\n");
  printf("** number of windows for occupied states: %d, unoccupied states: %d\n", wIdxv+1, wIdxc+1);
  
  /*debugging  
  printf("Ev: %lg %lg   Ec: %lg %lg\n",minEocc,maxEocc,minEunocc,maxEunocc);
  for(int i=0;i<maxnw;i++){
    for(int j=0;j<maxnw; j++){
      printf("at this window (%d,%d), minimum cost is: %lg\n",i+1,j+1,MINCOST[i][j]);
    }
  }
  */
  
  // now, we estimated the number of windows, and how to slice them to get the lowest computation.
  // report these values, and finish this routine

  int idx1 = minIdxv[wIdxv][wIdxc];
  int idx2 = minIdxc[wIdxv][wIdxc];

  nwocc = wIdxv+1;
  nwunocc = wIdxc+1;
  slice_occ = new int [nwocc];
  slice_unocc = new int [nwunocc];
  for(int i=0; i<nwocc; i++){
    slice_occ[i] = slice[wIdxv][idx1][i];
  }
  for(int i=0; i<nwunocc; i++){
    slice_unocc[i] = slice[wIdxc][idx2][i];
  }
  // total windows
  totnw = nwocc * nwunocc;
  printf("** total number of windows: %d\n",totnw);

      
}


void LAPLACE::set_bin(){
  
  bin_occ = new double [nwocc+1];
  bin_unocc = new double [nwunocc+1];

  bin_occ[0] = minEocc;
  bin_occ[nwocc] = maxEocc;
  bin_unocc[0] = minEunocc;
  bin_unocc[nwunocc] = maxEunocc;

  int totslice=0;
  for (int i=0; i<nwunocc; i++){
    totslice += slice_unocc[i];
  }
  double step_occ = (maxEocc-minEocc)/double(totslice);
  double step_unocc = (maxEunocc-minEunocc)/double(totslice);

  // set-up the bin
  for (int i=1; i<nwocc; i++){
    bin_occ[i] = bin_occ[i-1] + slice_occ[i-1]*step_occ;
  }
  for (int i=1; i<nwunocc; i++){
    bin_unocc[i] = bin_unocc[i-1] + slice_unocc[i-1]*step_unocc;
  }
}
    
    

// now we set a number of nodes at each windows
// number of windows for occupied and unoccupied states are given
void LAPLACE::set_numnodes(){
  // allocate array for Nnodes and opta
  Nnodes = new int *[nwocc];
  opta = new double *[nwocc];
  for (int i=0; i<nwocc ; i++){
    Nnodes[i] = new int [nwunocc];
    opta[i] = new double [nwunocc];
  }

  int maxnnode = 100; // maximum number of nodes to be used

  for (int i=0; i<nwocc; i++){
    for (int j=0; j<nwunocc; j++){

      double G; // gap
      double W; // bandwidth
      G = bin_unocc[j] - bin_occ[i+1];
      W = bin_unocc[j+1] - bin_occ[i];

      // search for the number of nodes
      // while I was testing this method with matlab, I realized that the maximum error always occurs
      // when Ec-Ev = gap or Ec-Ev = width, so I'll just check these two values to expedite the process
      for (int inode=1; inode<maxnnode+1; inode++){
	
	// find optimized a
	double a;
	func_opta(G,W,inode,a);

	// calculate percentage error
	double exact;
	double estimated;

	// 1. Ec-Ev = gap
	exact = 1/G;
	estimated = 0;
	for (int k=0; k<inode; k++){
	  estimated += w[inode-1][k] * exp( -1.0 * n[inode-1][k] * (G/a-1) ) ;
	}
	estimated *= 1.0/a;
	double errorG = abs(exact - estimated) * G * 100;
 
	// 2. Ec-Ev = bandwidth
	exact = 1/W;
	estimated = 0;
	for (int k=0; k<inode; k++){
	  estimated += w[inode-1][k] * exp( -1.0 * n[inode-1][k] * (W/a-1) ) ;
	}
	estimated *= 1.0/a;
	double errorW = abs(exact - estimated) * W * 100;

	if ( errorG < ptol &&  errorW < ptol ){
	  Nnodes[i][j] = inode;
 	  opta[i][j] = a;
	  break;
	}
      } // search for Nnodes ends
      
    }
  }

}


// function that optimize a to minimize the error by a given number of npts (number of G-L nodes)
// original script by Sohrab (func_opta.m)
void LAPLACE::func_opta(double gap, double width, int npts, double &aopt){

  double G = gap; // min(Ec-Ev)
  double W = width; // max(Ec-Ev)
  int ng = npts; // number of Gauss-Laguerre quadrature node
  
  
  // list of u: logarithmic grid from 10^-2 to 10^2 with 10,000 points
  int nmax = 10000;
  double ulist[nmax];
  double step;
  step = 4/double((nmax-1));

  // make ulist
  for (int i=0; i<nmax; i++){
    double p = -2+step*i;
    ulist[i] = pow(10,p);
  }

  
  double exact = 1; // exact value of the integrand
  double errdiff[nmax];
  for (int i=0; i<nmax; i++){
    // 1) integrand of u (=E/a)
    double sum_u = 0;
    for (int ig=0; ig<ng; ig++){
      sum_u += w[ng-1][ig]*exp(-(ulist[i]-1)*n[ng-1][ig]);
    }
    sum_u *= ulist[i];
    // 2) integrand of u*W/G
    double sum_uWG = 0;
    for (int ig=0; ig<ng; ig++){
      sum_uWG += w[ng-1][ig]*exp(-(ulist[i]*W/G-1)*n[ng-1][ig]);
    }
    sum_uWG *= ulist[i]*W/G;
    
    errdiff[i] = abs( sum_u - sum_uWG );
  }

  // find minimum error index
  int idx = 0;
  double minerrdiff = 1;
  for (int i=0; i<nmax; i++){
    // search only if ulist < 1
    if( ulist[i] < 1 ){
      if( errdiff[i] < minerrdiff ){
	minerrdiff = errdiff[i];
	idx = i;
      }
    }
  }

  // return aopt value
  aopt = G/ulist[idx];
  
}



void LAPLACE::set_noWindow(){

  // reset windowing values
  nwocc = 1;
  nwunocc = 1;
  totnw = 1;

  slice_occ = new int [nwocc];
  slice_unocc = new int [nwunocc];
  slice_occ[0] = 10;
  slice_unocc[0] = 10;

  printf("\n WARNING: Windowing not used\n");

}  
  
