#include "HFCalculator.h"
#include "main/cpaimd.h"
#include "charm++.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <fftw3.h>

//============================================================================
extern CkVec <CProxy_CP_State_RealSpacePlane> UrealSpacePlaneProxy;
extern CkGroupID mCastGrpId;
extern Config config;
extern CPcharmParaInfo simReadOnly;
extern bool HartreeFockOn;
extern CProxy_HFCalculator HFCalculatorProxy;

HFCalculator::HFCalculator() {
  CkPrintf("Memory usage on node: <%d> : <%.2lf> MB before construction\n",
    CkMyNode(), CmiMemoryUsage()/(1024.0*1024.0));
  mygridx = simReadOnly.sizeX;
  mygridy = simReadOnly.sizeY;
  mygridz = simReadOnly.sizeZ;
  oneStateSize = mygridx * mygridy * mygridz;
  mynode = CkMyNode();
  CkPrintf("HFCalculator created on: <node: %d>\n", mynode);
  totalnodes = CkNumNodes();
  totalStates = simReadOnly.nstates;
  regularSize = oneStateSize / totalnodes;
  remSize = oneStateSize % totalnodes;
  extendedSize = regularSize + remSize;
  mynumrows = regularSize;
  myrowstart = mynode * regularSize;
  if ((mynode + 1) == totalnodes) {
    mynumrows = extendedSize;
  }
  oneState = new double[oneStateSize];
  myouterproduct = new complex*[mynumrows];
  transposedouterproduct = new complex*[mynumrows];
  for(int i = 0 ; i < mynumrows ; i++) {
    myouterproduct[i] = new complex[oneStateSize];
    transposedouterproduct[i] = new complex[oneStateSize];
  }
#ifdef _HFCALCULATOR_VERBOSE_
  CkPrintf("totalnodes: <%d>, myrowstart: <%d>, mynumrows: <%d>\n",
    totalnodes, myrowstart, mynumrows);
  CkPrintf("Memory usage on node: <%d> : <%.2lf> MB after construction\n",
    CkMyNode(), CmiMemoryUsage()/(1024.0*1024.0));
  CkPrintf("HFCalculator: totalstates, gridxyz, oneStateSize (in MB): <%d,%d,%d,%d,%.2f>\n",
    totalStates, mygridx, mygridy, mygridz, (double) oneStateSize * sizeof(double)/(1024.0*1024.0));
#endif
}

void HFCalculator::doRowFFTs(complex **mymatrix, int mydir) {
  int N = oneStateSize;
  fftw_complex *in, *out;
  fftw_plan p;
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  if (mydir == -1) {
    p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  }
  else if (mydir == 1) {
    p = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
  }
  else {
    CkAbort("Invalid value of direction for doing ffts for Hartree-Fock phase 2.\n");
  }
  for (int nrow = 0 ; nrow < mynumrows ; nrow++) {
    for(int i = 0 ; i < N ; i++) {
      in[i][0] = mymatrix[nrow][i].re;
      in[i][1] = mymatrix[nrow][i].im;
    }
    fftw_execute(p);
    for(int i = 0 ; i < N ; i++) {
      mymatrix[nrow][i].re = out[i][0];
      mymatrix[nrow][i].im = out[i][1];
    }
  }
  fftw_destroy_plan(p);
  fftw_free(in);
  fftw_free(out);
}

void HFCalculator::fillOneState(HFInputMsg *msg) {
  int chunksize = mygridx * mygridy;
  CkAssert(chunksize == msg->inputSize);
  int offset = chunksize * (msg->chunknum);
  CmiMemcpy(&oneState[offset], msg->inputPsi, chunksize * sizeof(double));
  delete msg;
}

void HFCalculator::initializeOuterproduct() {
  for(int i = 0 ; i < mynumrows ; i++) {
    for(int j = 0 ; j < oneStateSize ; j++) {
      myouterproduct[i][j].re = 0.0;
      myouterproduct[i][j].im = 0.0;
      transposedouterproduct[i][j].re = 0.0;
      transposedouterproduct[i][j].im = 0.0;
    }
  }
}

void HFCalculator::aggregateOuterproduct() {
  for(int i = 0 ; i < mynumrows ; i++) {
    double firstTerm = oneState[i];
    for(int j = 0 ; j < oneStateSize ; j++) {
      myouterproduct[i][j].re += firstTerm * oneState[j];
    }
  }
}

#include "HartreeFock.def.h"
