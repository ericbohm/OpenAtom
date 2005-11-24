//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
/** \file util.h
 *
 */
//===================================================================================

#ifndef __PFFTUTIL_H__
#define __PFFTUTIL_H__
//===================================================================================

#include "converse.h"

//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
class CPcharmParaInfo; //forward decl
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <fftw.h>
#include <charm++.h>
#include <pairutil.h>

#define CAREFUL 1

#define ENERGY_EHART 0
#define ENERGY_ENL   1
#define ENERGY_EKE   2
#define ENERGY_EGGA  3
#define ENERGY_EEXC  4
#define ENERGY_EEXT  5
#define ENERGY_EWD   6
#define ENERGY_FMAG  7
#define ENERGY_FICTEKE 8

#define NUM_ENERGIES 9

#define MAX_CHAR_ARRAY_LENGTH 1024

//===================================================================================



//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
// Size or location in a regular 2D array
//===================================================================================
class size2d {
  int sz[2];
public:
  inline size2d() {sz[0]=sz[1]=0;}
  inline size2d(int ni,int nj) {sz[0]=ni; sz[1]=nj;}
  inline int ni(void) const {return sz[0];}
  inline int nj(void) const {return sz[1];}
  inline int nx(void) const {return sz[0];}
  inline int ny(void) const {return sz[1];}
  inline int operator[](int i) const {return sz[i];}
  inline int getIndex(int i,int j) const {
#if CAREFUL
    if (i<0 || i>=sz[0]) {CkAbort("Index out of bounds");}
    if (j<0 || j>=sz[1]) {CkAbort("Index out of bounds");}
#endif
    return i*sz[1] + j;
  }
  inline int getIndex(const size2d &at) const {
    return getIndex(at[0],at[1]);
  }
  inline int getTotal(void) const {return sz[0]*sz[1];}
  inline int getVolume(void) const {return getTotal();}
};
#if ! CMK_BLUEGENE_CHARM
PUPbytes(size2d);
#endif
//===================================================================================


//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
// Regular rectangular 2D array
//===================================================================================
template <class T>
class array2d {
  size2d size;

  //void free(void) {delete[] data; data=NULL;}

  inline void allocate(void) {
    //data=new T[size.getTotal()];
    data = (T*) malloc(size.getTotal() * sizeof(T)); 
  }
  
  inline void destroy(void) {
    //delete[] data; 

    free((void *) data);
    data=NULL;
  }


  inline array2d(const array2d &a);
  inline void operator=(const array2d &a);
  
public:

  T *data;  //probably not a good idea  JCHANGED
  inline array2d(void) {data=NULL;}
  inline array2d(const size2d &sz):size(sz){
    allocate();
    zero();
  }
  inline ~array2d() {destroy();}
  
  inline void zero(void) {
    int l=size.getTotal();
    //for (int i=0;i<l;i++) data[i]=T(0);

    memset(data, 0, l * sizeof(T));
  }

  inline const size2d &getSize(void) const {return size;}
  inline T *getData(void) const {return data;}

  inline T &operator[](const size2d &at) {return data[size.getIndex(at)];}
  inline const T &operator[](const size2d &at) const {return data[size.getIndex(at)];}

  inline T &operator()(int i,int j) {return data[size.getIndex(i,j)];}
  inline const T &operator()(int i,int j) const {return data[size.getIndex(i,j)];}
  
  void pup(PUP::er &p) {
	  p|size;
	  if (p.isUnpacking()) allocate();
	  p(data, size.getTotal());
  }

};
//===================================================================================

#include "../../include/RunDescriptor.h"

using namespace std;

//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
class Config {
 public:
	Config() { numSet = 0; numParams = 24;}
	int numParams;
	int numSet;
	int nstates;
	int gSpacePPC;
	int realSpacePPC;
	int rhoGPPC;
	int maxIter;
	int sGrainSize;
	int gSpaceNumChunks;
	int rhoGHelpers;
	int multicastDelayMS;
	int numMulticastMsgs;
	int numPartialReduction;
	int reductionDelayMS;
	int useCommlibMulticast;
	int useCommlib;
	int fftuseCommlib;
	int useGMulticast;
	int useGReduction;
	int checkForces;
	int atomIndex;
	int stateIndex;
	int planeIndex;
	int xIndex;
	int yIndex;
	int numFFTPoints;
        int numData;      // nPacked = number of non-zero g-space data pts
	int delayComputeZ;
	int pesPerState;
	int RpesPerState;
	int GpesPerState;
	int doublePack;
	int low_x_size;
	int high_x_size;
        int nchareG;
        int nchareRhoG;
	int inPlaceFFT;
	int conserveMemory;
	int prioFFTMsg; 
	int localSF;
	int delayCompStruct;
	int lbgspace;
	int lbpaircalc;
	int gspacesum;
        int numSfGrps;
        int numSfDups;
	int sfpriority;
	int rsfftpriority;
	int gsfftpriority;
	int rsifftpriority;
	int gsifftpriority;
	int lambdapriority;
	int psipriority;
	int rhorpriority;
	int rhogpriority;
	int priority;
	int fftprogresssplit;
        int stateOutputOn;
	int toleranceInterval;
        double gExpandFact;
        double gExpandFactRho;
	double displacement;
	char dataPath[MAX_CHAR_ARRAY_LENGTH];
	void print(char *fname_in);
	static void readConfig(const char *, Config &,int,int,int,int,int,int,int);
};
#if ! CMK_BLUEGENE_CHARM
PUPbytes(Config);
#endif
//===================================================================================



//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
// rho space run descriptors for spherical cutoff fft support
//===================================================================================
void make_rho_runs(CPcharmParaInfo *sim);

void get_rho_kvectors(double ecut4, double *hmati, int **kx_ret, int **ky_ret, 
                      int **kz_ret, int *nline_tot_ret, int *nPacked_ret,
                      int sizeX, int sizeY, int sizeZ);
//===================================================================================



//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
//   input/output routines, including states, run descriptors, and line balancing
//===================================================================================
void   readStateInfo(int &, int &, int &, int &, int &, int &,const char *, int );

void   readStateIntoRuns(int , complex *, CkVec<RunDescriptor> &, 
                         const char *,int ,int *,int *,int *,int *,int *,int *,
                         int **, int **, int);

void readState(int nPacked, complex *arrCP, const char *fromFile,
	       int ibinary_opt, int *nline_tot_ret,int *nplane_ret, int *kx, 
	       int *ky, int *kz, int *nx_ret, int *ny_ret, int *nz_ret,
               int *istrt_lgrp,int *iend_lgrp,int *npts_lgrp,int *nline_lgrp,
               int iget_decomp);

void create_line_decomp_descriptor(CPcharmParaInfo *sim);

void writePartState(int ,complex *,complex *,
                    int *,int *,int *,int ,int,int,int,char *,char *);
//===================================================================================

#endif //__PFFTUTIL_H__


