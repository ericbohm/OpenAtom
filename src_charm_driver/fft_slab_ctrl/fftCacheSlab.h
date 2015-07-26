//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
/** \file fftCacheSlab.h
 * Add type declarations for simulationConstants class (readonly vars) 
 * and once class for each type of object we currently have
 *(GStateXYslab, RealStateXZslab, RealRhoSlab, GRhoSlab, etc.) 
 */
//==============================================================================
#include "utility/util.h"
#include "uber/Uber.h"

#ifndef _fftcacheslab_h_
#define _fftcacheslab_h_

#include "RunDescriptor.h"

#ifdef _CP_USE_BLAS_GATHER_
include "src_mathlib/mathlib.h"
#endif

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
class GStateSlab {
  //==============================================================================
  // data 
  public:

    void pup(PUP::er &);
    int cp_min_opt;
    size2d planeSize;
    int numNonZeroPlanes;
    int numRuns, numLines;      // numLines=numRun/2=number of lines in collection
    int numPoints;              // number of non-zero pts in the collection
    int numFull;                // expanded data : nfftz*numLines
    int numFullNL;              // expanded data : nfftz*numLines
    int ees_nonlocal;
    bool fftReqd;               // flags to indicate whether this data set of pencils
    int mysizeX;                // size of the ffT in X dimension 
    int S_grainSize;            // PC decomposition
    int xdim, ydim, zdim;       // FFT sizes
    int ngridaNL,ngridbNL,ngridcNL;       // FFT sizes for non-local
    int iplane_ind;             // collection index
    int istate_ind;             // state index of this collection
    int ihave_kx0;              // plane zero is in the collection (stored consecutively)
    int ihave_g000;
    int ind_g000;
    int kx0_strt;               // starting pt
    int kx0_end;                // ending pt
    int nkx0,nkx0_uni,nkx0_red; // split kx=0 into unique and redundant parts
    int nkx0_zero;              // ncoef_true=numPoints-nkx0_red

    double eke_ret;            // kinetic energy
    double fictEke_ret;        // fictitious kinetic energy
    double ekeNhc_ret;         // NHC energies
    double potNHC_ret;
    double degfree;            // Degrees of freedom (ncoef_true+num_nhc-1)
    double degfreeNHC;         // Degrees of freedom (num_nhc-1)*len_nhc

    complex *packedPlaneData;   // Non-zero data pts [numPoints]
    complex *packedPlaneDataTemp; 
    complex *packedPlaneDataTemp1; 
    complex *packedPlaneDataTemp2; 
    complex *packedPlaneDataTemp3; 
    complex *packedPlaneDataScr; 
    complex *packedForceData; 
    complex *packedVelData; 
    complex *packedRedPsi;
    complex *packedRedPsiV;
    int     len_nhc_cp;
    int     num_nhc_cp;
    int     nck_nhc_cp;
    int     *istrNHC;
    int     *iendNHC;
    double  kTCP;
    double  tauNHCCP;
    double  ***xNHC;
    double  ***xNHCP;
    double  ***vNHC;
    double  ***fNHC;
    double  *degFreeSplt;
    double  *mNHC;
    double  *v0NHC,*a2NHC,*a4NHC;

    //==============================================================================
    // Constuctor, Destructor and utilities

    GStateSlab() {packedPlaneData=NULL; packedPlaneDataTemp=NULL; 
      packedForceData=NULL; packedPlaneDataScr=NULL; 
      packedVelData=NULL;}
    ~GStateSlab();

    void initNHC(int _len_nhc_cp, int _num_nhc_cp, int _nck_nhc_cp){
      nck_nhc_cp = _nck_nhc_cp; 
      num_nhc_cp = _num_nhc_cp;
      len_nhc_cp = _len_nhc_cp;
      xNHC       = new double **[nck_nhc_cp];
      xNHCP      = new double **[nck_nhc_cp];
      vNHC       = new double **[nck_nhc_cp];
      fNHC       = new double **[nck_nhc_cp];
      for(int k = 0;k<nck_nhc_cp;k++){
        xNHC[k]   = new double *[num_nhc_cp];
        xNHCP[k]  = new double *[num_nhc_cp];
        vNHC[k]   = new double *[num_nhc_cp];
        fNHC[k]   = new double *[num_nhc_cp];
        for(int i=0;i<num_nhc_cp;i++){
          xNHC[k][i]   = new double[len_nhc_cp];
          xNHCP[k][i]  = new double[len_nhc_cp];
          vNHC[k][i]   = new double[len_nhc_cp];
          fNHC[k][i]   = new double[len_nhc_cp];
        }//endfor
      }//endfor
      degFreeSplt = new double[nck_nhc_cp];
      istrNHC     = new int [nck_nhc_cp];
      iendNHC     = new int [nck_nhc_cp];
      mNHC  = new double[len_nhc_cp];
      v0NHC = new double[num_nhc_cp];
      a2NHC = new double[num_nhc_cp];
      a4NHC = new double[num_nhc_cp];
      for(int k=0;k<nck_nhc_cp;k++){
        for(int i=0;i<num_nhc_cp;i++){
          for(int j=0;j<len_nhc_cp;j++){
            xNHC[k][i][j]      = 0.0;
            xNHCP[k][i][j]     = 0.0;
            vNHC[k][i][j]      = 0.0;
            fNHC[k][i][j]      = 0.0;
          }}}//endfor
    }//end routine

    void destroyNHC(){
      for(int k=0;k<nck_nhc_cp;k++){
        for(int i=0;i<num_nhc_cp;i++){
          delete []xNHC[k][i];
          delete []xNHCP[k][i];
          delete []vNHC[k][i];
          delete []fNHC[k][i];
        }//endfor
        delete []xNHC[k];
        delete []xNHCP[k];
        delete []vNHC[k];
        delete []fNHC[k];
      }//endfor
      delete []xNHC;
      delete []xNHCP;
      delete []vNHC;
      delete []fNHC;
      delete []degFreeSplt;
      delete []istrNHC;
      delete []iendNHC;
      delete []mNHC;
      delete []v0NHC;
      delete []a2NHC;
      delete []a4NHC;
    }//end routine

    void setKRange(int , int *, int *, int *);
    void addForces(complex *,const int *);

    //-----------------------------------------------------------------------------
};
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
class RealStateSlab {
  public:
    complex *planeArr;
    double *planeArrR;   //planeArr cast to double
    int ngrid_a;
    int ngrid_b;
    int thisState;    
    int thisPlane;      
    int numPlanesToExpect;
    int nsize;
    int rsize;
    int size;
    double e_gga; 
    RealStateSlab() {}
    ~RealStateSlab(); 
    void zeroOutPlanes();
    void allocate();
    void destroy();
    void pup(PUP::er &);

    //-----------------------------------------------------------------------------
};
//==============================================================================

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
//    Holder classes for the plans : Allows many fft libaries to be used
//==============================================================================
typedef struct essl_work {
  int num;
  double *work1, *work2;
} ESSL_WORK;

typedef struct fftplanholder {
  int nchare;
  int nsplit;
  int option;             // 0= fftw, 1=essl
  int nfft;               // essl stuff
  int ostride;            // essl stuff
  int odist;              // essl stuff
  int isign;              // essl stuff
  int nwork1,    nwork2;  // essl stuff
  double scale;           // essl stuff
  int nval,nmax;          // essl stuff
  int *mapp;              // essl stuff
  ESSL_WORK *essl_work;   // essl stuff
  double *work1, *work2;  // essl stuff
  fftw_plan fftwPlan;     // fftw stuff
} FFTplanHolder;

typedef struct rfftplanholder {
  int nchare;
  int nsplit;
  int option;             // 0= fftw, 1=essl
  int nfft;               // essl stuff
  int ostride;            // essl stuff
  int odist;              // essl stuff
  int isign;              // essl stuff
  int nwork1,    nwork2;  // essl stuff
  double scale;           // essl stuff
  int nval,nmax;          // essl stuff
  int *mapp;              // essl stuff
  ESSL_WORK *essl_work;   // essl stuff
  double *work1, *work2;  // essl stuff
  rfftwnd_plan rfftwPlan; // fftw stuff
} RFFTplanHolder;
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
class FFTcache: public Group {
  public:
    //-----------------------------------------------------------
    // FFT-Sizes and Ees method Options
    int ngrida;
    int ngridb;
    int ngridc;
    int ngridaEext;    
    int ngridbEext;
    int ngridcEext;
    int ngridaNL;
    int ngridbNL;
    int ngridcNL;
    int ees_eext_on;
    int ees_NL_on;
    int nchareGState;
    int nchareRState;
    int nchareGNL;
    int nchareRNL;
    int nchareGRho;
    int nchareRRho;
    int nchareRRhoTot;
    int nchareGEext;
    int nchareREext;
    int nchareREextTot;
    int nsplitR;
    int nsplitG;
    int rhoRsubPlanes;

    int cacheMemFlag;
    char cacheMemName[1000];
    static CmiNodeLock fftw_plan_lock;
    //-----------------------------------------------------------
    // Generic plane temporaries used everywhere possible to avoid CmiMemcpys
    complex *tmpData; 
    double  *tmpDataR;

    //-----------------------------------------------------------
    // Da Plans : Use the holder class to allow other fft libs
    //            All plans are double pack plans

    RFFTplanHolder fwdXPlanState,   bwdXPlanState;    // state
    FFTplanHolder  fwdXPlanStateK,  bwdXPlanStateK;  // state for kpoints
    FFTplanHolder  fwdYPlanState,   bwdYPlanState;
    FFTplanHolder  fwdYPlanStateLower,  bwdYPlanStateLower; // state for kpoints     
    FFTplanHolder  fwdZPlanState,   bwdZPlanState;

    RFFTplanHolder fwdXPlanRho,  bwdXPlanRho;      // density 
    FFTplanHolder  fwdYPlanRho,  bwdYPlanRho;     
    FFTplanHolder  fwdYPlanRhoS, bwdYPlanRhoS;     // special subplane plan
    FFTplanHolder  fwdZPlanRho,  bwdZPlanRho;

    FFTplanHolder  fwdZPlanRhoHart;

    RFFTplanHolder fwdXPlanNL,  bwdXPlanNL;         // ees NL 
    FFTplanHolder  fwdXPlanNLK, bwdXPlanNLK;       // ees NL for kpoints
    FFTplanHolder  fwdYPlanNL,  bwdYPlanNL;
    FFTplanHolder  fwdYPlanNLLower, bwdYPlanNLLower; // for kpoints
    FFTplanHolder  fwdZPlanNL,  bwdZPlanNL;    


    RFFTplanHolder fwdXPlanEext,  bwdXPlanEext;    // ees Eext
    FFTplanHolder  fwdYPlanEext,  bwdYPlanEext;
    FFTplanHolder  fwdYPlanEextS, bwdYPlanEextS;   // special subplane plan
    FFTplanHolder  fwdZPlanEext,  bwdZPlanEext; 
    const UberCollection thisInstance;

    //-----------------------------------------------------------
    // The constructor 
    FFTcache(     int _ngrida, int _ngridb, int _ngridc, 
        int _ngridaEext, int _ngridbEext, int _ngridcEext, 
        int _ees_eext_on, int _ngridaNL, int _ngridbNL, int _ngridcNL, 
        int _ees_NL_on, int _nlines_max, int _nlines_max_rho,
        int _nchareGState, int _nchareRState,
        int _nchareGNL,    int _nchareRNL, 
        int _nchareGRho,   int _nchareRRho,  int _nchareRRhoTot,
        int _nchareGEext,  int _nchareREext, int _nchareREextTot,
        int  *numGState,   int  *numRXState, int *numRYState,int *numRYStateLower,
        int  *numGNL,      int  *numRXNL,    int *numRYNL,int *numRYNLLower,
        int  *numGRho,     int  *numRXRho,   int *numRYRho ,
        int  *numGEext,    int  *numRXEext,  int *numRYEext ,
        int _fftopt,       int _nsplitR,     int _nsplitG,
        int _rhoRsubPlanes, UberCollection _thisInstance);
    //-----------------------------------------------------------
    // cache control 
    void getCacheMem(const char *name){
      if(cacheMemFlag==1){
        CkPrintf("%s stealing from %s\n",cacheMemName,name);
        CkExit();
      }//endif
      cacheMemFlag = 1;
      strcpy(cacheMemName,name);
    }//end routine
    void freeCacheMem(const char *name){
      if(cacheMemFlag==0){
        CkPrintf("Bad cache memory free from %s\n",name);
        CkExit();
      }//endif
      cacheMemFlag = 0;
    }//end routine

    //-----------------------------------------------------------
    // Generic G-space expanders and contractors
    void expandGSpace(complex* data, complex *packedData, 
        RunDescriptor *runs, int numRuns, int numFull,
        int numPoints, int nfftz);
    void packGSpace(complex* data, complex *packedPlaneData, 
        RunDescriptor *runs, int numRuns, int numFull,
        int numPoints, int nfftz);

    //-----------------------------------------------------------
    // Density FFTs
    void doHartFFTGtoR_Gchare(complex *,complex *,int , int ,int , int , 
        RunDescriptor *, int ,int);
    void doRhoFFTRtoG_Gchare(complex *,complex *,int ,int ,int ,int ,RunDescriptor *, 
        int ,int ,int);
    void doRhoFFTGtoR_Gchare(complex *,complex *,int ,int ,int ,int ,RunDescriptor *, 
        int ,int ,int);
    void doRhoFFTRtoG_Rchare(complex *,double *,int , int ,int ,int);
    void doRhoFFTGtoR_Rchare(complex *,double *,int , int ,int ,int);

    void doRhoFFTRxToGx_Rchare(complex *,double *,int , int ,int ,int);
    void doRhoFFTRyToGy_Rchare(complex *,double *,int , int ,int ,int);
    void doRhoFFTGxToRx_Rchare(complex *,double *,int , int ,int ,int);
    void doRhoFFTGyToRy_Rchare(complex *,double *,int , int ,int ,int);

    //-----------------------------------------------------------
    // State FFTs
    void doStpFFTRtoG_Gchare(complex *,complex *,int, int ,int ,int, RunDescriptor *,int,int);
    void doStpFFTGtoR_Gchare(complex *,complex *,int, int ,int ,int, RunDescriptor *,int,int);
    void doStpFFTGtoR_Rchare(complex *,double *,int , int ,int ,int);
    void doStpFFTRtoG_Rchare(complex *,double *,int , int ,int ,int);

    //-----------------------------------------------------------
    // non-local fft
    void doNlFFTRtoG_Gchare(complex *,complex *,int, int ,int ,int, RunDescriptor *,int,int);
    void doNlFFTGtoR_Gchare(complex *,complex *,int, int ,int ,int, RunDescriptor *,int,int);
    void doNlFFTRtoG_Rchare(complex *,double *,int ,int ,int ,int);
    void doNlFFTGtoR_Rchare(complex *,double *,int ,int ,int ,int);

    //-----------------------------------------------------------
    // eext fft
    void doEextFFTRtoG_Gchare(complex *,int, int ,int ,int, RunDescriptor *,int,int);
    void doEextFFTGtoR_Gchare(complex *,complex *,int, int ,int ,int, RunDescriptor *,int,int);

    void doEextFFTRtoG_Rchare(complex *,double *,int ,int ,int ,int);
    void doEextFFTGtoR_Rchare(complex *,double *,int ,int ,int ,int);

    void doEextFFTRxToGx_Rchare(complex *,double *,int ,int ,int ,int);
    void doEextFFTRyToGy_Rchare(complex *,double *,int ,int ,int ,int);
    void doEextFFTGxToRx_Rchare(complex *,double *,int ,int ,int ,int);
    void doEextFFTGyToRy_Rchare(complex *,double *,int ,int ,int ,int);

    //-----------------------------------------------------------------------------
};
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
class GSlabInfo {
  public:
    CkIndex2D index;
    int numPoints;
    int numPlanes; 
    GSlabInfo() {}
    GSlabInfo(const GSlabInfo& info) { 
      index.x = info.index.x; index.y = info.index.y;
      numPoints = info.numPoints;
      numPlanes = info.numPlanes;
    }
    GSlabInfo(CkIndex2D idx, int _numPoints, int _numPlanes){
      index.x = idx.x; index.y = idx.x;
      numPoints = _numPoints;
      numPlanes = _numPlanes;
    }
};
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
// slab initialization helpers

void initGStateSlab(GStateSlab *gs, int sizeX, int sizeY, int sizeZ, int gSpaceUnits, 
    int realSpaceUnits, int s_grain,int iplane_ind,int istate_ind,
    int len_nhc_cp, int num_nhc_cp,int nck_nhc_cp);
void initRealStateSlab(RealStateSlab *rs, int ngrid_a, int ngrid_b, int ngrid_c,
    int gSpaceUnits, int realSpaceUnits, int stateIndex, int thisPlane);

//==============================================================================
// Eric's really cool BG/L progress callers

void fft_split(FFTplanHolder *fftplanholder, int howmany, 
    fftw_complex *in,  int istride, int idist, 
    fftw_complex *out, int ostride, int odist, int split, int index);

void rfftwnd_complex_to_real_split(RFFTplanHolder *rfftplanholder, int howmany, 
    fftw_complex *in, int istride, int idist, 
    fftw_real *out,   int ostride, int odist, int split, int index);

void rfftwnd_real_to_complex_split(RFFTplanHolder *rfftplanholder, int howmany, 
    fftw_real *in,     int istride, int idist, 
    fftw_complex *out, int ostride, int odist, int split, int index);

void initFFTholder  ( FFTplanHolder *,int *,int *,int *,double *,int *,int *,int *,int *,
    int ,int *,int *);
void initRCFFTholder(RFFTplanHolder *,int *,int *,int *,double *,int *,int *,int *,int *,
    int ,int *,int *);
void initCRFFTholder(RFFTplanHolder *,int *,int *,int *,double *,int *,int *,int *,int *,
    int ,int *,int *);
void make_essl_work_map(int ,int *,int *,int *,int *,int *, int *,int);

//==============================================================================
#endif
