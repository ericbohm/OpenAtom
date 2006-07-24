//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
/** \file fftCacheSlab.h
 * Add type declarations for simulationConstants class (readonly vars) 
 * and once class for each type of object we currently have
 *(GStateXYslab, RealStateXZslab, RealRhoSlab, GRhoSlab, etc.) 
 */
//==============================================================================

#ifndef _fftcacheslab_h_
#define _fftcacheslab_h_

#include "../../include/RunDescriptor.h"

#ifdef _CP_USE_BLAS_GATHER_
include "../../src_mathlib/mathlib.h"
#endif

#define LEN_NHC_CP 4;
#define NUM_NHC_CP 20;


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
class GStateSlab {
//==============================================================================
// data 
public:

   void pup(PUP::er &);
   int cp_min_opt;
   int numNonZeroPlanes;
   int numRuns, numLines;      // numLines=numRun/2=number of lines in collection
   int numPoints;              // number of non-zero pts in the collection
   int numFull;                // expanded data : nfftz*numLines
   int numFullNL;              // expanded data : nfftz*numLines
   int ees_nonlocal;
   size2d planeSize;           // size of the ffT in Y*Z dimension
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
   double degfree;            // Degrees of freedom (ncoef_true+num_nhc-1)
   double degfreeNHC;         // Degrees of freedom (num_nhc-1)*len_nhc
   double gammaNHC;           // Degrees of freedom degfree/(degfree+1.0)

   complex *packedPlaneData;   // Non-zero data pts [numPoints]
   complex *packedPlaneDataTemp; 
   complex *packedPlaneDataScr; 
   complex *packedForceData; 
   complex *packedVelData; 
   complex *packedRedPsi;
   int     len_nhc_cp;
   int     num_nhc_cp;
   double  kTCP;
   double  tauNHCCP;
   double  xNHC;
   double  mNHC;
   double **vNHC;
   double **vNHC_scr;
   double **fNHC;

//==============================================================================
// Constuctor, Destructor and utilities

   GStateSlab() {packedPlaneData=NULL; packedPlaneDataTemp=NULL; 
                 packedForceData=NULL; packedPlaneDataScr=NULL; 
		 xNHC=0.0;
                 packedVelData=NULL;}
   ~GStateSlab();

   void copyVNHC(){
     for(int i=0;i<num_nhc_cp;i++){
     for(int j=0;j<len_nhc_cp;j++){
       vNHC_scr[i][j] = vNHC[i][j];
     }}//endfor
   }//end routine

   void initNHC(){
     vNHC     = new double *[20];
     vNHC_scr = new double *[20];
     fNHC     = new double *[20];
     for(int i=0;i<20;i++){
       vNHC[i]    =new double[4];
       vNHC_scr[i]=new double[4];
       fNHC[i]    =new double[4];
     }//endfor
     for(int i=0;i<20;i++){
     for(int j=0;j<4;j++){
       vNHC[i][j]     = 0.0;
       vNHC_scr[i][j] = 0.0;
       fNHC[i][j]     = 0.0;
     }}//endfor
   }//end routine

   void destroyNHC(){
     for(int i=0;i<20;i++){
       delete []vNHC[i];
       delete []vNHC_scr[i];
       delete []fNHC[i];
     }//endfor
     delete []vNHC;
     delete []vNHC_scr;
     delete []fNHC;
   }//end routine

   complex* doFwFFT();
   void doBwFFT(complex*);
   void setKRange(int , int *, int *, int *);
   void addForces(const complex *,const int *);
   void expandGSpace(complex* data, int type, complex *packedPlaneData);
   void compressGSpace(const complex *points, int type);
   void printGSpace(int type);

//-----------------------------------------------------------------------------
  };
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
class RealStateSlab {
public:
	size2d planeSize;         // size of the state in Z*X 
	complex *planeArr;
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
class RhoRealSlab {
public:
   int sizeX, sizeY, sizeZ;

   double *Vks;            // energy/potential  : packed
   double *density;        // rho(r) : packed
   double *doFFTonThis;    // fft sized guy
   double *rhoIRX, *rhoIRY, *rhoIRZ, *gradientCorrection; // fft sized guys

   /* return values from rhoRSubroutine in subroutine.C */
   double exc_ret, muxc_ret, exc_gga_ret;

   /* used in the subroutines doDensitySum and doRhoRSubroutine */
   int size;
   int trueSize;
   int xdim, ydim, zdim;
   int startx, starty, startz; 
 
   RhoRealSlab() {
       sizeX= sizeY= sizeZ= size= trueSize=xdim= ydim= zdim= startx= starty= startz=0 ;
       Vks=NULL;
       density=NULL;
       doFFTonThis=NULL; 
       rhoIRX=NULL;
       rhoIRY=NULL;
       rhoIRZ=NULL;
       gradientCorrection=NULL; // fft sized guys

   }
   ~RhoRealSlab();
   void doFwFFTGtoR(int,double);
   void uPackAndScale(double *, double *,double );
   void pup(PUP::er &);

//-----------------------------------------------------------------------------
 };
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
class RhoGSlab {
 public:
	/* Stuff that is currently in CP_Rho_GSpacePlane's private variables */
	int sizeX, sizeY, sizeZ;
	int runsToBeSent;
	int numRuns;
	int numLines;
	int numFull;
	int numPoints;
	int nPacked;
	RunDescriptor *runs;

	complex *Rho;       // Is big enough to be expanded!
	complex *divRhoX;   // Is big enough to be expanded!
	complex *divRhoY;   // Is big enough to be expanded!
	complex *divRhoZ;   // Is big enough to be expanded!
	complex *packedRho;
        complex *packedVks;
        complex *Vks;       // Is big enough to be expanded!
	                    
 
        int *k_x, *k_y, *k_z;

	/* return values from rhoGSubroutine in subroutine.C */
	double ehart_ret, eext_ret, ewd_ret;
	
	/* used in the subroutine doRhoGSubroutine */
	int size;
	int xdim, ydim, zdim; 
	
	RhoGSlab() { //initialization paranoia
	    sizeX= sizeY= sizeZ=runsToBeSent= numRuns= numLines=numFull=0;
	    numPoints= nPacked=size=xdim=ydim= zdim=0; 
	    k_x=NULL;
	    k_y=NULL;
	    k_z=NULL;
	    Rho=NULL;       
	    divRhoX=NULL;   
	    divRhoY=NULL;   
	    divRhoZ=NULL;   
	    packedRho=NULL;
	    packedVks=NULL;
	    Vks=NULL;       
	    runs=NULL;
	}
	~RhoGSlab();
	void doBwFFTRtoG(int); 
	void doFwFFTGtoR(int,int); 
	void setKVectors(int *n);
	void compressGSpace(const complex *, int );
	void expandRhoGSpace(complex* , complex *);
        void divRhoGdot(double *,double );
        void createWhiteByrd(double *, double );
	void pup(PUP::er &p);

//-----------------------------------------------------------------------------
  };
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
class FFTcache: public Group {
 public:
        size2d planeSize;  // sizes for state/rho methods
        int ngridaEext;    // sizes for ees methods
        int ngridbEext;
        int ngridcEext;
        int ees_eext_on;
        int ngridaNL;
        int ngridbNL;
        int ngridcNL;
        int ees_NL_on;

	complex *fftData;  // temporary for state methods
        complex *tmpData;  // temporary for Ees methods

	fftw_plan    fwdZ1DdpPlan, bwdZ1DdpPlan;        // state and density 
	fftw_plan    fwdYPlan,     bwdYPlan;            // double pack plans
	rfftwnd_plan fwdX1DdpPlan, bwdX1DdpPlan;

	fftw_plan    fwdY1DdpPlanNL, bwdY1DdpPlanNL;    // ees NL 
	fftw_plan    fwdZPlanNL,     bwdZPlanNL;        // double pack plans
	rfftwnd_plan fwdX1DdpPlanNL, bwdX1DdpPlanNL;

	fftw_plan    fwdY1DdpPlanEext, bwdY1DdpPlanEext;// ees Eext
	fftw_plan    fwdZPlanEext,     bwdZPlanEext;    // double pack plans
	rfftwnd_plan fwdX1DdpPlanEext, bwdX1DdpPlanEext;

      //-----------------------------------------------------------
      // Generic puppies
	FFTcache(size2d planeSIZE, int , int , int , 
                 int , int , int , int , int );
        void expandGSpace(complex* data, complex *packedData, 
                          RunDescriptor *runs, int numRuns, int numFull,
 		          int numPoints, int nfftz);
        void packGSpace(complex* data, complex *packedPlaneData, 
                        RunDescriptor *runs, int numRuns, int numFull,
	                int numPoints, int nfftz);
      //-----------------------------------------------------------
      // State and Density FFTs
        complex *doGSRealFwFFT(complex *packedPlaneData, RunDescriptor *runs, 
                               int numRuns, int numLines,int numFull, int numPoints,
		               int nfftz, bool fftReqd);
	double *doRealFwFFT(complex *);
	void doRealBwFFT(const double *vks, complex *,int ,int);
	void doRhoRealtoRhoG(double *realArr);
      //-----------------------------------------------------------
      // non-local fft
        void doNlFftRtoG_Gchare(complex *,int, int ,int ,int, RunDescriptor *,int);
        void doNlFftGtoR_Gchare(complex *,int, int ,int ,int, RunDescriptor *,int);
        void doNlFFTRtoG_Rchare(complex *,double *,int ,int ,int );
        void doNlFFTGtoR_Rchare(complex *,double *,int ,int ,int );
      //-----------------------------------------------------------
      // eext fft
        void doEextFFTRtoG_Gchare(complex *,int, int ,int ,int, RunDescriptor *,int);
        void doEextFFTGtoR_Gchare(complex *,int, int ,int ,int, RunDescriptor *,int);
        void doEextFFTRtoG_Rchare(complex *,double *,int ,int ,int );
        void doEextFFTGtoR_Rchare(complex *,double *,int ,int ,int );

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

void initGStateSlab(GStateSlab *gs, int sizeX, size2d size, int gSpaceUnits, 
                    int realSpaceUnits, int s_grain,int iplane_ind,int istate_ind);
void initRealStateSlab(RealStateSlab *rs, size2d planeSize, int gSpaceUnits, 
                       int realSpaceUnits, int stateIndex, int thisPlane);
void initRhoRealSlab(RhoRealSlab *rho_rs, int xdim, int ydim, int zdim, 
                     int numRealSpace, int numRhoG, int myIndexX,int myIndexY);

//==============================================================================
// Eric's BG/L progress callers

void fft_split(fftw_plan plan, int howmany, fftw_complex *in, int istride,
	       int idist, fftw_complex *out, int ostride, int odist, int split);

void rfftwnd_complex_to_real_split(rfftwnd_plan plan, int howmany, 
				   fftw_complex *in, int istride,
				   int idist, fftw_real *out, int ostride, 
				   int odist, int split);

void rfftwnd_real_to_complex_split(rfftwnd_plan plan, int howmany, 
				   fftw_real *in, int istride,
				   int idist, fftw_complex *out, 
				   int ostride, int odist, int split);

//==============================================================================
#endif
