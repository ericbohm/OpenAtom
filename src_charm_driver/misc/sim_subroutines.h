//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
/** \file sim_subroutines.C
 * Add type declarations for simulationConstants class (readonly vars) 
 * and once class for each type of object we currently have
 *(GStateXYslab, RealStateXZslab, RealRhoSlab, GRhoSlab, etc.) 
 */
//==============================================================================

#ifndef _initsim_h_
#define _initsim_h_

#ifdef _CP_USE_BLAS_GATHER_
include "../../src_mathlib/mathlib.h"
#endif


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
class GStateSlab {

public:
   void pup(PUP::er &);
   int numNonZeroPlanes;
   int numRuns, numLines;      // numLines=numRun/2=number of lines in collection
   int numPoints;              // number of non-zero pts in the collection
   int numFull;                // expanded data : nfftz*numLines
   int mysizeX;                // size of the ffT in X dimension 
   size2d planeSize;           // size of the ffT in Y*Z dimension
   bool fftReqd;               // flags to indicate whether this data set of pencils
   int S_grainSize;            // PC decomposition
   int xdim, ydim, zdim;       // don't know what this is? doIntegrate?
   int iplane_ind;             // collection index
   int istate_ind;             // state index of this collection
   int ihave_kx0;              // plane zero is in the collection (stored consecutively)
   int kx0_strt;               // starting pt
   int kx0_end;                // ending pt

   double eke_ret;            // kinetic energy
   double fovlap_loc;         // overlap

   RunDescriptor *runs;        // information about the lines in the collection [numRuns]
   complex *packedPlaneData;   // Non-zero data pts [numPoints]
   complex *packedPlaneDataTemp; 
   complex *packedPlaneDataCG; 
   complex *packedForceData; 


   GStateSlab() {packedPlaneData=NULL; packedPlaneDataTemp=NULL; 
                 packedForceData=NULL; packedPlaneDataCG=NULL;}
   ~GStateSlab();
   complex* doFwFFT();
   void doBwFFT(complex*);
   void setKVectors(int *, int **, int **, int **);
   void addForces(const complex *,const int *);
   void expandGSpace(complex* data, int type, complex *packedPlaneData);
   void compressGSpace(const complex *points, int type);
   void printGSpace(int type);
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
};
//==============================================================================



//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
class RhoRealSlab {
public:
   complex *doFFTonThis;    // copy of CP_Rho_RealSpacePlane Density data to do FFT on
   int sizeX, sizeY, sizeZ;
   fftwnd_plan fft2dFwPlan, fft2dBwPlan;
   complex *Vks;           // copy of CP_Rho_RealSpacePlane energy/potential 
   complex *density;       // copy of CP_Rho_RealSpacePlane density 

   /* return values from rhoRSubroutine in subroutine.C */
   double exc_ret, muxc_ret, exc_gga_ret;

   /* used in the subroutines doDensitySum and doRhoRSubroutine */
   int size;
   int xdim, ydim, zdim;
   int startx, starty, startz; 
 
   RhoRealSlab() {}
   ~RhoRealSlab();
   void doFwFFT();
   complex *doBwFFT();
};
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
class RhoGSlab {
 public:
	/* Stuff that is currently in CP_Rho_GSpacePlane's private variables */
	int sizeX, sizeY, sizeZ;
	fftw_plan fft1dFwPlan, fft1dBwPlan;
	complex *chunk;        
	
	/* return values from rhoGSubroutine in subroutine.C */
	double ehart_ret, eext_ret, ewd_ret;
	
	/* used in the subroutine doRhoGSubroutine */
	int size;
	int xdim, ydim, zdim; 
	int startx, starty, startz;
	
	RhoGSlab() {}
	~RhoGSlab();
	void doBwFFT(int ); 
	void doFwFFT(); 
};
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
class FFTcache: public Group {
 public:
 	fftw_plan fwdZ1DPlan, bwdZ1DPlan;
	fftw_plan fwdZ1DdpPlan, bwdZ1DdpPlan;
	fftw_plan fwdX1DPlan, bwdX1DPlan;
	fftw_plan fwdYPlan;
	fftw_plan bwdYPlan;

	rfftwnd_plan fwdX1DdpPlan;
	rfftwnd_plan bwdX1DdpPlan;
	size2d planeSize;
	complex *fftData;
	FFTcache(size2d planeSIZE, int ArraySize);
	void doFFT(complex *);
        void expandGSpace(complex* data, complex *packedPlaneData, 
                            RunDescriptor *runs, int numRuns, int numFull,
		   int numPoints, int nfftz);
        complex *doGSRealFwFFT(complex *packedPlaneData, RunDescriptor *runs, 
                          int numRuns, int numLines,int numFull, int numPoints,
		       int nfftz, bool fftReqd);
	double* doRealFwFFT(complex *);
	void doRealBwFFT(const double *vks, complex *,int ,int);
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

void initGStateSlab(GStateSlab *gs, int sizeX, size2d size, int gSpaceUnits, 
                    int realSpaceUnits, int s_grain,int iplane_ind,int istate_ind);
void initRealStateSlab(RealStateSlab *rs, size2d planeSize, int gSpaceUnits, 
                       int realSpaceUnits, int stateIndex, int thisPlane);
void initRhoRealSlab(RhoRealSlab *rho_rs, int xdim, int ydim, int zdim, 
                     int numRealSpace, int numRhoG, int myIndex);
void initRhoGSlab(RhoGSlab *rho_gs, int xdim, int ydim, int zdim, int numRealSpace, 
                  int numRealSpaceDensity, int myIndex);
void configureCPcharmParaInfoAndAtoms(double ecut_cp, int cp_min_opt, int, size2d, 
                                      CPcharmParaInfo *, int&, Atom **, const char *);

//==============================================================================

#endif
//==============================================================================
