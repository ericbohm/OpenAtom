//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
#ifndef _PLANE_H_
#define _PLANE_H_

//============================================================================
#include "charm++.h"
#include "ckmulticast.h"
#include "fftlib.h"
#include "RTH.h"
#include "StructFactorCache.h"
#include "StructureFactor.h"
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class ProductMsg : public CkMcastBaseMsg, public CMessage_ProductMsg {
 public:
	int datalen, hops;
	double *data;
	int idx;
};
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class EnlCookieMsg : public CkMcastBaseMsg, public CMessage_EnlCookieMsg {
 public:
  int foo;
};
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class RSDummyResume: public CMessage_RSDummyResume {
};
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class TMsg: public CMessage_TMsg {
public:
	int datalen;
	complex *data;
};
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class GSIFFTMsg: public CMessage_GSIFFTMsg {
public:
	int size;
	int offset;
	complex *data;
/*	GSIFFTMsg(int _size, int _offset, complex *_data) : size(_size), offset(_offset) {
	    CmiMemcpy(data, _data, sizeof(complex)*size);
	}
*/
};
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class RSFFTMsg: public CMessage_RSFFTMsg {
public:
    int size; 
    int senderIndex;
    int numPlanes;
    complex *data;
};
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class CP_State_GSpacePlane: public CBase_CP_State_GSpacePlane {
 public:
	CP_State_GSpacePlane(int, size2d, int, int, int);
	CP_State_GSpacePlane(CkMigrateMessage *m);
	~CP_State_GSpacePlane(); 
	void initGSpace(int, RunDescriptor *, int, complex *,int,int,int);
	void doIFFT(GSIFFTMsg *);
	void doFFT();
	void startNewIter ();
	void sendPsi();
	void sendLambda();
	void makePCproxies();
	void acceptNewPsi(CkReductionMsg *msg);
	void acceptAllLambda(CkReductionMsg *msg);
        void psiCgOvlap(CkReductionMsg *msg);
	void acceptLambda(CkReductionMsg *msg);
	void getForcesAndIntegrate();
	void integrateModForce();
	void pup(PUP::er &);
	friend class CP_State_ParticlePlane;
	void isAtSync(int);
	void ResumeFromSync();
	void computeEnergies(void *param, void *msg);
	
	void computeEnergies(int p, double d);
        void computeCgOverlap();

        void run ();
        void resumeThread (PPDummyMsg *dmsg);
        void sendFFTData ();
        void resetIterState ();
        void doIFFT ();
	void readFile();
        bool integrateForces ();
	// for experimental barrier
	bool allDoneIFFT() {return allgdoneifft;}
	void gdoneIFFT(CkReductionMsg *msg);
	void releaseComputeZ();
        int first_step;

 private:
	bool allgdoneifft;
	bool doneDoingIFFT;
	bool initialized;
	int iteration;
	int count;
	GStateSlab gs; 
	bool flagsSent;
	int partialCount;
	int gSpaceNumPoints;
        int* k_x, *k_y, *k_z;
        double *coef_mass;
	int sendFFTDataSize;
	double ehart_total;
	double enl_total;
	double eke_total;
        double fmagPsi_total;
        double fovlap;
        double fovlap_old;
	double egga_total;
	double eexc_total;
	double eext_total;
	double ewd_total;
	int ecount, displace_count;
	double total_energy;
	int allEnergiesReceived;	
	int localState;
	int AllExpected;
	bool acceptedPsi;

	complex *ffttempdata; 
	int ireset;
	CProxy_CP_State_RealSpacePlane real_proxy;
	CProxySection_StructureFactor sfCompSectionProxy;
	CProxySection_PairCalculator lambdaproxy;
	CProxySection_PairCalculator psiproxy;
	CProxySection_PairCalculator psiproxyother;
	RTH_Runtime* run_thread;
	PairCalcID gpairCalcID1;
	PairCalcID gpairCalcID2;
};
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class CP_State_RealSpacePlane : public CBase_CP_State_RealSpacePlane {
 public:
	CP_State_RealSpacePlane(size2d, int, int);
	CP_State_RealSpacePlane(CkMigrateMessage *m) {};
	~CP_State_RealSpacePlane() { };
	void doFFT(RSFFTMsg *);
	void doFFT();
	void doProduct(ProductMsg *);
	void doProduct(int, const double *);
	void doProduct();
	void resumeProduct(RSDummyResume *msg);
	void run();
	void setNumPlanesToExpect(int num);
	void printData();
	void init(ProductMsg *);
	void doReduction();
	void ResumeFromSync();	
 private:
	int numPlanes;
	int count;
	double *vks;
	int size;
	bool initialized, flagsRecd;
	RealStateSlab rs;
	CkSectionInfo cookie;
	CProxy_CP_State_GSpacePlane gproxy;
	RTH_Runtime* run_thread;
	int sendFFTDataSize;
};
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class CP_Rho_RealSpacePlane : public CBase_CP_Rho_RealSpacePlane {
 public:	
	CP_Rho_RealSpacePlane(CkMigrateMessage *m){}
	CP_Rho_RealSpacePlane(int, size2d, int, int, bool, ComlibInstanceHandle);
	~CP_Rho_RealSpacePlane();
	void acceptDensity(CkReductionMsg *);
	void acceptDensity();
	void acceptDensity(int, double *, int);
	void acceptEnergyForSumming(int, complex *, int);
	void acceptEnergyForSumming();
	void doRingMulticast();
	void energyComputation();
	void startTranspose();
	void doneFFT();
	void run();
	void ResumeFromSync();
	void resumeThread(PPDummyMsg *dmsg);
	void foo(void);
	bool gotAllRhoEnergy;
 	bool doneDoingFFT;
 private:

	double volumeFactor;        
	double probScale;             
	int count, countFFTdata, numMcastSent;
	CProxySection_CP_State_RealSpacePlane realSpaceSectionProxy;
        //Comlib multicast proxy
        CProxySection_CP_State_RealSpacePlane realSpaceSectionCProxy;
	CProxy_CP_Rho_GSpacePlane rhoGProxy_com;
	RhoRealSlab rho_rs; 
	int cmid;
	int *pes;
	bool doneGradientCorrection, doneRhoGStuff;
	complex *rhoIRX, *rhoIRY, *rhoIRZ, *gradientCorrection;
	void doMulticast();
	void doneFFT(int);
	
	int vectorFFTCount;
	int id;
	RTH_Runtime* run_thread;
};
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class CP_Rho_GSpacePlane: public NormalSlabArray {
 public:
	CP_Rho_GSpacePlane(CkMigrateMessage *m) {}
	CP_Rho_GSpacePlane(int, size2d, int, int, bool, ComlibInstanceHandle);
	~CP_Rho_GSpacePlane();
	void acceptData(int, complex *, int, int, int);
	void acceptData();
	void recvProcessedPart(int, complex *, int);
	void recvProcessedPart();
	void run();
	void ResumeFromSync();
 private:
	int count, helperCount, helperWidth;
        CProxy_CP_Rho_RealSpacePlane rhoRealProxy_com;
	RhoGSlab rho_gs;
	AtomsGrp *atom;
	complex *rhoIGX, *rhoIGY, *rhoIGZ, *gradientCorrection;
	void doneIFFT(int);
	void doneIFFT();
	int vectorIFFTCount;
	int* k_x, *k_y, *k_z;
        void computeK(const int size,
		      const int startx, const int starty, const int startz,
		      const int xdim, const int ydim, const int zdim);
	RTH_Runtime* run_thread;
};
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class CP_Rho_GSpacePlaneHelper: public CBase_CP_Rho_GSpacePlaneHelper {
public:
	CP_Rho_GSpacePlaneHelper(CkMigrateMessage *m) {}
	CP_Rho_GSpacePlaneHelper(int, size2d, int);
	~CP_Rho_GSpacePlaneHelper();
	void recvRhoGPart(int, complex *);
private:
	RhoGSlab rho_gs;
	int *k_x, *k_y, *k_z;
	int helperWidth, helperPos;
	complex *temp;
};
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class CP_State_ParticlePlane: public CBase_CP_State_ParticlePlane {
 public:
	CP_State_ParticlePlane(CkMigrateMessage *m) {}
	CP_State_ParticlePlane(int, int, int, int,int ,int,int);
	~CP_State_ParticlePlane();
	void computeZ(PPDummyMsg *m);
	void setEnlCookie(EnlCookieMsg *m);
	void ResumeFromSync();
	void reduceZ(int, int, complex *,complex *,complex *,complex *);
	void sumEnergies(double);
	void getForces(int, int, complex *);
	void pup(PUP::er &);
	friend class CP_State_GSpacePlane;
	int *k_x, *k_y, *k_z;
	int gSpaceNumPoints;
        int natm_nl;
        int natm_nl_grp_max;
        int numSfGrps;
 private:
	void initKVectors(GStateSlab *);
	bool doneGettingForces;
	complex *myForces, *gspace;
	complex *zmatrixSum, *zmatrix;
        double enl;
        double enl_total;
	double totalEnergy;
	int *haveSFAtmGrp;
	int *count;
	int doneEnl;
	int doneForces;
	int zsize, energy_count;
	int sizeX, sizeY, sizeZ, gSpacePlanesPerChare;
	int reductionPlaneNum;
	// DY_NONLOCAL
        complex *zmatrix_fx,*zmatrix_fy,*zmatrix_fz;
        complex *zmatrixSum_fx,*zmatrixSum_fy,*zmatrixSum_fz;
        // DY_NONLOCAL 
	CkSectionInfo enlCookie; 
	CProxySection_CP_State_ParticlePlane particlePlaneENLProxy;
 public: 
};
//============================================================================

//============================================================================
#endif // #ifndef _PLANE_H_
//============================================================================


