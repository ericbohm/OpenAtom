//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \file CP_State_Plane.h
 *
 */
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
};
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class RhoGSFFTMsg: public CMessage_RhoGSFFTMsg {
public:
	int size;
	int offset;
	int iopt;
	complex *data;
};
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class RhoGHartMsg: public CMessage_RhoGHartMsg {
public:
	int size;
	complex *data;
};
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class RhoRSFFTMsg: public CMessage_RhoRSFFTMsg {
public:
    int size; 
    int senderIndex;
    int iopt;
    complex *data;
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
	void initGSpace(int, RunDescriptor *, int, complex *,int ,complex *,
                        int,int,int,int);
	void syncpsi();
	void requirePsiV();
	void doIFFT(GSIFFTMsg *);
	void doFFT();
	void startNewIter ();
	void sendPsi();
	void sendPsiV();
	bool weneedPsiV();
        void screenOutputPsi();
	void sendLambda();
	void makePCproxies();
	void acceptNewPsi(CkReductionMsg *msg);
	void acceptNewPsiV(CkReductionMsg *msg);
	void acceptAllLambda(CkReductionMsg *msg);
        void psiCgOvlap(CkReductionMsg *msg);
	void acceptLambda(CkReductionMsg *msg);
	void combineForcesGetEke();
	void integrateModForce();
        void writeStateDumpFile();
	void pup(PUP::er &);
	friend class CP_State_ParticlePlane;
	void isAtSync(int);
	void ResumeFromSync();
	void computeEnergies(void *param, void *msg);
	
	void computeEnergies(int p, double d);
        void computeCgOverlap();

        void run ();
        void resumePsiV (CkReductionMsg *msg);
        void resumeThread (PPDummyMsg *dmsg);
        void sendFFTData ();
        void doIFFT ();
	void readFile();
        bool completedExtExcNlForces ();
	// for experimental barrier
	bool allDoneIFFT() {return allgdoneifft;}
	void gdoneIFFT(CkReductionMsg *msg);
	void gdonePsiV(CkReductionMsg *msg);
	void releaseSFComputeZ();
        int first_step;
        int iwrite_now;
        void psiWriteComplete(CkReductionMsg *msg);
	int iteration;
 private:
	bool needPsiV;
	bool allgdoneifft;
	bool doneDoingIFFT;
	bool initialized;

        int istart_typ_cp;
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
	double fictEke_total;
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
	CP_Rho_RealSpacePlane(int, size2d, int, int, bool);
	~CP_Rho_RealSpacePlane();
	void acceptDensity(CkReductionMsg *);
	void acceptDensity();
	void acceptDensity(int, double *, int);
	void acceptEnergyForSumming(int, complex *, int);
	void acceptEnergyForSumming();
	void energyComputation();
	void doneFFT();
	void run();
	void ResumeFromSync();
	void fftRhoRtoRhoG();
        void acceptGradRhoVks(RhoRSFFTMsg *);
        void acceptHartVks(RhoRSFFTMsg *);
        void GradCorr();
        void whiteByrdFFT();
        void sendPartlyFFTtoRhoG(int );
        void acceptWhiteByrd(RhoRSFFTMsg *msg);
	void doMulticast();
 private:
	double FFTscale;        
	double volumeFactor;        
	double probScale;             
	int count, countFFTdata, numMcastSent;
	CProxySection_CP_State_RealSpacePlane realSpaceSectionProxy;
        int countGradVks[5], doneGradRhoVks;
        int countWhiteByrd;
	bool doneWhiteByrd;
	bool doneHartVks;
        //Comlib multicast proxy
        CProxySection_CP_State_RealSpacePlane realSpaceSectionCProxy;
	CProxy_CP_Rho_GSpacePlane rhoGProxy_com;
	CProxy_CP_Rho_GSpacePlane rhoGProxyIGX_com;
	CProxy_CP_Rho_GSpacePlane rhoGProxyIGY_com;
	CProxy_CP_Rho_GSpacePlane rhoGProxyIGZ_com;
	RhoRealSlab rho_rs; 
	int cmid;
	int *pes;
	RTH_Runtime* run_thread;
};
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class CP_Rho_GSpacePlane:  public CBase_CP_Rho_GSpacePlane {
 public:
	CP_Rho_GSpacePlane(CkMigrateMessage *m) {}
	CP_Rho_GSpacePlane(int, size2d, int, int, bool);
	~CP_Rho_GSpacePlane();
	void acceptData(RhoGSFFTMsg *msg);
	void acceptData();
	void recvProcessedPart(int, complex *, int);
	void recvProcessedPart();
	void run();
	void ResumeFromSync();
        void divRhoVksGspace();
        void RhoGSendRhoR(int );
        void acceptWhiteByrd(RhoGSFFTMsg *);
        void acceptWhiteByrd();
 private:
        int nPacked;
	int count;
        int countWhiteByrd[4];
        int doneWhiteByrd;
        CProxy_CP_Rho_RealSpacePlane rhoRealProxy0_com;
        CProxy_CP_Rho_RealSpacePlane rhoRealProxy1_com;
        CProxy_CP_Rho_RealSpacePlane rhoRealProxy2_com;
        CProxy_CP_Rho_RealSpacePlane rhoRealProxy3_com;
        CProxy_CP_Rho_RealSpacePlane rhoRealProxyByrd_com;
	RhoGSlab rho_gs;
	AtomsGrp *atom;
	int vectorIFFTCount;
	int *k_x, *k_y, *k_z;
	RTH_Runtime* run_thread;
};
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class CP_Rho_GHartExt:  public CBase_CP_Rho_GHartExt {
 public:
	CP_Rho_GHartExt(CkMigrateMessage *m) {}
	CP_Rho_GHartExt(size2d);
	~CP_Rho_GHartExt();
	void acceptData(RhoGHartMsg *msg);
	void HartExtVksG();
	void sendVks();
 private:
	RhoGSlab rho_gs;
	int iopt;
        CProxy_CP_Rho_RealSpacePlane rhoRealProxy_com;
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


