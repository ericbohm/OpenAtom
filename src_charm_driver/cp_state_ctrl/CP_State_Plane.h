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
#include "RTH.h"
#include "StructFactorCache.h"
#include "StructureFactor.h"
void getSplitDecomp(int *,int *,int *,int , int ,int );
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
class GHartDummyMsg: public CMessage_GHartDummyMsg {
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
class GSAtmMsg: public CMessage_GSAtmMsg {
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
        int senderIndex;
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
class RhoHartRSFFTMsg: public CMessage_RhoHartRSFFTMsg {
public:
    int size; 
    int senderBigIndex;
    int senderStrtLine;
    int iopt;
    complex *data;
};
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class GSRedPsiMsg: public CMessage_GSRedPsiMsg {
public:
    int size; 
    int senderIndex;
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
class GStateOutMsg: public CMessage_GStateOutMsg {
public:
    int size; 
    int senderIndex;
    complex *data;
    complex *vdata;
    int *k_x;
    int *k_y;
    int *k_z;
};
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class CP_State_GSpacePlane: public CBase_CP_State_GSpacePlane {
 public:
        int first_step; //control flags and functions reference by thread are public
        int iwrite_now;
	int iteration;
        int nrotation;
        int myatom_integrate_flag; // 0 after I launch, 1 after return of atoms
        int myenergy_reduc_flag;   // 0 after I launch eke, 1 after return of energy
        int exitFlag;
        int finishedRedPsi;
        int finishedCpIntegrate;
        int isuspend_energy;
        int isuspend_atms;

	friend class CP_State_ParticlePlane;
	CP_State_GSpacePlane(int, size2d, int, int, int);
	CP_State_GSpacePlane(CkMigrateMessage *m);
	~CP_State_GSpacePlane(); 
	void pup(PUP::er &);
	void initGSpace(int, RunDescriptor *, int, complex *,int ,complex *,
                        int,int,int,int);

        void launchAtoms();
	void syncpsi();
	void requirePsiV();
	void doFFT();
	void startNewIter ();
	void sendPsi();
	void sendPsiV();
        void screenOutputPsi();
	void sendLambda();
	void makePCproxies();
        void doneRedPsiIntegrate();
        void sendRedPsi();
	void combineForcesGetEke();
	void integrateModForce();
        void writeStateDumpFile();
	void isAtSync(int);
	void ResumeFromSync();
	bool weneedPsiV();
        void acceptNLForces ();
        bool doneNLForces();
	bool allDoneIFFT() {return allgdoneifft;}
	void doIFFT(GSIFFTMsg *);
        void acceptAtoms(GSAtmMsg *msg);
        void acceptEnergy(GSAtmMsg *msg);
	void gdoneIFFT(CkReductionMsg *msg);
	void gdonePsiV(CkReductionMsg *msg);
        void resumePsiV (CkReductionMsg *msg);
        void psiWriteComplete(CkReductionMsg *msg);
	void releaseSFComputeZ();
	void acceptNewPsi(CkReductionMsg *msg);
        void collectFileOutput(GStateOutMsg *msg);
	void acceptNewPsiV(CkReductionMsg *msg);
	void acceptAllLambda(CkReductionMsg *msg);
        void psiCgOvlap(CkReductionMsg *msg);
	void acceptLambda(CkReductionMsg *msg);
        void acceptRedPsi(GSRedPsiMsg *msg);  
        void computeCgOverlap();
        void run ();
        void sendFFTData ();
        void doIFFT ();
	void readFile();
	void computeEnergies(int p, double d);
	void startFFT(CkReductionMsg *msg);
	
 private:
	int ireset_cg;
        int numReset_cg;
        int istart_typ_cp;
        int numRecvRedPsi;
	int countIFFT;
        int countFileOut;
        int countRedPsi;
	int ecount;
	int countPsi;
	int countVPsi;
	int countLambda;
	int *countPsiO;
	int *countVPsiO;
	int *countLambdaO;
	int AllPsiExpected;
	int AllLambdaExpected;
        int itemp;
        int jtemp;
	bool needPsiV;
	bool allgdoneifft;
	bool doneDoingIFFT;
	bool initialized;
	bool acceptedPsi;
	bool acceptedVPsi;
        bool doneNewIter;
	bool acceptedLambda;
	double ehart_total;
	double enl_total;
	double eke_total;
	double fictEke_total;
        double fmagPsi_total;
        double fmagPsi_total_old;
        double fmagPsi_total0;
        double fovlap;
        double fovlap_old;
	double egga_total;
	double eexc_total;
	double eext_total;
	double ewd_total;
	double total_energy;
        double cpuTimeNow;
	int gSpaceNumPoints;
	GStateSlab gs; 
        double *coef_mass;
        int *k_x, *k_y, *k_z;
	int *tk_x,*tk_y,*tk_z;  // Temp memory for output (size could be 0)
        complex *tpsi;          // Temp memory for output (needs careful pup)
        complex *tvpsi;         // Temp memory for output
	complex *ffttempdataGrp;// temp memory from group for fft (unpuppable!)
	complex *ffttempdata;   // temp memory malloced and freed in scope for fft
                                // probably should be pupped (size could be 0)
	CProxy_CP_State_RealSpacePlane real_proxy;
	CProxySection_StructureFactor sfCompSectionProxy;
	CProxySection_PairCalculator *lambdaproxy;
	CProxySection_PairCalculator *lambdaproxyother;
	CProxySection_PairCalculator *psiproxy;
	CProxySection_PairCalculator *psiproxyother;
	PairCalcID gpairCalcID1;
	PairCalcID gpairCalcID2;

	RTH_Runtime* run_thread; // why is this private?
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
	void doReduction(double *);
	void ResumeFromSync();	
	void pup(PUP::er &);
 private:
	int numPlanes;
	int count;
	int sendFFTDataSize;
	int size;
	double *vks;
	bool initialized, flagsRecd;
	RealStateSlab rs;
	CkSectionInfo cookie;
	CProxy_CP_State_GSpacePlane gproxy;
	RTH_Runtime* run_thread;
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
        void acceptHartVks(RhoHartRSFFTMsg *);
        void GradCorr();
        void whiteByrdFFT();
        void sendPartlyFFTtoRhoG(int );
        void acceptWhiteByrd(RhoRSFFTMsg *msg);
	void doMulticast();
	void pup(PUP::er &);
        int cp_grad_corr_on;
	void isAtSync(int iter){AtSync();};
 private:
	double FFTscale;        
	double volumeFactor;        
	double probScale;             
	int count, countFFTdata, numMcastSent;
        int rhoGHelpers;
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
	void pup(PUP::er &p);
	void isAtSync(int iter){AtSync();};
        int cp_grad_corr_on;
 private:
        int nPacked;
	int count;
        int countWhiteByrd[4];
        int doneWhiteByrd;
        int rhoGHelpers;
        int *numSplit;
        int *istrtSplit;
        int *iendSplit;
        CProxy_CP_Rho_RealSpacePlane rhoRealProxy0_com;
        CProxy_CP_Rho_RealSpacePlane rhoRealProxy1_com;
        CProxy_CP_Rho_RealSpacePlane rhoRealProxy2_com;
        CProxy_CP_Rho_RealSpacePlane rhoRealProxy3_com;
        CProxy_CP_Rho_RealSpacePlane rhoRealProxyByrd_com;
	RhoGSlab rho_gs;
	int vectorIFFTCount;
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
	void pup(PUP::er &);
	void isAtSync(int iter){AtSync();};
 private:
	RhoGSlab rho_gs;
	int iopt;
        int iteration;
        int ind_x;       // This chares index=thisIndex.x.
        int ind_xdiv;    // This chare is a subcollection of rhog(ind_xdiv).
        int ind_xrem;    // The subcollection index  0<= ind_rem < rhoGHelpers.
        int rhoGHelpers; // The number of subcolletions of each rhog(ind_xdiv).
        int istrt_lines;  // start of my subdivion of lines in rhog()
        int iend_lines;   // end of my subdivion of lines in rhog()
        int numLines;    // Number of lines in my subdivision
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
	int calcReductionPlaneNum(int);
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
        complex *zmatrix_fx,*zmatrix_fy,*zmatrix_fz;
        complex *zmatrixSum_fx,*zmatrixSum_fy,*zmatrixSum_fz;
	CkSectionInfo enlCookie; 
	CProxySection_CP_State_ParticlePlane particlePlaneENLProxy;
 public: 
};
//============================================================================

//============================================================================
#endif // #ifndef _PLANE_H_
//============================================================================


