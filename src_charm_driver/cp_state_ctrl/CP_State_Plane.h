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

class NLDummyMsg: public CMessage_NLDummyMsg {
 public:
  int iteration;
};

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
class GSPPIFFTMsg: public CMessage_GSPPIFFTMsg {
public:
	int size;
	int offset;
        int iterNL;
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
        int offset;
        int iter;
	complex *data;
};
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class RhoRHartMsg: public CMessage_RhoRHartMsg {
public:
	int size;
        int senderIndex;
        int iopt;
        int iter;
	complex *data;
};
//============================================================================



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class NLFFTMsg: public CMessage_NLFFTMsg {
public:
	int size;
        int senderIndex;
        int step;
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
class RSPPFFTMsg: public CMessage_RSPPFFTMsg {
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
        int registrationFlag;
        int istate_ind;
        int iplane_ind;
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
	int numChunks;
        int ees_nonlocal;
        int cleanExitCalled;
        int isuspendNLForces;
	bool doneDoingIFFT;
	bool triggerNL;
	friend class CP_State_ParticlePlane;
	CP_State_GSpacePlane(int, size2d, int, int, int, int);
	CP_State_GSpacePlane(CkMigrateMessage *m);
	~CP_State_GSpacePlane(); 
	void pup(PUP::er &);
	void initGSpace(int, complex *,int ,complex *,
                        int,int,int,int,int,int,int);
        void startNLEes();
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
        void acceptNLForcesEes();
        bool doneNLForces();
	bool allDoneIFFT() {return allgdoneifft;}
	void doIFFT(GSIFFTMsg *);
        void acceptAtoms(GSAtmMsg *msg);
        void acceptEnergy(GSAtmMsg *msg);
	void gdoneIFFT(CkReductionMsg *msg);
	void gdonePsiV(CkReductionMsg *msg);
	void gdonePsi(CkReductionMsg *msg);
	bool allAcceptedPsiDone(){return(allAcceptedPsi);}
        void resumePsiV (CkReductionMsg *msg);
        void psiWriteComplete(CkReductionMsg *msg);
	void releaseSFComputeZ();
	void acceptNewPsi(CkReductionMsg *msg);
	void acceptNewPsi(partialResultMsg  *msg);
	void doNewPsi();
        void collectFileOutput(GStateOutMsg *msg);
	void acceptNewPsiV(CkReductionMsg *msg);
	void acceptNewPsiV(partialResultMsg *msg);
	void doNewPsiV();
	void acceptAllLambda(CkReductionMsg *msg);
        void psiCgOvlap(CkReductionMsg *msg);
	void acceptLambda(CkReductionMsg *msg);
	void acceptLambda(partialResultMsg *msg);
	void doLambda();
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
	bool initialized;
	bool allAcceptedPsi;
	bool acceptedPsi;
	bool acceptedVPsi;
	bool allAcceptedVPsi;
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
	int *tk_x,*tk_y,*tk_z;  // Temp memory for output (size could be 0)
        complex *tpsi;          // Temp memory for output (needs careful pup)
        complex *tvpsi;         // Temp memory for output
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
	CP_State_RealSpacePlane(size2d, int, int,int,int,int);
	CP_State_RealSpacePlane(CkMigrateMessage *m) {};
	~CP_State_RealSpacePlane() { };
	void doFFT(RSFFTMsg *);
	void doFFT();
	void doProduct(ProductMsg *);
	void doProductThenFFT();
        void sendFPsiToGSP();
	void run();
	void setNumPlanesToExpect(int num);
	void printData();
	void init(ProductMsg *);
	void doReduction();
	void ResumeFromSync();	
	void pup(PUP::er &);
 private:
        int ngrida;
        int ngridb;
        int ngridc;
	int count;
	int rsize;
	int csize;
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
        int cp_grad_corr_on;
	int ees_eext_on;
        int nplane_rho_x;
        int ngridcEext;
	int ngrida;
        int ngridb;
        int ngridc;
        int nptsExpnd;
        int npts;
	CP_Rho_RealSpacePlane(CkMigrateMessage *m){}
	CP_Rho_RealSpacePlane(int, size2d, bool,int,int);
       ~CP_Rho_RealSpacePlane();
	void run();
	void acceptDensity(CkReductionMsg *);
	void acceptDensity();
	void energyComputation();
	void fftRhoRtoRhoG();
        void acceptGradRhoVks(RhoRSFFTMsg *);
        void acceptHartVks(RhoHartRSFFTMsg *);
        void GradCorr();
        void whiteByrdFFT();
        void acceptWhiteByrd(RhoRSFFTMsg *msg);
        void sendPartlyFFTtoRhoG(int );
	void doMulticast();
	void isAtSync(int iter){AtSync();};
	void ResumeFromSync();
	void pup(PUP::er &);
 private:
        int rhoGHelpers;
        int countGradVks[5]; // number of collections that have arrived
        int countWhiteByrd;  // number of collections that have arrived
        int doneGradRhoVks; // count 1,2,3 = x,y,z all done
	bool doneWhiteByrd;
	bool doneHartVks;
	double FFTscale;        
	double volumeFactor;        
	double probScale;             
	RhoRealSlab rho_rs; 
        //Comlib multicast proxy
	CProxySection_CP_State_RealSpacePlane realSpaceSectionProxy;
        CProxySection_CP_State_RealSpacePlane realSpaceSectionCProxy;
	CProxy_CP_Rho_GSpacePlane rhoGProxy_com;
	CProxy_CP_Rho_GSpacePlane rhoGProxyIGX_com;
	CProxy_CP_Rho_GSpacePlane rhoGProxyIGY_com;
	CProxy_CP_Rho_GSpacePlane rhoGProxyIGZ_com;
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
	void run();
	void acceptRhoData(RhoGSFFTMsg *msg);
	void acceptRhoData(); // refine the name
	void ResumeFromSync();
        void divRhoVksGspace();
        void RhoGSendRhoR(int );
        void acceptWhiteByrd(RhoGSFFTMsg *);
        void acceptWhiteByrd();
	void pup(PUP::er &p);
	void isAtSync(int iter){AtSync();};
        int cp_grad_corr_on;
        int ees_eext_on;
        int ngridcEext;
 private:
        int nPacked;
	int count;
        int countWhiteByrd[4];
        int doneWhiteByrd;
        int rhoGHelpers;
        int *numSplit;
        int *istrtSplit;
        int *iendSplit;
	RhoGSlab rho_gs;
        CProxy_CP_Rho_RealSpacePlane rhoRealProxy0_com;
        CProxy_CP_Rho_RealSpacePlane rhoRealProxy1_com;
        CProxy_CP_Rho_RealSpacePlane rhoRealProxy2_com;
        CProxy_CP_Rho_RealSpacePlane rhoRealProxy3_com;
        CProxy_CP_Rho_RealSpacePlane rhoRealProxyByrd_com;
};
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class CP_Rho_RHartExt:  public CBase_CP_Rho_RHartExt {
 public:
        int registrationFlag;
        int launchFlag;
        int ngrida;
        int ngridb;
        int ngridc;
        int ees_eext_on;
        int natmTyp;
        int countFFT[2];
        int iteration;
        int iterAtmTyp;
        int csize;
        int nAtmTypRecv;
        complex *atmSFC;
        double  *atmSFR;
        complex *atmForcC;
        double  *atmForcR;

        complex *atmEwdSFC;
        double  *atmEwdSFR;
        complex *atmEwdForcC;
        double  *atmEwdForcR;
        CProxy_CP_Rho_GHartExt rhoGHartProxy_com;

	CP_Rho_RHartExt(CkMigrateMessage *m) {}
	CP_Rho_RHartExt(int , int , int , int , int );
	~CP_Rho_RHartExt();
	void pup(PUP::er &p);
        void startEextIter();
	void computeAtmSF();
	void FFTSFBck();
	void SendAtmSFRhoGHart();
	void recvAtmForcFromRhoGHart(RhoRHartMsg *msg);
	void FFTAtmForcFwd(int);
	void computeAtmForc(int);
        void registrationDone(CkReductionMsg *msg);
};

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class CP_Rho_GHartExt:  public CBase_CP_Rho_GHartExt {
 public:
	CP_Rho_GHartExt(CkMigrateMessage *m) {}
	CP_Rho_GHartExt(size2d,int , int , int , int ,int );
	~CP_Rho_GHartExt();
	void pup(PUP::er &);
	void acceptData(RhoGHartMsg *msg);
	void HartExtVksG();
	void FFTVks();
	void sendVks();
        void recvAtmSFFromRhoRHart(RhoGHartMsg *msg);
        void FFTEesBck();
        void getHartEextEes();
        void FFTEesFwd(int );
        void sendAtmSF(int );
	void isAtSync(int iter){AtSync();};
        int ngridaEext;
        int ngridbEext;
        int ngridcEext;
        int ees_eext_on;
        int natmTyp;
        int iterAtmTyp;
        int nsendAtmTyp;
        int numFullEext;
        int registrationFlag;
        int launchFlag;
        complex *atmSF;
        complex *atmSFtot;
        double ehart_ret;
        double eext_ret;
        double ewd_ret;
        void registrationDone(CkReductionMsg *msg);
 private:
	RhoGSlab rho_gs;
        int atmSFHere;
        int densityHere;
        int countEextFFT;
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
        CProxy_CP_Rho_RHartExt       rhoRHartProxy_com0;
        CProxy_CP_Rho_RHartExt       rhoRHartProxy_com1;
};
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class CP_State_ParticlePlane: public CBase_CP_State_ParticlePlane {
 public:
	CP_State_ParticlePlane(CkMigrateMessage *m) {}
	CP_State_ParticlePlane(int ,int ,int ,int ,int ,int ,int ,int ,int ,int ,
                               int ,int ,int ,int ,int );
	~CP_State_ParticlePlane();
	void pup(PUP::er &);
	void startNLEes(int);
	void lPrioStartNLEes(NLDummyMsg *m);
	void computeZ(PPDummyMsg *m);
	void setEnlCookie(EnlCookieMsg *m);
	void ResumeFromSync();
	void reduceZ(int, int, complex *,complex *,complex *,complex *);
	void getForces(int, int, complex *);

        void createNLEesFFTdata();
        void FFTNLEesFwd();
        void sendToEesRPP();
        void recvFromEesRPP(GSPPIFFTMsg *msg);
        void FFTNLEesBck();
        void computeNLEesForces();
        void registrationDone(CkReductionMsg *msg);

	friend class CP_State_GSpacePlane;
        int myChareG;
	int iteration;
        int iterNL;
        int numNLiter;
        int ees_nonlocal;
        int ngridaNL;
        int ngridbNL;
        int ngridcNL;
	int gSpaceNumPoints;
        int numLines;
        int numFullNL;
        int natm_nl;
        int natm_nl_grp_max;
        int numSfGrps;
	int nstates;
	int nchareG;
	int Gstates_per_pe;
        int countNLIFFT;
        int sendDone;
        int registrationFlag;
 private:
	int calcReductionPlaneNum(int);
	void initKVectors(GStateSlab *);
	bool doneGettingForces;
	complex *myForces, *gspace, *projPsiG;
	complex *zmatrixSum, *zmatrix;
        double *dyp_re,*dyp_im;
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
	CProxy_CP_State_RealParticlePlane realPP_proxy;
 public: 
};
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class CP_State_RealParticlePlane: public CBase_CP_State_RealParticlePlane {
 public:
  // Variables
   int ees_nonlocal;
   int nChareR;           // Real Space chares=# C-planes
   int nChareG;           // G Space chares
   int Rstates_per_pe;    // Real Space topomap variable
   int myPlane;           // Real space plane number

   int numIterNl;         // # of non-local iterations per time step
   int countZ;
   int countEnl;
   int count;             // fft communication counter
   int iterNL;            // Nl iteration counter
   int itime;             // time step counter;
   int recvBlock;
 
   int ngridA;            // FFT grid size along a
   int ngridB;            // FFT grid size along b
   int ngridC;            // FFT grid size along c
   int planeSize;         // expanded plane size for FFTing
   int planeSizeT;        // true plane size 
   int csize;             // complex variable size for FFT
   int zmatSizeMax;       // zmatrix size for projector
   int reductionPlaneNum; // Reduction Plane number
   int itimeRed;

   int registrationFlag;

   double cp_enl;         // Non-local energy
   double cp_enlTot;      // Reduced Non-local energy
   double *projPsiR;      // real/final form of projector (after gx,gy FFTs)
   double *zmat;          // Non-local matrix
   double *zmatScr;      // Non-local matrix
   complex *projPsiC;     // complex/intermediate form of projector (before gx,gy FFTs)

  //-----------
  // Proxies

   CProxySection_CP_State_RealParticlePlane rPlaneSectProxy; // Section Red proxy zmat
   CProxySection_CP_State_RealParticlePlane rPlaneENLProxy;  // Section Red proxy cp_enl
   CkSectionInfo rPlaneRedCookie;   // Section Red cookie for zmat
   CkSectionInfo rEnlCookie;        // Section Red cookie for cp_enl
   CProxy_CP_State_ParticlePlane gPP_proxy;

  //-----------
  // Functions
   CP_State_RealParticlePlane(CkMigrateMessage *m) {}
   CP_State_RealParticlePlane(int , int , int ,int , int ,int ,int,int);
  ~CP_State_RealParticlePlane();
   void pup(PUP::er &);
   void printEnlR(CkReductionMsg *m);
   void printEnlRSimp(double,int,int);
   void recvFromEesGPP(NLFFTMsg *);
   void FFTNLEesFwdR();
   void computeZmatEes();
   void recvZMatEes(CkReductionMsg *);
   void computeAtmForcEes(int, double *,int );
   void FFTNLEesBckR();
   void sendToEesGPP();
   void setPlaneRedCookie(EnlCookieMsg *);
   void setEnlCookie(EnlCookieMsg *);
   int calcReductionPlaneNum(int );
   void registrationDone(CkReductionMsg *msg);
   void recvZMatEesSimp(int , double *,int,int,int);
};
//============================================================================


//============================================================================
#endif // #ifndef _PLANE_H_
//============================================================================


