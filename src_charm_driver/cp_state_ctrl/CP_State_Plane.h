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
#include "debug_flags.h"
#include "ckmulticast.h"
#include "structure_factor/StructFactorCache.h"
#include "structure_factor/StructureFactor.h"
#include "RTH.h"
//#include "ckPairCalculator.h"
void getSplitDecomp(int *,int *,int *,int , int ,int );




//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class ProductMsg : public CkMcastBaseMsg, public CMessage_ProductMsg {
 public:
	int datalen, hops;
	int subplane;
	double *data;
	int idx;
};
//============================================================================
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class CompAtmForcMsg: public CkMcastBaseMsg, public CMessage_CompAtmForcMsg {
 public:
   int nZmat;
   int iterNL;
   double *zmat;
   void init(int _nzmat, double *_zmat, int _iterNL)
     {
       nZmat = _nzmat;
       iterNL = _iterNL;
       CmiMemcpy(zmat, _zmat, sizeof(double)* nZmat);
     }
   friend class CMessage_CompAtmForcMsg;
};
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class GHartDummyMsg: public CMessage_GHartDummyMsg {
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
class RhoGSFFTMsg: public CMessage_RhoGSFFTMsg {
public:
	int size;
	int offset;   
	int offsetGx; 
	int iopt;
        int num;
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
	int offsetGx; 
        int iter;
        int iopt;
        int num;
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
    int index;
    int senderBigIndex;
    int senderStrtLine;
    int iopt;
    complex *data;
};
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class RPPPFFTMsg: public CMessage_RPPPFFTMsg {
public:
    int size; 
    int senderIndex;
    int numPlanes;
    complex *data;
};
//============================================================================


class CP_State_RealSpacePlane : public CBase_CP_State_RealSpacePlane {
 public:
	CP_State_RealSpacePlane(int, int,int,int,int,int,int, UberCollection);
	CP_State_RealSpacePlane(CkMigrateMessage *m) {};
	~CP_State_RealSpacePlane() { if(cookie!=NULL) delete [] cookie; };
	void acceptFFT(RSFFTMsg *);
	void doFFT();
        void doVksFFT();
	void acceptProduct(ProductMsg *);
	void doProductThenFFT();
        void sendFPsiToGSP();
	void run();
	void setNumPlanesToExpect(int num);
	void printData();
	void init(ProductMsg *);
	void doReduction();
	void ResumeFromSync();	
	void pup(PUP::er &);
	void rdoneVks(CkReductionMsg *msg);
	bool allVksDone()
	  {
	    return vksDone;
	  }
 private:
	const UberCollection thisInstance;
	int forwardTimeKeep;
	int backwardTimeKeep;
        int iplane_ind;
        int iteration;
        int rhoRsubplanes;
        int ngrida;
        int ngridb;
        int ngridc;
	int count;
	int rsize;
	int csize;
        int countProduct;
	int numCookies;
        int istate;
	UberCollection RhoReductionDest;
	bool vksDone;
	RealStateSlab rs;
	CkSectionInfo *cookie;
	CProxy_CP_State_GSpacePlane gproxy;
	RTH_Runtime* run_thread;
};
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class CP_Rho_RealSpacePlane : public CBase_CP_Rho_RealSpacePlane {
 public:	
        int listSubFlag;
        int cp_grad_corr_on;
	int ees_eext_on;
        int nplane_rho_x;
        int ngridcEext;
	int ngrida;
        int ngridb;
        int ngridc;
        int iplane_ind;
        int myNgridb;      // I don't have all ngridb lines of x now
        int myNplane_rho;  // I don't have all nplane_rho_x lines of y now
        int nptsExpndA;
        int nptsExpndB;
        int nptsA;
        int nptsB;
        int myBoff;
        int myAoff;
        int countDebug;
        int countFFTRyToGy;
        int rhoRsubplanes;
        int myTime;
	bool doneRhoReal;
	bool doneRHart;
	int recvCountFromGRho;
	int recvCountFromGHartExt;
	CP_Rho_RealSpacePlane(CkMigrateMessage *m){}
	CP_Rho_RealSpacePlane(int, bool,int,int,int, UberCollection);
	void init();
       ~CP_Rho_RealSpacePlane();
	void pup(PUP::er &);
	void acceptDensity(CkReductionMsg *);
	void handleDensityReduction();
        void launchEextRNlG();
	void energyComputation();
	void fftRhoRtoRhoG();
        void launchNLRealFFT();
        void sendPartlyFFTRyToGy(int iopt);
        void acceptRhoGradVksRyToGy(RhoGSFFTMsg *msg);
        void fftRhoRyToGy(int iopt);
        void sendPartlyFFTtoRhoG(int );
        void acceptGradRhoVks(RhoRSFFTMsg *);
        void sendPartlyFFTGxToRx(int );
        void acceptRhoGradVksGxToRx(RhoGSFFTMsg *msg);
        void GradCorr();
        void whiteByrdFFT();
        void acceptWhiteByrd(RhoRSFFTMsg *msg);
        void addWhiteByrdVks();
        void acceptHartVks(RhoHartRSFFTMsg *);
        void addHartEextVks();
	void RHartReport();
        void doMulticastCheck();
	void doMulticast();
        void exitForDebugging();
	void isAtSync(int iter){AtSync();};
	void ResumeFromSync();
 	void sendPartlyFFTtoRhoGall();
        void acceptGradRhoVksAll(RhoRSFFTMsg *msg);

 private:
	const UberCollection thisInstance;
	int rhoKeeperId;
        int rhoGHelpers;
        int countGradVks[5]; // number of collections that have arrived
        int countIntRtoG[5]; // our internal transpose friends.
        int countIntGtoR[5]; // our internal transpose friends.
        int countWhiteByrd;  // number of collections that have arrived
	int countRHart;
	int countRHartValue;
        int doneGradRhoVks; // count 1,2,3 = x,y,z all done
	bool doneWhiteByrd;
	bool doneHartVks;
	double FFTscale;        
	double volumeFactor;        
	double probScale;             
	RhoRealSlab rho_rs; 
        //Comlib multicast proxies
	CProxySection_CP_State_RealSpacePlane *realSpaceSectionProxyA;
        CProxySection_CP_State_RealSpacePlane *realSpaceSectionCProxyA;
	CProxy_CP_Rho_GSpacePlane rhoGProxy_com;
	CProxy_CP_Rho_GSpacePlane rhoGProxyIGX_com;
	CProxy_CP_Rho_GSpacePlane rhoGProxyIGY_com;
	CProxy_CP_Rho_GSpacePlane rhoGProxyIGZ_com;
	int redCount;
	CkReductionMsg *RedMsg;
	RTH_Runtime* run_thread;
};
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class CP_Rho_GSpacePlane:  public CBase_CP_Rho_GSpacePlane {
 public:
	CP_Rho_GSpacePlane(CkMigrateMessage *m) {}
	CP_Rho_GSpacePlane(int, int, int, bool, UberCollection);
	~CP_Rho_GSpacePlane();
	void run();
	void init();
	void acceptRhoData(RhoGSFFTMsg *msg);
	void doRhoFFT(); // refine the name
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
        int rhoRsubplanes;
        int countDebug;
        void exitForDebugging();
        void acceptWhiteByrdAll(RhoGSFFTMsg *msg);
        void RhoGSendRhoRall(); 
	void launchNlG();
 private:
	const UberCollection thisInstance;
	CProxySection_GSpaceDriver nlsectproxy;
	int myTime;
	int recvCountFromRRho;
        int nPacked;
        int count_stuff;
	int count;
        int countWhiteByrd[4];
        int doneWhiteByrd;
        int rhoGHelpers;
        int iplane_ind;
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
	const UberCollection thisInstance;
        int listSubFlag;
        int nplane_rho_x;
        int rhoRsubplanes;
        int registrationFlag;
        int launchFlag;
        int ngrida;
        int ngridb;
        int ngridc;
        int iplane_ind;
        int ees_eext_on;
        int natmTyp;
        int countFFT[2];
        int countIntRtoG;
        int countIntGtoR[2];
        int iteration;
        int iterAtmTyp;
        int csize;
        int nAtmTypRecv;
        int csizeInt;
        int myNgridb;
        int myBoff;
        int nptsB;
        int nptsExpndB;
        int myNplane_rho;
        int myAoff;
        int nptsA;
        int nptsExpndA;
        int countDebug;
	int recvCountFromGHartExt;
	int nchareHartAtmT;
	int natmTypTot;
	int atmTypoff;

        complex *atmSFC;
        double  *atmSFR;
        complex *atmForcC;
        double  *atmForcR;

        complex *atmEwdSFC;
        double  *atmEwdSFR;
        complex *atmEwdForcC;
        double  *atmEwdForcR;

        complex *atmSFCint;
        double  *atmSFRint;
        complex *atmForcCint;
        double  *atmForcRint;

        complex *atmEwdSFCint;
        double  *atmEwdSFRint;
        complex *atmEwdForcCint;
        double  *atmEwdForcRint;

        CProxy_CP_Rho_GHartExt rhoGHartProxy_com;
	CP_Rho_RHartExt(CkMigrateMessage *m) {}
	CP_Rho_RHartExt(int , int , int , int , int, UberCollection );
	~CP_Rho_RHartExt();
	void init();
	void pup(PUP::er &p);
        void startEextIter();
	void computeAtmSF();
        void registrationDone(CkReductionMsg *msg);
        void fftAtmSfRtoG();
        void sendAtmSfRyToGy();
        void recvAtmSfRyToGy(RhoGHartMsg *msg);
	void sendAtmSfRhoGHart();
        void recvAtmForcFromRhoGHart(RhoRHartMsg *msg);
	void fftAtmForcGtoR(int flagEwd);
        void sendAtmForcGxToRx(int iopt);
        void recvAtmForcGxToRx(RhoGHartMsg *msg);
	void computeAtmForc(int);
        void exitForDebugging();
};
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class CP_Rho_GHartExt:  public CBase_CP_Rho_GHartExt {
 public:
	CP_Rho_GHartExt(CkMigrateMessage *m) {}
	CP_Rho_GHartExt(int , int , int , int ,int, UberCollection );
	void init();
	~CP_Rho_GHartExt();
	void pup(PUP::er &);
	void acceptData(RhoGHartMsg *msg);
	void HartExtVksG();
	void FFTVks();
	void sendVks();
	void acceptVks(int size, complex * inVks);
	void acceptAtmSFTot(int size, complex * inAtm);
        void recvAtmSFFromRhoRHart(RhoGHartMsg *msg);
        void FFTEesBck();
        void getHartEextEes();
        void FFTEesFwd(int );
        void sendAtmSF(int );
	void isAtSync(int iter){AtSync();};
        int rhoRsubplanes;
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
        int CountDebug;
        int iperd;
        complex *atmSF;
        complex *atmSFtot;
        double ehart_ret;
        double eext_ret;
        double ewd_ret;
        void registrationDone(CkReductionMsg *msg);
        void exitForDebugging();
 private:
	const UberCollection thisInstance;
	complex *atmSFtotRecv;
	complex *VksRecv;
	int countAtmSFtot;
	int countVksTot;
	int nchareHartAtmT;
	int natmTypTot;
	int atmTypoff;
	int recvCountFromRHartExt;
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
	int **index_pack_tran;
};
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================





class CP_State_RealParticlePlane: public CBase_CP_State_RealParticlePlane {
 public:
  // Variables
  const UberCollection thisInstance;
   int ees_nonlocal;
   int nChareR;           // Real Space chares=# C-planes
   int nChareG;           // G Space chares
   int Rstates_per_pe;    // Real Space topomap variable
   int myPlane;           // Real space plane number

   int rhoRTime;
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

   bool launchFFT; 
   bool fftDataDone;
   bool rPlaneRedCookieSet;
   bool planeReductionReady;
   bool rEnlCookieSet;
   bool energyReductionReady;

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
#ifdef _CP_GS_DEBUG_COMPARE_VKS_
  complex *savedprojpsiC;
  complex *savedProjpsiCScr;
  double *savedProjpsiRScr;
  double *savedzmat;
  double **savedmn;
  double **saveddmn_x;
  double **saveddmn_y;
  double **saveddmn_z;
  int **savedigrid;

#endif
  //-----------
  // Functions
   CP_State_RealParticlePlane(CkMigrateMessage *m) {}
   CP_State_RealParticlePlane(int , int , int ,int , int ,int ,int,int, UberCollection);
   void init();
  ~CP_State_RealParticlePlane();
   void launchFFTControl(int );
   void pup(PUP::er &);
   void printEnlR(CkReductionMsg *m);
   void printEnlRSimp(double,int,int);
   void recvFromEesGPP(NLFFTMsg *);
   void FFTNLEesFwdR();
   void computeZmatEes();
   void planeReduction();
   void recvZMatEes(CkReductionMsg *);
   void computeAtmForcEes(CompAtmForcMsg *msg);
   void energyReduction();
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
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
//============================================================================
/**
 * The Large Sparse RhoG is where we interface with NAMD in QMMM for
 * the large grid
*/
class CP_LargeSP_RhoGSpacePlane: public CBase_CP_LargeSP_RhoGSpacePlane {
 public:
  const UberCollection thisInstance;
  // Functions
   CP_LargeSP_RhoGSpacePlane(CkMigrateMessage *m) {}
   CP_LargeSP_RhoGSpacePlane(UberCollection);
   void init();
  ~CP_LargeSP_RhoGSpacePlane();
   void pup(PUP::er &);
   void acceptMDSg();
   void acceptLSPRhoR();

};
//============================================================================

//============================================================================  
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
//============================================================================
/**
 * The Large Sparse RhoR is where we interpolate dense RhoR onto the large
 * grid for QMMM
*/
class CP_LargeSP_RhoRealSpacePlane: public CBase_CP_LargeSP_RhoRealSpacePlane {
 public:
  const UberCollection thisInstance;
  // constructor Functions
   CP_LargeSP_RhoRealSpacePlane(CkMigrateMessage *m) {}
   CP_LargeSP_RhoRealSpacePlane(UberCollection);
   // entry methods
   void init();
   void acceptLSPRhoG();
   void acceptRhoR();

   // sequential stuff
  ~CP_LargeSP_RhoRealSpacePlane();
   void pup(PUP::er &);

};
//============================================================================

//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class CP_VanderWaalsR:public CBase_CP_VanderWaalsR{
 public:	
	CP_VanderWaalsR(CkMigrateMessage *m){}
	CP_VanderWaalsR(UberCollection);
	void init();
       ~CP_VanderWaalsR();
	void pup(PUP::er &);
 private:
	const UberCollection thisInstance;
	RTH_Runtime* run_thread;
};
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class CP_VanderWaalsG:  public CBase_CP_VanderWaalsG {
 public:
	CP_VanderWaalsG(CkMigrateMessage *m) {}
	CP_VanderWaalsG(UberCollection);
	~CP_VanderWaalsG();
	void init();
	void pup(PUP::er &p);
 private:
	const UberCollection thisInstance;
};
//============================================================================


//============================================================================
#endif // #ifndef _PLANE_H_
//============================================================================


