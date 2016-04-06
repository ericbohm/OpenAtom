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
#include "fft_types.h"
//#include "ckPairCalculator.h"
extern bool HartreeFockOn;
void getSplitDecomp(int *,int *,int *,int , int ,int );

#define MY_X 0
#define MY_Y 1
#define MY_Z 2
#define MY_A MY_Z
#define MY_B MY_Y
#define MY_C MY_X

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//RAZ:  Added message for sending Dn Density to Up Instance
//      Used in sendRhoDnToRhoUp();
//      See also addition to cpaim.ci file
class RhoRDnMsg : public CMessage_RhoRDnMsg{
  public:
  int datalen;
  int time;
  double *data;
  int idx;
  };

class VksHartMsg : public CMessage_VksHartMsg{
 public:
  int datalen;
  int time;
  double *data;
  int idx;
};
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class InitDensity : public CkMcastBaseMsg, public CMessage_InitDensity {
  public:
    int grid_offset_b, grid_num_b;
    int pencil_offset_x, pencil_offset_y;
};

class VksMsg : public CkMcastBaseMsg, public CMessage_VksMsg {
  public:
    int pencil_offset_y;
    int myspin;
    double *data;
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
class FFT_Done_Msg: public CkMcastBaseMsg, public CMessage_FFT_Done_Msg {
};
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class RhoGHartMsg: public CMessage_RhoGHartMsg {
  public:
    int size;
    int mySpinIndex;
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

/** @addtogroup RealSpaceState
 *  @{
 *
 * \brief Chare Array 2D chare array nplanex X nstates.  Handles
 * electronic structure in real space. The points of plane-wave
 * pseudo-potential are cut along the x-dimension for finer
 * parallelization.
 */
class CP_State_RealSpacePlane : public CBase_CP_State_RealSpacePlane {
  public:
    CP_State_RealSpacePlane_SDAG_CODE;

    CP_State_RealSpacePlane(int, int,int,int,int,int,int, UberCollection);
    CP_State_RealSpacePlane(CkMigrateMessage *m) {};
    ~CP_State_RealSpacePlane() {
      if(cookie!=NULL) delete [] cookie;
    };

    void unpackFFT(RSFFTMsg *);
    void doFFT();
    void doVksFFT();
    void unpackVks(VksMsg *);
    void doProductThenFFT();
    void sendFPsiToGSP();
    void setNumPlanesToExpect(int num);
    void printData();
    void init(InitDensity *);
    void doReduction();
    void ResumeFromSync();
    void pup(PUP::er &);
  private:
    const UberCollection thisInstance;
    int forwardTimeKeep;
    int backwardTimeKeep;
    int iplane_ind;
    int ibead_ind,kpoint_ind, itemper_ind;
    int iteration;
    int ngrida;
    int ngridb;
    int ngridc;
    int count;
    int rsize;
    int csize;
    int countProduct;
    int numCookies;
    int istate;
    //RAZ: Added spin index:
    int mySpinIndex;

    int *grid_offset_b, *grid_num_b;
    double* hartree1;
    int rho_rpencil_offset_x, rho_rpencil_num_y;
    UberCollection RhoReductionDest;
    RealStateSlab rs;
    CkSectionInfo *cookie;
    CProxy_CP_State_GSpacePlane gproxy;
};
//============================================================================
/*@}*/

/** @addtogroup Density
  @{
 */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class CP_Rho_RealSpacePlane : public CBase_CP_Rho_RealSpacePlane {
  public:
    int myGrid_length[3];
    int myGrid_start[3], myGrid_end[3];
    int myGrid_size;
    int fft_xoffset, fft_yoffset, fft_zoffset, fft_hartoffset;
    int myTime;

    CP_Rho_RealSpacePlane(CkMigrateMessage *m) { }
    CP_Rho_RealSpacePlane(int, UberCollection);
    ~CP_Rho_RealSpacePlane();
    void pup(PUP::er &);

    void init();
    void acceptDensity(CkReductionMsg *);
    void handleDensityReduction();


    //RAZ: Added Spin dn declarations:
    void handleDensityReductionDn();  
    void sendRhoDnToRhoUp();          
    void acceptDensityDn(RhoRDnMsg *);

    
    void launchEextRNlG();
    void energyComputation();
    void launchNLRealFFT();
    void whiteByrdFFT();
    void acceptWhiteByrd();
    void doMulticastCheck();
    void doMulticast();
    void exitForDebugging();
    void isAtSync(int iter) {
      AtSync();
    };
    void ResumeFromSync() { }
    void scaleData(double *scaledData, double scaleFac);
    void acceptHartVks();
    void RHartReport();
    void acceptGradRhoVks();

    //RAZ:  Added vks to vks_dn routine:
    void sendVksHartToVksDn();
    void acceptVksHartDn(VksHartMsg *);
    void fftRhoRtoRhoG();
    void launchNlG();
    void launchEextR();

  private:
    const UberCollection thisInstance;
    //RAZ:  added spin vars here:
    int mySpinIndex;
    int cp_lsda;
    bool doneRhoUp;
    bool doneRhoDn;
    int rhoKeeperId;
    int doneGradRhoVks; // count 1,2,3 = x,y,z all done
    bool doneWhiteByrd;
    bool doneHartVks;
    bool doneRHart;
    int countRHart;
    int countRHartValue;
    double FFTscale;
    double volumeFactor;
    double probScale;
    double exc_ret, muxc_ret, exc_gga_ret;  //energy

    int *redCount, num_redn_complete;
    CkReductionMsg **RedMsg;
    CProxySection_CP_State_RealSpacePlane **realSpaceSectionProxyA;

    //data
    double *Vks;
    double *VksDn;
    double *density;
    double *densityDn;
    double *rhoIRX,*rhoIRY,*rhoIRZ;
    double *rhoIRXDn,*rhoIRYDn,*rhoIRZDn;
    double *VksHart;
    // complex pointers to the same memory as the corresponding double array
    complex *VksC;
    complex *VksDnC;
    complex *densityC;
    complex *densityDnC;
    complex *rhoIRXC,*rhoIRYC,*rhoIRZC;
    complex *rhoIRXDnC,*rhoIRYDnC,*rhoIRZDnC;
    complex *VksHartC;

};
//============================================================================
/*@}*/

/** @addtogroup Density
  @{
 */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class CP_Rho_GSpacePlane:  public CBase_CP_Rho_GSpacePlane {
  public:
    CP_Rho_GSpacePlane(CkMigrateMessage *m) {}
    CP_Rho_GSpacePlane(UberCollection);
    ~CP_Rho_GSpacePlane();
    void ResumeFromSync() { }
    void pup(PUP::er &p);
    void isAtSync(int iter){
      AtSync();
    };

    void run();
    void init();
    void acceptRhoData();
    void divRhoVksGspace();
    void acceptWhiteByrd();
    void exitForDebugging();
    void launchNlG();

    //RAZ: added spin dn send to Vks routine 
    void sendRhoGDntoHartVks();

    std::vector< gridPoint > * myPoints;
    int myGrid_length[3];
    int myGrid_start[3], myGrid_end[3];
    int myGrid_size, numPoints;
    int fft_xoffset, fft_yoffset, fft_zoffset;

  private:
    const UberCollection thisInstance;
    CProxySection_GSpaceDriver nlsectproxy;
    int myTime;
    int doneWhiteByrd;

    complex *divRhoX;
    complex *divRhoY;
    complex *divRhoZ;

    //RAZ:  added spin vars here:
    int mySpinIndex;
    int cp_lsda;


    /* return values from rhoGSubroutine in subroutine.C */
    double ehart_ret, eext_ret, ewd_ret;
};
//============================================================================
/*@}*/

/** @addtogroup Density
  @{
 */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class CP_Rho_RHartExt:  public CBase_CP_Rho_RHartExt {
  public:
    const UberCollection thisInstance;
    int registrationFlag, launchFlag;
    int natmTyp, atmTypoff;

    int iteration;
    int iterAtmTyp;
    int nAtmTypRecv;
    int countDebug;
    int myGrid_length[3];
    int myGrid_start[3], myGrid_end[3];
    int myGrid_size, numPoints;
    int fft_atmSFOffset;
    int fft_atmSFTotOffset;
    int eesRHart_index;


    complex *atmSFC;
    double  *atmSFR;
    complex *atmForcC;
    double  *atmForcR;

    complex *atmEwdSFC;
    double  *atmEwdSFR;
    complex *atmEwdForcC;
    double  *atmEwdForcR;

    CP_Rho_RHartExt(CkMigrateMessage *m) {}
    CP_Rho_RHartExt( UberCollection );
    ~CP_Rho_RHartExt();
    void init();
    void pup(PUP::er &p);
    void startEextIter();
    void computeAtmSF();
    void registrationDone();
    void recvAtmForcFromRhoGHart();
    void recvAtmForcTotFromRhoGHart();
    void computeAtmForc(int);
    void exitForDebugging();
    void doneAtmSF_FFT();
    void doneAtmSF_Multicast(FFT_Done_Msg*);
    void doneAtmSFTot_FFT();
    void doneAtmSFTot_Multicast(FFT_Done_Msg*);

};
//============================================================================
/*@}*/

/** @addtogroup Density
  @{
 */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class CP_Rho_GHartExt:  public CBase_CP_Rho_GHartExt {
  public:
    CP_Rho_GHartExt(CkMigrateMessage *m) {}
    CP_Rho_GHartExt(UberCollection );
    ~CP_Rho_GHartExt();
    void pup(PUP::er &);
    void isAtSync(int iter){
      AtSync();
    };
    void init();

    void acceptData(RhoGHartMsg *msg);
    void HartExtVksG();
    void FFTVks();
    void exitForDebugging();

    void acceptVks(RhoGHartMsg*);
    void acceptAtmSFTot(RhoGHartMsg*);
    void recvAtmSFFromRhoRHart();
    void FFTEesBck();
    void getHartEextEes();
    void FFTEesFwd(int );
    void registrationDone();
    void doneAtmSF_FFT();
    void doneAtmSF_Multicast(FFT_Done_Msg*);

    
    void operateOnData();

    std::vector< gridPoint > * myPoints;
    //RAZ:  added spin vars here:
    int numAcceptDensity;
    int mySpinIndex;
    int cp_lsda;

    int myGrid_length[3];
    int myGrid_start[3], myGrid_end[3];
    int myGrid_size, numPoints;
    int fft_hartoffset;
    int densityHere;
    int iteration;

    //ees method
    int registrationFlag;
    int launchFlag;
    int atmSFHere;
    int eesGHart_index;
    int iterAtmTyp;
    int nsendAtmTyp;
    int CountDebug;
    int natmTyp;
    int atmTypoff;
    std::vector< gridPoint > * myPoints_ext;
    int myGrid_length_ext[3];
    int myGrid_start_ext[3], myGrid_end_ext[3];
    int myGrid_size_ext, numPoints_ext;
    int fft_atmSFOffset;
    int fft_atmSFTotOffset;

    //dataValues
    double ehart_ret;
    double eext_ret;
    double ewd_ret;
    const UberCollection thisInstance;
    int countAtmSFtot;
    int countVksTot;
    complex *Rho, *Vks;
    complex *atmSF, *atmSFtot, *atmSFtotRecv;
    complex *VksRecv;
};
//============================================================================
/*@}*/

/** @addtogroup Particle
  @{
 */

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
    int ibead_ind,kpoint_ind, itemper_ind;

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
    int state0ReductionPlaneNum;
    int itimeRed;

    int registrationFlag;

    bool doOutput;
    bool launchFFT; 
    bool fftDataDone;
    bool planeRedSectionComplete;
    bool enlSectionComplete;
    bool initDone;

    double cp_enl;         // Non-local energy
    double cp_enlTot;      // Reduced Non-local energy
    double *projPsiR;      // real/final form of projector (after gx,gy FFTs)
    complex *projPsiC;     // complex/intermediate form of projector (before gx,gy FFTs)
    double *zmat;          // Non-local matrix for doublePack
    double *zmatScr;       // Non-local matrix for doublePack
    complex *zmatC;        // Non-local matrix complex for not doublePack
    complex *zmatScrC;     // Non-local matrix complex for not doublePack

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
    void planeRedSectDone(CkReductionMsg *m);
    void enlSectDone(CkReductionMsg *m);
    void initComplete();
    void launchFFTControl(int );
    void pup(PUP::er &);
    void printEnlR(CkReductionMsg *m);
    void printEnlRSimp(double,int,int);
    void recvFromEesGPP(NLFFTMsg *);
    void FFTNLEesFwdR();
    void computeZmatEes();
    void recvZMatEes(CkReductionMsg *);
    void computeAtmForcEes(CompAtmForcMsg *msg);
    void FFTNLEesBckR();
    void sendToEesGPP();
    void setPlaneRedCookie(EnlCookieMsg *);
    void setEnlCookie(EnlCookieMsg *);
    int calcReductionPlaneNum(int );
    void registrationDone();
    void recvZMatEesSimp(int , double *,int,int,int);
};
//============================================================================
/*@}*/

/** @addtogroup LargeSparse
  @{
 */

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
/*@}*/

/** @addtogroup LargeSparse
  @{
 */

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

/*@}*/


//============================================================================
#endif // #ifndef _PLANE_H_
//============================================================================


