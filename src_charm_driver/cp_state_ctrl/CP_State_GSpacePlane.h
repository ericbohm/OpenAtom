#include "debug_flags.h"
#include "charm++.h" /// @note: Needed for ckcomplex.h ?
#include "ckcomplex.h"
#include "gStatePlane.decl.h"
#include "uber/Uber.h"
#include "fft_slab_ctrl/fftCacheSlab.h"
#include "paircalc/pairCalculator.h"
struct EnergyStruct; /// @warning: Forward declarations of structs seem to choke the ci parser. It doesnt recognize the keyword struct.
#include "paircalc/ckPairCalculator.h"
#include "cpaimd.decl.h"

#ifndef CP_STATE_GSPACE_PLANE_H
#define CP_STATE_GSPACE_PLANE_H

/// Forward declarations 
class GStateSlab;
class RDMApair_GSP_PC;
template <class tokenType> class RDMASetupConfirmationMsg;

 
class GSIFFTMsg: public CMessage_GSIFFTMsg 
{
    public:
        int size;
        int offset;
        complex *data;
};




class GSRedPsiMsg: public CMessage_GSRedPsiMsg 
{
    public:
        int size;
        int senderIndex;
        complex *data;
};




class GStateOutMsg: public CMessage_GStateOutMsg 
{
    public:
        int size;
        int senderIndex;
        complex *data;
        complex *vdata;
        int *k_x;
        int *k_y;
        int *k_z;
};




class RSFFTMsg: public CMessage_RSFFTMsg {
public:
    int size; 
    int senderIndex;
    int numPlanes;
    complex *data;
};




class CP_State_GSpacePlane: public CBase_CP_State_GSpacePlane 
{
    public:
        friend class CP_State_ParticlePlane;
        // ----------- Flags and counters used in the GSpace driver code -----------
        /// My state index
		int istate_ind;
		/// My plane index
		int iplane_ind;
		///
		int iteration;
		///
		int iRecvRedPsi;
		///
		int cleanExitCalled;
		///
		bool doneDoingIFFT;

        int halfStepEvolve;
        int redPlane;
        int registrationFlag;
        int nrotation;
        int exitFlag;
        int iRecvRedPsiV;
        int iSentRedPsi;
        int iSentRedPsiV;
        int finishedCpIntegrate;
        int numRecvRedPsi;
        int iterRotation;
        double ake_old;
        
        bool acceptedVPsi;
        CP_State_GSpacePlane( int, int, int, int,int,int, UberCollection);
        CP_State_GSpacePlane(CkMigrateMessage *m);
        ~CP_State_GSpacePlane();
        /// Gets called from the PairCalc data receivers to confirm the setup of an RDMA link
        void completeRDMAhandshake(RDMASetupConfirmationMsg<RDMApair_GSP_PC> *msg);
        void pup(PUP::er &);
        void initGSpace(int, complex *,int ,complex *,int,int,int,int,int,int,int);
        void launchAtoms();
        void launchOrthoT();
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
        /// @entry This is used to receive data from all the corresponding RealSpacePlanes, upon which the inverse FFTs are triggered
        void acceptIFFT(GSIFFTMsg *);
        void doIFFT();
        void acceptNewPsi(CkReductionMsg *msg);
        void acceptNewPsi(partialResultMsg  *msg);
        /// @RTH resume
        void doNewPsi();
        void collectFileOutput(GStateOutMsg *msg);
        void acceptNewPsiV(CkReductionMsg *msg);
        void acceptNewPsiV(partialResultMsg *msg);
        /// @RTH resume
        void doNewPsiV();
        /// @warning: Deprecated! @todo: Deprecated
        void acceptAllLambda(CkReductionMsg *msg);
        /// @entry @RTH resume
        void psiCgOvlap(CkReductionMsg *msg);
        void acceptLambda(CkReductionMsg *msg);
        void acceptLambda(partialResultMsg *msg);
        /// @RTH resume
        void doLambda();
        /// @RTH resume
        void acceptRedPsi(GSRedPsiMsg *msg);
        
        void computeCgOverlap();
        void sendFFTData ();
        
        void readFile();
        void computeEnergies(int p, double d);
        void startFFT(CkReductionMsg *msg);
        void sendRedPsiV();
        void acceptRedPsiV(GSRedPsiMsg *msg);
        void doneRedPsiVIntegrate();
        
    private:
        const UberCollection thisInstance;
        int gotHandles;
        int forwardTimeKeep;
        int backwardTimeKeep;
        int ireset_cg;
        int numReset_cg;
        int istart_typ_cp;
        int countIFFT;
        int countFileOut;
        int countRedPsi;
        int countRedPsiV;
        int ecount;
        int countPsi;
        int countVPsi;
        int countLambda;
        int *countPsiO;
        int *countVPsiO;
        int *countLambdaO;
        int AllPsiExpected;
        int AllLambdaExpected;
        /// The number of symmetric and asymmetric PCs that communicate with me
        int numRDMAlinksSymm, numRDMAlinksAsymm;
        int itemp;
        int jtemp;
        bool initialized;
        bool acceptedPsi;
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
        CProxySection_PairCalculator *lambdaproxy;
        CProxySection_PairCalculator *lambdaproxyother;
        CProxySection_PairCalculator *psiproxy;
        CProxySection_PairCalculator *psiproxyother;
        PairCalcID gpairCalcID1;
        PairCalcID gpairCalcID2;
        #ifdef _CP_GS_DEBUG_COMPARE_VKS_
            complex *savedvksBf;
            complex *savedforceBf;
        #endif
        #ifdef  _CP_GS_DEBUG_COMPARE_PSI_
        // place to keep data loaded from files for comparison
            complex *savedpsiBfp;
            complex *savedpsiBf;
            complex *savedpsiAf;
            complex *savedlambdaBf;
            complex *savedlambdaAf;
        #endif
};

#endif // CP_STATE_GSPACE_PLANE_H


