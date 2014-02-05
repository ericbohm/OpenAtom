#include "debug_flags.h"
#include "charm++.h"
#include "ckcomplex.h"
#include "gStatePlane.decl.h"
#include "pcCommManager.h"
#include "uber/Uber.h"
#include "fft_slab_ctrl/fftCacheSlab.h"
struct EnergyStruct; /// @warning: Forward declarations of structs seem to choke the ci parser. It doesnt recognize the keyword struct.
#include "paircalc/ckPairCalculator.h"
#include "cpaimd.decl.h"
#include "utility/matrix2file.h"

using namespace cp::gspace; ///< @note: Should be temporary until GSpace chares live within namespace gspace

#ifndef CP_STATE_GSPACE_PLANE_H
#define CP_STATE_GSPACE_PLANE_H
/** @addtogroup GSpaceState
    @{
*/

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
    int senderJndex;
    int senderKndex;
    int numPlanes;
    complex *data;
};



/**
 * @brief 2D chare array [\ref nchareG][\ref nstates] Handles the
 * electronic structure in Fourier space (referred to as GSpace). Due
 * to spherical cutoff sparsity, GSpace is rearranged into
 * approximately equal size chunks \re nchareg. The number of chunks
 * is a free parameter selected at runtime.
 * A 2D chare array [\ref nchareG][\ref nstates] holding the
 * electronic structure in Fourier space (referred to as GSpace). Due
 * to spherical cutoff sparsity, GSpace is rearranged into
 * approximately equal size chunks \re nchareg. The number of chunks
 * is a free parameter selected at runtime.
 *
 * This is a description of the 'life' of a CP_State_GSpacePlane
 * object.  At the beginning of the program, the constructor
 * CP_State_GSpacePlane() is called, to initialize the
 * CP_State_GSpacePlane array. The GSpaceSlab within the
 * CP_State_GSpacePlane is initialized using the initGSpace(...)
 * method.
 *
 * To start off an iteration of program execution, the method doFFT()
 * is called. As part of this method, the CP_State_GSpacePlane object
 * does a forward 1-D fft on its slab of data, and sends off this data
 * to the next stage in the computational loop. After this, the
 * CP_State_GSpacePlane is idle, waiting for a message to trigger some
 * computation/communication.
 *
 * The idle period of the CP_State_GSpacePlane is terminated by a called to 
 * the method doIFFT(). In this method, the CP_State_GSpacePlane receives the 
 * partially processed data from the CP_State_RealSpacePlanes and performs 
 * 1-D inverse FFT on its slab of data. Then the calculation of the "S" 
 * matrix is started by  calling the sendPsi() method
 *
 * After the forces are calculated using the inverse FFT, a check is done 
 * to see if the forces from the particle calculations are ready. If so, the 
 * forces are added up. The sendPsi() method is not called until the forces
 * from the particle calculations and the forces from the quantum computation
 * are ready.
 *
 * The object is idle until the corrected g-space data from orthonormalization
 * is received through the acceptNewPsi() method.
 */

class CP_State_GSpacePlane: public CBase_CP_State_GSpacePlane 
{
    public:
		CP_State_GSpacePlane_SDAG_CODE
        friend class CP_State_ParticlePlane;
        // ----------- Flags and counters used in the GSpace driver code -----------
        /// My state index
		int istate_ind;
		/// My plane index
		int iplane_ind;
                int ibead_ind,kpoint_ind, itemper_ind,ispin_ind;
		///
		int iteration;
		///
		int iRecvRedPsi;
		///
		int cleanExitCalled;
		///
		bool doneDoingIFFT;
        /// A proxy for the my ortho chare array so I can interact with it
        CProxy_Ortho myOrtho;

        int halfStepEvolve;
        int redPlane;
        int registrationFlag;
        int nrotation;
        int exitFlag;
        int exitFlagMin;
        int outputFlag;
        int iRecvRedPsiV;
        int iSentRedPsi;
        int iSentRedPsiV;
        int finishedCpIntegrate;
        int numRecvRedPsi;
        int iterRotation;
	int myBeadIndex;
	int myKptIndex;
	int myTemperIndex;
	int mySpinIndex;

        double ake_old;
        
        bool acceptedVPsi;
        CP_State_GSpacePlane( int, int, int, int,int, int, UberCollection);
        CP_State_GSpacePlane(CkMigrateMessage *m);
        ~CP_State_GSpacePlane();
        /// Gets called from the PairCalc data receivers to confirm the setup of an RDMA link
        void completeRDMAhandshake(RDMASetupConfirmationMsg<RDMApair_GSP_PC> *msg);
        void pup(PUP::er &);
        void createSymPCmap(const int boxSize, PeListFactory getPeList, UberCollection thisInstance);
        void createAsymPCmap(const int boxSize, PeListFactory getPeList, UberCollection thisInstance);
        void acceptPairCalcAIDs(pcSetupMsg *msg);
        void initGSpace(int, complex *,int ,complex *,int,int,int,int,int,int,int);
        void launchAtoms();
        void launchOrthoT();
        void doFFT();
        void startNewIter ();
        void sendPsi();
        void sendPsiV();
        void screenOutputPsi(int);
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
        void unpackNewPsi(CkReductionMsg *msg);
        void unpackNewPsi(partialResultMsg  *msg);
        /// @RTH resume
        void doNewPsi();
        void collectFileOutput(GStateOutMsg *msg);
        void unpackNewPsiV(CkReductionMsg *msg);
        void unpackNewPsiV(partialResultMsg *msg);
        /// @RTH resume
        void doNewPsiV();
        /// @entry @RTH resume
        void psiCgOvlap(CkReductionMsg *msg);
        void unpackLambda(CkReductionMsg *msg);
        void unpackLambda(partialResultMsg *msg);
        /// @RTH resume
        void doLambda();
        /// @RTH resume
        void unpackRedPsi(GSRedPsiMsg *msg);
        
	void initBeadCookie(ICCookieMsg *m);
	void minimizeSync(ICCookieMsg *m);
        void computeCgOverlap();
        void sendFFTData ();
	void setExitFlag();
        void readFile();
        void computeEnergies(int p, double d);
        void startFFT(CkReductionMsg *msg);
        void sendRedPsiV();
        void unpackRedPsiV(GSRedPsiMsg *msg);
        void doneRedPsiVIntegrate();
        void screenPrintWallTimes();
	void acceptNewTemperature(double temp);
        const UberCollection thisInstance;        
	CkSectionInfo beadCookie;
    private:
	FILE *temperScreenFile;
        double *wallTimeArr;//only used on [0,0]
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
        /// Manages communication with the symmetric paircalc array
        PCCommManager symmPCmgr;
        /// Manages communication with the asymmetric paircalc array
        PCCommManager asymmPCmgr;
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
#if USE_PERSISTENT
        PersistentHandle   *fftHandler;
        void setupFFTPersistent();
#endif
};
/*@}*/
#endif // CP_STATE_GSPACE_PLANE_H


