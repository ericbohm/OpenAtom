#include "debug_flags.h"
#include "paircalc/pcConfig.h"

#ifndef PC_COMM_MANAGER_H
#define PC_COMM_MANAGER_H

// Forward declarations
class PairCalcID;
class CProxySection_PairCalculator;
struct ckcomplex;
typedef ckcomplex complex;
/// Forward declaration of the handshake token
struct RDMApair_GSP_PC;

namespace cp {
    /// Namespace aliases
    namespace pc = paircalc;

    namespace gspace {

///
class PCCommManager
{
    public:
        /// Constructor
        PCCommManager(const pc::pcConfig &_cfg): pcCfg(_cfg) {}
        PCCommManager() {} ///< @warning: Just to appease charm migration constructors. pffouggh...
        /// Creates multicast trees to the appropriate PC chare array sections used in the symmetric / asymmetric loops
        static void createPCarray(const pc::pcConfig pcCfg, PairCalcID* pcid, CkGroupID *mapid);
        void makeLeftTree(PairCalcID* pid, int myS, int myZ);
        /// Creates a multicast tree that includes the PC chare arrays used in the asymmetric loop
        void makeRightTree(PairCalcID* pid, int myS, int myZ);
        /// Starts the forward path work (Psi, Lambda and PsiV cases) by multicasting an entry method call to the appropriate PC chare array section
        void sendLeftData (PairCalcID* aid, int n, complex* ptr, int myS, int myZ, bool psiV);
        /// Starts the forward path work (along with startPairCalcLeft()) in the asymmetric (Lambda) case
        void sendRightData(PairCalcID* aid, int n, complex* ptr, int myS, int myZ, bool psiV);

        //@{
        /// Initialize an array section that is used to reduce the results from the PCs back to the GSP chares
        CProxySection_PairCalculator makeOneResultSection_asym(PairCalcID* pcid, int state, int plane, int chunk);
        CProxySection_PairCalculator makeOneResultSection_asym_column(PairCalcID* pcid, int state, int plane, int chunk);
        CProxySection_PairCalculator makeOneResultSection_sym1(PairCalcID* pcid, int state, int plane, int chunk);
        CProxySection_PairCalculator makeOneResultSection_sym2(PairCalcID* pcid, int state, int plane, int chunk);
        //@}

    private:
        /// Multicasts the left matrix data to the PC section
        void sendLeftDataMcast (PairCalcID* aid, int n, complex* ptr, int myS, int myZ, bool psiV);
        /// Multicasts the right matrix data to the PC section
        void sendRightDataMcast(PairCalcID* aid, int n, complex* ptr, int myS, int myZ, bool psiV);
        /// Sends left matrix data via RDMA
        void sendLeftDataRDMA  (PairCalcID* aid, int n, complex* ptr, int myS, int myZ, bool psiV);
        /// Sends right matrix data via RDMA
        void sendRightDataRDMA (PairCalcID* aid, int n, complex* ptr, int myS, int myZ, bool psiV);
        /// Send RDMA setup requests to all the destination PC chares that will be getting left data
        void sendLeftRDMARequest (PairCalcID *pid, RDMApair_GSP_PC idTkn, int totalsize, CkCallback cb);
        /// Send RDMA setup requests to all the destination PC chares that will be getting right data
        void sendRightRDMARequest(PairCalcID *pid, RDMApair_GSP_PC idTkn, int totalsize, CkCallback cb);
        /// Send out a dummy mcast to prod CkMulticast into setting up the result reduction trees etc
        void setResultProxy(CProxySection_PairCalculator *sectProxy,int state, int GrainSize,  CkGroupID mCastGrpId, bool lbsync, CkCallback synccb);


        /// Input configurations for the paircalcs
        cp::paircalc::pcConfig pcCfg;
};




inline void PCCommManager::sendLeftData(PairCalcID* pcid, int n, complex* ptr, int myS, int myPlane, bool psiV)
{
    #ifdef PC_USE_RDMA
        sendLeftDataRDMA(pcid,n,ptr,myS,myPlane,psiV);
    #else
		sendLeftDataMcast(pcid,n,ptr,myS,myPlane,psiV);
    #endif
}




inline void PCCommManager::sendRightData(PairCalcID* pcid, int n, complex* ptr, int myS, int myPlane, bool psiV)
{
    #ifdef PC_USE_RDMA
        sendRightDataRDMA(pcid,n,ptr,myS,myPlane,psiV);
    #else
		sendRightDataMcast(pcid,n,ptr,myS,myPlane,psiV);
    #endif
}

    } // end namespace gspace
} // end namespace cp
#endif // PC_COMM_MANAGER_H
