#include "debug_flags.h"
#include "paircalc/pcConfig.h"
#include "paircalc/pcFwdDeclarations.h"
#include "inputDataHandler.decl.h"

#ifndef PC_COMM_MANAGER_H
#define PC_COMM_MANAGER_H

// Forward declarations
class CP_State_GSpacePlane;
struct ckcomplex;
typedef ckcomplex complex;

namespace cp {
    namespace gspace {

///
class PCCommManager
{
    friend class ::CP_State_GSpacePlane; ///< @note: Temporary until paircalc startup moves completely to GSpace

    public:
        /// Constructor
        PCCommManager(const pc::pcConfig &_cfg);
        PCCommManager() {} ///< @warning: Just to appease charm migration constructors. pffouggh...
        /// Create a paircalc array using info in the supplied pcConfig object. Originally createPairCalculator()
        void createPCarray(PairCalcID* pcid, CkGroupID *mapid);
        /// Creates multicast trees to the appropriate PC chare array sections used in the symmetric / asymmetric loops
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
        void setResultProxy(CProxySection_PairCalculator *sectProxy,int state, int GrainSize,  bool lbsync, CkCallback synccb);


        /// Input configurations for the paircalcs
        cp::paircalc::pcConfig pcCfg;
		/// The array ID of the PC chare array instance I am managing
		CkArrayID pcAID;
		/// The array ID of the PC's input handler chare array
		CkArrayID ipHandlerAID;
        /// The group ID of the multicast manager that will handle the multicasts to the PC array
        CkGroupID mCastMgrGID;

        /** Array section which receives left matrix block data
         *
         * symm instance: post-diagonal chares on row 's' that get data from this GSpace[s,p] chare
         * asymm instance: all chares on row 's' that get data from this GSpace[s,p] chare
         */
        CProxySection_InputDataHandler<CollatorType,CollatorType> *sectionGettingLeft;

        /** Array section which receives right matrix block data
         *
         * symm instance: pre-diagonal chares on column 's' that get data from this GSpace[s,p] chare
         * asymm instance: all chares on column 's' that get data from this GSpace[s,p] chare
         */
        CProxySection_InputDataHandler<CollatorType,CollatorType> *sectionGettingRight;

        /// A proxy to the PC input handler chare array
        CProxy_InputDataHandler<CollatorType,CollatorType> handlerProxy;
        /// A list of PC array elements which expect left matrix data from owning GSpace chare
        CkVec <CkArrayIndex4D> listGettingLeft;
        /// A list of PC array elements which expect right matrix data from owning GSpace chare
        CkVec <CkArrayIndex4D> listGettingRight;
        /// RDMA handles for each PC chare's input data handler that will receive data from the owner of this object (a GSpace[s,p] chare)
        CkVec<rdmaHandleType> leftDestinationHandles, rightDestinationHandles;
        /// True if a proxy for the destination PC array section including a (portion of a) row exists
        bool existsLproxy;
        /// True if a proxy for the destination PC array section including a (portion of a) column exists
        bool existsRproxy;
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
