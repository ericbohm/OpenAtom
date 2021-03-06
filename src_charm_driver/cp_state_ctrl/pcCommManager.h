#include "debug_flags.h"
#include "paircalc/pcConfig.h"
#include "paircalc/pcInstanceIDs.h"
#include "paircalc/pcFwdDeclarations.h"

#ifndef PC_COMM_MANAGER_H
#define PC_COMM_MANAGER_H

// Forward declarations
class CP_State_GSpacePlane;
struct ckcomplex;
typedef ckcomplex complex;
/** @addtogroup GSpaceState
  @{
 */

namespace cp {
  namespace gspace {

    /**
     * Manages communication with a single paircalc array
     *
     * It needs a handle to an already created paircalc instance and the config settings for that instance
     * Given this info, and the index of its owner (GSpace) chare, it creates appropriate mcast/redn sections
     * and manages these communications. GSpace simply delegates the data sends to the appropriate comm manager
     * instances in every iteration.
     *
     * @todo: This class also manages the RDMA setup requests, although it could also assume the responsibility
     * for processing the ack msg that arrives in GSpace to complete an RDMA handshake.
     */
    class PCCommManager
    {
      friend class ::CP_State_GSpacePlane; ///< @note: Temporary until paircalc startup moves completely to GSpace

      public:
      /// Constructor
      PCCommManager(const CkIndex2D gspaceIdx, const pc::pcConfig &_cfg, const pc::InstanceIDs _pcHandle);
      PCCommManager() {} ///< @warning: Just to appease charm migration constructors. pffouggh...
      /// Starts the forward path work (Psi, Lambda and PsiV cases) by multicasting an entry method call to the appropriate PC chare array section
      void sendLeftData (int numPoints, complex* ptr, bool psiV);
      /// Starts the forward path work (along with startPairCalcLeft()) in the asymmetric (Lambda) case
      void sendRightData(int numPoints, complex* ptr, bool psiV);

      //@{
      /// Initialize an array section that is used to reduce the results from the PCs back to the GSP chares
      CProxySection_PairCalculator makeOneResultSection_asym(int chunk);
      CProxySection_PairCalculator makeOneResultSection_asym_column(int chunk);
      CProxySection_PairCalculator makeOneResultSection_sym1(int chunk);
      CProxySection_PairCalculator makeOneResultSection_sym2(int chunk);
      //@}


      private:
      /// Creates multicast trees to the appropriate PC chare array sections used in the symmetric / asymmetric loops
      void makeLeftTree();
      /// Creates a multicast tree that includes the PC chare arrays used in the asymmetric loop
      void makeRightTree();
      /// Multicasts the left matrix data to the PC section
      void sendLeftDataMcast (int numPoints, complex* ptr, bool psiV);
      /// Multicasts the right matrix data to the PC section
      void sendRightDataMcast(int numPoints, complex* ptr, bool psiV);
      /// Sends left matrix data via RDMA
      void sendLeftDataRDMA  (int numPoints, complex* ptr, bool psiV);
      /// Sends right matrix data via RDMA
      void sendRightDataRDMA (int numPoints, complex* ptr, bool psiV);
      /// Send RDMA setup requests to all the destination PC chares that will be getting left data
      void sendLeftRDMARequest (RDMApair_GSP_PC idTkn, int totalsize, CkCallback cb);
      /// Send RDMA setup requests to all the destination PC chares that will be getting right data
      void sendRightRDMARequest(RDMApair_GSP_PC idTkn, int totalsize, CkCallback cb);
      /// Send out a dummy mcast to prod CkMulticast into setting up the result reduction trees etc
      void setResultProxy(CProxySection_PairCalculator *sectProxy, bool lbsync, CkCallback synccb);


      /// The array index of the owner GSpace chare
      CkIndex2D gspaceIndex;
      /// Input configurations for the paircalcs
      cp::paircalc::pcConfig pcCfg;
      /// Handles to the paircalc array and related entities that I will be managing comm with
      cp::paircalc::InstanceIDs pcHandle;

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




    inline void PCCommManager::sendLeftData(int n, complex* ptr, bool psiV)
    {
#ifdef PC_USE_RDMA
      sendLeftDataRDMA(n,ptr,psiV);
#else
      sendLeftDataMcast(n,ptr,psiV);
#endif
    }




    inline void PCCommManager::sendRightData(int n, complex* ptr, bool psiV)
    {
#ifdef PC_USE_RDMA
      sendRightDataRDMA(n,ptr,psiV);
#else
      sendRightDataMcast(n,ptr,psiV);
#endif
    }

  } // end namespace gspace
} // end namespace cp
/*@}*/
#endif // PC_COMM_MANAGER_H

