
/** \file pairCalculator.h
 *
 */
 
#ifndef _pairCalculator_h
#define _pairCalculator_h
#include "ckPairCalculator.h"
#include "pcConfig.h"


/* delegated paircalc proxies perform like fermented dung on BG/L */
#ifdef CMK_BLUEGENEL
#define _PAIRCALC_DO_NOT_DELEGATE_ 1
#endif
// Do not use comlib for multicasts within paircalc
#define _PC_COMMLIB_MULTI_ 0

//============================================================================



/// A place to keep the section proxies for the reduction
class PairCalcID 
{
	public:
		/// The array ID of this PC chare array instance
		CkArrayID Aid;
		/// The array ID of the PC's input handler chare array
		CkArrayID ipHandlerID;
		///@note: This doesnt seem to be getting used. There are no references to Gid
		CkGroupID Gid;
		int GrainSize;
		int numChunks;
		int nstates;
		//@{
		///@todo: (RV) These are repeated in here and in PC config data. Understand how they are necessary in the message
		bool Symmetric;
		bool useComlib;
		bool useDirectSend;
		bool isDoublePacked;
		bool conserveMemory;
		bool lbpaircalc;
		//@}
		/// True if a proxy for the destination PC array section including a (portion of a) row exists
		bool existsLproxy;
		/// True if a proxy for the destination PC array section including a (portion of a) column exists
		bool existsRproxy;

		CkVec <CkGroupID> mCastGrpId;
		int priority;


		/** Array section which receives left matrix block data from the owner of this object (a Gspace chare)
		 * Symmetric loop : Includes the post-diagonal chares on row 's' that get data from this GSpace[s,p] chare
		 * Asymmetric loop: Includes all the chares on row 's' that get data from this GSpace[s,p] chare
		 */
		CProxySection_InputDataHandler<CollatorType,CollatorType> *sectionGettingLeft;
		/** Array section which receives right matrix block data from the owner of this object (a Gspace chare)
		 * Symmetric loop : Includes the pre-diagonal chares on column 's' that get data from this GSpace[s,p] chare
		 * Asymmetric loop: Includes all the chares on column 's' that get data from this GSpace[s,p] chare
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

		#ifdef _CP_SUBSTEP_TIMING_
		CkCallback beginTimerCB;
		CkCallback endTimerCB;
		int forwardTimerID;
		int backwardTimerID;
		#endif



		PairCalcID() {
		    sectionGettingLeft=NULL;
		    sectionGettingRight=NULL;
		}



		~PairCalcID() {
		  if(existsLproxy)
		  	delete [] sectionGettingLeft;
		  if(existsRproxy)
		  	delete [] sectionGettingRight;
		}
		


		void Init(CkArrayID aid, CkArrayID handlerID, int grain, int _numChunks, int s, bool sym, bool _useComlib,  bool _dp, bool _conserveMemory, bool _lbpaircalc, int _priority,  bool _useDirectSend) {
		  Aid = aid;
		  ipHandlerID = handlerID;
		  handlerProxy = CProxy_InputDataHandler<CollatorType,CollatorType> (handlerID);
		  GrainSize = grain;
		  numChunks = _numChunks;
		  nstates = s;
		  Symmetric = sym;
		  useComlib = _useComlib;
		  useDirectSend = _useDirectSend;
		  conserveMemory = _conserveMemory;
		  existsRproxy=false;
		  existsLproxy=false;
		  isDoublePacked = _dp;
		  lbpaircalc=_lbpaircalc;
		  priority=_priority;
		}



void resetProxy()
{
    CkAbort("need to adjust for having plane instance of multicastmgr");
    CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(mCastGrpId[0]).ckLocalBranch();       
    for(int chunk=0;chunk<numChunks;chunk++)
      {

        if(useComlib && _PC_COMMLIB_MULTI_)
          {
#ifdef USE_COMLIB
            /*if(existsRproxy)
            	ComlibResetSectionProxy(&sectionGettingRight[chunk]);
            if(existsLproxy)
            	ComlibResetSectionProxy(&sectionGettingLeft[chunk]);*/
#endif
          }
        else
          {
            if(existsRproxy)
            	mcastGrp->resetSection(sectionGettingRight[chunk]);
            if(existsLproxy)
            	mcastGrp->resetSection(sectionGettingLeft[chunk]);
          }
      }
}



PairCalcID &operator=(const PairCalcID& pid) {
  Aid=pid.Aid;
  ipHandlerID = pid.ipHandlerID;
  Gid=pid.Gid;    
  GrainSize=pid.GrainSize;
  numChunks=pid.numChunks;
  nstates=pid.nstates;
  Symmetric=pid.Symmetric;
  useComlib=pid.useComlib;
  useDirectSend=pid.useDirectSend;
  isDoublePacked=pid.isDoublePacked;
  conserveMemory=pid.conserveMemory;
  lbpaircalc=pid.lbpaircalc;
  existsLproxy=pid.existsLproxy;
  existsRproxy=pid.existsRproxy;
  priority=pid.priority;
  mCastGrpId=pid.mCastGrpId;
#ifdef _CP_SUBSTEP_TIMING_
    forwardTimerID=pid.forwardTimerID;
    backwardTimerID=pid.backwardTimerID;
    beginTimerCB=pid.beginTimerCB;
    endTimerCB=pid.endTimerCB;
#endif
    // everyone has to make their own proxies
    return *this;
  }



  void pup(PUP::er &p) {
    p|Aid;
    p|ipHandlerID;
    p|Gid;
    p|GrainSize;
    p|numChunks;
    p|nstates;
    p|Symmetric;
    p|useComlib;
    p|useDirectSend;
    p|isDoublePacked;
    p|conserveMemory;
    p|lbpaircalc;
    p|existsLproxy;
    p|existsRproxy;
    p|mCastGrpId;
    p|priority;
#ifdef _CP_SUBSTEP_TIMING_
    p|forwardTimerID;
    p|backwardTimerID;
    p|beginTimerCB;
    p|endTimerCB;
#endif
    if(p.isUnpacking())
      {
	if(existsLproxy)
	    sectionGettingLeft=new CProxySection_InputDataHandler<CollatorType,CollatorType>[numChunks];
	if(existsRproxy)
	    sectionGettingRight=new CProxySection_InputDataHandler<CollatorType,CollatorType>[numChunks];
      }
    if(existsLproxy)
      {
	if(useDirectSend)
	  p|handlerProxy;
	PUParray(p,sectionGettingLeft,numChunks);
	if(useDirectSend)
	  p|listGettingLeft;
      }
    if(existsRproxy)
      {
	PUParray(p,sectionGettingRight,numChunks);
	if(useDirectSend) {
	  p|listGettingRight; }
      }
  }

};

/// Creates the PC chare array. Called separately for the symm / asymm instances
void createPairCalculator(const cp::paircalc::pcConfig pcCfg, PairCalcID* aid, int flag, CkGroupID *mapid, int priority, CkVec <CkGroupID> mCastGrpId);


/// Forward declaration of the handshake token
struct RDMApair_GSP_PC;
/// Send out RDMA setup requests to all the destination PC chares that will be getting left data 
void sendLeftRDMARequest (PairCalcID *pid, RDMApair_GSP_PC idTkn, int totalsize, CkCallback cb);
/// Send out RDMA setup requests to all the destination PC chares that will be getting right data 
void sendRightRDMARequest(PairCalcID *pid, RDMApair_GSP_PC idTkn, int totalsize, CkCallback cb);
/// 
void isAtSyncPairCalc(PairCalcID* pcid);


void setResultProxy(CProxySection_PairCalculator *sectProxy,int state, int GrainSize,  CkGroupID mCastGrpId, bool lbsync, CkCallback synccb);

//@{
/// Matrix read/write utils
void dumpMatrixDouble(const char *infilename, double *matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric);
void loadMatrixDouble(const char *infilename, double *matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric);
void dumpMatrix2DDouble(const char *infilename, double **matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric);
void loadMatrix2DDouble(const char *infilename, double **matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric);
void dumpMatrix2DInt(const char *infilename, int **matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric);
void loadMatrix2DInt(const char *infilename, int **matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric);
//@}

//@{
///
bool reorder_elem_list(CkArrayIndexMax *elems, int numelems, int newstart);
bool reorder_elem_list_4D(CkArrayIndex4D *elems, int numelems, int newstart);
bool reorder_elem_list_max(CkArrayIndexMax *elems, int numelems, int newstart);
//@}

#endif
