
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
		/// True if a proxy for the destination PC array section including a (portion of a) row exists
		bool existsLproxy;
		/// True if a proxy for the destination PC array section including a (portion of a) column exists
		bool existsRproxy;

		CkGroupID mCastGrpId;


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
		


		void Init(CkArrayID aid, CkArrayID handlerID) {
		  handlerProxy = CProxy_InputDataHandler<CollatorType,CollatorType> (handlerID);
		  existsRproxy=false;
		  existsLproxy=false;
		}


PairCalcID &operator=(const PairCalcID& pid) {
  existsLproxy=pid.existsLproxy;
  existsRproxy=pid.existsRproxy;
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
    p|existsLproxy;
    p|existsRproxy;
    p|mCastGrpId;
#ifdef _CP_SUBSTEP_TIMING_
    p|forwardTimerID;
    p|backwardTimerID;
    p|beginTimerCB;
    p|endTimerCB;
#endif
  }

};


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
