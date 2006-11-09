
/** \file pairCalculator.h
 *
 */
 
#ifndef _pairCalculator_h
#define _pairCalculator_h
#include "ckPairCalculator.h"


/* delegated paircalc proxies perform like fermented dung on BG/L */
/*#ifdef CMK_VERSION_BLUEGENE
#define _PAIRCALC_DO_NOT_DELEGATE_ 1
#endif
*/
/* a place to keep the section proxies for the reduction */

class PairCalcID {
 public:
  CkArrayID Aid;
  CkGroupID Gid;
  int GrainSize;
  int numChunks;
  int nstates; 
  bool Symmetric;
  bool useComlib;
  bool useEtoM;
  bool useDirectSend;
  bool isDoublePacked;
  bool conserveMemory;
  bool lbpaircalc;
  bool existsLproxy;
  bool existsLNotFromproxy;
  bool existsRproxy;
  CkVec <CkGroupID> mCastGrpId;
  CkGroupID orthomCastGrpId;
  CkGroupID orthoRedGrpId;
  int priority;
  CProxySection_PairCalculator proxySym;
  CProxySection_PairCalculator proxyAsym;
  CProxySection_PairCalculator *proxyLFrom;
  CProxySection_PairCalculator *proxyLNotFrom;
  CProxySection_PairCalculator *proxyRNotFrom;
  CProxy_PairCalculator cproxy;
  CkVec <CkArrayIndex4D> listLFrom;
  CkVec <CkArrayIndex4D> listLNotFrom;
  CkVec <CkArrayIndex4D> listRNotFrom;
  PairCalcID() {}
  ~PairCalcID() {
    if(existsLproxy)
      delete [] proxyLFrom;
    if(existsLNotFromproxy)
      delete [] proxyLNotFrom;
    if(existsRproxy)
      delete [] proxyRNotFrom;
  }

  void Init(CkArrayID aid, int grain, int _numChunks, int s, bool sym, bool _useComlib,  bool _dp, bool _conserveMemory, bool _lbpaircalc, int _priority, bool _useEtoM, bool _useDirectSend) {
    Aid = aid;
    GrainSize = grain;
    numChunks = _numChunks;
    nstates = s;
    Symmetric = sym;
    useComlib = _useComlib;
    useEtoM = _useEtoM;
    useDirectSend = _useDirectSend;
    conserveMemory = _conserveMemory;
    existsRproxy=false;
    existsLproxy=false;
    existsLNotFromproxy=false;
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
	      if(existsLNotFromproxy)
		ComlibResetSectionProxy(&proxyLNotFrom[chunk]);
	      if(existsRproxy)
		ComlibResetSectionProxy(&proxyRNotFrom[chunk]);
	      if(existsLproxy)
		ComlibResetSectionProxy(&proxyLFrom[chunk]);
	    }
	  else
	    {
	      if(existsRproxy)
		{
		  mcastGrp->resetSection(proxyRNotFrom[chunk]);
		}
	      if(existsLproxy)
		{
		  mcastGrp->resetSection(proxyLFrom[chunk]);
		}
	      if(existsLNotFromproxy)
		{
		  mcastGrp->resetSection(proxyLNotFrom[chunk]);
		}
	    }
	}
    }
  PairCalcID &operator=(const PairCalcID& pid) {
    Aid=pid.Aid;
    Gid=pid.Gid;    
    GrainSize=pid.GrainSize;
    numChunks=pid.numChunks;
    nstates=pid.nstates;
    Symmetric=pid.Symmetric;
    useComlib=pid.useComlib;
    useEtoM=pid.useEtoM;
    useDirectSend=pid.useDirectSend;
    isDoublePacked=pid.isDoublePacked;
    conserveMemory=pid.conserveMemory;
    lbpaircalc=pid.lbpaircalc;
    existsLproxy=pid.existsLproxy;
    existsLNotFromproxy=pid.existsLNotFromproxy;
    existsRproxy=pid.existsRproxy;
    priority=pid.priority;
    mCastGrpId=pid.mCastGrpId;
    orthomCastGrpId=pid.orthomCastGrpId;
    orthoRedGrpId=pid.orthoRedGrpId;
    cproxy=pid.cproxy;
    // everyone has to make their own proxies
    return *this;
  }

  void pup(PUP::er &p) {
    p|Aid;
    p|Gid;
    p|GrainSize;
    p|numChunks;
    p|nstates;
    p|Symmetric;
    p|useComlib;
    p|useEtoM;
    p|useDirectSend;
    p|isDoublePacked;
    p|conserveMemory;
    p|lbpaircalc;
    p|existsLproxy;
    p|existsLNotFromproxy;
    p|existsRproxy;
    p|mCastGrpId;
    p|orthomCastGrpId;
    p|orthoRedGrpId;
    p|priority;
    if(p.isUnpacking())
      {
	if(existsLproxy)
	  {
	    proxyLFrom=new CProxySection_PairCalculator[numChunks];
	  }
	if(existsLNotFromproxy)
	  {
	    proxyLNotFrom=new CProxySection_PairCalculator[numChunks];
	  }
	if(existsRproxy)
	  {
	    proxyRNotFrom=new CProxySection_PairCalculator[numChunks];
	  }
      }
    if(existsLproxy)
      {
	if(useEtoM)
	  p|cproxy;
	PUParray(p,proxyLFrom,numChunks);
	if(useEtoM)
	  p|listLFrom;
      }
    if(existsLNotFromproxy)
      {
	PUParray(p,proxyLNotFrom,numChunks);
	if(useEtoM)
	  p|listLNotFrom;
      }
    if(existsRproxy)
      {
	PUParray(p,proxyRNotFrom,numChunks);
	if(useEtoM)
	  p|listRNotFrom;
      }
  }

};

void createPairCalculator(bool sym, int w, int grainSize, int numZ, int* z,  CkCallback cb, PairCalcID* aid, int ep, int ep2, CkArrayID cbid, int flag, CkGroupID *mapid, int flag_dp, bool conserveMemory, bool lbpaircalc, int priority, CkVec <CkGroupID> mCastGrpId, CkGroupID orthomcastgrpid, CkGroupID orthoredgrpid, int numChunks, int orthoGrainSize, int usePairEtoM, bool collectTiles, bool streamBWout, bool delayBWSend, int streamFW, bool useDirectSend, bool gSpaceSum, int gpriority, bool phantomSym, bool useBWBarrier, int gemmSplitFWk, int gemmSplitFWm, int gemmSplitBW);

void startPairCalcLeft(PairCalcID* aid, int n, complex* ptr, int myS, int myZ, bool psiV);

void startPairCalcRight(PairCalcID* aid, int n, complex* ptr, int myS, int myZ);
void makeLeftTree(PairCalcID* pid, int myS, int myZ);

void makeRightTree(PairCalcID* pid, int myS, int myZ);

extern "C" void finishPairCalcSection(int n, double *ptr, PairCalcID *pcid, int orthoX, int orthoY, int actionType, int priority);

extern "C" void finishPairCalcSection2( int n, double *ptr1, double *ptr2, PairCalcID *pcid, int orthoX, int orthoY, int actionType, int priority);

//CProxySection_PairCalculator initOneRedSect( int numZ, int* z, int blkSize,  PairCalcID* pcid, CkCallback cb, int s1, int s2, int o1, int o2, int ograin, bool phantom, bool direct, bool commlib);
void initOneRedSect( int numZ, int* z, int blkSize,  PairCalcID* pcid, CkCallback cb, int s1, int s2, int o1, int o2, int ograin, bool phantom, bool direct, bool commlib);

//void startPairCalcLeftAndFinish(PairCalcID* pcid, int n, complex* ptr, int myS, int myZ);

//void startPairCalcRightAndFinish(PairCalcID* pcid, int n, complex* ptr, int myS, int myZ);

void isAtSyncPairCalc(PairCalcID* pcid);

/* These are the classic no multicast version for comparison and debugging */
void startPairCalcLeftSlow(PairCalcID* aid, int n, complex* ptr, int myS, int myZ);

void startPairCalcRightSlow(PairCalcID* aid, int n, complex* ptr, int myS, int myZ);

CProxySection_PairCalculator makeOneResultSection_asym(PairCalcID* pcid, int state, int plane, int chunk);
CProxySection_PairCalculator makeOneResultSection_asym_column(PairCalcID* pcid, int state, int plane, int chunk);
CProxySection_PairCalculator makeOneResultSection_sym1(PairCalcID* pcid, int state, int plane, int chunk);
CProxySection_PairCalculator makeOneResultSection_sym2(PairCalcID* pcid, int state, int plane, int chunk);
void setGredProxy(CProxySection_PairCalculator *sectProxy, CkGroupID mCastGrpId, CkCallback cb, bool lbsync, CkCallback synccb, int orthoX, int orthoY);
void setResultProxy(CProxySection_PairCalculator *sectProxy,int state, int GrainSize,  CkGroupID mCastGrpId, bool lbsync, CkCallback synccb);

void dumpMatrixDouble(const char *infilename, double *matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric);
bool reorder_elem_list(CkArrayIndexMax *elems, int numelems, int newstart);
#endif
