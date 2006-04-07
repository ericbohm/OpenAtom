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
  int S; 
  bool Symmetric;
  bool useComlib;
  bool isDoublePacked;
  bool conserveMemory;
  bool lbpaircalc;
  bool existsLproxy;
  bool existsLNotFromproxy;
  bool existsRproxy;
  CkGroupID mCastGrpId;
  int priority;
  CProxySection_PairCalculator *proxyLFrom;
  CProxySection_PairCalculator *proxyLNotFrom;
  CProxySection_PairCalculator *proxyRNotFrom;

  PairCalcID() {}
  ~PairCalcID() {}

  void Init(CkArrayID aid, int grain, int _numChunks, int s, bool sym, bool _useComlib,  bool _dp, bool _conserveMemory, bool _lbpaircalc, CkGroupID _mCastGrpId, int _priority) {
    Aid = aid;
    GrainSize = grain;
    numChunks = _numChunks;
    S = s;
    Symmetric = sym;
    useComlib = _useComlib;
    conserveMemory = _conserveMemory;
    existsRproxy=false;
    existsLproxy=false;
    existsLNotFromproxy=false;
    isDoublePacked = _dp;
    lbpaircalc=_lbpaircalc;
    mCastGrpId=_mCastGrpId;
    priority=_priority;
  }
  void resetProxy()
    {
      CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch();       
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
    S=pid.S;
    Symmetric=pid.Symmetric;
    useComlib=pid.useComlib;
    isDoublePacked=pid.isDoublePacked;
    conserveMemory=pid.conserveMemory;
    lbpaircalc=pid.lbpaircalc;
    existsLproxy=pid.existsLproxy;
    existsLNotFromproxy=pid.existsLNotFromproxy;
    existsRproxy=pid.existsRproxy;
    priority=pid.priority;
    mCastGrpId=pid.mCastGrpId;
    for(int chunk;chunk<numChunks;chunk++)
      {
	proxyLFrom[chunk]=pid.proxyLFrom[chunk];
	proxyLNotFrom[chunk]=pid.proxyLNotFrom[chunk];
	proxyRNotFrom[chunk]=pid.proxyRNotFrom[chunk];
      }
    return *this;
  }

  void pup(PUP::er &p) {
    p|Aid;
    p|Gid;
    p|GrainSize;
    p|numChunks;
    p|S;
    p|Symmetric;
    p|useComlib;
    p|isDoublePacked;
    p|conserveMemory;
    p|lbpaircalc;
    p|existsLproxy;
    p|existsLNotFromproxy;
    p|existsRproxy;
    p|mCastGrpId;
    p|priority;
    if(existsLproxy)
      PUParray(p,proxyLFrom,numChunks);
    if(existsLNotFromproxy)
      PUParray(p,proxyLNotFrom,numChunks);
    if(existsRproxy)
      PUParray(p,proxyRNotFrom,numChunks);
  }

};

void createPairCalculator(bool sym, int w, int grainSize, int numZ, int* z,  CkCallback cb, PairCalcID* aid, int ep, int ep2, CkArrayID cbid, int flag, CkGroupID *mapid, int flag_dp, bool conserveMemory, bool lbpaircalc, int priority, CkGroupID mCastGrpId, int numChunks);

void startPairCalcLeft(PairCalcID* aid, int n, complex* ptr, int myS, int myZ, bool psiV);

void startPairCalcRight(PairCalcID* aid, int n, complex* ptr, int myS, int myZ);
void makeLeftTree(PairCalcID* pid, int myS, int myZ);

void makeRightTree(PairCalcID* pid, int myS, int myZ);

extern "C" void finishPairCalcSection(int n, double *ptr,CProxySection_PairCalculator sectionProxy, int actionType);

extern "C" void finishPairCalcSection2( int n, double *ptr1, double *ptr2,CProxySection_PairCalculator sectionProxy, int actionType);

CProxySection_PairCalculator initOneRedSect( int numZ, int* z, int blkSize,  PairCalcID* pcid, CkCallback cb, int s1, int s2);

void startPairCalcLeftAndFinish(PairCalcID* pcid, int n, complex* ptr, int myS, int myZ);

void startPairCalcRightAndFinish(PairCalcID* pcid, int n, complex* ptr, int myS, int myZ);

void isAtSyncPairCalc(PairCalcID* pcid);

/* These are the classic no multicast version for comparison and debugging */
void startPairCalcLeftSlow(PairCalcID* aid, int n, complex* ptr, int myS, int myZ);

void startPairCalcRightSlow(PairCalcID* aid, int n, complex* ptr, int myS, int myZ);

CProxySection_PairCalculator makeOneResultSection_asym(PairCalcID* pcid, int state, int plane, int chunk);
CProxySection_PairCalculator makeOneResultSection_asym_column(PairCalcID* pcid, int state, int plane, int chunk);
CProxySection_PairCalculator makeOneResultSection_sym1(PairCalcID* pcid, int state, int plane, int chunk);
CProxySection_PairCalculator makeOneResultSection_sym2(PairCalcID* pcid, int state, int plane, int chunk);
void setGredProxy(CProxySection_PairCalculator *sectProxy, CkGroupID mCastGrpId, CkCallback cb, bool lbsync, CkCallback synccb);
void setResultProxy(CProxySection_PairCalculator *sectProxy,int state, int GrainSize,  CkGroupID mCastGrpId, bool lbsync, CkCallback synccb);

void dumpMatrixDouble(const char *infilename, double *matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric);

#endif
