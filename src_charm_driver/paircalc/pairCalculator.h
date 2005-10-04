
#ifndef _pairCalculator_h
#define _pairCalculator_h
#include "ckPairCalculator.h"


/* delegated paircalc proxies perform like fermented dung on BG/L */
#ifdef CMK_VERSION_BLUEGENE
#define _PAIRCALC_DO_NOT_DELEGATE_ 1
#endif

/* a place to keep the section proxies for the reduction */

class PairCalcID {
 public:
  CkArrayID Aid;
  CkGroupID Gid;
  int GrainSize;
  int BlkSize;
  int S; 
  bool Symmetric;

  ComlibInstanceHandle cinst;
  ComlibInstanceHandle minst;
  ComlibInstanceHandle rinst;

  bool useComlib;
  bool isDoublePacked;
  bool conserveMemory;
  bool lbpaircalc;
  bool existsLproxy;
  bool existsLNotFromproxy;
  bool existsRproxy;
  int priority;
  int resultcount;
  CProxySection_PairCalculator proxyLFrom;
  CProxySection_PairCalculator proxyLNotFrom;
  CProxySection_PairCalculator proxyRNotFrom;
  CkVec <CProxySection_PairCalculator> *resultProxies;
  CkGroupID mCastGrpId;
  PairCalcID() {}
  ~PairCalcID() {}

  void Init(CkArrayID aid, CkGroupID gid, int grain, int blk, int s, bool sym, bool _useComlib,  ComlibInstanceHandle h, bool _dp, bool _conserveMemory, bool _lbpaircalc, CkGroupID _mCastGrpId, int _priority) {
    resultcount=0;
    Aid = aid;
    Gid = gid;
    GrainSize = grain;
    BlkSize = blk;
    S = s;
    Symmetric = sym;
    useComlib = _useComlib;
    conserveMemory = _conserveMemory;
    cinst = h;
    existsRproxy=false;
    existsLproxy=false;
    existsLNotFromproxy=false;
    isDoublePacked = _dp;
    lbpaircalc=_lbpaircalc;
    mCastGrpId=_mCastGrpId;
    priority=_priority;
    resultProxies=NULL;
  }
  void resetProxy()
    {
      CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch();       

      if(useComlib && _PC_COMMLIB_MULTI_)
	{
	  if(existsLNotFromproxy)
	    ComlibResetSectionProxy(&proxyLNotFrom);
	  if(existsRproxy)
	    ComlibResetSectionProxy(&proxyRNotFrom);
	  if(existsLproxy)
	    ComlibResetSectionProxy(&proxyLFrom);
	}
      else
	{
	  if(existsRproxy)
	    {
	      mcastGrp->resetSection(proxyRNotFrom);
	    }
	  if(existsLproxy)
	    {
	      mcastGrp->resetSection(proxyLFrom);
	    }
	  if(existsLNotFromproxy)
	    {
	      mcastGrp->resetSection(proxyLNotFrom);
	    }
	}
      for(int i=0;i<resultProxies->size();i++)
	mcastGrp->resetSection((*resultProxies)[i]);
    }
  PairCalcID &operator=(const PairCalcID& pid) {
    Aid=pid.Aid;
    Gid=pid.Gid;
    GrainSize=pid.GrainSize;
    BlkSize=pid.BlkSize;
    S=pid.S;
    Symmetric=pid.Symmetric;
    cinst=pid.cinst;
    minst=pid.minst;
    rinst=pid.rinst;
    useComlib=pid.useComlib;
    isDoublePacked=pid.isDoublePacked;
    conserveMemory=pid.conserveMemory;
    lbpaircalc=pid.lbpaircalc;
    existsLproxy=pid.existsLproxy;
    existsLNotFromproxy=pid.existsLNotFromproxy;
    existsRproxy=pid.existsRproxy;
    proxyLFrom=pid.proxyLFrom;
    proxyLNotFrom=pid.proxyLNotFrom;
    proxyRNotFrom=pid.proxyRNotFrom;
    mCastGrpId=pid.mCastGrpId;
    priority=pid.priority;
    resultcount=pid.resultcount;
    if(lbpaircalc && resultcount>0)
      {
	resultProxies= new CkVec <CProxySection_PairCalculator> (resultcount);
	memcpy(resultProxies->getVec(),pid.resultProxies->getVec(),resultcount*sizeof(CProxySection_PairCalculator));
      }
    else
      {
	resultProxies=NULL;
      }
    return *this;
  }

  void pup(PUP::er &p) {
    p|Aid;
    p|Gid;
    p|GrainSize;
    p|BlkSize;
    p|S;
    p|Symmetric;
    p|cinst;
    p|minst;
    p|rinst;
    p|useComlib;
    p|isDoublePacked;
    p|conserveMemory;
    p|lbpaircalc;
    p|existsLproxy;
    p|existsLNotFromproxy;
    p|existsRproxy;
    p|proxyLFrom;
    p|proxyLNotFrom;
    p|proxyRNotFrom;
    p|mCastGrpId;
    p|priority;
    p|resultcount;
    if(p.isUnpacking())
      {
	if(lbpaircalc)
	  resultProxies= new CkVec <CProxySection_PairCalculator> (resultcount);
      }
    if(lbpaircalc)
      p(resultProxies->getVec(),resultcount);

  }

};

void createPairCalculator(bool sym, int w, int grainSize, int numZ, int* z, int op1, FuncType f1, int op2, FuncType f2, CkCallback cb, PairCalcID* aid, int ep, CkArrayID cbid, int flag, CkGroupID *gid, int flag_dp, bool conserveMemory, bool lbpaircalc, int priority, CkGroupID mCastGrpId);

void startPairCalcLeft(PairCalcID* aid, int n, complex* ptr, int myS, int myZ);

void startPairCalcRight(PairCalcID* aid, int n, complex* ptr, int myS, int myZ);
void makeLeftTree(PairCalcID* pid, int myS, int myZ);

void makeRightTree(PairCalcID* pid, int myS, int myZ);

extern "C" void finishPairCalc(PairCalcID* aid, int n, double *ptr);

extern "C" void finishPairCalc2(PairCalcID* pcid, int n, double *ptr1, double *ptr2);

extern "C" void finishPairCalcSection(int n, double *ptr,CProxySection_PairCalculator sectionProxy);

extern "C" void finishPairCalcSection2( int n, double *ptr1, double *ptr2,CProxySection_PairCalculator sectionProxy);

CProxySection_PairCalculator initOneRedSect( int numZ, int* z, int blkSize,  PairCalcID* pcid, CkCallback cb, int s1, int s2, int c);

void startPairCalcLeftAndFinish(PairCalcID* pcid, int n, complex* ptr, int myS, int myZ);

void startPairCalcRightAndFinish(PairCalcID* pcid, int n, complex* ptr, int myS, int myZ);

void isAtSyncPairCalc(PairCalcID* pcid);

/* These are the classic no multicast version for comparison and debugging */
void startPairCalcLeftSlow(PairCalcID* aid, int n, complex* ptr, int myS, int myZ);

void startPairCalcRightSlow(PairCalcID* aid, int n, complex* ptr, int myS, int myZ);

CProxySection_PairCalculator makeOneResultSection_asym(PairCalcID* pcid, int state, int plane);
CProxySection_PairCalculator makeOneResultSection_sym1(PairCalcID* pcid, int state, int plane);
CProxySection_PairCalculator makeOneResultSection_sym2(PairCalcID* pcid, int state, int plane);
#endif
