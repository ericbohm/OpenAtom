//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \file cpaimd.h
 * Some basic data structures and the array map classes are defined
 * here. These maps are used to map the array elements to the correct
 * processors.  Please read ../doc/README for a detailed
 * documentation on the structure of the program.
 */
//============================================================================
//#define _NAN_CHECK_

#include "debug_flags.h"
#ifndef _CPAIMD_H
#define _CPAIMD_H
//#define MAP_DEBUG 1
#include "CPcharmParaInfoGrp.h"
#include "load_balance/PeList.h"
#include "uber/Uber.h"
#include "ckhashtable.h"

#ifdef CMK_BALANCED_INJECTION_API
#include "ckBIconfig.h"
#endif

#include "PlatformSpecific.decl.h"

class PlatformSpecific : public CBase_PlatformSpecific
{
 public:
  PlatformSpecific(){}
  void reset_BI();
};

#define USE_INT_MAP
#include "load_balance/IntMap.h"
#include "load_balance/MapTable.h"

//#define BARRIER_CP_GSPACE_PSI 1

#define LOAD_BALANCE_STEP 100000000

#define PRE_BALANCE_STEP 2

#define FIRST_BALANCE_STEP 100000000

#if CMK_TRACE_ENABLED
#define TRACE_ON_STEP 4
#define TRACE_OFF_STEP 6
#endif
#ifdef CMK_BLUEGENEP
#define HPM_ON_STEP 4
#define HPM_OFF_STEP 5
#endif
#ifndef CmiMemcpy
#define CmiMemcpy(dest, src, size) CmiMemcpy((dest), (src), (size))
#endif

extern bool fakeTorus;
extern CkVec <int> PIBImaptable;
extern CkVec <MapType1> EnergyCommMgrImaptable;
extern CkVec <MapType1> AtomImaptable;
extern CkVec <MapType2> GSImaptable;
extern CkVec <MapType2> RSImaptable;
extern CkVec <MapType2> RPPImaptable;
extern CkVec <MapType1> RhoGSImaptable;
extern CkVec <MapType2> RhoRSImaptable;
extern CkVec <MapType2> RhoGHartImaptable;
extern CkVec <MapType3> RhoRHartImaptable;
extern CkVec < CkVec <MapType2> > RhoYPencilImaptable;
extern CkVec < MapType2 > RhoHartYPencilImaptable;
extern CkVec < CkVec <MapType2> > AtmSFYPencilImaptable;
extern CkVec < CkVec <int> > UberPes;

extern CkHashtableT <intdual, int> GSmaptable;
extern CkHashtableT <intdual, int> RSmaptable;
extern CkHashtableT <intdual, int> RPPmaptable;
extern CkHashtableT <intdual, int> RhoGSmaptable;
extern CkHashtableT <intdual, int> RhoRSmaptable;
extern CkHashtableT <intdual, int> RhoGHartmaptable;
extern CkHashtableT <inttriple, int> RhoRHartmaptable;
extern CkHashtableT <intdual, int> AsymScalcmaptable;
extern CkHashtableT <intdual, int> SymScalcmaptable;
CkReductionMsg *sumFastDouble(int nMsg, CkReductionMsg **msgs);
void fastAdd (double *a, double *b, int nelem);
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \brief The class which creates the main chare.
 *
 *
 *
 */
//============================================================================

class main : public Chare {
  public:
    main(CkMigrateMessage *m) { }
    main(CkArgMsg *);
    ~main();
};
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \brief Helper class for map hashtables copied from femrefine.C.
 *
 */
//============================================================================




//============================================================================
/** \brief Base Class used for maptable based proc maps
 *
 *
 */
class CkArrayMapTable1 : public CkArrayMap
{
  // trivial 1d here for usage consistency
  public:
    MapType1 *maptable;
    UberCollection thisInstance;

    CkArrayMapTable1() {}
    inline int procNum(int, const CkArrayIndex &iIndex){
      int *index=(int *) iIndex.data();
      int proc;
      proc=maptable->get(index[0]);
      CkAssert(proc >= 0);
      return(proc);

    }
    void pup(PUP::er &p)
    {
      CkArrayMap::pup(p);
      p|thisInstance;
    }

    ~CkArrayMapTable1(){}
};


class CkArrayMapTable2 : public CkArrayMap
{
  public:
    MapType2 *maptable;
    UberCollection thisInstance;

    CkArrayMapTable2() {}
    inline int procNum(int, const CkArrayIndex &iIndex){
      int *index = (int *) iIndex.data();
      int proc;

      proc=maptable->get(index[0],index[1]);
      CkAssert(proc >= 0);
      return(proc);

    }
    void pup(PUP::er &p)
    {
      CkArrayMap::pup(p);
      p|thisInstance;
    }
    ~CkArrayMapTable2(){}

};

class CkArrayMapTable3 : public CkArrayMap
{
  public:
    MapType3 *maptable;
    UberCollection thisInstance;

    CkArrayMapTable3() {}
    inline int procNum(int, const CkArrayIndex &iIndex){
      int *index = (int *) iIndex.data();
      int proc;
      proc = maptable->get(index[0], index[1], index[2]);
      CkAssert(proc >= 0);
      return(proc);
    }
    void pup(PUP::er &p)
    {
      CkArrayMap::pup(p);
      p|thisInstance;
    }
    ~CkArrayMapTable3(){}

};

class CkArrayMapTable4 : public CkArrayMap
{
  public:
    MapType4 *maptable;
    UberCollection thisInstance;

    CkArrayMapTable4() {}
    inline int procNum(int, const CkArrayIndex &iIndex){
      int proc;
      short *sindex = (short *) iIndex.data();
      proc = maptable->get(sindex[0], sindex[1], sindex[2], sindex[3]);
      CkAssert(proc >= 0);
      return(proc);
    }
    void pup(PUP::er &p)
    {
      CkArrayMap::pup(p);
      p|thisInstance;
    }
    ~CkArrayMapTable4(){}

};

class EnergyCommMgrMap : public CkArrayMap {

 public:
  MapType1 *maptable;
  UberCollection thisInstance;

 EnergyCommMgrMap(UberCollection _instance ) : thisInstance(_instance)
  {
    maptable= &EnergyCommMgrImaptable[thisInstance.getPO()];
    for(int element=0; element< maptable->getmax(); element++)
      {
	maptable->set(element,-1);
      }
    for(int element=0; element< UberPes[thisInstance.proxyOffset].length(); element++)
      {
	maptable->set(UberPes[thisInstance.proxyOffset][element],UberPes[thisInstance.proxyOffset][element]);
	//	CkPrintf("{%d}[%d] EnergyCommMgrMap set ith %d key %d to pe=%d\n",thisInstance.proxyOffset, CkMyPe(), element,  UberPes[thisInstance.proxyOffset][element], UberPes[thisInstance.proxyOffset][element]);
      }  
  }
  inline int procNum(int, const CkArrayIndex &iIndex){
    int *index=(int *) iIndex.data();
    int proc;
    proc = maptable->get(index[0]);
    proc = (proc >=0) ? proc : 0;
    //    CkPrintf("{%d}[%d] procNum index %d on pe=%d\n",thisInstance.proxyOffset, CkMyPe(), index[0],proc);
    CkAssert(proc >= 0);
    return(proc);
  }

  void pup(PUP::er &p)
  {
    CkArrayMap::pup(p);
    p|thisInstance;
  }


  void populateInitial(int arrayHdl,int numInitial, void *msg,CkArrMgr *mgr)
  {
    CkAbort("hail hail, populateInitial was actually used!\n");
    CkPrintf("pop init pe=%d\n",CkMyPe());
    if (numInitial==0) return; //No initial elements requested
    if(maptable->get(CkMyPe())>=0)
      {
	CkPrintf("pop init pe=%d\n",CkMyPe());
	mgr->insertInitial(CkArrayIndex1D(CkMyPe()),CkCopyMsg(&msg));
      }
    CkPrintf("pop init pe=%d\n",CkMyPe());
  }

};

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \brief Class used for instantiation of G-space group objects.
 *
 *
 */
//============================================================================


class AtomComputeMap : public CkArrayMapTable1 {
  public:
    AtomComputeMap(UberCollection _instance)
    {
      thisInstance=_instance;
      maptable= &AtomImaptable[thisInstance.getPO()];
    }

    ~AtomComputeMap()
    {
    }

    void pup(PUP::er &p)
    {
      CkArrayMap::pup(p);
      maptable= &AtomImaptable[thisInstance.getPO()];
    }
    inline int procNum(int, const CkArrayIndex &iIndex){
      int *index=(int *) iIndex.data();
      int proc;
      proc = maptable->get(index[0]);
      CkAssert(proc >= 0);
      return(proc);
    }

};

class GSMap: public CkArrayMapTable2 {

  public:
    GSMap(UberCollection _instance)
    {
      thisInstance=_instance;
      maptable = &GSImaptable[thisInstance.getPO()];
      if(CkMyPe()) {
        if(maptable == NULL) CkAbort("maptable does not exist in GSMap\n");
      }
    }
    void pup(PUP::er &p)
    {
      CkArrayMapTable2::pup(p);
      maptable= &GSImaptable[thisInstance.getPO()];
    }
    inline int procNum(int, const CkArrayIndex &iIndex){
      int *index = (int *) iIndex.data();
      int proc;
      proc = maptable->get(index[0],index[1]);
      CkAssert(proc >= 0);
      return(proc);
    }

    ~GSMap(){ }
};
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \brief Class used for instantiation of real-space group objects.
 *
 *
 */
//============================================================================

class RSMap: public CkArrayMapTable2 {

  public:
    RSMap(UberCollection _instance)
    {
      thisInstance = _instance;
      maptable= &RSImaptable[thisInstance.getPO()];
    }
    void pup(PUP::er &p)
    {
      CkArrayMapTable2::pup(p);
      maptable= &RSImaptable[thisInstance.getPO()];
    }
    inline int procNum(int, const CkArrayIndex &iIndex){
      int *index = (int *) iIndex.data();
      int proc;
      proc = maptable->get(index[0],index[1]);
      CkAssert(proc >= 0);
      return(proc);
    }

    ~RSMap(){ }
};
//============================================================================

class RPPMap: public CkArrayMapTable2 {

  public:
    RPPMap(UberCollection _instance)
    {
      thisInstance=_instance;
      maptable= &RPPImaptable[thisInstance.getPO()];
    }
    void pup(PUP::er &p)
    {
      CkArrayMapTable2::pup(p);
      maptable= &RPPImaptable[thisInstance.getPO()];
    }
    ~RPPMap(){ }

    inline int procNum(int, const CkArrayIndex &iIndex){
      int *index = (int *) iIndex.data();
      int proc;
      proc = maptable->get(index[0],index[1]);
      CkAssert(proc >= 0);
      return(proc);
    }

};
//============================================================================


//============================================================================
/**
 * provide procnum mapping for RhoR
 */
class RhoRSMap : public CkArrayMapTable2 {
  public:
    RhoRSMap(UberCollection _instance)
    {
      thisInstance=_instance;
      maptable= &RhoRSImaptable[thisInstance.getPO()];
    }

    ~RhoRSMap() {
    }

    void pup(PUP::er &p)
    {
      CkArrayMapTable2::pup(p);
      maptable = &RhoRSImaptable[thisInstance.getPO()];
    }

    inline int procNum(int, const CkArrayIndex &iIndex){
      int *index=(int *) iIndex.data();
      int proc;
      proc = maptable->get(index[0],index[1]);
      CkAssert(proc >= 0);
      return(proc);
    }
};


/**
 * provide procnum mapping for RhoG
 */
class RhoGSMap : public CkArrayMapTable1 {
  public:
    int dimindex;
    RhoGSMap(UberCollection _instance, int _index=0): dimindex(_index) {
      thisInstance = _instance;
      maptable= &RhoGSImaptable[thisInstance.getPO()];
    }
    ~RhoGSMap() {
    }

    inline int procNum(int, const CkArrayIndex &iIndex){
      int *index = (int *) iIndex.data();
      int proc;
      proc = maptable->get(index[dimindex]);
      CkAssert(proc >= 0);
      return(proc);
    }

    void pup(PUP::er &p)
    {
      CkArrayMapTable1::pup(p);
      maptable = &RhoGSImaptable[thisInstance.getPO()];
      p|dimindex;
    }
};

class RhoGHartMap : public CkArrayMapTable2 {
  public:
    int dimindex1, dimindex2, fixIndex;
    RhoGHartMap(UberCollection _instance, int _index1, int _index2, int
        _fixIndex) : dimindex1(_index1), dimindex2(_index2), fixIndex(_fixIndex) {
      thisInstance = _instance;
      maptable = &RhoGHartImaptable[thisInstance.getPO()];
    }

    ~RhoGHartMap() { }

    void pup(PUP::er &p)
    {
      CkArrayMapTable2::pup(p);
      maptable= &RhoGHartImaptable[thisInstance.getPO()];
      p|dimindex1;
      p|dimindex2;
      p|fixIndex;
    }
    inline int procNum(int, const CkArrayIndex &iIndex){
      int *index=(int *) iIndex.data();
      int proc;
      if(fixIndex == -1) {
        proc = maptable->get(index[dimindex1], index[dimindex2]);
      } else {
        proc = maptable->get(index[dimindex1], fixIndex);
      }
      CkAssert(proc>=0);
      return(proc);
    }
};

class RhoYPencilMap : public CkArrayMapTable2 {
 public:
    int ffttype, fftoffset;

    RhoYPencilMap(UberCollection _instance, int _ffttype, int _fftoffset)
    :  ffttype(_ffttype), fftoffset(_fftoffset) {
      thisInstance = _instance;
      if(ffttype == 0) {
	maptable= &RhoYPencilImaptable[fftoffset][thisInstance.getPO()];
      } else if(ffttype == 1) {
        maptable= &RhoHartYPencilImaptable[thisInstance.getPO()];
      } else if(ffttype == 2) {
        maptable= &AtmSFYPencilImaptable[fftoffset][thisInstance.getPO()];
      } else {
        CkAbort("RhoYPencilMap constructed with unknown ffttype\n");
      }
    }

    ~RhoYPencilMap() { }

    void pup(PUP::er &p)
    {
      CkArrayMapTable2::pup(p);
      p|ffttype;
      p|fftoffset;
      if(ffttype==0) {
	maptable= &RhoYPencilImaptable[fftoffset][thisInstance.getPO()];
      } else if(ffttype == 1) {
        maptable= &RhoHartYPencilImaptable[thisInstance.getPO()];
      } else if(ffttype == 2) {
        maptable= &AtmSFYPencilImaptable[fftoffset][thisInstance.getPO()];
      } else {
        CkAbort("RhoYPencilMap pupped with unknown ffttype\n");
      }
    }

    inline int procNum(int, const CkArrayIndex &iIndex){
      int *index = (int *) iIndex.data();
      int proc;
      proc = maptable->get(index[0],index[2]);
      CkAssert(proc >= 0);
      return(proc);
    }
};

class RhoRHartMap : public CkArrayMapTable3 {
  public:
    int dimindex1, dimindex2, fixIndex;
    RhoRHartMap(UberCollection _instance, int _index1, int _index2, int
        _fixIndex) : dimindex1(_index1), dimindex2(_index2), fixIndex(_fixIndex) {
      thisInstance=_instance;
      maptable= &RhoRHartImaptable[thisInstance.getPO()];
    }

    ~RhoRHartMap() { }

    void pup(PUP::er &p)
    {
      CkArrayMapTable3::pup(p);
      maptable= &RhoRHartImaptable[thisInstance.getPO()];
      p|dimindex1;
      p|dimindex2;
      p|fixIndex;
    }
    inline int procNum(int, const CkArrayIndex &iIndex){
      int *index=(int *) iIndex.data();
      int proc;
      if(fixIndex == -1) {
        proc=maptable->get(index[0],index[1],index[2]);
      } else {
        proc=maptable->get(index[dimindex1],index[dimindex2],fixIndex);
      }
      CkAssert(proc >= 0);
      return(proc);
    }
};


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================


/******** Optimization flags ********/
//#define _CP_USE_BLAS_GATHER_

/******** User Trace Events *********/
#define DoFFTContribute_     1100
#define doRealFwFFT_         1110
#define doRealBwFFT_         1120
#define GspaceFwFFT_         1130
#define GspaceBwFFT_         1140
#define RhoRtoGFFT_          1150

#define PostByrdfwFFTGtoR_   1152
#define BwFFTRtoG_           1153
#define ByrdanddoFwFFTGtoR_  1154
#define HartExcVksG_         1160
#define divRhoVksGspace_     1161
#define AcceptStructFact_    1162
#define eesHartExcG_         1163
#define eesEwaldG_           1164
#define eesAtmForcR_         1165
#define eesAtmBspline_       1166
#define eesZmatR_            1167
#define eesEnergyAtmForcR_   1168
#define eesProjG_            1169
#define fwFFTGtoR0_          1170
#define fwFFTGtoRnot0_       1171
#define doNlFFTGtoR_         1172
#define doNlFFTRtoG_         1173
#define eesPsiForcGspace_    1174
#define GradCorrGGA_         1180
#define WhiteByrdFFTX_       1181
#define WhiteByrdFFTY_       1182
#define WhiteByrdFFTZ_       1183
#define enlMatrixCalc_       1300
#define enlAtmForcCalc_      1301
#define enlForcCalc_         1302
#define doEextFFTRtoG_       1303
#define doEextFFTGtoR_       1304
#define doEextFFTGxtoRx_     1305
#define doEextFFTRytoGy_     1306
#define doRhoFFTRytoGy_      1307
#define doRhoFFTGxtoRx_      1308
#define OrthoDGEMM1_         1401
#define OrthoDGEMM2_         1402
//200-300 reserved for paircalculator
#define IntegrateModForces_  1000
#define Scalcmap_            2000
#define GHartAtmForcCopy_    3000
#define GHartAtmForcSend_    4000
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/*
 * Initialization routines to initialize the GSpace, Real Space,
 * Particle space planes, CP_Rho_GSpacePlane & CP_Rho_RealSpacePlane arrays
 and the S_Calculators
 */
//============================================================================
class size2d; //forward decl to shup the compiler
namespace cp { namespace paircalc { class pcConfig; } }
namespace pc = cp::paircalc;

void init_SF_non_EES(int natm_nl, int natm_nl_grp_max, int numSfGrps,
		      CPcharmParaInfo *sim, UberCollection thisInstance);
void build_all_maps(CPcharmParaInfo *sim, UberCollection thisInstance);
void build_uber_maps(CPcharmParaInfo *sim, UberCollection thisInstance);
void lst_sort_clean(int , int *, int *);
void init_PIBeads(CPcharmParaInfo *sim, UberCollection thisInstance);

void init_state_chares(int natm_nl,int natm_nl_grp_max,int numSfGrps,
    int doublePack, CPcharmParaInfo *sim, UberCollection thisInstance);

void init_eesNL_chares(int natm_nl,int natm_nl_grp_max,
    int doublePack, PeList *exclusion, CPcharmParaInfo *sim, UberCollection thisInstance);
int init_rho_chares(CPcharmParaInfo*, UberCollection thisInstance);
void control_physics_to_driver(UberCollection thisInstance, CPcharmParaInfo *sim);
void get_grp_params(int natm_nl, int numSfGrps, int indexSfGrp, int planeIndex,
    int *n_ret, int *istrt_ret, int *iend_ret);
int atmGrpMap(int istart, int nsend, int listsize, int *listpe, int AtmGrp,
    int dup, int planeIndex);
int gsprocNum(CPcharmParaInfo *sim,int state, int plane, int numInst);
void create_Rho_fft_numbers(int ,int ,int , int, int, int, int *,int *,int *,int *, int *);
void setTraceUserEvents();
void computeMapOffsets();

//============================================================================


//============================================================================

// stuff to be include before the decl or else
#include "Atoms.h"
#include "energy.h"
#include "paircalc/ckPairCalculator.h"
#include "cpaimd.decl.h"
void paircalcstartup(pc::pcConfig *cfgSymmPC, pc::pcConfig *cfgAsymmPC, CPcharmParaInfo *sim, int doublePack);
void orthostartup( cp::ortho::orthoConfig *orthoCfg, pc::pcConfig *cfgSymmPC, pc::pcConfig *cfgAsymmPC, CPcharmParaInfo *sim, PeListFactory *peList4PCmapping);

#endif

