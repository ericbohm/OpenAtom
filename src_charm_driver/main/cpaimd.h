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

#ifndef _CPAIMD_H
#define _CPAIMD_H
//#define MAP_DEBUG 1
#include "EachToManyMulticastStrategy.h"
#include "RingMulticastStrategy.h"
#include "StreamingStrategy.h"
#include "pairCalculator.h"
#include "ckhashtable.h"
#include "PeList.h"
#define USE_INT_MAP
#ifndef USE_INT_MAP
class IntMap2
{
 public:
      void pup(PUP::er &p)
      {
      }
};

typedef IntMap2 IntMap4;
#else
#include "IntMap.h"
typedef IntMap4 MapType4;
typedef IntMap3 MapType3;
class MapType2 : public IntMap2on2 {
 public:
  int getCentroid (int);
};

#endif

#include "MapTable.h"

#ifdef CMK_VERSION_BLUEGENE
//#include "builtins.h"
#endif


//#define GPSI_BARRIER 1

#define LOAD_BALANCE_STEP 100000000

#define PRE_BALANCE_STEP 2

#define FIRST_BALANCE_STEP 100000000

#ifndef CMK_OPTIMIZE
#define TRACE_ON_STEP 16
#define TRACE_OFF_STEP 20
#endif

#ifndef CmiMemcpy
#define CmiMemcpy(dest, src, size) CmiMemcpy((dest), (src), (size))
#endif

extern MapType2 GSImaptable;
extern MapType2 RSImaptable;
extern MapType2 RPPImaptable;
extern MapType2 RhoGSImaptable;
extern MapType2 RhoRSImaptable;
extern MapType2 RhoGHartImaptable;
extern MapType3 RhoRHartImaptable;
extern MapType4 AsymScalcImaptable;
extern MapType4 SymScalcImaptable;

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
    void doneInit(CkReductionMsg *msg);
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
class CkArrayMapTable2 : public CkArrayMap
{
 public:
  MapType2 *maptable;


  CkArrayMapTable2() {}
  inline int procNum(int, const CkArrayIndex &iIndex){
    int *index=(int *) iIndex.data();
#ifdef USE_INT_MAP
	return(maptable->get(index[0],index[1]));
#else
	return(maptable->get(intdual(index[0],index[1])));
#endif
  }
  void pup(PUP::er &p)
    {
      CkArrayMap::pup(p);
    }
  ~CkArrayMapTable2(){}
  
};

class CkArrayMapTable3 : public CkArrayMap
{
 public:
  MapType3 *maptable;


  CkArrayMapTable3() {}
  inline int procNum(int, const CkArrayIndex &iIndex){
    int *index=(int *) iIndex.data();
#ifdef USE_INT_MAP
	return(maptable->get(index[0],index[1],index[2]));
#else
	return(maptable->get(inttriple(index[0],index[1],index[2])));
#endif
  }
  void pup(PUP::er &p)
    {
      CkArrayMap::pup(p);
    }
  ~CkArrayMapTable3(){}
  
};

class CkArrayMapTable4 : public CkArrayMap
{
 public:
  MapType4 *maptable;


  CkArrayMapTable4() {}
  inline int procNum(int, const CkArrayIndex &iIndex){
#ifdef USE_INT_MAP
	short *sindex=(short *) iIndex.data();
	return(maptable->get(sindex[0], sindex[1], sindex[2], sindex[3]));
#else
	int *index=(int *) iIndex.data();
	return(maptable->get(intdual(index[0], index[1])));
#endif
  }
  void pup(PUP::er &p)
    {
      CkArrayMap::pup(p);
    }
  ~CkArrayMapTable4(){}
  
};

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \brief Class used for instantiation of G-space group objects.
 *
 *
 */
//============================================================================


 
class GSMap: public CkArrayMapTable2 {

 public:
  GSMap()
      { 
#ifdef USE_INT_MAP	
	maptable= &GSImaptable;
#else
	maptable= &GSmaptable;
#endif
      }
  void pup(PUP::er &p)
	{
	    CkArrayMapTable2::pup(p);
#ifdef USE_INT_MAP	
	    maptable= &GSImaptable;
#else
	    maptable= &GSmaptable;
#endif
	}
  inline int procNum(int, const CkArrayIndex &iIndex){
    int *index=(int *) iIndex.data();
#ifdef USE_INT_MAP
	return(maptable->get(index[0],index[1]));
#else
	return(maptable->get(intdual(index[0],index[1])));
#endif
  }

  //  int procNum(int, const CkArrayIndex &);
  ~GSMap(){
  }
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
  RSMap()
      { 
#ifdef USE_INT_MAP
	maptable= &RSImaptable;
#else
	maptable= &RSmaptable;
#endif
      }
  void pup(PUP::er &p)
	{
	    CkArrayMapTable2::pup(p);
#ifdef USE_INT_MAP
	    maptable= &RSImaptable;
#else
	    maptable= &RSmaptable;
#endif
	}
  //  int procNum(int, const CkArrayIndex &);
  inline int procNum(int, const CkArrayIndex &iIndex){
    int *index=(int *) iIndex.data();
#ifdef USE_INT_MAP
	return(maptable->get(index[0],index[1]));
#else
	return(maptable->get(intdual(index[0],index[1])));
#endif
  }

  ~RSMap(){
  }
};
//============================================================================

class RPPMap: public CkArrayMapTable2 {

 public:
  RPPMap()
      { 
#ifdef USE_INT_MAP
	maptable= &RPPImaptable;
#else
	maptable= &RPPmaptable;
#endif
      }
  void pup(PUP::er &p)
	{
	    CkArrayMapTable2::pup(p);
#ifdef USE_INT_MAP
	    maptable= &RPPImaptable;
#else
	    maptable= &RPPmaptable;
#endif

	}
  //  int procNum(int, const CkArrayIndex &);
  ~RPPMap(){

  }
  inline int procNum(int, const CkArrayIndex &iIndex){
    int *index=(int *) iIndex.data();
#ifdef USE_INT_MAP
	return(maptable->get(index[0],index[1]));
#else
	return(maptable->get(intdual(index[0],index[1])));
#endif
  }

};
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//============================================================================
/** \brief Class used for instantiation of pair-calculator group objects.
 *
 *
 */
//============================================================================

class SCalcMap : public CkArrayMapTable4 {
  bool symmetric;
 public:
    SCalcMap(bool _symmetric): symmetric(_symmetric){
#ifdef USE_INT_MAP
	if(symmetric)
	  maptable= &SymScalcImaptable;
	else
	  maptable= &AsymScalcImaptable;
#else
	if(symmetric)
	  maptable= &SymScalcmaptable;
	else
	  maptable= &AsymScalcmaptable;
#endif
    }
    void pup(PUP::er &p)
	{
	    CkArrayMapTable4::pup(p);
	    p|symmetric;
#ifdef USE_INT_MAP
	    if(symmetric)
	      maptable= &SymScalcImaptable;
	    else
	      maptable= &AsymScalcImaptable;
#else
	    if(symmetric)
	      maptable= &SymScalcmaptable;
	    else
	      maptable= &AsymScalcmaptable;
#endif
	}
    //  int procNum(int, const CkArrayIndex &);
  inline int procNum(int, const CkArrayIndex &iIndex){
#ifdef USE_INT_MAP
	short *sindex=(short *) iIndex.data();
	return(maptable->get(sindex[0], sindex[1], sindex[2], sindex[3]));
#else
	int *index=(int *) iIndex.data();
	return(maptable->get(intdual(index[0], index[1])));
#endif
  }

};
//============================================================================
/**
 * provide procnum mapping for RhoR
 */
class RhoRSMap : public CkArrayMapTable2 {
  public:
    int nchareRhoR;
    RhoRSMap()
    {
#ifdef USE_INT_MAP
      maptable= &RhoRSImaptable;
#else
      maptable= &RhoRSmaptable;
#endif
    }
    
    ~RhoRSMap() {
    }
  
    void pup(PUP::er &p)
      {
	CkArrayMapTable2::pup(p);
#ifdef USE_INT_MAP
	maptable= &RhoRSImaptable;
#else
	maptable= &RhoRSmaptable;
#endif
      }
    
    //    int procNum(int arrayHdl, const CkArrayIndex &idx);
  inline int procNum(int, const CkArrayIndex &iIndex){
    int *index=(int *) iIndex.data();
#ifdef USE_INT_MAP
	return(maptable->get(index[0],index[1]));
#else
	return(maptable->get(intdual(index[0],index[1])));
#endif
  }

};


/**
 * provide procnum mapping for RhoG
 */
class RhoGSMap : public CkArrayMapTable2 {
  public:
    RhoGSMap()
    {
#ifdef USE_INT_MAP
      maptable= &RhoGSImaptable;
#else
      maptable= &RhoGSmaptable;
#endif
    }
    
    ~RhoGSMap() {
    }
    
    //    int procNum(int arrayHdl, const CkArrayIndex &idx);
  inline int procNum(int, const CkArrayIndex &iIndex){
    int *index=(int *) iIndex.data();
#ifdef USE_INT_MAP
	return(maptable->get(index[0],index[1]));
#else
	return(maptable->get(intdual(index[0],index[1])));
#endif
  }
    
    void pup(PUP::er &p)
      {
	CkArrayMapTable2::pup(p);
      }
};


class RhoGHartMap : public CkArrayMapTable2 {
  public:
  RhoGHartMap()
  {
#ifdef USE_INT_MAP
    maptable= &RhoGHartImaptable;
#else
    maptable= &RhoGHartmaptable;
#endif
  }
  
  ~RhoGHartMap()
  {
  }
  
  //  int procNum(int arrayHdl, const CkArrayIndex &idx);
    
  void pup(PUP::er &p)
      {
	CkArrayMapTable2::pup(p);
#ifdef USE_INT_MAP
	maptable= &RhoGHartImaptable;
#else
	maptable= &RhoGHartmaptable;
#endif
      }
  inline int procNum(int, const CkArrayIndex &iIndex){
    int *index=(int *) iIndex.data();
#ifdef USE_INT_MAP
	return(maptable->get(index[0],index[1]));
#else
	return(maptable->get(intdual(index[0],index[1])));
#endif
  }


};

class RhoRHartMap : public CkArrayMapTable3 {
  public:
  RhoRHartMap()
  {
#ifdef USE_INT_MAP
    maptable= &RhoRHartImaptable;
#else
    maptable= &RhoRHartmaptable;
#endif
  }
  
  ~RhoRHartMap()
  {
  }
  
  //  int procNum(int arrayHdl, const CkArrayIndex &idx);
    
  void pup(PUP::er &p)
      {
	CkArrayMapTable3::pup(p);
#ifdef USE_INT_MAP
	maptable= &RhoRHartImaptable;
#else
	maptable= &RhoRHartmaptable;
#endif
      }
  inline int procNum(int, const CkArrayIndex &iIndex){
    int *index=(int *) iIndex.data();
#ifdef USE_INT_MAP
	return(maptable->get(index[0],index[1],index[2]));
#else
	return(maptable->get(inttriple(index[0],index[1],index[2])));
#endif
  }

};


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//============================================================================
/**
 * \brief class CPcharmParaInfoGrp.
 *
 *
 */
//============================================================================

#include "../../include/CPcharmParaInfo.h"
class CPcharmParaInfoGrp: public Group {
 public:
    CPcharmParaInfoGrp(CkMigrateMessage *m) {}
    CPcharmParaInfoGrp(CPcharmParaInfo &s);
    ~CPcharmParaInfoGrp();
    CPcharmParaInfo *cpcharmParaInfo;
};
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================


/******** Optimization flags ********/
 //#define _CP_USE_BLAS_GATHER_

/******** User Trace Events *********/
#define DoFFTContribute_     100
#define doRealFwFFT_         110
#define doRealBwFFT_         120
#define GspaceFwFFT_         130
#define GspaceBwFFT_         140
#define RhoRtoGFFT_          150

#define PostByrdfwFFTGtoR_   152
#define BwFFTRtoG_           153
#define ByrdanddoFwFFTGtoR_  154
#define HartExcVksG_         160
#define divRhoVksGspace_     161
#define AcceptStructFact_    162
#define eesHartExcG_         163
#define eesEwaldG_           164
#define eesAtmForcR_         165
#define eesAtmBspline_       166
#define eesZmatR_            167
#define eesEnergyAtmForcR_   168
#define eesProjG_            169
#define fwFFTGtoR0_          170
#define fwFFTGtoRnot0_       171
#define doNlFFTGtoR_         172
#define doNlFFTRtoG_         173
#define eesPsiForcGspace_    174
#define GradCorrGGA_         180
#define WhiteByrdFFTX_       181
#define WhiteByrdFFTY_       182
#define WhiteByrdFFTZ_       183
#define enlMatrixCalc_       300
#define enlAtmForcCalc_      301
#define enlForcCalc_         302
#define doEextFFTRtoG_       303
#define doEextFFTGtoR_       304
#define doEextFFTGxtoRx_     305
#define doEextFFTRytoGy_     306
#define doRhoFFTRytoGy_      307
#define doRhoFFTGxtoRx_      308
#define OrthoDGEMM1_         401
#define OrthoDGEMM2_         402
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
void init_pair_calculators(int nstates, int indexSize, int *indexZ, int doublePack, 
                           CPcharmParaInfo *sim, int boxSize);
void init_ortho_chares(int, int, int *);

void init_commlib_strategies(int, int,int);
void lst_sort_clean(int , int *, int *);
void init_state_chares(size2d,int,int,int,int,CPcharmParaInfo *);
void init_eesNL_chares(size2d sizeYZ, int natm_nl,int natm_nl_grp_max,
                       int doublePack, PeList *exclusion, CPcharmParaInfo *sim);
void init_rho_chares(size2d ,CPcharmParaInfo*);
void control_physics_to_driver();
void get_grp_params(int natm_nl, int numSfGrps, int indexSfGrp, int planeIndex,
		    int *n_ret, int *istrt_ret, int *iend_ret);
int atmGrpMap(int istart, int nsend, int listsize, int *listpe, int AtmGrp, 
              int dup, int planeIndex);
int gsprocNum(CPcharmParaInfo *sim,int state, int plane);
bool findCuboid(int &x, int &y, int &z, int maxX, int maxY, int maxZ, int volume, int &order, int vn);
void create_Rho_fft_numbers(int ,int ,int , int, int, int, int *,int *,int *,int *, int *);

//============================================================================


//============================================================================

// stuff to be include before the decl or else
#include "../../include/Atoms.h"
#include "energy.h"
#include "cpaimd.decl.h"

#endif
