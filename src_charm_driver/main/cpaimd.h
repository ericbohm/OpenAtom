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

#ifndef _CPAIMD_H
#define _CPAIMD_H
#define USE_TOPOMAP 1
#include "EachToManyMulticastStrategy.h"
#include "RingMulticastStrategy.h"
#include "StreamingStrategy.h"
#include "pairCalculator.h"
#include "ckhashtable.h"
#include "PeList.h"
#include "MapTable.h"


//#define GPSI_BARRIER 1

#define LOAD_BALANCE_STEP 10

#define PRE_BALANCE_STEP 2

#define FIRST_BALANCE_STEP 10

#ifndef CMK_OPTIMIZE
#define TRACE_ON_STEP 4
#define TRACE_OFF_STEP 7
#endif

#ifndef CmiMemcpy
#define CmiMemcpy(dest, src, size) memcpy((dest), (src), (size))
#endif

extern CkHashtableT <intdual, int> GSmaptable;
extern CkHashtableT <intdual, int> RSmaptable;
extern CkHashtableT <intdual, int> RhoGSmaptable;
extern CkHashtableT <intdual, int> RhoRSmaptable;
extern CkHashtableT <intdual, int> RhoGHartmaptable;
extern CkHashtableT <intdual, int> AsymScalcmaptable;
extern CkHashtableT <intdual, int> SymScalcmaptable;


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
class CkArrayMapTable : public CkArrayMap
{
 public:
  CkHashtableT<intdual, int> *maptable;


  CkArrayMapTable() {}
  int procNum(int, const CkArrayIndex &){CkAbort("do not call the base procnum");return(0);}
  void pup(PUP::er &p)
    {
      CkArrayMap::pup(p);
    }
  ~CkArrayMapTable(){}
  
};

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \brief Class used for instantiation of G-space group objects.
 *
 *
 */
//============================================================================


 
class GSMap: public CkArrayMapTable {

 public:
  GSMap()
      { 
	maptable= &GSmaptable;
      }
  void pup(PUP::er &p)
	{
	    CkArrayMapTable::pup(p);
	    maptable= &GSmaptable;
	}
  int procNum(int, const CkArrayIndex &);
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

class RSMap: public CkArrayMapTable {

 public:
  RSMap()
      { 
	maptable= &RSmaptable;
      }
  void pup(PUP::er &p)
	{
	    CkArrayMapTable::pup(p);
	    maptable= &RSmaptable;
	}
  int procNum(int, const CkArrayIndex &);
  ~RSMap(){
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

class SCalcMap : public CkArrayMapTable {
  bool symmetric;
 public:
    SCalcMap(bool _symmetric): symmetric(_symmetric){
	if(symmetric)
	  maptable= &SymScalcmaptable;
	else
	  maptable= &AsymScalcmaptable;
    }
    void pup(PUP::er &p)
	{
	    CkArrayMapTable::pup(p);
	    p|symmetric;
	    if(symmetric)
	      maptable= &SymScalcmaptable;
	    else
	      maptable= &AsymScalcmaptable;
	}
  int procNum(int, const CkArrayIndex &);
};
//============================================================================
/**
 * provide procnum mapping for RhoR
 */
class RhoRSMap : public CkArrayMapTable {
  public:
    int nchareRhoR;
    RhoRSMap()
    {
      maptable= &RhoRSmaptable;
    }
    
    ~RhoRSMap() {
    }
  
    void pup(PUP::er &p)
      {
	CkArrayMapTable::pup(p);
	maptable= &RhoRSmaptable;
      }
    
    int procNum(int arrayHdl, const CkArrayIndex &idx);

};


/**
 * provide procnum mapping for RhoG
 */
class RhoGSMap : public CkArrayMapTable {
  public:
    RhoGSMap()
    {
      maptable= &RhoGSmaptable;
    }
    
    ~RhoGSMap() {
    }
    
    int procNum(int arrayHdl, const CkArrayIndex &idx);
    
    void pup(PUP::er &p)
      {
	CkArrayMapTable::pup(p);
      }
};


class RhoGHartMap : public CkArrayMapTable {
  public:
  RhoGHartMap()
  {
    maptable= &RhoGHartmaptable;
  }
  
  ~RhoGHartMap()
  {
  }
  
  int procNum(int arrayHdl, const CkArrayIndex &idx);
    
  void pup(PUP::er &p)
      {
	CkArrayMapTable::pup(p);
	maptable= &RhoGHartmaptable;
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
#define fwFFTGtoR0_          170
#define fwFFTGtoRnot0_       171
#define GradCorrGGA_         180
#define WhiteByrdFFTX_       181
#define WhiteByrdFFTY_       182
#define WhiteByrdFFTZ_       183
#define OrthoDGEMM1_         401
#define OrthoDGEMM2_         402
//200-300 reserved for paircalculator
#define IntegrateModForces_  1000
#define Scalcmap_            2000
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
void init_pair_calculators(int nstates, int indexSize, int *indexZ, int doublePack, CPcharmParaInfo *sim);
void init_ortho_chares(int, int, int *);

void init_commlib_strategies(int, int);
void lst_sort_clean(int , int *, int *);
void init_state_chares(size2d,int,int,int,int,CPcharmParaInfo *);
void init_rho_chares(size2d ,CPcharmParaInfo*);
void control_physics_to_driver();
void get_grp_params(int natm_nl, int numSfGrps, int indexSfGrp, int planeIndex,
		    int *n_ret, int *istrt_ret, int *iend_ret);
int atmGrpMap(int istart, int nsend, int listsize, int *listpe, int AtmGrp, 
              int dup, int planeIndex);
int gsprocNum(CPcharmParaInfo *sim,int state, int plane);

//============================================================================

// stuff to be include before the decl or else
#include "../../include/Atoms.h"
#include "energy.h"
#include "cpaimd.decl.h"

#endif
