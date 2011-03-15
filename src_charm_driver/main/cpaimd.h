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
#include "EachToManyMulticastStrategy.h"
#include "RingMulticastStrategy.h"
#include "StreamingStrategy.h"
#include "ckhashtable.h"

#define CORES_PER_NODE 4

#undef OLD_COMMLIB 
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
#include "load_balance/IntMap.h"
#endif

#include "load_balance/MapTable.h"

#ifdef CMK_BLUEGENEL
//#include "builtins.h"
#endif


//#define BARRIER_CP_GSPACE_PSI 1

#define LOAD_BALANCE_STEP 100000000

#define PRE_BALANCE_STEP 2

#define FIRST_BALANCE_STEP 100000000

#if CMK_TRACE_ENABLED
#define TRACE_ON_STEP 4
#define TRACE_OFF_STEP 7
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
extern CkVec <MapType2> GSImaptable;
extern CkVec <MapType2> RSImaptable;
extern CkVec <MapType2> RPPImaptable;
extern CkVec <MapType2> VdWGSImaptable;
extern CkVec <MapType3> VdWRSImaptable;
extern CkVec <MapType2> RhoGSImaptable;
extern CkVec <MapType2> RhoRSImaptable;
extern CkVec <MapType2> RhoGHartImaptable;
extern CkVec <MapType3> RhoRHartImaptable;

extern CkHashtableT <intdual, int> GSmaptable;
extern CkHashtableT <intdual, int> RSmaptable;
extern CkHashtableT <intdual, int> RPPmaptable;
extern CkHashtableT <intdual, int> VdWGSmaptable;
extern CkHashtableT <inttriple, int> VdWRSmaptable;
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
class CkArrayMapTable2 : public CkArrayMap
{
 public:
  MapType2 *maptable;
  UberCollection thisInstance;

  CkArrayMapTable2() {}
  inline int procNum(int, const CkArrayIndex &iIndex){
    int *index=(int *) iIndex.data();
    int proc;
    
#ifdef USE_INT_MAP
    proc=maptable->get(index[0],index[1]);
#else
    proc=maptable->get(intdual(index[0],index[1]));
#endif
    CkAssert(proc>=0);
    if(numPes!=CkNumPes())
      return(proc%CkNumPes());
    else
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
    int *index=(int *) iIndex.data();
    int proc;
    
#ifdef USE_INT_MAP
    proc=maptable->get(index[0],index[1],index[2]);
#else
    proc=maptable->get(inttriple(index[0],index[1],index[2]));
#endif
    CkAssert(proc>=0);
    if(numPes!=CkNumPes())
      return(proc%CkNumPes());
    else
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
    
#ifdef USE_INT_MAP
	short *sindex=(short *) iIndex.data();
	proc=maptable->get(sindex[0], sindex[1], sindex[2], sindex[3]);
#else
	int *index=(int *) iIndex.data();
	proc=maptable->get(intdual(index[0], index[1]));
#endif
    CkAssert(proc>=0);
    if(numPes!=CkNumPes())
      return(proc%CkNumPes());
    else
      return(proc);

  }
  void pup(PUP::er &p)
    {
      CkArrayMap::pup(p);
      p|thisInstance;
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
  GSMap(UberCollection _instance) 
      { 
	thisInstance=_instance;
#ifdef USE_INT_MAP
	maptable = &GSImaptable[thisInstance.getPO()];
	if(CkMyPe()) {
	  if(maptable == NULL)
	    CkAbort("hey 2");
	}
#else
	maptable = &GSmaptable;
#endif
      }
  void pup(PUP::er &p)
	{
	    CkArrayMapTable2::pup(p);
#ifdef USE_INT_MAP	
	    maptable= &GSImaptable[thisInstance.getPO()];
#else
	    maptable= &GSmaptable;
#endif
	}
  inline int procNum(int, const CkArrayIndex &iIndex){
    int *index=(int *) iIndex.data();
    int proc;
#ifdef USE_INT_MAP
    if(maptable == NULL) {
      CkPrintf("hey %d\n", CkMyPe());
      CkAbort("hey");
    }
    proc=maptable->get(index[0],index[1]);
#else
    proc=maptable->get(intdual(index[0],index[1]));
#endif
    CkAssert(proc>=0);
    if(numPes!=CkNumPes())
      return(proc%CkNumPes());
    else
      return(proc);
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
  RSMap(UberCollection _instance)
      { 
	thisInstance=_instance;
#ifdef USE_INT_MAP
	maptable= &RSImaptable[thisInstance.getPO()];
#else
	maptable= &RSmaptable;
#endif
      }
  void pup(PUP::er &p)
	{
	    CkArrayMapTable2::pup(p);
#ifdef USE_INT_MAP
	    maptable= &RSImaptable[thisInstance.getPO()];
#else
	    maptable= &RSmaptable;
#endif
	}
  //  int procNum(int, const CkArrayIndex &);
  inline int procNum(int, const CkArrayIndex &iIndex){
    int *index=(int *) iIndex.data();
    int proc;
    
#ifdef USE_INT_MAP
    proc=maptable->get(index[0],index[1]);
#else
    proc=maptable->get(intdual(index[0],index[1]));
#endif
    CkAssert(proc>=0);
    if(numPes!=CkNumPes())
      return(proc%CkNumPes());
    else
      return(proc);

  }

  ~RSMap(){
  }
};
//============================================================================

//Maps even sized stripes of a 2D array with strides across nodes
class NodeMap2DArray: public CkArrayMap {

public:
	//Need the size of one dimension, and the total size, and the offset of
	//the core within a node to map to e.g. an offset of 0 will map in a node
	//starting at core 0, an offset of 1 will start with core 1.  If more cores
	//than CORES_PER_NODE-offset are required, it will wrap around mapping
	//within a node (i.e. if offset = 2, cores per node = 4, and 5 chares per
	//node, the chares are mapped to cores 2, 3, 0, 1, 2)

	NodeMap2DArray(int _x_size, int _size, int offset){

		//Check offset
		if(offset<0)
			CkAbort("NodeMap2DArray: received negative offset!\n");

		x_size = _x_size;
		int size = _size;

#ifdef CRAYDEBUG
		CkPrintf("NodeMap2DArray using %d CORES_PER_NODE\n",CORES_PER_NODE);
#endif

		int cores_per_node = CORES_PER_NODE;
		int node_count = CkNumPes()/cores_per_node;

		//If ran into case where node_count <=0 due to single processor run or incorrect CORES_PER_NODE defined
		if(node_count<=0){
			CkPrintf("NodeMap2DArray: Num PEs / cores_per_node <=0, degrading to single node case\n");
			node_count=1;
			cores_per_node=1;
			offset=0;
		}

#ifdef CRAYDEBUG
		CkPrintf("NodeMap2DArray using %d Nodes\n",node_count);
#endif

		int chares_per_node = size/node_count;

		//number of nodes that need an extra chare if size/node_count has remainder
		int big_nodes = size%node_count;

		map = new int[size];

		//map the IDs of chares that belong on a node with an extra chare
		int count=0;
		for(int i=0;i<big_nodes;i++){
			for(int j=0;j<chares_per_node+1;j++){
				map[count]=i;
				count++;
			}
		}

		//map the remaining chares to the remaining cores
		for(int i=big_nodes;i<node_count;i++){
			for(int j=0;j<chares_per_node;j++){
				map[count]=i;
				count++;
			}
		}

		//At this point, the map array contains the node to which each chare
		//must be mapped.  The node needs to be converted to the exact core.
		//If an offset was defined, that needs to be taken into account.
		int core_num=0;
		int current_node;
		int previous_node=-1;

		for(int i=0;i<size;i++){
			current_node = map[i];

			if(current_node!=previous_node)
				core_num=offset;

			map[i] = (core_num%cores_per_node)+current_node*cores_per_node;
#ifdef CRAYDEBUG
			if(map[i]<0 || map[i]>=CkNumPes())
			{
				CkPrintf("attempted to map chare %d to proc %d\n",i,map[i]);
				CkAbort("NodeMap2DArray: mapped to out of bounds processor!");
			}
#endif
			previous_node = current_node;
			core_num++;
		}
	}

	//If you don't specify an offset, assumes offset=0
	NodeMap2DArray(int _x_size, int _size){
#ifdef CRAYDEBUG
		CkPrintf("NodeMap2DArray called WITHOUT an offset\n");
#endif
			NodeMap2DArray(_x_size, _size, 0);
	}

	//  int procNum(int, const CkArrayIndex &);
	inline int procNum(int, const CkArrayIndex &iIndex){
		int *index=(int *) iIndex.data();

		int proc=map[index[0]*x_size+index[1] ];

		CkAssert(proc>=0);
		return(proc);
	}

	~NodeMap2DArray(){
		delete[] map;
	}

private:
	int x_size;
	int *map;
};

class BlockMap2DArray: public CkArrayMap {

public:
	//Need the size of one dimension, and the total size
	BlockMap2DArray(int _x_size, int _size){
		x_size = _x_size;
		size = _size;
	}

	//  int procNum(int, const CkArrayIndex &);
	inline int procNum(int, const CkArrayIndex &iIndex){
		int *index=(int *) iIndex.data();

		int proc=(float)(index[0]*x_size+index[1])/size*CkNumPes();
		CkAssert(proc>=0);
		return(proc);
	}

	~BlockMap2DArray(){}

private:
	int x_size;
	int size;
};

class RPPMap: public CkArrayMapTable2 {

 public:
  RPPMap(UberCollection _instance)
      { 
	thisInstance=_instance;
#ifdef USE_INT_MAP
	maptable= &RPPImaptable[thisInstance.getPO()];
#else
	maptable= &RPPmaptable;
#endif
      }
  void pup(PUP::er &p)
	{
	    CkArrayMapTable2::pup(p);
#ifdef USE_INT_MAP
	    maptable= &RPPImaptable[thisInstance.getPO()];
#else
	    maptable= &RPPmaptable;
#endif

	}
  //  int procNum(int, const CkArrayIndex &);
  ~RPPMap(){

  }
  inline int procNum(int, const CkArrayIndex &iIndex){
    int *index=(int *) iIndex.data();
    int proc;
    
#ifdef USE_INT_MAP
    proc=maptable->get(index[0],index[1]);
#else
    proc=maptable->get(intdual(index[0],index[1]));
#endif
    CkAssert(proc>=0);
    if(numPes!=CkNumPes())
      return(proc%CkNumPes());
    else
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
    int nchareRhoR;
    RhoRSMap(UberCollection _instance)
    {
	thisInstance=_instance;
#ifdef USE_INT_MAP
      maptable= &RhoRSImaptable[thisInstance.getPO()];
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
	maptable= &RhoRSImaptable[thisInstance.getPO()];
#else
	maptable= &RhoRSmaptable;
#endif
      }
    
    //    int procNum(int arrayHdl, const CkArrayIndex &idx);
  inline int procNum(int, const CkArrayIndex &iIndex){
    int *index=(int *) iIndex.data();
    int proc;
    
#ifdef USE_INT_MAP
    proc=maptable->get(index[0],index[1]);
#else
    proc=maptable->get(intdual(index[0],index[1]));
#endif
    CkAssert(proc>=0);
    if(numPes!=CkNumPes())
      return(proc%CkNumPes());
    else
      return(proc);
  }

};


/**
 * provide procnum mapping for RhoG
 */
class RhoGSMap : public CkArrayMapTable2 {
  public:
  RhoGSMap(UberCollection _instance)
    {
	thisInstance=_instance;
#ifdef USE_INT_MAP
      maptable= &RhoGSImaptable[thisInstance.getPO()];
#else
      maptable= &RhoGSmaptable;
#endif
    }
    
    ~RhoGSMap() {
    }
    
    //    int procNum(int arrayHdl, const CkArrayIndex &idx);
  inline int procNum(int, const CkArrayIndex &iIndex){
    int *index=(int *) iIndex.data();
    int proc;
    
#ifdef USE_INT_MAP
    proc=maptable->get(index[0],index[1]);
#else
    proc=maptable->get(intdual(index[0],index[1]));
#endif
    CkAssert(proc>=0);
    if(numPes!=CkNumPes())
      return(proc%CkNumPes());
    else
      return(proc);

  }
    
    void pup(PUP::er &p)
      {
	CkArrayMapTable2::pup(p);
      }
};


class RhoGHartMap : public CkArrayMapTable2 {
  public:
  RhoGHartMap(UberCollection _instance)
  {
	thisInstance=_instance;
#ifdef USE_INT_MAP
    maptable= &RhoGHartImaptable[thisInstance.getPO()];
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
	maptable= &RhoGHartImaptable[thisInstance.getPO()];
#else
	maptable= &RhoGHartmaptable;
#endif
      }
  inline int procNum(int, const CkArrayIndex &iIndex){
    int *index=(int *) iIndex.data();
    int proc;
#ifdef USE_INT_MAP
    proc=maptable->get(index[0],index[1]);
#else
    proc=maptable->get(intdual(index[0],index[1]));
#endif
    CkAssert(proc>=0);
    if(numPes!=CkNumPes())
      return(proc%CkNumPes());
    else
      return(proc);

  }


};

class RhoRHartMap : public CkArrayMapTable3 {
  public:
  RhoRHartMap(UberCollection _instance)
  {
    thisInstance=_instance;
#ifdef USE_INT_MAP
    maptable= &RhoRHartImaptable[thisInstance.getPO()];
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
	maptable= &RhoRHartImaptable[thisInstance.getPO()];
#else
	maptable= &RhoRHartmaptable;
#endif
      }
  inline int procNum(int, const CkArrayIndex &iIndex){
    int *index=(int *) iIndex.data();
    int proc;
    
#ifdef USE_INT_MAP
    proc=maptable->get(index[0],index[1],index[2]);
#else
    proc=maptable->get(inttriple(index[0],index[1],index[2]));
#endif
    CkAssert(proc>=0);
    if(numPes!=CkNumPes())
      return(proc%CkNumPes());
    else
      return(proc);
  }

};

//============================================================================
/**
 * provide procnum mapping for RhoR
 */
class VdWRSMap : public CkArrayMapTable3 {
  public:
    int nchareRhoR;
    VdWRSMap(UberCollection _instance)
    {
	thisInstance=_instance;
#ifdef USE_INT_MAP
      maptable= &VdWRSImaptable[thisInstance.getPO()];
#else
      maptable= &VdWRSmaptable;
#endif
    }
    
    ~VdWRSMap() {
    }
  
    void pup(PUP::er &p)
      {
	CkArrayMapTable3::pup(p);
#ifdef USE_INT_MAP
	maptable= &VdWRSImaptable[thisInstance.getPO()];
#else
	maptable= &VdWRSmaptable;
#endif
      }
    
    //    int procNum(int arrayHdl, const CkArrayIndex &idx);
  inline int procNum(int, const CkArrayIndex &iIndex){
    int *index=(int *) iIndex.data();
    int proc;
    
#ifdef USE_INT_MAP
    proc=maptable->get(index[0],index[1],index[2]);
#else
    proc=maptable->get(inttriple(index[0],index[1],index[2]));
#endif
    CkAssert(proc>=0);
    if(numPes!=CkNumPes())
      return(proc%CkNumPes());
    else
      return(proc);
  }

};


/**
 * provide procnum mapping for VdWG
 */
class VdWGSMap : public CkArrayMapTable2 {
  public:
  VdWGSMap(UberCollection _instance)
    {
	thisInstance=_instance;
#ifdef USE_INT_MAP
      maptable= &VdWGSImaptable[thisInstance.getPO()];
#else
      maptable= &VdWGSmaptable;
#endif
    }
    
    ~VdWGSMap() {
    }
    
    //    int procNum(int arrayHdl, const CkArrayIndex &idx);
  inline int procNum(int, const CkArrayIndex &iIndex){
    int *index=(int *) iIndex.data();
    int proc;
    
#ifdef USE_INT_MAP
    proc=maptable->get(index[0],index[1]);
#else
    proc=maptable->get(intdual(index[0],index[1]));
#endif
    CkAssert(proc>=0);
    if(numPes!=CkNumPes())
      return(proc%CkNumPes());
    else
      return(proc);

  }
    
    void pup(PUP::er &p)
      {
	CkArrayMapTable2::pup(p);
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

void init_commlib_strategies(int, int,int, UberCollection thisInstance);
void lst_sort_clean(int , int *, int *);
void init_PIBeads(CPcharmParaInfo *sim, UberCollection thisInstance);

void init_state_chares(int natm_nl,int natm_nl_grp_max,int numSfGrps,
                       int doublePack, CPcharmParaInfo *sim, UberCollection thisInstance);

void init_eesNL_chares(int natm_nl,int natm_nl_grp_max,
                       int doublePack, PeList *exclusion, CPcharmParaInfo *sim, UberCollection thisInstance);
int init_rho_chares(CPcharmParaInfo*, UberCollection thisInstance);
void init_VdW_chares(CPcharmParaInfo*, UberCollection thisInstance);
void control_physics_to_driver(UberCollection thisInstance);
void get_grp_params(int natm_nl, int numSfGrps, int indexSfGrp, int planeIndex,
		    int *n_ret, int *istrt_ret, int *iend_ret);
int atmGrpMap(int istart, int nsend, int listsize, int *listpe, int AtmGrp, 
              int dup, int planeIndex);
int gsprocNum(CPcharmParaInfo *sim,int state, int plane, int numInst);
bool findCuboid(int &x, int &y, int &z, int &order, int maxX, int maxY, int maxZ, int maxT, int volume, int vn);
void create_Rho_fft_numbers(int ,int ,int , int, int, int, int *,int *,int *,int *, int *);
void setTraceUserEvents();
//============================================================================


//============================================================================

// stuff to be include before the decl or else
#include "Atoms.h"
#include "energy.h"
#include "paircalc/ckPairCalculator.h"
#include "cpaimd.decl.h"

#endif

