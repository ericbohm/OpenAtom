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
	//than cores_per_node-offset are required, it will wrap around mapping
	//within a node (i.e. if offset = 2, cores per node = 4, and 5 chares per
	//node, the chares are mapped to cores 2, 3, 0, 1, 2)

	//the input num_cores_to_use (per node) and core_offset can be used together to exclude certain cores in a node

	//Definitions:
	//cores_per_node: total number of cores per node (pass in TopoManager's getProcsPerNode() for CrayXT5)
	//offset: core offset within the given set of cores and nodes
	//core_offset: used to exclude cores from a mapping when used with num_nodes_to_use
	//             i.e. core_offset=1 and num_cores_to_use_per_node=cores_per_node-1 would exclude core 0 of every node
	//node_offset: used to exclude nodes when used with num_nodes_to_use
	//             i.e. node_offset=num_nodes/2 and num_nodes_to_use=num_nodes/2 excludes the first half of the nodes
	//stripe: value of 0 stripes 2D arrays across Y dimension, value of 1 stripes across X dimension.  Irrelevant for 1D or 3D arrays

	//big_nodes: if size/num_nodes has a remainder, the number of nodes that have size/num_nodes+1 chares
	//big_cores: number of cores on all big nodes (big_nodes*cores_per_node)
	//small_cores: number of cores on all small nodes ((num_nodes-big_nodes)*cores_per_node)
	//chares_on_big_nodes: the total number of chares assigned to big nodes
	//chares_on_small_nodes: the total number of chares assigned to small nodes

	NodeMap2DArray(int _offset, int _num_cores_to_use_per_node, int _core_offset, int _num_nodes_to_use, int _node_offset, int _stripe){

		total_num_nodes = CmiNumPhysicalNodes();
		total_cores_per_node = CkNumPes()/total_num_nodes;

#ifdef CRAYDEBUG
		CkPrintf("NODEMAP INFO: detected %d cores per node and %d nodes\n",total_cores_per_node,total_num_nodes);
#endif

		stripe = _stripe;

		//If number of cores to use per node is out of bounds, set it to the value passed into _cores_per_node
		if(_num_cores_to_use_per_node > total_cores_per_node || _num_cores_to_use_per_node <= 0){
			CkPrintf("NODEMAP WARNING: num_cores_to_use_per_node = %d is out of range, setting to detected cores per node\n", _num_cores_to_use_per_node, total_cores_per_node);
			cores_per_node = total_cores_per_node;
		}
		else
			cores_per_node = _num_cores_to_use_per_node;

		//If number of nodes to use per node is out of bounds, set it to the value determined by the API
		if(_num_nodes_to_use <= 0 || _num_nodes_to_use > total_num_nodes){
			CkPrintf("NODEMAP WARNING: num_nodes_to_use = %d is out of range, setting to detected value of nodes: %d\n", _num_nodes_to_use, total_num_nodes);
			num_nodes = total_num_nodes;
		}
		else
			num_nodes = _num_nodes_to_use;

		//Check offset
		if(_offset < 0) {
			CkPrintf("NODEMAP WARNING: offset = %d is negative, setting to zero\n", _offset);
			offset = 0;
		}
		else
			offset = _offset;

		//if core_offset out of bounds, mod it to fit in bounds
		if(_core_offset < 0 || _core_offset >= total_cores_per_node) {
			CkPrintf("NODEMAP WARNING: core_offset = %d is out of bounds, setting it to core_offset MOD total_cores_per_node = %d\n", _core_offset, _core_offset % total_cores_per_node);
			core_offset = _core_offset % total_cores_per_node;
		}
		else
			core_offset = _core_offset;

		//if node_offset out of bounds, mod it to fit in bounds
		if(_node_offset < 0 || _node_offset >= total_num_nodes) {
			CkPrintf("NODEMAP WARNING: node_offset = %d is out of bounds, setting it to node_offset MOD total_num_nodes = %d\n", _node_offset, _node_offset % total_num_nodes);
			node_offset = _node_offset % total_num_nodes;
		}
		else
			node_offset = _node_offset;

#ifdef CRAYDEBUG
		CkPrintf("NODEMAP using %d Nodes\n",num_nodes);
#endif
	}

	int registerArray(CkArrayIndexMax& numElements,CkArrayID aid) {

		int dim = numElements.dimension;

		size = numElements.getCombinedCount();
		if(size <= 0)
			CkAbort("NODEMAP received size <= 0\n");

#ifdef CRAYDEBUG
		CkPrintf("==========NODEMAP ============== size = %d, dim = %d\n",size, dim);
#endif

		if(dim == 2 || dim == 3) {
			x_size = numElements.data()[0];
			if(x_size <= 0)
				CkAbort("NODEMAP received x_size <= 0\n");

			y_size = numElements.data()[1];
			if(x_size <= 0)
				CkAbort("NODEMAP received y_size <= 0\n");
		}

		if(dim > 3)
			CkAbort("NodeMap cannot handle Chare arrays with more than 3 dimensions\n");

#ifdef CRAYDEBUG
		CkPrintf("==========NODEMAP ============== x_size = %d, y_size = %d\n",x_size,y_size);
#endif

		chares_per_node = size/num_nodes;

		//number of nodes that need an extra chare if size/num_nodes has remainder, a big node
		//big nodes have (chares_per_node+1) number of chares each
		big_nodes = size%num_nodes;

		//total number of chares that go on a big node
		chares_on_big_nodes = (chares_per_node+1)*big_nodes;

		//small nodes hold chares_per_node number of chares each
		chares_on_small_nodes = size-chares_on_big_nodes;

		//Number of cores on all big nodes and small nodes
		big_cores = big_nodes*total_cores_per_node;
		small_cores = (num_nodes-big_nodes)*total_cores_per_node;

		return 0;
	}

	//  int procNum(int, const CkArrayIndex &);
	inline int procNum(int, const CkArrayIndex &iIndex){
		int dim = iIndex.dimension;
		int *index=(int *) iIndex.data();

		//index for striping across lower dimensions of chare array
		int chare_num;

		if(dim == 1)
			chare_num = index[0];
		else if(dim == 2) {
			if(stripe == 0)
				chare_num = index[0]*y_size+index[1]; //stripe across y-dimension
			else
				chare_num = index[0]+x_size*index[1]; //stripe across x-dimension
		}
		else if(dim == 3)
			chare_num = index[0]+x_size*index[1]+x_size*y_size*index[2];
		else
			CkAbort("NodeMap cannot handle Chare arrays with more than 3 dimensions\n");

		return getProc(chare_num);
	}

	~NodeMap2DArray(){}

	inline int getProc(int chare_num) {
		//Calculate what proc and node this would go to for a block mapping scheme spread across all processors
		//If this chare belongs on a big node:
		int block_proc;
		if(chare_num < chares_on_big_nodes)
			block_proc = (float)chare_num/chares_on_big_nodes*big_cores;
		//otherwise it belongs on a small node, note that small nodes start after big cores, so they are offset by big_cores
		else
			block_proc = ((float)chare_num-chares_on_big_nodes)/chares_on_small_nodes*small_cores + big_cores;

		//node assignment of a chare
		int my_node = block_proc/total_cores_per_node;

		//Calculate the chare_num index relative to nodes
		int node_index;
		//If this chare belongs on a big node
		if(chare_num < chares_on_big_nodes)
			node_index = chare_num % (chares_per_node + 1);
		else
			node_index = chare_num % chares_per_node;

		//Calculate which proc within a node this chare is assigned to including all offsets
		int node_proc = ( ( (node_index + offset) % cores_per_node) + core_offset) % total_cores_per_node;

		//calculate final proc index
		int proc = node_proc + ( (my_node + node_offset) % total_num_nodes) * total_cores_per_node;

#ifdef CRAYDEBUG
		CkPrintf("chare_num %d mapped to %d globally\n",chare_num,proc);
#endif

		CkAssert(proc>=0 && proc<CkNumPes());
		return proc;
	}

	int x_size;
	int y_size;
	int size;
	int offset;
	int stripe;
	int core_offset;
	int cores_per_node;
	int total_cores_per_node;
	int total_num_nodes;
	int node_offset;
	int num_nodes;
	int chares_per_node;
	int big_nodes;
	int big_cores;
	int small_cores;
	int chares_on_big_nodes;
	int chares_on_small_nodes;
};

class BlockMapPCArray: public CkArrayMap {

public:
	//Need the size of one dimension, and the total size
	BlockMapPCArray(int _numPlanes, int _numStates, int _numChunks, int _grainSize, bool _symm){

		numPlanes = _numPlanes;
		numStates = _numStates;
		numChunks = _numChunks;
		grainSize = _grainSize;
		actualStateSize = numStates/_grainSize;
		size = _numPlanes*actualStateSize*actualStateSize*_numChunks;

		map = new int[size];

		int count = 0;
		int pcMaxStateDimIndex = (_numStates / _grainSize - 1) * _grainSize;

		//The nesting of these loops determines what is adjacent
		//Currently chunks, then planes, then s2, then s1
		//meaning chunks are closest, and s1 states are spread apart
		//effectively an entire s1 state is grouped together when linearized
		for (int s1 = 0; s1 <= pcMaxStateDimIndex; s1 += _grainSize)
		{
			// Make the symmetric array a triangular prism of chares only if phantoms are not needed
			int s2start = _symm ? s1 : 0;
			for (int s2 = s2start; s2 <= pcMaxStateDimIndex; s2 += _grainSize)
			{
				for(int numX = 0; numX < _numPlanes; numX ++)
				{
					for (int c = 0; c < _numChunks; c++)
					{
						int index = s1/grainSize*_numChunks*_numPlanes*actualStateSize+s2/grainSize*_numChunks*_numPlanes+numX*_numChunks+c;
#ifdef CRAYDEBUG
						CkPrintf("Accessing map[%d] with map size = %d\n",index, size);
						CkPrintf("Accessing s1=%d, s2=%d, numX=%d, c=%d\n",s1,s2,numX,c);
						CkPrintf("actualStateSize=%d, numStates=%d, numPlanes=%d, numChunks=%d\n",actualStateSize,numStates,numPlanes,numChunks);
#endif
						map[index] = count++;
					}
				}
			}
		}
	}

	//  int procNum(int, const CkArrayIndex &);
	inline int procNum(int, const CkArrayIndex &iIndex){
		int dim = iIndex.dimension;
		short *index=(short *) iIndex.data();

		//index for striping across lower dimensions of chare array
		int chare_num;

		if(dim == 4)
			//Chunks adjacent, planes adjacent, state2 adjacent, then state1 adjacent
			chare_num = index[1]/grainSize*numChunks*numPlanes*actualStateSize+index[2]/grainSize*numChunks*numPlanes+index[0]*numChunks+index[3];
		else
			CkAbort("BlockMapPC cannot handle Chare arrays != 4 dimensions - this mapping scheme is made for PairCalcs\n");

		int proc=(float)map[chare_num]/size*CkNumPes();
		CkAssert(proc>=0);
		return(proc);
	}

	~BlockMapPCArray(){}

private:
	int numPlanes;
	int numStates;
	int numChunks;
	int grainSize;
	int actualStateSize;
	int size;
	int* map;
};

class NodeMapPCArray: public NodeMap2DArray {
public:
	NodeMapPCArray(int _numPlanes, int _numStates, int _numChunks, int _grainSize, bool _symm, int _num_cores_to_use_per_node, int _core_offset, int _num_nodes_to_use, int _node_offset)
	: NodeMap2DArray(0, _num_cores_to_use_per_node, _core_offset, _num_nodes_to_use, _node_offset, 0) {

		numPlanes = _numPlanes;
		numStates = _numStates;
		numChunks = _numChunks;
		grainSize = _grainSize;
		actualStateSize = numStates/_grainSize;
		size = _numPlanes*actualStateSize*actualStateSize*_numChunks;

		map = new int[size];

		int count = 0;
		int pcMaxStateDimIndex = (_numStates / _grainSize - 1) * _grainSize;

		//The nesting of these loops determines what is adjacent
		//Currently chunks, then planes, then s2, then s1
		//meaning chunks are closest, and s1 states are spread apart
		//effectively an entire s1 state is grouped together when linearized
		for (int s1 = 0; s1 <= pcMaxStateDimIndex; s1 += _grainSize)
		{
			// Make the symmetric array a triangular prism of chares only if phantoms are not needed
			int s2start = _symm ? s1 : 0;
			for (int s2 = s2start; s2 <= pcMaxStateDimIndex; s2 += _grainSize)
			{
				for(int numX = 0; numX < _numPlanes; numX ++)
				{
					for (int c = 0; c < _numChunks; c++)
					{
						int index = s1/grainSize*_numChunks*_numPlanes*actualStateSize+s2/grainSize*_numChunks*_numPlanes+numX*_numChunks+c;
#ifdef CRAYDEBUG
						CkPrintf("Accessing map[%d] with map size = %d\n",index, size);
						CkPrintf("Accessing s1=%d, s2=%d, numX=%d, c=%d\n",s1,s2,numX,c);
						CkPrintf("actualStateSize=%d, numStates=%d, numPlanes=%d, numChunks=%d\n",actualStateSize,numStates,numPlanes,numChunks);
#endif
						map[index] = count++;
					}
				}
			}
		}

		chares_per_node = size/num_nodes;

		//number of nodes that need an extra chare if size/num_nodes has remainder, a big node
		//big nodes have (chares_per_node+1) number of chares each
		big_nodes = size%num_nodes;

		//total number of chares that go on a big node
		chares_on_big_nodes = (chares_per_node+1)*big_nodes;

		//small nodes hold chares_per_node number of chares each
		chares_on_small_nodes = size-chares_on_big_nodes;

		//Number of cores on all big nodes and small nodes
		big_cores = big_nodes*total_cores_per_node;
		small_cores = (num_nodes-big_nodes)*total_cores_per_node;
	}

	int registerArray(CkArrayIndexMax& numElements,CkArrayID aid) {}

	//  int procNum(int, const CkArrayIndex &);
	inline int procNum(int, const CkArrayIndex &iIndex){
		int dim = iIndex.dimension;
		short *index=(short *) iIndex.data();

		//index for striping across lower dimensions of chare array
		int chare_num;

		if(dim == 4)
			//Chunks adjacent, planes adjacent, state2 adjacent, then state1 adjacent
			chare_num = index[1]/grainSize*numChunks*numPlanes*actualStateSize+index[2]/grainSize*numChunks*numPlanes+index[0]*numChunks+index[3];
		else
			CkAbort("NodeMapPC cannot handle Chare arrays != 4 dimensions - this mapping scheme is made for PairCalcs\n");

		return getProc(map[chare_num]);
	}

private:
	int numPlanes;
	int numStates;
	int numChunks;
	int grainSize;
	int actualStateSize;
	int *map;

};

class BlockMapOrthoArray: public CkArrayMap {

public:
	//Need the size of one dimension, and the total size
	BlockMapOrthoArray(int _numStates, int _grainSize){

		numStates = _numStates;
		grainSize = _grainSize;
		actualStateSize = numStates/_grainSize;
		size = actualStateSize*actualStateSize;

		map = new int[size];

		int count = 0;
		int maxorthostateindex=(numStates/grainSize-1) * grainSize;

		//The nesting of these loops determines what is adjacent
		//s2, then s1 are placed adjacent
		//effectively an entire s1 state is grouped together when linearized
		//this linearization matches the PC mapping
		for (int s1 = 0; s1 <= maxorthostateindex; s1 += grainSize)
		{
			for (int s2 = 0; s2 <= maxorthostateindex; s2 += grainSize)
			{
				int index = s1/grainSize*actualStateSize+s2/grainSize;
				map[index] = count++;
			}
		}
	}

	//  int procNum(int, const CkArrayIndex &);
	inline int procNum(int, const CkArrayIndex &iIndex){
		int dim = iIndex.dimension;
		int *index=(int *) iIndex.data();

		//index for striping across lower dimensions of chare array
		int chare_num;

		if(dim == 2)
			//state2 adjacent, then state1 adjacent
			chare_num = index[0]/grainSize*actualStateSize+index[1]/grainSize;
		else
			CkAbort("BlockMapOrtho cannot handle Chare arrays != 2 dimensions - this mapping scheme is made for Orthos\n");

		int proc=(float)map[chare_num]/size*CkNumPes();
		CkAssert(proc>=0);
		return(proc);
	}

	~BlockMapOrthoArray(){}

private:
	int numStates;
	int grainSize;
	int actualStateSize;
	int size;
	int* map;
};

class NodeMapOrthoArray: public NodeMap2DArray {
public:
	NodeMapOrthoArray(int _numStates, int _grainSize, int _num_cores_to_use_per_node, int _core_offset, int _num_nodes_to_use, int _node_offset)
	: NodeMap2DArray(0, _num_cores_to_use_per_node, _core_offset, _num_nodes_to_use, _node_offset, 0) {

		numStates = _numStates;
		grainSize = _grainSize;
		actualStateSize = numStates/_grainSize;
		size = actualStateSize*actualStateSize;

		map = new int[size];

		int count = 0;
		int maxorthostateindex=(numStates/grainSize-1) * grainSize;

		//The nesting of these loops determines what is adjacent
		//s2, then s1 are placed adjacent
		//effectively an entire s1 state is grouped together when linearized
		//this linearization matches the PC mapping
		for (int s1 = 0; s1 <= maxorthostateindex; s1 += grainSize)
		{
			for (int s2 = 0; s2 <= maxorthostateindex; s2 += grainSize)
			{
				int index = s1/grainSize*actualStateSize+s2/grainSize;
				map[index] = count++;
			}
		}

		chares_per_node = size/num_nodes;

		//number of nodes that need an extra chare if size/num_nodes has remainder, a big node
		//big nodes have (chares_per_node+1) number of chares each
		big_nodes = size%num_nodes;

		//total number of chares that go on a big node
		chares_on_big_nodes = (chares_per_node+1)*big_nodes;

		//small nodes hold chares_per_node number of chares each
		chares_on_small_nodes = size-chares_on_big_nodes;

		//Number of cores on all big nodes and small nodes
		big_cores = big_nodes*total_cores_per_node;
		small_cores = (num_nodes-big_nodes)*total_cores_per_node;
	}

	int registerArray(CkArrayIndexMax& numElements,CkArrayID aid) {}

	//  int procNum(int, const CkArrayIndex &);
	inline int procNum(int, const CkArrayIndex &iIndex){
		int dim = iIndex.dimension;
		int *index=(int *) iIndex.data();

		//index for striping across lower dimensions of chare array
		int chare_num;

		if(dim == 2)
			//State 2 is adjacent, then state 1 adjacent to match PairCalc mapping
			chare_num = index[0]/grainSize*actualStateSize+index[1]/grainSize;
		else
			CkAbort("NodeMapOrtho cannot handle Chare arrays != 2 dimensions - this mapping scheme is made for Orthos\n");

		return getProc(map[chare_num]);
	}

private:
	int numStates;
	int grainSize;
	int actualStateSize;
	int *map;
};

class BlockMap2DArray: public CkArrayMap {

public:
	//Need the size of one dimension, and the total size
	BlockMap2DArray(int _y_size, int _size){
		y_size = _y_size;
		size = _size;
	}

	//  int procNum(int, const CkArrayIndex &);
	inline int procNum(int, const CkArrayIndex &iIndex){
		int *index=(int *) iIndex.data();

		int proc=(float)(index[0]*y_size+index[1])/size*CkNumPes();
		CkAssert(proc>=0);
		return(proc);
	}

	~BlockMap2DArray(){}

private:
	int y_size;
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

