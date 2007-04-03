/** \file MapTable.h
 *  Author: Eric J Bohm
 *  Date Created: June 4th, 2006
 *
 *  Given all necessary inputs, the MapTable creates a ckhashtable
 *  or an int array which maps ckarrayindices to processors.
 *
 *  Subclasses are made for each ckarray which needs different input
 *  in its partitioning scheme.
 *
 *  Maptables are sequential objects which can be used in a readonly
 *  global or local context.  Essentially just a factory for creating
 *  CkHashtable <intdual, int> maps for use by CkArrayMap::procnum
 *  functions.
 */
#ifndef _MAPTABLE_H_
#define _MAPTABLE_H_
#include "../../include/debug_flags.h"

class inttriple {
 private:
    int x, y, z;
 public:
    inttriple(){x=y=z=0;}

    inttriple(int _x,int _y,int _z) : x(_x), y(_y), z(_z) {}
    void pup(PUP::er &p)
      {
	  p|x;
	  p|y;
	  p|z;
      }
    inline int getx(){return x;};
    inline int gety(){return y;};
    inline int getz(){return z;};
    // silenty assumes that X is the heavy hitter in keyspace
    // safe for our purposes
    inline CkHashCode hash() const {
	return (CkHashCode)((x<<16)|(y<<8)|z);
    }
    static CkHashCode staticHash(const void *k,size_t){
	return ((inttriple *)k)->hash();
    }
    inline int compare(inttriple &t) const{
	return (t.getx() == x && t.gety() == y && t.getz() ==z);
    }
    static int staticCompare(const void *a,const void *b,size_t){
	return ((inttriple *)a)->compare((*(inttriple *)b));
    }
   
};

class intdual {
 private:
    int x, y;
 public:
    intdual(){x=y=0;}

    intdual(int _x,int _y) : x(_x), y(_y){}
    void pup(PUP::er &p)
      {
	  p|x;
	  p|y;
      }
    inline int getx(){return x;};
    inline int gety(){return y;};
    inline CkHashCode hash() const {
	return (CkHashCode)((x<<10)|y);
    }
    static CkHashCode staticHash(const void *k,size_t){
	return ((intdual *)k)->hash();
    }
    inline int compare(intdual &t) const{
	return (t.getx() == x && t.gety() == y);
    }
    static int staticCompare(const void *a,const void *b,size_t){
	return ((intdual *)a)->compare((*(intdual *)b));
    }
   
};


#ifndef USE_INT_MAP
//typedef CkHashtableT <intdual, int > MapType2;
typedef CkHashtableT <intdual, int > MapType4;
typedef CkHashtableT <inttriple, int > MapType3;
#else
#endif


/**
 * Abstract base class.
 *  
 */
class MapTable2 
{
 public:
  MapType2 *maptable;
  PeList *availprocs;
  void dump()
    {
#ifndef USE_INT_MAP
      CkHashtableIterator *it=maptable->iterator();
      it->seekStart();
      CkPrintf("Map dump\n");
      intdual *key;
      while(it->hasNext())
	{
	  it->next((void **) &key);
	  int proc =maptable->get(key[0]);
#ifdef _MAP_VERBOSE_
	  CkPrintf("%d %d %d\n", key[0].getx(), key[0].gety(),proc);
#endif
	}
      delete it;
#else
      maptable->dump();
#endif
    }

 protected:
   CkVec <intdual> *reverseMap;
   
  MapTable2()
    {
      availprocs=NULL;
      maptable=NULL;
      reverseMap=NULL;
    }
  ~MapTable2()
    {
      if(reverseMap!=NULL)
	delete [] reverseMap;
    }
  void makeReverseMap();
  
  /**
   * return ckvec containing the  reverse map of  all elements on
   *  given proc 
   */
  inline CkVec <intdual>  ProcByArrIndex(int proc) 
    { 
      if(reverseMap==NULL)
	makeReverseMap();
      return(reverseMap[proc]); 
    }
  
  //! return processor at topological center of this list
  int getCentroid();
  
};

class MapTable3
{
 public:
  MapType3 *maptable;
  PeList *availprocs;
  void dump()
    {
#ifndef USE_INT_MAP
      CkHashtableIterator *it=maptable->iterator();
      it->seekStart();
      CkPrintf("Map dump\n");
      inttriple *key;
      while(it->hasNext())
	{
	  it->next((void **) &key);
	  int proc =maptable->get(key[0]);
#ifdef _MAP_VERBOSE_
	  CkPrintf("%d %d %d %d\n", key[0].getx(), key[0].gety(), key[0].getz(),proc);
#endif
	}
      delete it;
#else
      maptable->dump();
#endif
    }

 protected:
   CkVec <inttriple> *reverseMap;
   
  MapTable3()
    {
      availprocs=NULL;
      maptable=NULL;
      reverseMap=NULL;
    }
  ~MapTable3()
    {
      if(reverseMap!=NULL)
	delete [] reverseMap;
    }
  void makeReverseMap();
  
  /**
   * return ckvec containing the  reverse map of  all elements on
   *  given proc 
   */
  inline CkVec <inttriple>  ProcByArrIndex(int proc) 
    { 
      if(reverseMap==NULL)
	makeReverseMap();
      return(reverseMap[proc]); 
    }
  
  //! return processor at topological center of this list
  int getCentroid();
  
};


/**
 * Abstract base class.
 *  
 */
class MapTable4 
{
 public:
  MapType4 *maptable;
  PeList *availprocs;
  void dump()
    {
#ifndef USE_INT_MAP
      CkHashtableIterator *it=maptable->iterator();
      it->seekStart();
      CkPrintf("Map dump\n");
      intdual *key;
      while(it->hasNext())
	{
	  it->next((void **) &key);
	  int proc =maptable->get(key[0]);
#ifdef _MAP_VERBOSE_
	  CkPrintf("%d %d %d\n", key[0].getx(), key[0].gety(),proc);
#endif
	}
      delete it;
#else
      maptable->dump();
#endif
    }

 protected:
   CkVec <intdual> *reverseMap;

  MapTable4()
    {
      availprocs=NULL;
      maptable=NULL;
      reverseMap=NULL;
    }
  ~MapTable4()
    {
      if(reverseMap!=NULL)
	delete [] reverseMap;
    }
  void makeReverseMap();
  
  /**
   * return ckvec containing the  reverse map of  all elements on
   *  given proc 
   */
  inline CkVec <intdual>  ProcByArrIndex(int proc) 
    { 
      if(reverseMap==NULL)
	makeReverseMap();
      return(reverseMap[proc]); 
    }


};

PeList *subListPlane(int plane, int nstates, MapType2 *smap);
PeList *subListState(int state, int nplanes, MapType2 *smap);
PeList *subListState2(int state1, int state2, int nplanes, int numChunks, MapType4 *smap);


class GSMapTable : public MapTable2
{
 public:

  int nchareG;
  int nstates;
  int Gstates_per_pe;
  
  double state_load;
  int planes_per_pe;

  GSMapTable(MapType2  *_map, PeList *_availprocs, int _nchareG,
	       int _nstates, int _Gstates_per_pe, bool useCuboidMap);

  GSMapTable()
    {
    }
    
};


class SCalcMapTable : public MapTable4
{

  int max_states, nchareG, grainsize;
  CmiBool symmetric;
  int scalc_per_plane;
  int planes_per_pe;
  int numChunksAsym;
  int numChunksSym;
  double totalload;
    
 public:
    SCalcMapTable(MapType4  *_map, PeList *_availprocs, int _nstates, 
	     int _nchareG,  int gs, CmiBool _flag, int _scalc_per_plane,
	     int _planes_per_pe, int _numChunksA, int _numChunksS, MapType2  *_gmap, bool useCuboidMap, bool useCentroid, int boxSize);

  void dump()
    {
#ifndef USE_INT_MAP
      CkHashtableIterator *it=maptable->iterator();
      it->seekStart();
      CkPrintf("Map dump\n");
      intdual *key;
      while(it->hasNext())
	{
	  it->next((void **) &key);
	  int proc =maptable->get(key[0]);
	  short *four=(short*) key;
#ifdef _MAP_VERBOSE_
	  CkPrintf("%d %d %d %d %d\n", four[0], four[1], four[2], four[3],proc);
#endif
	}
      delete it;
#else
      maptable->dump();
#endif
    }
  void sortByCentroid(PeList *avail, int plane, int stateX, int stateY, int grainsize, MapType2 *gsmap);
  
};

class RSMapTable  : public MapTable2
{
 public:
  int nstates;
  int sizeZ;
  int Rstates_per_pe;
  RSMapTable(MapType2  *_map, PeList *_availprocs,
	int _nstates, int _sizeZ, int _Rstates_per_pe, bool useCuboid, MapType2 *gsmap, int nchareG);
  RSMapTable(){}
};


class RPPMapTable  : public MapTable2
{
 public:
  int nstates;
  int sizeZNL;
  int Rstates_per_pe;
  RPPMapTable(MapType2  *_map, PeList *_availprocs,
	      PeList *exclude, 	int _nstates, int _sizeZNL, 
	      int _Rstates_per_pe, int boxSize, bool useCuboidMap, 
	      int nchareG, MapType2 *ppmap) ;
  RPPMapTable(){}
};


class OrthoMapTable : public MapTable2
{
  public:
    int nstates;
    int orthoGrainSize;

    OrthoMapTable(MapType2 *_map, PeList *_availprocs, int _nstates, int _orthograinsize, MapType4 *scalcmap, int nplanes, int numChunks, int sGrainSize, PeList *);
    OrthoMapTable() { }
    void sortByCentroid(PeList *avail, int nplanes, int state1, int state2, int numChunks, MapType4 *smap);
};

class OrthoHelperMapTable : public MapTable2
{
  public:
    int nstates;
    int orthoGrainSize;

    OrthoHelperMapTable(MapType2 *_map, int _nstates, int _orthograinsize, MapType2 *omap, PeList*, PeList*);
    OrthoHelperMapTable() { }
};

class RhoRSMapTable  : public MapTable2
{
 public:
  int nchareRhoR;
  int rhoRsubplanes;
  RhoRSMapTable(MapType2  *_map, PeList *_availprocs,
	int _nchareRhoR, int _rhoRsubplanes, int maxstates, bool useCentroid, MapType2 *rsmap, PeList *exclude);
  void sortByCentroid(PeList *avail, int plane, int nstates, MapType2 *rsmap);
  RhoRSMapTable(){}
};


class RhoRHartMapTable  : public MapTable3
{
 public:
  int nchareRhoRHart;
  RhoRHartMapTable(MapType3  *_map, PeList *_availprocs,
	int _nchareRhoRHart,  int rhoRsubplanes, int _nchareHartAtmT,
		   PeList *exclude);
  RhoRHartMapTable(){}
};

class RhoGHartMapTable  : public MapTable2
{
 public:
  int nchareRhoGHart;
  RhoGHartMapTable(MapType2  *_map, PeList *_availprocs,
	int _nchareRhoGHart, int _nchareHartAtmT, int useCentroid, 
		   MapType3 *rhartmap,   PeList *exclude);
  RhoGHartMapTable(){}
};

class RhoGSMapTable  : public MapTable2
{
 public:
  int nchareRhoG;
  RhoGSMapTable(MapType2  *_map, PeList *_availprocs,
	int _nchareRhoG,  bool useCentroid, MapType2 *rhorsmap, PeList *exclude);
  RhoGSMapTable(){}
};


#endif
