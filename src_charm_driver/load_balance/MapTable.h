/** \file MapTable.h
 *  Author: Eric J Bohm
 *  Date Created: June 4th, 2006
 *
 *  Given all necessary inputs, the MapTable creates a ckhashtable
 *  which maps ckarrayindices to processors.
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
//#define  _VERBOSE_MAP_OFF_
#define  _VERBOSE_MAP_

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
typedef CkHashtableT <intdual, int > MapType2;
typedef CkHashtableT <intdual, int > MapType4;
#else
#endif


/**
 * Abstract base class.
 *  
 */
class MapTable 
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
#ifdef _VERBOSE_MAP_
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

  MapTable()
    {
      availprocs=NULL;
      maptable=NULL;
      reverseMap=NULL;
    }
  ~MapTable()
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
#ifdef _VERBOSE_MAP_
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

class GSMapTable : public MapTable
{
 public:

  int nchareG;
  int nstates;
  int Gstates_per_pe;
  
  double *lines_per_chareG;
  double *pts_per_chareG;
  double state_load;
  int planes_per_pe;

  GSMapTable(MapType2  *_map, PeList *_availprocs, int _nchareG,
	       double *_lines_per_chareG, double *_pts_per_chareG, int _nstates,  
	       int _Gstates_per_pe, bool useCuboidMap);

  GSMapTable()
    {
    }
    
};

class SCalcMapTable : public MapTable4
{

  int max_states,nchareG, grainsize;
  CmiBool symmetric;
  int max_planes;
  double *lines_per_chareG;
  double *pts_per_chareG;
  int scalc_per_plane;
  int planes_per_pe;
  int numChunksAsym;
  int numChunksSym;
  double totalload;
    
 public:
    SCalcMapTable(MapType4  *_map, PeList *_availprocs, int _nstates, 
	     int _nchareG,  int gs, CmiBool _flag, int _nplanes,
	     double *_lines_per_chareG, double *_pts_per_chareG, int _scalc_per_plane,
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
#ifdef _VERBOSE_MAP_
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

class RSMapTable  : public MapTable
{
 public:
  int nstates;
  int sizeZ;
  int Rstates_per_pe;
  RSMapTable(MapType2  *_map, PeList *_availprocs,
	int _nstates, int _sizeZ, int _Rstates_per_pe, bool useCuboid, MapType2 *gsmap, int nchareG);
  RSMapTable(){}
};


class RPPMapTable  : public MapTable
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


class RhoRSMapTable  : public MapTable
{
 public:
  int nchareRhoR;
  int rhoRsubplanes;
  RhoRSMapTable(MapType2  *_map, PeList *_availprocs,
	int _nchareRhoR, int _rhoRsubplanes, int maxstates, bool useCentroid, MapType2 *rsmap, PeList *exclude);
  void sortByCentroid(PeList *avail, int plane, int nstates, MapType2 *rsmap);
  RhoRSMapTable(){}
};


class RhoRHartMapTable  : public MapTable
{
 public:
  int nchareRhoRHart;
  RhoRHartMapTable(MapType2  *_map, PeList *_availprocs,
	int _nchareRhoRHart, int rhoRsubplanes, PeList *exclude);
  RhoRHartMapTable(){}
};

class RhoGHartMapTable  : public MapTable
{
 public:
  int nchareRhoGHart;
  RhoGHartMapTable(MapType2  *_map, PeList *_availprocs,
	int _nchareRhoGHart, PeList *exclude);
  RhoGHartMapTable(){}
};

class RhoGSMapTable  : public MapTable
{
 public:
  int nchareRhoG;
  RhoGSMapTable(MapType2  *_map, PeList *_availprocs,
	int _nchareRhoG, PeList *exclude);
  RhoGSMapTable(){}
};


#endif
