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
	return (CkHashCode)((x<<10)+y);
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
/**
 * Abstract base class.
 *  
 */
class MapTable 
{
 public:
  CkHashtableT <intdual, int > *maptable;
  PeList *availprocs;
  void dump()
    {
      CkHashtableIterator *it=maptable->iterator();
      it->seekStart();
      CkPrintf("Map dump\n");
      intdual *key;
      while(it->hasNext())
	{
	  it->next((void **) &key);
	  int proc =maptable->get(key[0]);
	  CkPrintf("%d %d %d\n", key[0].getx(), key[0].gety(),proc);
	}
      delete it;
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

  GSMapTable(CkHashtableT <intdual, int > *_map, PeList *_availprocs, int _nchareG,
	       double *_lines_per_chareG, double *_pts_per_chareG, int _nstates,  
	       int _Gstates_per_pe);

  GSMapTable()
    {
    }
    
};

class SCalcMapTable : public MapTable
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
    SCalcMapTable(CkHashtableT <intdual, int> *_map, PeList *_availprocs, int _nstates, 
	     int _nchareG,  int gs, CmiBool _flag, int _nplanes,
	     double *_lines_per_chareG, double *_pts_per_chareG, int _scalc_per_plane,
	     int _planes_per_pe, int _numChunksA, int _numChunksS);

  void dump()
    {
      CkHashtableIterator *it=maptable->iterator();
      it->seekStart();
      CkPrintf("Map dump\n");
      intdual *key;
      while(it->hasNext())
	{
	  it->next((void **) &key);
	  int proc =maptable->get(key[0]);
	  short *four=(short*) key;
	  CkPrintf("%d %d %d %d %d\n", four[0], four[1], four[2], four[3],proc);
	}
      delete it;
    }

};

class RSMapTable  : public MapTable
{
 public:
  int nstates;
  int sizeY;
  int Rstates_per_pe;
  RSMapTable(CkHashtableT <intdual, int > *_map, PeList *_availprocs,
	int _nstates, int _sizeY, int _Rstates_per_pe) ;
  RSMapTable(){}
};

class RhoRSMapTable  : public MapTable
{
 public:
  int nchareRhoR;
  RhoRSMapTable(CkHashtableT <intdual, int > *_map, PeList *_availprocs,
	int _nchareRhoR);
  RhoRSMapTable(){}
};

class RhoGHartMapTable  : public MapTable
{
 public:
  int nchareRhoGHart;
  RhoGHartMapTable(CkHashtableT <intdual, int > *_map, PeList *_availprocs,
	int _nchareRhoGHart);
  RhoGHartMapTable(){}
};

class RhoGSMapTable  : public MapTable
{
 public:
  int nchareRhoG;
  RhoGSMapTable(CkHashtableT <intdual, int > *_map, PeList *_availprocs,
	int _nchareRhoG);
  RhoGSMapTable(){}
};


