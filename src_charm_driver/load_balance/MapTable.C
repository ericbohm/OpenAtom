/*****************************************************************************
 * $Source$
 * $Author$
 * $Date$
 * $Revision$
 *****************************************************************************/

/** \file MapTable.C
 *
 */

#include "charm++.h"
#include "PeList.h"
#include "utility/MapFile.h"
#define USE_INT_MAP 1
#include "IntMap.h"
#include "MapTable.h"
#define DO_SORT 0
#include "fft_charm.h"
#include "include/CPcharmParaInfo.h"
extern TopoManager *topoMgr;
extern Config config;
extern CPcharmParaInfo                          simReadOnly;
#define TIMER_SET(a)  

PeList *subListPlanes(int start_plane, int end_plane, int nstates, MapType2 *smap);

//============================================================================
int MapType1::getCentroid(int torusMap) {
  int matchPe, dims[10];
  getCentroid(torusMap, dims);
  if(config.simpleTopo) {
    TopoManager_getPeRank(&matchPe, dims);
  } else {
    matchPe = dims[0];
  }
  return matchPe;
}

void MapType1::getCentroid(int torusMap, int *rdims) {
  int points;
  if(config.simpleTopo) {
    int ndims;
    int sum[10], dim[10];
    points=0;
    TopoManager_getDimCount(&ndims);
    for(int di = 0; di < ndims; di++) {
      sum[di] = 0;
    }
    for(int i=0;i<getXmax();i++)
    {
      CkAssert(get(i)>=0);
      TopoManager_getPeCoordinates(get(i), dim);
      for(int di = 0; di < ndims; di++) {
        sum[di] += dim[di];
      }
      points++;
    }
    for(int di = 0; di < ndims; di++) {
      rdims[di] = sum[di]/points;
    }
    rdims[ndims] = 0;
  } else {
    int sum=0;
    points=0;
    for(int i=0;i<getXmax();i++)
    {
      CkAssert(get(i)>=0);
      sum+=get(i);
      points++;
    }
    rdims[0]=sum/points;
  }
}

int MapType2::getCentroid(int torusMap) {
  int matchPe, dims[10];
  getCentroid(torusMap, dims);
  if(config.simpleTopo) {
    TopoManager_getPeRank(&matchPe, dims);
  } else {
    matchPe = dims[0];
  }
  return matchPe;
}

void MapType2::getCentroid(int torusMap, int *rdims) {
  int points;
  if(config.simpleTopo) {
    int ndims;
    int sum[10], dim[10];
    points=0;
    TopoManager_getDimCount(&ndims);
    for(int di = 0; di < ndims; di++) {
      sum[di] = 0;
    }
    for(int i=0;i<getXmax();i++)
      for(int j=0;j<getYmax();j++) {
        CkAssert(get(i,j)>=0);
        TopoManager_getPeCoordinates(get(i,j), dim);
        for(int di = 0; di < ndims; di++) {
          sum[di] += dim[di];
        }
        points++;
      }
    for(int di = 0; di < ndims; di++) {
      rdims[di] = sum[di] / points;
    }
    rdims[ndims] = 0;
  }
  else {
    int sum=0;
    points=0;
    for(int i=0;i<getXmax();i++)
      for(int j=0;j<getYmax();j++) {
        CkAssert(get(i,j)>=0);
        sum+=get(i,j);
        points++;
      }
    rdims[0] = sum/points;
  }
}

int MapType3::getCentroid(int torusMap) {
  int matchPe, dims[10];
  getCentroid(torusMap, dims);
  if(config.simpleTopo) {
    TopoManager_getPeRank(&matchPe, dims);
  } else {
    matchPe = dims[0];
  }
  return matchPe;
}

void MapType3::getCentroid(int torusMap, int *rdims) {
  int points;
  if(config.simpleTopo) {
    int ndims;
    int sum[10], dim[10];
    points=0;
    TopoManager_getDimCount(&ndims);
    for(int di = 0; di < ndims; di++) {
      sum[di] = 0;
    }
    for(int i=0;i<getXmax();i++)
      for(int j=0;j<getYmax();j++) {
        for(int k=0;k<getZmax();k++) {
          CkAssert(get(i,j,k)>=0);
          TopoManager_getPeCoordinates(get(i,j,k), dim);
          for(int di = 0; di < ndims; di++) {
            sum[di] += dim[di];
          }
          points++;
        }
      }
    for(int di = 0; di < ndims; di++) {
      rdims[di] = sum[di] / points;
    }
    rdims[ndims] = 0;
  }
  else {
    int sum=0;
    points=0;
    for(int i=0;i<getXmax();i++)
      for(int j=0;j<getYmax();j++){
        for(int k=0;k<getZmax();k++){
          CkAssert(get(i,j,k)>=0);
          sum+=get(i,j,k);
          points++;
        }
      }
    rdims[0] = sum/points;
  }
}

void IntMap1::translate(IntMap1 *fromMap, int offsetX, int offsetY, int offsetZ, bool torus )
{

  keyXmax=fromMap->keyXmax;
  CkAssert(keyXmax>0);
  Map= new int[keyXmax];
  for(int xind=0; xind<keyXmax; xind++)
    set(xind,(fromMap->get(xind)+offsetX)%config.numPes);
};


void IntMap2on2::translate(IntMap2on2 *fromMap, int offsetX, int offsetY, int offsetZ, bool torus )
{
  keyXmax=fromMap->keyXmax;
  keyYmax=fromMap->keyYmax;
  CkAssert(keyXmax>0);
  CkAssert(keyYmax>0);
  CkAssert(keyXmax<10000000);
  CkAssert(keyYmax<10000000);
  Map= new int*[keyXmax];

  int *mapbuf=new int[keyXmax*keyYmax];
  for(int x=0;x<keyXmax;x++)
  {
    Map[x]  =  mapbuf +  keyYmax * x;
    memset(Map[x],-1,keyYmax*sizeof(int));
  }

  for(int xind=0; xind<keyXmax; xind++)
    for(int yind=0; yind<keyYmax; yind++) {
      set(xind, yind,(fromMap->get(xind, yind)+offsetX)%config.numPes);
    }
};

void IntMap3::translate(IntMap3 *fromMap, int offsetX, int offsetY, int offsetZ, bool torus )
{

  keyXmax=fromMap->keyXmax;
  keyYmax=fromMap->keyYmax;
  keyZmax=fromMap->keyZmax;
  CkAssert(keyXmax>0);
  CkAssert(keyYmax>0);
  CkAssert(keyZmax>0);
  Map=new int**[keyXmax];
  int **mappointbuf = new int*[keyXmax*keyYmax];
  int *mapbuf= new int[keyXmax*keyYmax*keyZmax];
  for(int x=0;x<keyXmax;x++)
  {
    Map[x]= mappointbuf + (x*keyYmax);
    for(int y=0;y<keyYmax;y++)
    {
      Map[x][y]= mapbuf + (x*keyYmax+y)*keyZmax;
      memset(Map[x][y],-1,keyZmax*sizeof(int));
    }
  }
  for(int xind=0; xind<keyXmax; xind++)
    for(int yind=0; yind<keyYmax; yind++) {
      for(int zind=0; zind<keyZmax; zind++) {
        set(xind, yind, zind, (fromMap->get(xind, yind,zind)+offsetX)%config.numPes);
      }
    }
};

void IntMap4::translate(IntMap4 *fromMap, int offsetX, int offsetY, int offsetZ, bool torus )
{
  keyWmax=fromMap->keyWmax;
  keyXmax=fromMap->keyXmax;
  keyYmax=fromMap->keyYmax;
  keyZmax=fromMap->keyZmax;
  keyStep=fromMap->keyStep;
  CkAssert(keyWmax>0);
  CkAssert(keyXmax>0);
  CkAssert(keyYmax>0);
  CkAssert(keyZmax>0);

  Map=new int***[keyWmax];
  int ***mappointpointbuf = new int**[keyWmax*keyXmax];
  int **mappointbuf = new int*[keyWmax*keyXmax*keyYmax];
  int *mapbuf= new int[keyWmax*keyXmax*keyYmax*keyZmax];
  for(int w=0;w<keyWmax;w++)
  {
    Map[w]=   mappointpointbuf + (w*keyXmax);
    for(int x=0;x<keyXmax;x++)
    {
      Map[w][x]= mappointbuf + (w*keyXmax+x)*keyYmax;
      for(int y=0;y<keyYmax;y++)
        Map[w][x][y]= mapbuf + ((w*keyXmax+x)*keyYmax+y)*keyZmax;
    }
  }
  for(int wind=0; wind<keyWmax; wind++)
    for(int xind=0; xind<keyXmax; xind++){
      for(int yind=0; yind<keyYmax; yind++) {
        for(int zind=0; zind<keyZmax; zind++) {
          set(wind, xind*keyStep, yind*keyStep, zind, (fromMap->get(wind, xind*keyStep, yind*keyStep,zind)+offsetX)%config.numPes);
        }
      }
    }

  // The stepTable is translate-invariant. Simply generate the same values
  if (keyXmax > 0 && keyStep > 0)
  {
    if (stepTable) delete [] stepTable;
    stepTable= new int [keyXmax*keyStep];
    for(int s=0; s<keyXmax*keyStep; s++)
      stepTable[s] = s/keyStep;
  }
}

AtomMapTable::AtomMapTable(MapType1 *_tomap, PeList *availprocs, int numInst,
    int _nchareAtoms): nchareAtoms(_nchareAtoms)
{
  maptable = _tomap;
  for(int element=0; element<nchareAtoms; element++)
  {
    int dest =  availprocs->findNext();
    maptable->set(element, dest);
  }  
}


FFTPencilMapTable::FFTPencilMapTable(MapType3  *_tomap, PeList *availprocs,
				     int *_dim, PeList *exclude)
{
  maptable = _tomap;
  // find the unused dimension(s)
  // optimized versions of this will probably use these
  /** future extensions should consider extensions which use the
   * placement of the input and output for the pencils.  Similarly,
   * the transpose from x -> y -> z, and back, will benefit from
   * topology aware mapping on some networks.
   */
  numFlat = 0;
  for(int i = 0; i < 3; ++i) {
    dim[i] = _dim[i];
    if(dim[i] == 1) {
      flat[i] = true;
      ++numFlat;
    } else {
      flat[i] = false;
    }
  }

  for(int d1 = 0; d1 < dim[0]; ++d1)
    for(int d2 = 0; d2 < dim[0]; ++d2)
      for(int d3 = 0; d3 < dim[0]; ++d3) {
        maptable->set(d1,d2,d3,availprocs->findNext());
      }
}

RhoYPencilMapTable::RhoYPencilMapTable(MapType2  *_tomap, PeList *_availprocs,
    int _nchareInterX, int  _nchareInterZ, bool useCentroid, MapType2 *rhorsmap,
    PeList *exclude, int offset) : nchareInter_x(_nchareInterX),
    nchareInter_z(_nchareInterZ) {
  maptable = _tomap;
  availprocs = _availprocs;
  int numchares = nchareInter_z * nchareInter_x;
  if(availprocs->count() < numchares)
    availprocs->reset();

  PeList *avail= new PeList(*availprocs);
  avail->deleteList(*exclude, 0, 0);
  PeList *pencillist;
  if(avail->count() > numchares)
  {
    if(config.simpleTopoCentroid) {
      int dims[10];
      rhorsmap->getCentroid(config.torusMap, dims);
      avail->sortSource(dims, 0);
    }
    pencillist = avail->distributeAcrossPelist(numchares, offset);
    printf("Pencil map: offset %d, avail %d, used %d\n", offset, 
      avail->size(), numchares);
  }
  else
  {
    if(config.simpleTopoCentroid) {
      // get centroid of rsmap  use it to sort the avail list
      int dims[10];
      rhorsmap->getCentroid(config.torusMap, dims);
      availprocs->sortSource(dims, 1);
    }
    pencillist = availprocs->distributeAcrossPelist(numchares, offset);
    printf("Pencil map: offset %d, avail %d, used %d\n", offset, 
      availprocs->size(), numchares);
  }
  delete avail;

  int destpe = pencillist->findNext();
  for(int x = 0; x < nchareInter_x; x++)
  {
    for(int z = 0; z < nchareInter_z; z++)
      {
	maptable->set(x, z, destpe);
	exclude->checkAndAdd(destpe);
	destpe=pencillist->findNext();
      }
  }
  delete pencillist;
#ifdef _MAP_DEBUG_
  CkPrintf("RhoGSMap created on processor %d\n", CkMyPe());
  dump();
#endif
}

GSMapTable::GSMapTable(MapType2 *_frommap, MapType2 *_tomap, PeList *_availprocs, 
    int _nchareG, int _nstates, int _Gstates_per_pe, bool useCuboidMap, int numInst) :
  nchareG(_nchareG), nstates(_nstates), Gstates_per_pe(_Gstates_per_pe)
{
  maptable = _tomap;
  availprocs = _availprocs;

  /** The first instance creates the map and the other instances just use the
   *  map with a translation
   */
  state_load = 0.0;
  int l, m, pl, pm, srem, rem, i=0;

  /** The first thing is to find the size of blocks (chunks) of GSpace chares
   *  which will be put on each processor. This depends on the number of
   *  states to be given to each processor (input by the user)
   *
   *			  <- m ->
   *		<--- nchareG --->
   *		*****************  ^  ^
   *		|		|  |  l
   *		|		|  |  _
   *		|		|
   *		|		| nstates
   *		|		|		Y
   *		|		|  |
   *		|		|  |		|
   *		|		|  |		|
   *		*****************  -		|________ X
   *
   */
  l = Gstates_per_pe;		// no of states in one chunk
  pl = nstates / l;		// no of procs on y axis
  if(nstates % l == 0)
    srem = 0;			// remainder states
  else
  {
    srem = nstates % pl;
  }
  pm = availprocs->count() / pl;		// no of procs on x axis

  if(!config.simpleTopo && pm == 0) {
    CkPrintf("Choose a larger Gstates_per_pe than %d such that { no. of processors [%d] / (no. of states [%d] / Gstates_per_pe [%d]) } is > 0 \n", 
        l, availprocs->count(), nstates, l);
    CkAssert(availprocs->count() / (nstates/l) > 0);
  }
  m = nchareG / pm;		// no of planes in one chunk
  rem = nchareG % pm;		// remainder of planes left to be mapped

  planes_per_pe=m;

  /*if(CkMyPe()==0) 
    {
    CkPrintf("nstates %d nchareG %d Pes %d\n", nstates, nchareG, availprocs->count());
    CkPrintf("l %d, m %d pl %d pm %d srem %d rem %d\n", l, m, pl, pm, srem, rem);
    }*/

  // Initialize pelist
  int srcpe=0;
  if(config.simpleTopo)
  {
    int *Pecount= new int[config.numPes];
    bzero(Pecount, config.numPes *sizeof(int));

    int procsPerPlane = config.numPesPerInstance/nchareG;
    int cubeGstates_per_pe = nstates/procsPerPlane;
    int charesperpe = nchareG*nstates/config.numPesPerInstance;
    int cubesrem = nstates%procsPerPlane;
    if(cubesrem)
      cubeGstates_per_pe++;
    if(cubesrem)
      charesperpe++;

    CkPrintf("procsPerPlane %d Gstates_per_pe %d remainder %d\n", procsPerPlane, cubeGstates_per_pe, cubesrem);
    for(int plane=0; plane<nchareG; plane++)
    {
      // slice us off our plane's processors
      // planeProcs only needs deletion of current index - hence set is not created
      PeList *planeProcs=new PeList(0, 1, *availprocs, plane*procsPerPlane, procsPerPlane);
      int destpe=planeProcs->findNext();
      int workingGsPerPe=cubeGstates_per_pe;
      bool unallocateRem= (cubesrem) ? true: false;
      // non power of two systems need some exclusion logic
      for(int state=0;state<nstates;state+=workingGsPerPe)
      {
        if(unallocateRem)
          if(state>=cubesrem)
          {

            workingGsPerPe--;
            unallocateRem=false;
            //CkPrintf("rem %d complete at state %d gsperpenow %d\n",cubesrem, state, workingGsPerPe);
          }

        // we should block these better
        for(int stateperpe=0;(stateperpe<workingGsPerPe)&&((state+stateperpe)<nstates);stateperpe++)
        {
#ifdef USE_INT_MAP
          maptable->set(state+stateperpe, plane, destpe);
#else
          maptable->put(intdual(state+stateperpe, plane))=destpe;
#endif
          if(cubesrem)
          {
            Pecount[destpe]++;
            if(((stateperpe+1<workingGsPerPe)&&((state+stateperpe+1)<nstates)) || state+workingGsPerPe<nstates)
            {
              // we will need another proc from this list
              destpe=planeProcs->findNext();
              if(destpe < 0)
              {
                CkPrintf("GSMap exceeding count on plane %d state %d after pe %d "
                    "workingGsPerPe %d cubesrem %d\n",state,plane,destpe, workingGsPerPe, 
                    cubesrem);
              }
            }
          }
        }
        if(!cubesrem)
        {
          destpe=planeProcs->findNext();
        }
      }
      delete planeProcs;
    }
    delete [] Pecount;
  }
  else
  {
    int destpe=availprocs->findNext();

    int orig_l=l;
    for(int ychunk=0; ychunk<nchareG; ychunk=ychunk+m)
    {
      if(ychunk==(pm-rem)*m)
        m=m+1;
      l=orig_l;
      for(int xchunk=0; xchunk<nstates; xchunk=xchunk+l)
      {
        if(xchunk==(pl-srem)*l)
          l=l+1;
        for(int state=xchunk; state<xchunk+l && state<nstates; state++)
        {
          for(int plane=ychunk; plane<ychunk+m && plane<nchareG; plane++)
          {
#ifdef USE_INT_MAP
            maptable->set(state, plane,destpe);
#else
            maptable->put(intdual(state, plane))=destpe;
#endif
          }
        }
        srcpe=destpe;
        destpe=availprocs->findNext();
      }
    }
  }
#ifdef _MAP_DEBUG_
  CkPrintf("GSMap created on processor %d\n", CkMyPe());
  dump();
  int size[2] = {128, 12};
  MapFile *mf = new MapFile("GSMap", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
  mf->dumpMap(maptable);
#endif
}

SCalcMapTable::SCalcMapTable(MapType4  *_map, PeList *_availprocs, 
    int _nstates, int _nchareG,  int _grainsize, 
    bool _flag, int _scalc_per_plane,
    int _planes_per_pe, 
    int _numChunksA, 
    int _numChunksS, 
    MapType2  *gsmap, bool useCuboidMap, 
    bool useCentroid, int boxSize):

  max_states(_nstates), nchareG(_nchareG),  
  grainsize(_grainsize), symmetric(_flag), 
  scalc_per_plane(_scalc_per_plane), planes_per_pe(_planes_per_pe), 
  numChunksAsym(_numChunksA), numChunksSym(_numChunksS)
{ 

  int scobjs_per_pe, rem;
  int count=0, procno=0;
  int intidx[2];
  int lesser_scalc = 0;
  maptable=_map;
  availprocs=_availprocs;
  availprocs->reset();
  int maxstateindex=max_states/grainsize*grainsize;
  if(planes_per_pe == 0)
    CkAbort("Choose a larger nChareG to avoid this crash\n");
  if(symmetric)
  {
    for(int i=1; i<=max_states/grainsize; i++)
      lesser_scalc += i;
    scobjs_per_pe = lesser_scalc*nchareG*numChunksSym/availprocs->count();
    rem = lesser_scalc*nchareG*numChunksSym % availprocs->count();
    if(rem!=0)
      scobjs_per_pe+=1;
#ifdef _MAP_DEBUG_
    CkPrintf(" lesser_scalc %d *nchareG %d *numChunksSym %d %% availprocs->count() %d = rem %d and scobjs_per_pe is %d\n", lesser_scalc,scalc_per_plane,nchareG,numChunksAsym , availprocs->count(),rem, scobjs_per_pe);
#endif

    int srcpe=0,destpe=availprocs->findNext();

      for(int pchunk=0; pchunk<nchareG; pchunk=pchunk+planes_per_pe)
        for(int newdim=0; newdim<numChunksSym; newdim++)
          for(int xchunk=0; xchunk<maxstateindex; xchunk=xchunk+grainsize)
            for(int ychunk=xchunk; ychunk<maxstateindex; ychunk=ychunk+grainsize)
              //for(int newdim=0; newdim<numChunksSym; newdim++)
              for(int plane=pchunk; plane<pchunk+planes_per_pe && plane<nchareG; plane++)
              {
                CkArrayIndex4D idx4d(plane, xchunk, ychunk, newdim);
                CmiMemcpy(intidx,idx4d.index,2*sizeof(int));  // our 4 shorts are now 2 ints

                if(count<scobjs_per_pe)
                {
                  CkAssert(destpe<config.numPes);
                  CkAssert(destpe>=0);
#ifdef USE_INT_MAP 
                  maptable->set(plane, xchunk, ychunk, newdim,destpe);
#else
                  maptable->put(intdual(intidx[0], intidx[1]))=destpe;
#endif

                  count++;
                }
                else
                {
                  // new partition
                  procno++;
                  srcpe=destpe;
                  destpe=availprocs->findNext();

                  if(rem!=0 &&scobjs_per_pe>1)
                    if(procno==rem)
                      scobjs_per_pe-=1;
                  CkAssert(destpe<config.numPes);
#ifdef USE_INT_MAP
                  maptable->set(plane, xchunk, ychunk, newdim,destpe);
#else
                  maptable->put(intdual(intidx[0], intidx[1]))=destpe;
#endif
                  count=0;
                  count++;
                }

              }
#ifdef _MAP_DEBUG_
    CkPrintf("Symmetric SCalcMap created on processor %d\n", CkMyPe());
    dump();
#endif
  }
  else
  {
    scobjs_per_pe = scalc_per_plane*nchareG*numChunksAsym/availprocs->count();
    rem = scalc_per_plane*nchareG*numChunksAsym % availprocs->count();
    if(rem!=0)
      scobjs_per_pe+=1;

    int srcpe=0,destpe=0;
    destpe=availprocs->findNext();
      for(int pchunk=0; pchunk<nchareG; pchunk=pchunk+planes_per_pe)
        for(int newdim=0; newdim<numChunksAsym; newdim++)
          for(int xchunk=0; xchunk<maxstateindex; xchunk=xchunk+grainsize)
            for(int ychunk=0; ychunk<maxstateindex; ychunk=ychunk+grainsize)
              for(int plane=pchunk; plane<pchunk+planes_per_pe && plane<nchareG; plane++)
              {
                CkArrayIndex4D idx4d(plane, xchunk, ychunk, newdim);
                CmiMemcpy(intidx,idx4d.index,2*sizeof(int));  // our 4 shorts are now 2 ints

                if(count<scobjs_per_pe)
                {
#ifdef USE_INT_MAP		      
                  maptable->set(plane, xchunk, ychunk, newdim, destpe);
#else
                  maptable->put(intdual(intidx[0], intidx[1]))=destpe;
#endif
                  count++;
                }
                else
                {  // new partition
                  procno++;
                  srcpe=destpe;
                  destpe=availprocs->findNext();
                  if(rem!=0)
                    if(procno==rem)
                      scobjs_per_pe-=1;
#ifdef USE_INT_MAP		      
                maptable->set(plane, xchunk, ychunk, newdim,destpe);
#else
                maptable->put(intdual(intidx[0], intidx[1]))=destpe;
#endif
                count=0;
                count++;
              }
            }
#ifdef _MAP_DEBUG_
    CkPrintf("Asymmetric SCalcMap created on processor %d\n", CkMyPe());
    dump();
#endif
  }
}

/** helper function */
PeList * rebuildExclusion(int *Pecount, int rsobjs_per_pe)
{
  //needs set, use list for memory, size 0
  PeList *exclusionList=new PeList(1, 1, 0);
  for(int exc=0; exc<config.numPes; exc++)
  {
    if(Pecount[exc]>=rsobjs_per_pe)
    {
      exclusionList->checkAndAdd(exc);
    }
  }
  exclusionList->reset();
  return(exclusionList);
}

RSMapTable::RSMapTable(MapType2  *_frommap, MapType2 *_tomap, PeList *_availprocs,
    int _nstates, int _sizeZ, int _Rstates_per_pe, bool useCuboidMap, MapType2 *gsmap, 
    int nchareG, int numInst):
  nstates(_nstates), sizeZ(_sizeZ), Rstates_per_pe(_Rstates_per_pe)
{
  CkAssert(numInst == 0);
  maptable = _tomap;
  availprocs = _availprocs;
  availprocs->reset();

  double globalTime = 0, sortTime = 0, otherTime1 = 0, otherTime2 = 0, 
         otherTime3 = 0, otherTime4 = 0, otherTime5 = 0;
  double startTime;

  int l, m, pl, pm, srem, rem, i=0, rsobjs_per_pe;
  int *Pecount= new int [config.numPes];

  bzero(Pecount, config.numPes*sizeof(int));

  rsobjs_per_pe = nstates*sizeZ/config.numPesPerInstance;
  l = Rstates_per_pe;		// no of states in one chunk
  pl = nstates / l;
  if(nstates % l == 0)
    srem = 0;
  else
  {
    srem = nstates % pl;
  }
  pm = availprocs->count() / pl;

  if(pm == 0) {
    CkPrintf("Choose a larger Rstates_per_pe than %d such that { no. of processors [%d] / (no. of states [%d] / Rstates_per_pe [%d]) } is > 0 \n",
        l, availprocs->count(), nstates, l);
    CkAssert(availprocs->count() / (nstates/l) > 0);
  }

  m = sizeZ / pm;
  rem = sizeZ % pm;

  int srcpe=0;
  int destpe;

  if(config.simpleTopo)
  {
    int srem = (nstates*sizeZ) % config.numPesPerInstance;

    // exclusion mapping has the sad side effect of increasing the
    // number of exclusions with the state and plane number
    // until you end up increasing the cap too high
    // Topo mapping is less important than even distribution
    // so if you go over the cap, use the master list instead of the
    // state box.

    // this has the effect of creating an imbalance

    PeList *myavail=new PeList(1, 0, *availprocs);
    PeList *exclusionList = new PeList(1, 1, 0);
    for(int state=0; state < nstates; state++)
    {

      int srsobjs_per_pe=rsobjs_per_pe;
      TIMER_SET(startTime = CmiWallTimer();)
      PeList *thisStateBox = subListState(state, nchareG, gsmap);
      int samplePE = *(thisStateBox->list.begin());
      TIMER_SET(globalTime += (CmiWallTimer() - startTime);)
        bool useExclude = true;
      if(useExclude)
      {
        TIMER_SET(otherTime1 -= CmiWallTimer();)
        thisStateBox->deleteList(*exclusionList, 0, 1);
        TIMER_SET(otherTime1 += CmiWallTimer();)
      }
      if(thisStateBox->count() <= 0)
      {
        delete thisStateBox;
        TIMER_SET(otherTime1 -= CmiWallTimer();)
        myavail->addListorVector(0);
        thisStateBox=new PeList(1, 1, *myavail);
        TIMER_SET(otherTime1 += CmiWallTimer();)
        TIMER_SET(sortTime -= CmiWallTimer();)
#if DO_SORT
        thisStateBox->sortSource(samplePE, 1);
#endif
        TIMER_SET(sortTime += CmiWallTimer();)
      }
      if(thisStateBox->count()==0)
      {
        useExclude=false;
        delete thisStateBox;
        TIMER_SET(startTime = CmiWallTimer();)
        thisStateBox = subListState(state, nchareG, gsmap);
        TIMER_SET(globalTime += (CmiWallTimer() - startTime);)
      }

      for(int plane=0; plane < sizeZ; plane++)
      {
        if(thisStateBox->count()<=0)
        {
          if(myavail->size()<=0)
          {
            srsobjs_per_pe++;
            delete thisStateBox;
            delete myavail;
            TIMER_SET(otherTime1 -= CmiWallTimer();)
            myavail= new PeList(1, 0, *availprocs);
            TIMER_SET(otherTime1 += CmiWallTimer();)
            if(exclusionList!=NULL)
              delete exclusionList;
            TIMER_SET(startTime = CmiWallTimer();)
            exclusionList=rebuildExclusion(Pecount, srsobjs_per_pe);
            TIMER_SET(globalTime += (CmiWallTimer() - startTime);)
            TIMER_SET(otherTime5 -= CmiWallTimer();)
            myavail->deleteList(*exclusionList, 0, 0);
            TIMER_SET(otherTime5 += CmiWallTimer();)
            TIMER_SET(otherTime1 -= CmiWallTimer();)
            thisStateBox = new PeList(1, 1, *myavail);
            TIMER_SET(otherTime1 += CmiWallTimer();)
          }
          else
          {
            delete thisStateBox;
            TIMER_SET(otherTime1 -= CmiWallTimer();)
            myavail->addListorVector(0);
            thisStateBox=new PeList(1, 1, *myavail);
            TIMER_SET(otherTime1 += CmiWallTimer();)
          }
          TIMER_SET(sortTime -= CmiWallTimer();)
#if DO_SORT
          thisStateBox->sortSource(samplePE, 1);
#endif
          TIMER_SET(sortTime += CmiWallTimer();)
        }
        if(useExclude && thisStateBox->count()<=0)
        { // surrender
          useExclude=false;
          delete thisStateBox;
          CkAbort("RS map hopeless please examine configuration\n");
        }

        destpe = thisStateBox->findNext();

#ifdef USE_INT_MAP
        maptable->set(state, plane, destpe);
#else						
        maptable->put(intdual(state, plane))= destpe;
#endif
        Pecount[destpe]++;
        if(Pecount[destpe]>=srsobjs_per_pe)
        {
          TIMER_SET(otherTime2 -= CmiWallTimer();)
          exclusionList->checkAndAdd(destpe);
          TIMER_SET(otherTime2 += CmiWallTimer();)
          TIMER_SET(otherTime3 -= CmiWallTimer();)
          thisStateBox->deleteCurrent();
          TIMER_SET(otherTime3 += CmiWallTimer();)
          TIMER_SET(otherTime4 -= CmiWallTimer();)
          PeList one(0, 0, 0);
          one.addOne(destpe);
          myavail->deleteList(one, 1, 0);
          TIMER_SET(otherTime4 += CmiWallTimer();)
        }
      }
      delete thisStateBox;
    }
    if(exclusionList!=NULL)
      delete exclusionList;
    if(myavail != NULL)
      delete myavail;
  } else {

    // this remainder scheme is odd, creates imbalance.
    // doesn't use all processors
    destpe=availprocs->findNext();
    int orig_l=l;
    for(int ychunk=0; ychunk<sizeZ; ychunk=ychunk+m)
    {
      if(ychunk==(pm-rem)*m)
        m=m+1;
      l=orig_l;
      for(int xchunk=0; xchunk<nstates; xchunk=xchunk+l)
      {
        if(xchunk==(pl-srem)*l)
          l=l+1;
        for(int state=xchunk; state<xchunk+l && state<nstates; state++)
        {
          for(int plane=ychunk; plane<ychunk+m && plane<sizeZ; plane++)
          {
#ifdef USE_INT_MAP
            maptable->set(state, plane, destpe);
#else						
            maptable->put(intdual(state, plane))= destpe;
#endif
          }
        }
        destpe=availprocs->findNext();
      }
    }
  }
  delete [] Pecount;
#ifdef _MAP_DEBUG_
  CkPrintf("RSMap created on processor %d\n", CkMyPe());
  dump();
#endif
  CkPrintf("Time in state creation %.2lf s, sorting %.2lf, othertime1 %.2lf, othertime2 %.2lf "
      "othertime3 %.2lf othertime4 %.2lf othertime5 %.2lf\n", globalTime, sortTime, otherTime1,
      otherTime2, otherTime3, otherTime4, otherTime5);
}

RPPMapTable::RPPMapTable(MapType2  *_map, 
    PeList *_availprocs, PeList *exclusion,
    int _nstates, int _sizeZNL, int _Rstates_per_pe,
    int boxSize, bool usePPmap, int nchareG,
    MapType2 *pp_map) :
  nstates(_nstates), sizeZNL(_sizeZNL),
  Rstates_per_pe(_Rstates_per_pe) {

  int states_per_pe=Rstates_per_pe;
  int totalChares=nstates*sizeZNL;
  maptable=_map;
  availprocs=_availprocs;
  bool useExclusion=true;
  PeList *RPPlist=availprocs;
  int chares_per_pe=totalChares/config.numPesPerInstance;
  CkPrintf("CharesPerPe %d states_per_pe %d\n",chares_per_pe, states_per_pe);
  if(exclusion==NULL || exclusion->count()==0 || config.numPesPerInstance <=exclusion->count() )
    useExclusion=false;
  int afterExclusion=availprocs->count();
  if(useExclusion)
    afterExclusion=availprocs->count() - exclusion->count();
  if(useExclusion && afterExclusion > chares_per_pe*sizeZNL)
  { // we can fit the exclusion without blinking
    CkPrintf("RPP using density exclusion to avoid %d processors\n",exclusion->count());
    RPPlist->deleteList(*exclusion, 0, 0);
  }
  else
  {// so an rstates_per_pe chosen for realstate might be too big
    if(useExclusion && afterExclusion > (sizeZNL+nstates)*2)
    { // set states_per_pe to fit the exclusion

      states_per_pe=(int) ((float) (totalChares/afterExclusion))*0.75;
      CkPrintf("RPP adjusting states per pe from %d to %d to use density exclusion to stay within %d processors\n",Rstates_per_pe, states_per_pe, exclusion->count());
      RPPlist->deleteList(*exclusion, 0, 0);
    }
    else
    {
      useExclusion=false;
      CkPrintf("RPP with %d chares ignoring density exclusion which left %d out of %d processors\n",totalChares, exclusion->count(), availprocs->count());
    }
  }

  chares_per_pe=totalChares/RPPlist->count();
  // CkPrintf("nstates %d sizeZNL %d Pes %d\n", nstates, sizeZNL, RPPlist->count());	
  int srcpe=0;
  if(config.simpleTopo)
  {
    /*  
     * RPP(s,*) <-> PP(s,*)  
     */
    bool neednewexc=false;
    // make all the maps
    // subtract the exclusion from each 
    // keep the smallest count as the max charesperpe
    //      PeList **maps= new PeList* [nstates];
    // this code is too memory hoggy
    int maxcharesperpe=states_per_pe;
    int *usedPes= new int[config.numPes];
    bool *useExclude= new bool[nstates];
    bzero(usedPes, config.numPes * sizeof(int));
    for(int state=0; state < nstates ; state++)
    {
      // have variable number of exclusions per list
      // need to figure out what the max should be
      // for the border cases
      //maps[state]= subListState( state, nchareG, pp_map);
      PeList *state_map=subListState( state, nchareG, pp_map);
      int cur_count = state_map->size();
      CkAssert(state_map->size()>0);
      if(useExclusion)
      {
        state_map->deleteList(*exclusion, 1, 0);
        if(state_map->size()==0 || sizeZNL/state_map->size() > states_per_pe)
        { //not enough for exclusion
          useExclude[state]=false;
        }
        else
        {
          cur_count = state_map->size();
          useExclude[state]=true;
        }
      }
      else
      {
        useExclude[state]=false;
      }
      int thischaresperpe=sizeZNL/cur_count + 1;

      maxcharesperpe=(thischaresperpe>maxcharesperpe) ? thischaresperpe : maxcharesperpe;

      // if(states_per_pe>maxcharesperpe)
      //	    maxcharesperpe=states_per_pe;
      delete state_map;
    }
    int totcharesperpe=sizeZNL*nstates/config.numPesPerInstance+1;
    maxcharesperpe=(maxcharesperpe>totcharesperpe) ? maxcharesperpe : totcharesperpe;
    //      maxcharesperpe++;
    PeList *usedbyRPP = new PeList(1, 1, 0);
    CkPrintf("RPP maxcharesperpe is %d\n",maxcharesperpe);
    int origmaxcharesperpe=maxcharesperpe;
    PeList *excludedBigmap=new PeList(1, 0, 0);
    for(int state=0; state < nstates ; state++)
    {
      PeList *state_map=subListState( state, nchareG, pp_map);
      if(useExclude[state])
      {
        state_map->deleteList(*exclusion, 0, 1);
      }
      while(state_map->count()<=0 && maxcharesperpe<=sizeZNL*nstates)
      {
        // ditch topo scheme for overflow
        // use the RPPlist
        //CkPrintf("State %d  Ran out of procs in RPP centroid using full RPPlist\n",state);
        delete state_map;
        if(!neednewexc && excludedBigmap->size() != 0)
        {
          excludedBigmap->addListorVector(0);
          state_map= new  PeList(1, 1, *excludedBigmap);
        }
        else
        {
          state_map = new PeList(1, 1, *RPPlist);
          if(usedbyRPP!=NULL){
            state_map->deleteList(*usedbyRPP, 0, 1);
            if(excludedBigmap!=NULL)
              delete excludedBigmap;
            excludedBigmap=new PeList(1, 0, *state_map);
          }
        }
        if(state_map->count()<=0)
        {
          maxcharesperpe++;
          //CkPrintf("plane %d  Ran out of procs in RPP centroid using full RPPlist and bumping maxcharesperpe to %d\n",state, maxcharesperpe);
          if(usedbyRPP!=NULL)
            delete usedbyRPP;
          usedbyRPP=rebuildExclusion(usedPes, maxcharesperpe);
          delete state_map;
          state_map = new PeList(1, 1, *RPPlist);
          if(usedbyRPP!=NULL && usedbyRPP->count()>0){
            state_map->deleteList(*usedbyRPP, 0, 1);
          }
          if(excludedBigmap!=NULL)
            delete excludedBigmap;
          excludedBigmap= new PeList(1, 0, *state_map);
          neednewexc=false;
        }
        // ADD SORT here
        //state_map->sortSource(pp_map->get(state,0));
      }
      for(int plane=0; plane < sizeZNL; plane++)
      {
        while(state_map->count()<=0 && maxcharesperpe<=sizeZNL*nstates)
        {
          //CkPrintf("State %d  Ran out of procs in RPP centroid using full RPPlist\n",state);
          delete state_map;
          if(!neednewexc && excludedBigmap!=NULL && excludedBigmap->size()>0)
          {
            excludedBigmap->addListorVector(0);
            state_map= new PeList(1, 1, *excludedBigmap);
          }
          else
          {
            state_map = new PeList(1, 1, *RPPlist);
            if(usedbyRPP->size() != 0){
              state_map->deleteList(*usedbyRPP, 0, 1);
              if(excludedBigmap!=NULL)
                delete excludedBigmap;
              if(state_map->count()==0)
              { // man we're totally dry here.
                // this should be handled in the next block
                excludedBigmap=new PeList(1, 0, 0);
              }
              else
              {
                excludedBigmap= new PeList(1, 0, *state_map);
                neednewexc=false;
              }
            }
            if(state_map->count()<=0)
            {
              maxcharesperpe++;
              //CkPrintf("plane %d  Ran out of procs in RPP centroid using full RPPlist and bumping maxcharesperpe to %d\n",state, maxcharesperpe);

              if(usedbyRPP!=NULL)
                delete usedbyRPP;
              usedbyRPP=rebuildExclusion(usedPes, maxcharesperpe);
              delete state_map;
              state_map = new PeList(1, 1, *RPPlist);
              if(usedbyRPP!=NULL)
              {
                state_map->deleteList(*usedbyRPP, 0, 1);
              }
              if(state_map->count()==0)
              { // man we're totally dry here.
                // reboot
                //CkPrintf("plane %d  Ran out of procs in RPP centroid using full RPPlist after bumping maxcharesperpe to %d, clearing used list, resetting maxcharesperpe to  %d\n",state, maxcharesperpe, origmaxcharesperpe);

                maxcharesperpe=origmaxcharesperpe;
                bzero(usedPes, config.numPes * sizeof(int));
                if(usedbyRPP!=NULL)
                  delete usedbyRPP;
                usedbyRPP=rebuildExclusion(usedPes, maxcharesperpe);
                delete state_map;
                state_map = new PeList(1, 1, *RPPlist);
                if(usedbyRPP!=NULL)
                {
                  state_map->deleteList(*usedbyRPP, 0, 1);
                }
              }

              if(excludedBigmap!=NULL)
                delete excludedBigmap;
              excludedBigmap= new PeList(1, 0, *state_map);
              neednewexc=false;
            }
          }
#if DO_SORT
          state_map->sortSource(pp_map->get(state,0), 1);
#endif
        }
        int destpe=state_map->findNext(); 
#ifdef USE_INT_MAP
        maptable->set(state, plane, destpe);
#else
        maptable->put(intdual(state, plane))=destpe;
#endif
        usedPes[destpe]++;
        if(usedPes[destpe]>maxcharesperpe)
        {
          neednewexc=true;
          usedbyRPP->checkAndAdd(destpe);
          exclusion->checkAndAdd(destpe);
          PeList one(0,0,0);
          one.addOne(destpe);
          if(excludedBigmap->size() != 0) {
            neednewexc=false;
            excludedBigmap->deleteList(one, 1, 0);
          }
          state_map->deleteCurrent();
        }
      }
      delete state_map;
    }
    delete usedbyRPP;
    delete [] usedPes;
    delete [] useExclude;
    if(excludedBigmap!=NULL)
      delete excludedBigmap;
  }
  else
  {
    int destpe=RPPlist->findNext();
    // place by planes up to chares_per_pe to give us a balanced layout
    int charesOnThisPe=0;
    int rem = totalChares % chares_per_pe;
    int starting_cpp=chares_per_pe;
    //      if(rem) --chares_per_pe;
    CkPrintf("chares_per_pe %d rem %d\n",chares_per_pe, rem);
    for(int state=0;  state < nstates; state++)
    {
      for(int plane=0; plane < sizeZNL; plane++)
      {
#ifdef USE_INT_MAP
        maptable->set(state, plane, destpe);
#else
        maptable->put(intdual(state, plane))=destpe;
#endif
        if(++charesOnThisPe==chares_per_pe)
        {
          destpe=RPPlist->findNext();
          charesOnThisPe=0;
        }
      }
      if(rem && RPPlist->count()<2*nstates)
      {
        chares_per_pe=starting_cpp+1;
      }

    }
  }
#ifdef _MAP_DEBUG_
  CkPrintf("RPPMap created on processor %d\n", CkMyPe());
  dump();
#endif
}

OrthoMapTable::OrthoMapTable(MapType2 *_map, PeList *_availprocs, int _nstates, int _orthograinsize, MapType4 *scalcmap, int nplanes, int numChunks, int sGrainSize, PeList *exclusionList): nstates(_nstates), orthoGrainSize(_orthograinsize) 
{
  maptable = _map;
  availprocs = _availprocs;
  int oobjs_per_pe;
  int srcpe = 0, destpe;
  bool useSublist=true;
  bool useExclude=true;
  int maxorthoindex=(nstates/orthoGrainSize-1);
  int northo=(maxorthoindex+1)*(maxorthoindex+1);
  oobjs_per_pe = northo/(config.numPes);
  int *Pecount= new int [config.numPes];
  bzero(Pecount, config.numPesPerInstance*sizeof(int)); 
  int s1 = 0, s2 = 0;
  int maxpcstateindex=(nstates/sGrainSize-1)*sGrainSize;
  int maxorthostateindex=(nstates/orthoGrainSize-1)*orthoGrainSize;
  for(int state1 = 0; state1 <= maxorthostateindex; state1 += orthoGrainSize)
    for(int state2 = 0; state2 <= maxorthostateindex; state2 += orthoGrainSize)
    {
      s1 = (state1/sGrainSize);
      s1 = s1 * sGrainSize;
      s1 = (s1>maxpcstateindex) ? maxpcstateindex :s1;
      s2 = (state2/sGrainSize);
      s2 = s2 * sGrainSize;
      s2 = (s2>maxpcstateindex) ? maxpcstateindex :s2;
      PeList *thisStateBox=NULL;
      if(useSublist)
      {
        thisStateBox = subListState2(s1, s2, nplanes, numChunks, scalcmap);
        useExclude = true;
        thisStateBox->deleteList(*exclusionList, 0, 1);
      }
      else
      {
        thisStateBox = availprocs;
      }
      if(thisStateBox->size() == 0)
      {  // the sublist scheme failed 
        if(thisStateBox!=availprocs)
          delete thisStateBox;
        thisStateBox = availprocs;
      }

      if(thisStateBox->size() == 0)
      {
        //CkPrintf("Ortho %d %d ignoring exclusion\n", state1, state2);
        if(thisStateBox!=availprocs)
          delete thisStateBox;
        thisStateBox = subListState2(s1, s2, nplanes, numChunks, scalcmap);
        useExclude = false;
      }
      thisStateBox->addListorVector(1);
      destpe=minDistCentroid(thisStateBox, nplanes, s1, s2, numChunks, scalcmap);
      int os1 = (state1/orthoGrainSize);
      os1 = (os1>maxorthoindex) ? maxorthoindex :os1;
      int os2=(state2/orthoGrainSize);
      os2 = (os2>maxorthoindex) ? maxorthoindex :os2;
#ifdef USE_INT_MAP
      maptable->set(os1, os2, destpe);
#else
      maptable->put(intdual(os1, os2))=destpe;
#endif
      Pecount[destpe]++;	
      if(Pecount[destpe]>=oobjs_per_pe)
      {
        exclusionList->checkAndAdd(destpe);
        if(useExclude) {
          PeList one(0, 0, 0);
          one.addOne(destpe);
          availprocs->deleteList(one, 1, 0);
        }
      }
      if(thisStateBox!=availprocs)
        delete thisStateBox;
    }
    availprocs->addListorVector(1);
    delete [] Pecount;
}

OrthoHelperMapTable::OrthoHelperMapTable(MapType2 *_map, int _nstates, int _orthograinsize, MapType2 *omap, PeList *_avail, PeList *exclude): nstates(_nstates), orthoGrainSize(_orthograinsize) 
{
  maptable = _map;
  int destpe = 0;
  // map orthohelper near but not on ortho by removing all ortho procs
  availprocs= new PeList(1, 1, *_avail);
  availprocs->deleteList(*exclude, 0, 1);
  //  CkPrintf("exlude for helpers is \n");
  //  exclude->dump();
  bool useExclude=true;
  if(availprocs->size() < ((nstates/orthoGrainSize)*(nstates/orthoGrainSize)))
  {
    CkPrintf("There aren't enough processors to effectively parallelize orthoHelpers, you should disable orthoHelpers!!!\n");
    delete availprocs;
    availprocs= new PeList(1, 1, *_avail);
    useExclude=false;
  }
  int maxorthoindex=(nstates/orthoGrainSize-1);
  for(int state1 = 0; state1 <= maxorthoindex; state1++)
    for(int state2 = 0; state2 <= maxorthoindex; state2++)
    {

#ifdef USE_INT_MAP
      destpe=availprocs->minDist(omap->get(state1, state2));
      maptable->set(state1, state2, destpe);
#else
      destpe=availprocs->minDist(omap->get(intdual(state1, state2)));
      maptable->put(intdual(state1, state2))=destpe;
#endif
      if(useExclude)
      {
        availprocs->deleteCurrent();
      }
    }
  delete availprocs;
}


RhoRSMapTable::RhoRSMapTable(MapType2  *_map, PeList *_availprocs,
    int _nchareRhoR_x, int _nchareRhoR_y, int max_states, bool useCentroid,
    MapType2 *rsmap, PeList *exclusionList): nchareRhoR_x(_nchareRhoR_x),
    nchareRhoR_y(_nchareRhoR_y) {
  maptable = _map;
  availprocs = _availprocs;
  int srcpe=0;
  int numChares = nchareRhoR_x * nchareRhoR_y;
  int pesused=0;

  int destpe;
  if(config.simpleTopo)
  {
    for(int pencil_x=0; pencil_x<nchareRhoR_x; pencil_x++)
    {
      int start_plane, end_plane;
      FFT_START(start_plane, pencil_x, nchareRhoR_x, simReadOnly.sizeZ);
      FFT_END(end_plane, pencil_x, nchareRhoR_x, simReadOnly.sizeZ);
      PeList *initPlaneBox = subListPlanes(start_plane, end_plane, max_states, rsmap);
      bool useExclude=true;
      initPlaneBox->deleteList(*exclusionList, 0, 1);
      if(initPlaneBox->count() < nchareRhoR_y) { CkAbort("Ran out of PEs to map pencils onto");}
      sortByCentroid(initPlaneBox, pencil_x, max_states, rsmap);
      PeList *thisPlaneBox= initPlaneBox->distributeAcrossPelist(nchareRhoR_y);
      destpe=thisPlaneBox->findNext();
      for(int pencil_y=0 ; pencil_y<nchareRhoR_y ; pencil_y++)
      {
        pesused++;
        maptable->set(pencil_x, pencil_y, destpe);
        exclusionList->checkAndAdd(destpe);
        destpe=thisPlaneBox->findNext();
      }
      delete thisPlaneBox;
      delete initPlaneBox;
    }
  }
  else
  {
    int nprocs=0, objs=0;
    PeList *initPlaneBox = new PeList(*availprocs);
    initPlaneBox->deleteList(*exclusionList, 0, 1);
    if(initPlaneBox->count()==0) { CkAbort("Ran out of PEs to map pencils onto");}
    PeList *thisPlaneBox = initPlaneBox->distributeAcrossPelist(numChares);
    destpe = thisPlaneBox->findNext();
    for(int pencil_x=0; pencil_x<nchareRhoR_x; pencil_x++)
    {
      pesused++;
      for(int pencil_y=0; pencil_y<nchareRhoR_y; pencil_y++)
      {
        maptable->set(pencil_x, pencil_y, destpe);
	destpe=thisPlaneBox->findNext();
        exclusionList->checkAndAdd(destpe);
        destpe = thisPlaneBox->findNext();
      }
    }
    delete thisPlaneBox;
    delete initPlaneBox;
  }
  CkPrintf("Built RhoRS Map [%d, %d] on %d processors\n",nchareRhoR_x,nchareRhoR_y, pesused );
#ifdef _MAP_DEBUG_
  CkPrintf("RhoRSMap created on processor %d\n", CkMyPe());
  dump();
#endif

}
/**
 * RhoG and RhoR are mutually all to all. They should be mapped such
 * that they are relatively near but exclusive.  Meaning they should
 * share no processors if there are enough to go around.  Given that
 * RhoR is already placed.  We can take the centroid of that map. Then
 * sort the available list based on distance to that centroid.  Giving
 * us the cloud of available processors around the RhoRS processors as
 * our preferred placement.
 */

RhoGSMapTable::RhoGSMapTable(MapType1  *_map, PeList *_availprocs, int
_nchareRhoG,  bool useCentroid, MapType2 *rhorsmap, PeList *exclusionList):
nchareRhoG(_nchareRhoG)
{
  maptable=_map;
  availprocs=_availprocs;
  int rgsobjs_per_pe, rem;

  PeList *avail= new PeList(*availprocs);
  avail->deleteList(*exclusionList, 0, 0);

  PeList *rhoglist;
  if(avail->count() > nchareRhoG)
  {
    if(config.simpleTopoCentroid) {
      int dims[10];
      rhorsmap->getCentroid(config.torusMap, dims);
      avail->sortSource(dims, 0);
    }
    rhoglist = avail->distributeAcrossPelist(nchareRhoG);
  } else {
    if(config.simpleTopoCentroid) {
      // get centroid of rsmap  use it to sort the avail list
      int dims[10];
      rhorsmap->getCentroid(config.torusMap, dims);
      availprocs->sortSource(dims, 1);
    }
    rhoglist = availprocs->distributeAcrossPelist(nchareRhoG);
  }
  delete avail;

  int destpe=rhoglist->findNext();
  for(int chunk=0; chunk<nchareRhoG; chunk++)
  {
    maptable->set(chunk, destpe);
    exclusionList->checkAndAdd(destpe);
    destpe = rhoglist->findNext();
  }
  delete rhoglist;
#ifdef _MAP_DEBUG_
  CkPrintf("[%d] RhoGSMap created on pes %d\n", CkMyPe(), nchareRhoG);
  dump();
#endif
}

RhoRHartMapTable::RhoRHartMapTable(MapType3  *_map, PeList *_availprocs,
    int _nchareRhoRHart_x, int _nchareRhoRHart_y, int nchareHartAtmT,PeList *exclude ):
    nchareRhoRHart_x(_nchareRhoRHart_x), nchareRhoRHart_y(_nchareRhoRHart_y) {
  maptable=_map;
  availprocs=_availprocs;
  int rrsobjs_per_pe, rem;
  int numChares=nchareRhoRHart_x*nchareRhoRHart_y*nchareHartAtmT;
  PeList *avail= new PeList(1, 0, *availprocs);
  avail->deleteList(*exclude, 0, 0);
  if(avail->count()>numChares)
  {
    // try an exclusion
    CkPrintf("RhoRHart excluding %d from avail %d\n",exclude->count(), availprocs->count());
  }
  else
  {
    CkPrintf("Cannot use exclusion in rhoRhart\n");
    delete avail;
    avail = new PeList(1, 0, *availprocs);
  }

  if(avail->count()==1)
  {
    rrsobjs_per_pe = numChares;
    rem=0;
  }
  else
  {
    rrsobjs_per_pe= numChares/(avail->count());
    rem = numChares % (avail->count());
    if(numChares<avail->count())
    {
      rrsobjs_per_pe=1;
      rem=0;
    }
    if(rem!=0)
      rrsobjs_per_pe += 1;
  }

  PeList *finalList = avail->distributeAcrossPelist(numChares);
  delete avail;
  int destpe = finalList->findNext();
  if(CkMyPe()==0)
    CkPrintf("nchareRhoRHart_x %d nchareRhoRHart_y %d rrsobjs_per_pe %d rem %d\n", 
    nchareRhoRHart_x, nchareRhoRHart_y, rrsobjs_per_pe, rem);

  int nprocs=0, objs=0;

  for(int atmtype=0; atmtype< nchareHartAtmT;atmtype++)
  {
    for(int chunk=0; chunk<nchareRhoRHart_x; chunk++)
    {
      for(int subplane=0; subplane<nchareRhoRHart_y; subplane++)
      {
        if(rem!=0)
          if(nprocs==rem)
            rrsobjs_per_pe -= 1;
        maptable->set(chunk, subplane, atmtype,destpe);
        objs++;
        exclude->checkAndAdd(destpe);
        if(objs>=rrsobjs_per_pe)
        {
          destpe=finalList->findNext();
          objs=0;
          nprocs++;
        }
      }
    }
  }
  delete finalList;
#ifdef _MAP_DEBUG_
  CkPrintf("RhoRHartMap created on processor %d\n", CkMyPe());
  dump();
#endif
}

RhoGHartMapTable::RhoGHartMapTable(MapType2  *_map, PeList *_availprocs,
    int _nchareRhoGHart, int nchareHartAtmT, int useCentroid, MapType2 *rhorsmap,
    MapType3 *rhartmap, PeList *exclude): nchareRhoGHart(_nchareRhoGHart)
{
  maptable=_map;
  PeList *availprocs=_availprocs;
  PeList *avail= new PeList(*availprocs);
  avail->deleteList(*exclude, 0, 0);
  int numchares=nchareRhoGHart*nchareHartAtmT;
  bool excluded=true;
  PeList *rhohartlist=NULL;
  if(avail->size()>numchares)
  {
    // try an exclusion
    if(config.simpleTopoCentroid)
    {
      // get centroid of rhartmap  use it to sort the avail list
      if(rhorsmap != NULL) {
        int dims[10];
        rhorsmap->getCentroid(config.torusMap, dims);
        avail->sortSource(dims, 0);
      } else {
        int dims[10];
        rhartmap->getCentroid(config.torusMap, dims);
        avail->sortSource(dims, 0);
      }
    }

    CkPrintf("RhoGHart excluding %d from avail %d\n",exclude->count(), availprocs->count());
    rhohartlist = avail->distributeAcrossPelist(numchares);
    printf("RhoGHart map: avail %d, used %d,%d\n", avail->size(), 
      numchares, rhohartlist->size());
  }
  else
  {
    if(config.simpleTopoCentroid)
    {
      // get centroid of rhartmap  use it to sort the avail list
      if(rhorsmap != NULL) {
        int dims[10];
        rhorsmap->getCentroid(config.torusMap, dims);
        availprocs->sortSource(dims, 1);
      } else {
        int dims[10];
        rhartmap->getCentroid(config.torusMap, dims);
        printf("Centroid is %d %d %d %d %d\n", dims[0], dims[1],
        dims[2], dims[3], dims[4]);
        availprocs->sortSource(dims, 1);
      }
    }

    excluded=false;
    CkPrintf("Cannot use exclusion in rhoghart\n");
    rhohartlist=availprocs->distributeAcrossPelist(numchares);
    printf("RhoGHart map: avail %d, used %d,%d\n", availprocs->size(), 
      numchares, rhohartlist->size());
  }
  delete avail;
  int destpe=rhohartlist->findNext();
  for(int atmtype=0; atmtype< nchareHartAtmT;atmtype++)
  {
    for(int chunk=0; chunk<nchareRhoGHart; chunk++)
    {
      maptable->set(chunk, atmtype, destpe);
      exclude->checkAndAdd(destpe);
      destpe=rhohartlist->findNext();
    }
  }
  delete rhohartlist;
#ifdef _MAP_DEBUG_
  CkPrintf("RhoGHartMap created on processor %d\n", CkMyPe());
  dump();
#endif

}

/* We receive a copy of the RhoRS and then act independently with the
   bulk of our communication with VdWGS. Exclusion mapping will try to
   put us on different processors from RhoRS. Our connection to RS is
   more tenuous, so centroid mapping to RS is not as compelling. */
VdWRSMapTable::VdWRSMapTable(MapType3  *_map, PeList *_availprocs, int _nchareRhoR, 
    int _rhoRsubplanes, int _nchareVdW, int max_states, PeList *exclude): 
  nchareRhoR(_nchareRhoR), rhoRsubplanes(_rhoRsubplanes) ,nchareVdW(_nchareVdW) {
    maptable=_map;
    availprocs=_availprocs;
    int rvdwobjs_per_pe, rem;
    int srcpe=0;
    int numChares=nchareRhoR*rhoRsubplanes* nchareVdW;
    int pesused=0;
    if(availprocs->count()==1)
    {
      rvdwobjs_per_pe= numChares;
      rem=0;
    }
    else
    {
      rvdwobjs_per_pe= numChares/(availprocs->count());
      rem = numChares % (availprocs->count());
      if(numChares<availprocs->count())
      {
        rem=0;
        rvdwobjs_per_pe=1;
      }
      if(rem!=0)
        rvdwobjs_per_pe += 1;
    }
    int destpe;
    int nprocs=0, objs=0;
    destpe=availprocs->findNext();
    for(int chunk=0; chunk<nchareRhoR; chunk++)
    {
      for(int subplane=0; subplane<rhoRsubplanes; subplane++)
      {
        for(int vdw=0; vdw<nchareVdW; vdw++)
        {
          if(rem!=0)
            if(nprocs==rem)
              rvdwobjs_per_pe -= 1;
#ifdef USE_INT_MAP
          maptable->set(chunk, subplane, vdw, destpe);
#else
          maptable->put(inttriple(chunk, subplane,vdw))=destpe;
#endif
          objs++;
          if(objs==rvdwobjs_per_pe)
          {
            destpe=availprocs->findNext();
            objs=0;
            nprocs++;
          }
        }
      }
    }
    CkPrintf("Built VdWRS Map [%d, %d, %d]\n",nchareRhoR,rhoRsubplanes, nchareVdW); 
#ifdef _MAP_DEBUG_
    CkPrintf("VdWRSMap created on processor %d\n", CkMyPe());
    dump();
#endif
  }

/**
 * VdWG and VdWR are mutually all to all. They should be mapped such
 * that they are relatively near but exclusive.  Meaning they should
 * share no processors if there are enough to go around.

 * Given that VdWR is already placed we could map relative to it. But
 * the cross communication within each of R and G could dominate so
 * we'll start with a simple map and optimize it when we see how this
 * new beast works.

 */

VdWGSMapTable::VdWGSMapTable(MapType2  *_map, PeList *_availprocs, int _nchareRhoG,  int _nchareVdW, PeList *exclude): nchareRhoG(_nchareRhoG), nchareVdW(_nchareVdW)
{
  maptable=_map;
  availprocs=_availprocs;
  int rgsobjs_per_pe, rem;

  availprocs->reset();
  PeList *avail= new PeList(1, 0, *availprocs);
  avail->deleteList(*exclude, 1, 0);
  int numVdWObjs=nchareRhoG*nchareVdW;
  if(avail->size()>numVdWObjs)
  {
    // try an exclusion
    CkPrintf("RhoG excluding %d from avail %d\n",exclude->count(), availprocs->count());
    availprocs->deleteList(*exclude, 1, 0);
  }
  else
  {
    CkPrintf("cannot use exclusion in rhog\n");
    availprocs->reset();
  }
  delete avail;

  if(availprocs->count()==1)
  {
    rgsobjs_per_pe= numVdWObjs;
    rem=0;
  }
  else
  {
    rgsobjs_per_pe= numVdWObjs/(availprocs->count());
    rem = numVdWObjs % (availprocs->count());
    if(numVdWObjs<availprocs->count())
    {
      rem=0;
      rgsobjs_per_pe=1;
    }
    if(rem!=0)
      rgsobjs_per_pe += 1;
  }

  int destpe=availprocs->findNext();

  int thischunk=0;
  for(int i=0; i<nchareRhoG; i++)
  {
    for(int j=0; j<nchareVdW; j++)
    {

#ifdef USE_INT_MAP
      maptable->set(i, j, destpe);
#else
      maptable->put(intdual(i, j))=destpe;
#endif
      thischunk++;
      if(thischunk>rgsobjs_per_pe)
      {
        exclude->checkAndAdd(destpe);
        destpe=availprocs->findNext();
      }
    }
  }
#ifdef _MAP_DEBUG_
  CkPrintf("VdWGSMap created on processor %d\n", CkMyPe());
  dump();
#endif
}

PeList *subListPlane(int plane, int nstates, MapType2 *smap)
{
  //no set, use list, zero elems
  PeList *thisPlane = new PeList(1, 1, 0);
  for(int state=0; state<nstates; state++)
  {
    thisPlane->checkAndAdd(smap->get(state,plane));
  }
  thisPlane->reset();
  return(thisPlane);
}

PeList *subListState(int state, int nplanes, MapType2 *smap)
{
  //      CkPrintf("in sublist state %d\n",state);
  PeList *thisState= new PeList(1, 1, 0);
  for(int plane=0; plane < nplanes; plane++)
  {
    int pe=smap->get(state, plane);
    thisState->checkAndAdd(pe);
  }
  thisState->reset();
  return(thisState);
}

PeList *subListState2(int state1, int state2, int nplanes, int numChunks, MapType4 *smap)
{
  PeList *thisState = new PeList(1, 1, 0);
  for(int plane=0; plane<nplanes; plane++)
    for(int chunk=0; chunk<numChunks; chunk++)
    {
      int pe = smap->get(plane, state1, state2, chunk);
      thisState->checkAndAdd(pe);
    }
  thisState->reset();
  return(thisState);
}

PeList *subListPlanes(int start_plane, int end_plane, int nstates, MapType2 *smap)
{

  PeList *thisBox = new PeList(1, 1, 0);
  for(int plane=start_plane; plane < end_plane; plane++)
  {
    for(int state=0; state<nstates; state++)
    {
      thisBox->checkAndAdd(smap->get(state,plane));
    }
  }
  thisBox->reset();
  return(thisBox);
}

void RhoRSMapTable::sortByCentroid(PeList *avail, int plane, int nstates, MapType2 *rsmap)
{
  int points=0, bestPe;
  if(config.simpleTopo) {
    int ndims;
    int sum[10], dim[10];
    TopoManager_getDimCount(&ndims);
    for(int i = 0; i < ndims; i++) {
      sum[i] = 0;
    }
    for(int state=0;state<nstates;state++)
    {
      TopoManager_getPeCoordinates(rsmap->get(state, plane), dim);
      for(int i = 0; i < ndims; i++) {
        sum[i] += dim[i];
      }
      points++;
    }
    for(int i = 0; i < ndims; i++) {
      sum[i] /= points;
    }
    sum[ndims] = 0;
    avail->sortSource(sum, 1);
    avail->reset();
  } else {
    int sumPe=0;
    for(int state=0;state<nstates;state++)
    {
      sumPe+=rsmap->get(state,plane);
      points++;
    }
    bestPe=sumPe/points;
    avail->sortSource(bestPe, 1);
    avail->reset();
  }
}

void SCalcMapTable::sortByCentroid(PeList *avail, int plane, int stateX, int stateY, int grainsize, MapType2 *gsmap)
{
  int points=0, bestPe;

  if(config.simpleTopo) {
    int ndims;
    int sum[10], dim[10];
    TopoManager_getDimCount(&ndims);
    for(int i = 0; i < ndims; i++) {
      sum[i] = 0;
    }
    for(int state=stateX;(state<stateX+grainsize)&&(state<config.nstates);state++)
    {
      TopoManager_getPeCoordinates(gsmap->get(state, plane), dim);
      for(int i = 0; i < ndims; i++) {
        sum[i] += dim[i];
      }
      points++;
    }
    for(int state=stateY;(state<stateY+grainsize)&&(state<config.nstates);state++)
    {
      TopoManager_getPeCoordinates(gsmap->get(state, plane), dim);
      for(int i = 0; i < ndims; i++) {
        sum[i] += dim[i];
      }
      points++;
    }
    for(int i = 0; i < ndims; i++) {
      sum[i] /= points;
    }
    sum[ndims] = 0;
    avail->sortSource(sum, 1);
    avail->reset();
  }
  else {
    int sumPe=0;
    for(int state=stateX;(state<stateX+grainsize)&&(state<config.nstates);state++)
    {
      sumPe+=gsmap->get(state,plane);
      points++;
    }
    for(int state=stateY;(state<stateY+grainsize)&&(state<config.nstates);state++)
    {
      sumPe+=gsmap->get(state,plane);
      points++;
    }
    bestPe=sumPe/points;
    avail->sortSource(bestPe, 1);
    avail->reset();
  }
}

void OrthoMapTable::sortByCentroid(PeList *avail, int nplanes, int state1, int state2, int numChunks, MapType4 *smap)
{
  int points=0, bestPe;
  if(config.simpleTopo) {
    int ndims;
    int sum[10], dim[10];
    TopoManager_getDimCount(&ndims);
    for(int i = 0; i < ndims; i++) {
      sum[i] = 0;
    }
    for(int plane=0; plane<nplanes; plane++)
      for(int chunk=0; chunk<numChunks; chunk++)
      {  
        TopoManager_getPeCoordinates(smap->get(plane, state1, state2, chunk), dim);
        for(int i = 0; i < ndims; i++) {
          sum[i] += dim[i];
        }
        points++;
      }
    for(int i = 0; i < ndims; i++) {
      sum[i] /= points;
    }
    sum[ndims] = 0;
    avail->sortSource(sum, 1);
    avail->reset();
  } else {
    int sumPe = 0;
    for(int plane=0; plane<nplanes; plane++)
      for(int chunk=0; chunk<numChunks; chunk++)
      {
        sumPe += smap->get(plane, state1, state2, chunk);
        points++;
      }
    bestPe = sumPe/points;
    avail->sortSource(bestPe, 1);
    avail->reset();
  }
}

int OrthoMapTable::minDistCentroid(PeList *avail, int nplanes, int state1, int state2, int numChunks, MapType4 *smap)
{
  int points=0, bestPe;
  if(config.simpleTopo) {
    int ndims;
    int sum[10], dim[10];
    TopoManager_getDimCount(&ndims);
    for(int i = 0; i < ndims; i++) {
      sum[i] = 0;
    }
    for(int plane=0; plane<nplanes; plane++)
      for(int chunk=0; chunk<numChunks; chunk++)
      {    
        TopoManager_getPeCoordinates(smap->get(plane, state1, state2, chunk), dim);
        for(int i = 0; i < ndims; i++) {
          sum[i] += dim[i];
        }
        points++;
      }
    for(int i = 0; i < ndims; i++) {
      sum[i] /= points;
    }
    sum[ndims] = 0;
    return(avail->minDist(sum));
  }
  else {
    int sumPe = 0;
    for(int plane=0; plane<nplanes; plane++)
      for(int chunk=0; chunk<numChunks; chunk++)
      {
        sumPe += smap->get(plane, state1, state2, chunk);
        points++;
      }
    bestPe = sumPe/points;
    return(avail->minDist(bestPe));
  }
}

