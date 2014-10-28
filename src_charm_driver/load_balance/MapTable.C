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
extern TopoManager *topoMgr;
extern Config config;

//============================================================================
int MapType1::getCentroid(int torusMap) {
  int points, bestPe;
  if(torusMap==1) {
    int sumX=0;
    int sumY=0;
    int sumZ=0;
    int X=0,Y=0,Z=0,T=0;
    points=0;
    for(int i=0;i<getXmax();i++)
    {
      CkAssert(get(i)>=0);
      topoMgr->rankToCoordinates(get(i), X, Y, Z, T);
      sumX+=X;
      sumY+=Y;
      sumZ+=Z;
      points++;
    }
    int avgX=sumX/points;
    int avgY=sumY/points;
    int avgZ=sumZ/points;
    bestPe=topoMgr->coordinatesToRank(avgX, avgY, avgZ, 0);
  }
  else {
    int sum=0;
    points=0;
    for(int i=0;i<getXmax();i++)
    {
      CkAssert(get(i)>=0);
      sum+=get(i);
      points++;
    }
    bestPe=sum/points;
  }
  return(bestPe);
}

int MapType2::getCentroid(int torusMap) {
  int points, bestPe;
  if(torusMap==1) {
    int sumX=0;
    int sumY=0;
    int sumZ=0;
    int X=0,Y=0,Z=0,T=0;
    points=0;
    for(int i=0;i<getXmax();i++)
      for(int j=0;j<getYmax();j++) {
        CkAssert(get(i,j)>=0);
        topoMgr->rankToCoordinates(get(i, j), X, Y, Z, T);
        sumX+=X;
        sumY+=Y;
        sumZ+=Z;
        points++;
      }
    int avgX=sumX/points;
    int avgY=sumY/points;
    int avgZ=sumZ/points;
    bestPe=topoMgr->coordinatesToRank(avgX, avgY, avgZ, 0);
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
    bestPe=sum/points;
  }
  return(bestPe);
}

int MapType3::getCentroid(int torusMap) {
  int points, bestPe;
  if(torusMap==1) {
    int sumX=0;
    int sumY=0;
    int sumZ=0;
    int X=0,Y=0,Z=0,T=0;
    points=0;
    for(int i=0;i<getXmax();i++)
      for(int j=0;j<getYmax();j++) {
        for(int k=0;k<getZmax();k++) {
          CkAssert(get(i,j,k)>=0);
          topoMgr->rankToCoordinates(get(i, j, k), X, Y, Z, T);
          sumX+=X;
          sumY+=Y;
          sumZ+=Z;
          points++;
        }
      }
    int avgX=sumX/points;
    int avgY=sumY/points;
    int avgZ=sumZ/points;
    bestPe=topoMgr->coordinatesToRank(avgX, avgY, avgZ, 0);
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
    bestPe=sum/points;
  }
  return(bestPe);
}

void IntMap1::translate(IntMap1 *fromMap, int offsetX, int offsetY, int offsetZ, bool torus ) 
{

  keyXmax=fromMap->keyXmax;
  CkAssert(keyXmax>0);
  Map= new int[keyXmax];
  if(torus)
  {
    int x, y, z, t, destpe;
    for(int xind=0; xind<keyXmax; xind++)
    {
      topoMgr->rankToCoordinates(fromMap->get(xind), x, y, z, t);
      int newx=(x+offsetX)%topoMgr->getDimNX();
      int newy=(y+offsetY)%topoMgr->getDimNY();
      int newz=(z+offsetZ)%topoMgr->getDimNZ();
      destpe =  topoMgr->coordinatesToRank(newx, newy, newz, t);
      set(xind, destpe);
    }
  }
  else
  {
    for(int xind=0; xind<keyXmax; xind++)
      set(xind,fromMap->get(xind)+offsetX);
  }
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

  if(torus)
  {
    int x, y, z, t, destpe;
    for(int xind=0; xind<keyXmax; xind++)
      for(int yind=0; yind<keyYmax; yind++) {
        topoMgr->rankToCoordinates(fromMap->get(xind, yind), x, y, z, t);
        int newx=(x+offsetX)%topoMgr->getDimNX();
        int newy=(y+offsetY)%topoMgr->getDimNY();
        int newz=(z+offsetZ)%topoMgr->getDimNZ();
        destpe =  topoMgr->coordinatesToRank(newx, newy, newz, t);
        set(xind, yind, destpe);
      }
  }
  else
  {
    for(int xind=0; xind<keyXmax; xind++)
      for(int yind=0; yind<keyYmax; yind++) {
        set(xind, yind,fromMap->get(xind, yind)+offsetX);
      }
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
  if(torus)
  {
    int x, y, z, t, destpe;
    for(int xind=0; xind<keyXmax; xind++)
      for(int yind=0; yind<keyYmax; yind++) {
        for(int zind=0; zind<keyZmax; zind++) {
          topoMgr->rankToCoordinates(fromMap->get(xind, yind,zind), x, y, z, t);
          int newx=(x+offsetX)%topoMgr->getDimNX();
          int newy=(y+offsetY)%topoMgr->getDimNY();
          int newz=(z+offsetZ)%topoMgr->getDimNZ();
          destpe =  topoMgr->coordinatesToRank(newx, newy, newz, t);
          set(xind, yind, zind, destpe);
        }
      }
  }
  else
  {
    for(int xind=0; xind<keyXmax; xind++)
      for(int yind=0; yind<keyYmax; yind++) {
        for(int zind=0; zind<keyZmax; zind++) {
          set(xind, yind, zind, fromMap->get(xind, yind,zind)+offsetX);
        }
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
  if(torus)
  {
    int x, y, z, t, destpe;
    for(int wind=0; wind<keyWmax; wind++)
      for(int xind=0; xind<keyXmax; xind++){
        for(int yind=0; yind<keyYmax; yind++) {
          for(int zind=0; zind<keyZmax; zind++) {
            topoMgr->rankToCoordinates(fromMap->get(wind,xind*keyStep, yind*keyStep,zind), x, y, z, t);
            int newx=(x+offsetX)%topoMgr->getDimNX();
            int newy=(y+offsetY)%topoMgr->getDimNY();
            int newz=(z+offsetZ)%topoMgr->getDimNZ();
            destpe =  topoMgr->coordinatesToRank(newx, newy, newz, t);
            set(wind, xind*keyStep, yind*keyStep, zind, destpe);
          }
        }
      }
  }
  else
  {
    for(int wind=0; wind<keyWmax; wind++)
      for(int xind=0; xind<keyXmax; xind++){
        for(int yind=0; yind<keyYmax; yind++) {
          for(int zind=0; zind<keyZmax; zind++) {
            set(wind, xind*keyStep, yind*keyStep, zind, fromMap->get(wind, xind*keyStep, yind*keyStep,zind)+offsetX);
          }
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
};


AtomMapTable::AtomMapTable(MapType1 *_tomap, PeList *availprocs, int numInst,
    int _nchareAtoms): nchareAtoms(_nchareAtoms)
{
  maptable = _tomap;
  if(numInst == 0) // first instance does the hard stuff
  {  
    for(int element=0; element<nchareAtoms; element++)
    {
      if(availprocs->count()==0)
        availprocs->reset();
      maptable->set(element,availprocs->findNext());
    }
  }
  else
  { //
    CkAbort("this should only called on the first instance");
  }
}


GSMapTable::GSMapTable(MapType2 *_frommap, MapType2 *_tomap, PeList *_availprocs, 
    int _nchareG, int _nstates, int _Gstates_per_pe, bool useCuboidMap, int numInst,
    int offsetX, int offsetY, int offsetZ):
  nchareG(_nchareG), nstates(_nstates), Gstates_per_pe(_Gstates_per_pe)
{
  reverseMap = NULL;
  maptable = _tomap;
  availprocs = _availprocs;

  /** The first instance creates the map and the other instances just use the
   *  map with a translation
   */
  if(numInst == 0) {
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
      while(pow(2.0, (double)i) < pl)
        i++;
      // make it same as the nearest smaller power of 2
      pl = (int) pow(2.0, (double)(i-1));
      srem = nstates % pl;
    }
    pm = availprocs->count() / pl;		// no of procs on x axis

    if(!useCuboidMap && pm == 0) {
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
    if(useCuboidMap)
    {
      int *Pecount= new int[config.numPesPerInstance];
      bzero(Pecount, config.numPesPerInstance *sizeof(int));
      PeList *exclusionList = NULL;

      // here we require that nchareG tiles
      if(config.numPesPerInstance % nchareG != 0)
      {
        CkPrintf("To use CuboidMap nchareG %d should be set as a factor of numprocs %d using gExpandFact\n", nchareG, config.numPesPerInstance);
        CkExit();
      }
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
        PeList *planeProcs=new PeList(availprocs, plane*procsPerPlane, procsPerPlane);
        if(planeProcs->count()==0)
          planeProcs->reset();
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
              if(Pecount[destpe]>=charesperpe)
              {
                if(exclusionList==NULL)
                {
                  exclusionList=new PeList(1);
                  exclusionList->TheList[0]=destpe;
                  exclusionList->sortIdx[0]=0;
                }
                else
                {
                  exclusionList->mergeOne(destpe);
                  PeList one(1);
                  one.TheList[0]=destpe;
                  *planeProcs - one;
                  planeProcs->reindex();
                  planeProcs->reset();
                }
              }
              if(((stateperpe+1<workingGsPerPe)&&((state+stateperpe+1)<nstates)) || state+workingGsPerPe<nstates)
              {
                // we will need another proc from this list
                destpe=planeProcs->findNext();
                if(planeProcs->count()==0)
                {
                  planeProcs->reset();
                  if(planeProcs->count()==0)
                    CkPrintf("GSMap exceeding count on plane %d state %d after pe %d workingGsPerPe %d cubesrem %d\n",state,plane,destpe, workingGsPerPe, cubesrem);
                }
              }
            }
          }
          if(!cubesrem)
          {
            if(planeProcs->count()==0)
              planeProcs->reset();
            destpe=planeProcs->findNext();
          }
          /*if(availprocs->count()==0)
            {
            CkPrintf("GSMap created on processor %d\n", CkMyPe());
            dump();
            CkPrintf("GSMap ran out of nodes on plane %d state %d\n", plane, state);
            CkExit();
            }
           */
        }
        delete planeProcs;
      }
      delete [] Pecount;
      if(exclusionList!=NULL)
        delete exclusionList;
    }
    else
    {
      int destpe=availprocs->findNext();

      // foreach statechunk 
      //         foreach state in chunk
      //              map it
      //         new pe
      // done
      //
      //else old way 
      for(int ychunk=0; ychunk<nchareG; ychunk=ychunk+m)
      {
        if(ychunk==(pm-rem)*m)
          m=m+1;
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
          if(availprocs->count()==0)
            availprocs->reset();
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
  } else { // not instance 0
    if(config.torusMap)
    {
      int x, y, z, t, destpe;
      CkPrintf("{%d} GS using offsets X=%d Y=%d Z=%d\n",numInst, offsetX, offsetY, offsetZ);
      for(int state=0; state<nstates; state++)
        for(int plane=0; plane<nchareG; plane++) {
          topoMgr->rankToCoordinates(_frommap->get(state, plane), x, y, z, t);
          int newx=(x+offsetX)%topoMgr->getDimNX();
          int newy=(y+offsetY)%topoMgr->getDimNY();
          int newz=(z+offsetZ)%topoMgr->getDimNZ();
          destpe =  topoMgr->coordinatesToRank(newx, newy, newz, t);
          maptable->set(state, plane, destpe);
        }
    }
    else
    {
      CkPrintf("WARNING: using co-mapping for instances because I'm too lazy to partition the processors in the non topo case\n");
      for(int state=0; state<nstates; state++)
        for(int plane=0; plane<nchareG; plane++) {
          maptable->set(state,plane,_frommap->get(state, plane));
        }
    }
  }
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
  reverseMap=NULL;
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
    if(availprocs->count()==0)
      availprocs->reset();

    //if(CkMyPe()==0) CkPrintf("scobjs_per_pe %d grainsize %d nchareG %d scalc_per_plane %d planes_per_pe %d numChunks %d rem %d\n", scobjs_per_pe, grainsize, nchareG, scalc_per_plane, planes_per_pe, numChunksSym, rem);
    if(useCuboidMap|| useCentroid)
    { // in the cuboid map case we place all planes box by box
      for(int plane=0; plane<nchareG; plane++)
      { // could restrict list to the gs box here

        PeList *thisPlaneBox;
        if(useCentroid)
        {

          thisPlaneBox= subListPlane(plane, max_states, gsmap);
          thisPlaneBox->reset();
        }
        else
        {

          //		  CkPrintf("plane %d making pelist from boxSize %d\n",plane,boxSize);
          thisPlaneBox= new PeList(availprocs, plane*boxSize, boxSize);
        }

        //	      PeList *thisPlaneBox= availprocs;
        //	      thisPlaneBox->dump();
        if(!useCentroid)
        {
          destpe=thisPlaneBox->findNext();
          if(thisPlaneBox->count()==0)
            thisPlaneBox->reset();

        }
        for(int xchunk=0; xchunk<maxstateindex; xchunk=xchunk+grainsize)
          for(int ychunk=xchunk; ychunk<maxstateindex; ychunk=ychunk+grainsize)
          { // could find centroid here
            if(useCentroid)
            {
              thisPlaneBox->trimUsed();
              //			CkPrintf("trimmed list for %d %d %d\n", plane, xchunk, ychunk);
              //			thisPlaneBox->dump();
              sortByCentroid(thisPlaneBox, plane, xchunk, ychunk, grainsize, gsmap);
              //			CkPrintf("sorted by centroid\n");
              //			thisPlaneBox->dump();
              destpe=thisPlaneBox->findNext();
              if(thisPlaneBox->count()==0)
                thisPlaneBox->reset();
              count=0;
            }
            for(int newdim=0; newdim<numChunksSym; newdim++)
            {
              if(count<scobjs_per_pe)
              {
                // nothing to see here
              }
              else
              {  // new partition
                procno++;
                if(thisPlaneBox->count()==0)
                  thisPlaneBox->reset();

                destpe=thisPlaneBox->findNext();
                if(thisPlaneBox->count()==0)
                  thisPlaneBox->reset();

                if(rem!=0)
                  if(procno==rem)
                    scobjs_per_pe-=1;
                count=0;

              }
              //			CkPrintf("%d %d %d %d mapped to %d\n",plane,xchunk,ychunk,newdim, destpe);
              CkAssert(destpe<config.numPes);
#ifdef USE_INT_MAP
              maptable->set(plane, xchunk, ychunk, newdim,destpe);
#else
              CkArrayIndex4D idx4d(plane, xchunk, ychunk, newdim);
              CmiMemcpy(intidx,idx4d.index,2*sizeof(int));
              maptable->put(intdual(intidx[0], intidx[1]))=destpe;
#endif

              count++;

            }
          }
        delete thisPlaneBox;
      }

    }
    else
    {
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
                  //if(CkMyPe()==0) CkPrintf("plane %d x %d y %d newdim %d = proc %d\n", plane, xchunk, ychunk, newdim, fp.next[0]*x*y+fp.next[1]*x+fp.next[2]);
                  CkAssert(destpe<config.numPes);
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
                  if(availprocs->count()==0)
                    availprocs->reset();
                  //			  availprocs->sortSource(srcpe);
                  destpe=availprocs->findNext();
                  if(availprocs->count()==0)
                    availprocs->reset();

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
    } // else not cuboid
  }
  else
  {
    scobjs_per_pe = scalc_per_plane*nchareG*numChunksAsym/availprocs->count();
    rem = scalc_per_plane*nchareG*numChunksAsym % availprocs->count();


    if(rem!=0)
      scobjs_per_pe+=1;
#ifdef _MAP_DEBUG_
    CkPrintf(" scalc_per_plane %d *nchareG %d *numChunksAsym %d  availprocs->count() %d = rem %d and scobjs_per_pe is %d boxSize %d\n", scalc_per_plane,nchareG,numChunksAsym , availprocs->count(),rem, scobjs_per_pe, boxSize);
#endif
    //if(CkMyPe()==0) CkPrintf("scobjs_per_pe %d grainsize %d nchareG %d scalc_per_plane %d planes_per_pe %d numChunksAsym %d rem %d\n", scobjs_per_pe, grainsize, nchareG, scalc_per_plane, planes_per_pe, numChunksAsym, rem);
    int srcpe=0,destpe=0;
    if(!useCentroid && !useCuboidMap)
      destpe=availprocs->findNext();
    if(availprocs->count()==0)
      availprocs->reset();
    if(useCuboidMap|| useCentroid)
    { // in the cuboid map case we place all planes box by box
      for(int plane=0; plane<nchareG; plane++)
      { // could restrict list to the gs box here

        PeList *thisPlaneBox;
        if(useCentroid)
        {

          thisPlaneBox= subListPlane(plane, max_states, gsmap);
          //		  CkPrintf("plane %d making pelist with %d pes from plane's pes\n",plane, thisPlaneBox->count());
          //		  thisPlaneBox->dump();
        }
        else
        {
          //		  CkPrintf("plane %d making pelist from boxSize %d\n",plane,boxSize);
          thisPlaneBox= new PeList(availprocs, plane*boxSize, boxSize);
        }

        //	      PeList *thisPlaneBox= availprocs;
        //	      thisPlaneBox->dump();
        if(!useCentroid)
        {
          destpe=thisPlaneBox->findNext();
          if(thisPlaneBox->count()==0)
            thisPlaneBox->reset();
        }

        for(int xchunk=0; xchunk<maxstateindex; xchunk=xchunk+grainsize)
          for(int ychunk=0; ychunk<maxstateindex; ychunk=ychunk+grainsize)
          { // could find centroid here
            if(useCentroid)
            {
              thisPlaneBox->trimUsed();
              //			CkPrintf("plane %d xchunk %d ychunk %d trim\n",plane,xchunk,ychunk);
              //			thisPlaneBox->dump();
              sortByCentroid(thisPlaneBox, plane, xchunk, ychunk, grainsize, gsmap);
              //			CkPrintf("plane %d xchunk %d ychunk %d sortbycentroid\n",plane,xchunk,ychunk);
              //			thisPlaneBox->dump();
              destpe=thisPlaneBox->findNext();
              if(thisPlaneBox->count()==0)
                thisPlaneBox->reset();
              count=0;
            }
            for(int newdim=0; newdim<numChunksAsym; newdim++)
            {

              if(count<scobjs_per_pe)
              {
                // nothing to see here
              }
              else
              {  // new partition
                procno++;
                destpe=thisPlaneBox->findNext();
                if(thisPlaneBox->count()==0)
                  thisPlaneBox->reset();

                if(rem!=0&& scobjs_per_pe>1)
                  if(procno==rem)
                    scobjs_per_pe-=1;
                count=0;

              }
              //			CkPrintf("%d %d %d %d mapped to %d\n",plane,xchunk,ychunk,newdim, destpe);
#ifdef USE_INT_MAP
              maptable->set(plane, xchunk, ychunk, newdim,destpe);
#else
              CkArrayIndex4D idx4d(plane, xchunk, ychunk, newdim);
              CmiMemcpy(intidx,idx4d.index,2*sizeof(int));
              maptable->put(intdual(intidx[0], intidx[1]))=destpe;
#endif

              count++;

            }
          }
        delete thisPlaneBox;
      }

    }
    else
    {
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
                  //if(CkMyPe()==0) CkPrintf("plane %d x %d y %d newdim %d= proc %d\n", plane, xchunk, ychunk, newdim, assign[0]*x*y+assign[1]*x+assign[2]);
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
                  if(availprocs->count()==0)
                    availprocs->reset();
                  //			  availprocs->sortSource(srcpe);
                  destpe=availprocs->findNext();
                  if(availprocs->count()==0)
                    availprocs->reset();

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
    }
#ifdef _MAP_DEBUG_
    CkPrintf("Asymmetric SCalcMap created on processor %d\n", CkMyPe());
    dump();
#endif
  }

}

/** helper function
 */
PeList * rebuildExclusion(int *Pecount, int rsobjs_per_pe)
{
  PeList *exclusionList=NULL;
  for(int exc=0; exc<config.numPesPerInstance; exc++)
  {
    if(Pecount[exc]>=rsobjs_per_pe)
    {
      if(exclusionList==NULL)
      {
        exclusionList=new PeList(1);
        exclusionList->TheList[0]=exc;
        exclusionList->sortIdx[0]=0;
      }
      else
      {
        exclusionList->mergeOne(exc);
      }
    }
  }
  return(exclusionList);
}


RSMapTable::RSMapTable(MapType2  *_frommap, MapType2 *_tomap, PeList *_availprocs,
    int _nstates, int _sizeZ, int _Rstates_per_pe, bool useCuboidMap, MapType2 *gsmap, 
    int nchareG, int numInst, int offsetX, int offsetY, int offsetZ):
  nstates(_nstates), sizeZ(_sizeZ),
  Rstates_per_pe(_Rstates_per_pe)
{
  reverseMap = NULL;
  maptable = _tomap;
  availprocs = _availprocs;
  availprocs->reset();

  if(numInst == 0) {
    int l, m, pl, pm, srem, rem, i=0, rsobjs_per_pe;
    int *Pecount= new int [config.numPes];

    bzero(Pecount, config.numPes*sizeof(int));

    rsobjs_per_pe = nstates*sizeZ/config.numPesPerInstance;
    if(config.useStrictCuboid) 
      rsobjs_per_pe++; // you'll need the wiggle room
    l = Rstates_per_pe;		// no of states in one chunk
    pl = nstates / l;
    if(nstates % l == 0)
      srem = 0;
    else
    {
      while(pow(2.0, (double)i) < pl)
        i++;
      pl = (int) pow(2.0, (double)(i-1));		// make it same as the nearest smaller power of 2
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

    //CkPrintf("nstates %d sizeZ %d Pes %d rsobjs_per_pe %d\n", nstates, sizeZ, availprocs->count(), rsobjs_per_pe);	
    //CkPrintf("l %d, m %d pl %d pm %d srem %d rem %d\n", l, m,
    //	   pl, pm, srem, rem);
    int srcpe=0;
    int destpe=availprocs->findNext();

    if(availprocs->count()==0)
      availprocs->reset();	
    if(useCuboidMap)
    {
      int srem = (nstates*sizeZ) % config.numPesPerInstance;
      //      if(srem)
      //	rsobjs_per_pe++;

      // exclusion mapping has the sad side effect of increasing the
      // number of exclusions with the state and plane number
      // until you end up increasing the cap too high
      // Topo mapping is less important than even distribution
      // so if you go over the cap, use the master list instead of the
      // state box.

      // this has the effect of creating an imbalance

      PeList *myavail=new PeList(*availprocs);
      PeList *exclusionList = NULL;
      for(int state=0; state < nstates; state++)
      {

        int srsobjs_per_pe=rsobjs_per_pe;
        PeList *thisStateBox = subListState(state, nchareG, gsmap);
        int samplePe=thisStateBox->TheList[0];
        //  CkPrintf("RS state %d box has %d procs rsobjs_per_pe %d\n",state,thisStateBox->count(), rsobjs_per_pe);
        //	  bool useExclude=false;
        bool useExclude=true;
        if(exclusionList!=NULL && useExclude)
        {
          *thisStateBox - *exclusionList;
          thisStateBox->reindex();
          thisStateBox->reset();
        }
        //CkPrintf("RS state %d box has %d procs after exclusion\n",state,thisStateBox->count());
        // CkPrintf("RS state %d pe list after exclude \n",state);
        //	  thisStateBox->dump();
        if(thisStateBox->count()<=0)
        {
          if(config.useStrictCuboid)
          {
            // use old scheme of bumping srsobjs_per_pe 
            while(thisStateBox->count()<=0 && srsobjs_per_pe<=sizeZ*nstates)
            {
              CkPrintf("State %d  Ran out of procs in RS centroid map increasing rs objects per proc to %d\n",state,srsobjs_per_pe);
              srsobjs_per_pe++;
              if(exclusionList!=NULL)
                delete exclusionList;
              exclusionList = rebuildExclusion(Pecount, srsobjs_per_pe);
              delete thisStateBox;
              thisStateBox = subListState(state, nchareG, gsmap);
              if(exclusionList!=NULL){
                *thisStateBox - *exclusionList;
                thisStateBox->reindex();
                thisStateBox->reset();
              }

            }
          }
          else

          {
            CkPrintf("State %d  Ran out of procs in RS centroid map scheme, spilling over to master list\n",state,srsobjs_per_pe);
            //	      *myavail - *exclusionList;
            delete thisStateBox;
            thisStateBox=new PeList(*myavail);
            thisStateBox->reindex();
            thisStateBox->reset();
            thisStateBox->sortSource(samplePe);
          }
        }
        if(thisStateBox->count()==0)
        {
          useExclude=false;
          CkPrintf("RS state %d excluded to nil, ignoring exclusions\n",state);
          delete thisStateBox;
          thisStateBox = subListState(state, nchareG, gsmap);
        }

        // sort by centroid (not necessary)

        for(int plane=0; plane < sizeZ; plane++)
        {
          if(thisStateBox->count()==0)
            thisStateBox->reset();
          if(thisStateBox->count()<=0)
          {
            if(config.useStrictCuboid)
            {
              while(thisStateBox->count()<=0 && srsobjs_per_pe<=sizeZ*nstates)
              {
                srsobjs_per_pe++;
                CkPrintf("State %d Plane %d Ran out of procs in RS centroid map increasing rs objects per proc to %d\n",state,plane,srsobjs_per_pe);
                if(exclusionList!=NULL)
                  delete exclusionList;
                exclusionList=rebuildExclusion(Pecount, srsobjs_per_pe);
                delete thisStateBox;
                thisStateBox = subListState(state, nchareG, gsmap);
                if(exclusionList!=NULL && useExclude)
                {
                  *thisStateBox - *exclusionList;
                  thisStateBox->reset();
                  thisStateBox->reindex();
                }
              }
            }
            else
            {
              CkPrintf("State %d  Ran out of procs in RS centroid map scheme, spilling over to master list\n",state,srsobjs_per_pe);
              if(myavail->count()<=0)
              {
                srsobjs_per_pe++;
                CkPrintf("State %d  Ran out of procs in master, bumping srsobjs_per_pe\n",state,srsobjs_per_pe);
                delete thisStateBox;
                delete myavail;
                myavail= new PeList(*availprocs);
                if(exclusionList!=NULL)
                  delete exclusionList;

                exclusionList=rebuildExclusion(Pecount, srsobjs_per_pe);
                if(exclusionList!=NULL)
                {
                  *myavail - *exclusionList;
                }
                thisStateBox = new PeList(*myavail);
              }
              else
              {
                delete thisStateBox;
                thisStateBox=new PeList(*myavail);
              }
              thisStateBox->reindex();
              thisStateBox->reset();
              thisStateBox->sortSource(samplePe);
            }
          }
          if(useExclude && thisStateBox->count()<=0)
          { // surrender
            useExclude=false;
            delete thisStateBox;
            CkAbort("RS cuboid map hopeless please examine configuration\n");

          }

          destpe = thisStateBox->findNext();

#ifdef USE_INT_MAP
          maptable->set(state, plane, destpe);
#else						
          maptable->put(intdual(state, plane))= destpe;
#endif
          // if(CkMyPe()==0) CkPrintf("%d %d [%d]\n", state, plane, destpe);
          // CkAssert(destpe < config.numPesPerInstance);
          Pecount[destpe]++;
          if(Pecount[destpe]>=srsobjs_per_pe)
          {
            if(exclusionList==NULL)
            {
              exclusionList=new PeList(1);
              exclusionList->TheList[0]=destpe;
              exclusionList->sortIdx[0]=0;
            }
            else
            {
              exclusionList->mergeOne(destpe);
              PeList one(1);
              one.TheList[0]=destpe;
              *thisStateBox - one;
              *myavail - one;
              thisStateBox->reindex();
              thisStateBox->reset();

            }
          }
        }
        delete thisStateBox;
      }
      if(exclusionList!=NULL)
        delete exclusionList;
      if(myavail != NULL)
        delete myavail;
    }
    else
    {

      // this remainder scheme is odd, creates imbalance.
      // doesn't use all processors
      //destpe=availprocs->findNext();
      for(int ychunk=0; ychunk<sizeZ-rem; ychunk=ychunk+m)
      {
        /*if(ychunk==(pm-rem)*m)
          m=m+1;*/
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
          if(availprocs->count()==0)
            availprocs->reset();	

        }
      }
      if(rem!=0)
        for(int state=0; state<nstates; state++)
        {
          for(int plane=sizeZ-rem; plane<sizeZ; plane++)
          {
            if(availprocs->count()==0)
              availprocs->reset();
            destpe=availprocs->findNext();
#ifdef USE_INT_MAP
            maptable->set(state, plane, destpe);
#else						
            maptable->put(intdual(state, plane)) = destpe;
#endif
          }
        }
    }
    delete [] Pecount;
#ifdef _MAP_DEBUG_

    CkPrintf("RSMap created on processor %d\n", CkMyPe());
    dump();
#endif
  } else { // not instance 0
    if(config.torusMap==1)
    {
      int x, y, z, t, destpe;
      CkPrintf("{%d} RS using offsets X=%d Y=%d Z=%d\n",numInst, offsetX, offsetY, offsetZ);
      for(int state=0; state<nstates; state++)
        for(int plane=0; plane<sizeZ; plane++) {
          topoMgr->rankToCoordinates(_frommap->get(state, plane), x, y, z, t);
          int newx=(x+offsetX)%topoMgr->getDimNX();
          int newy=(y+offsetY)%topoMgr->getDimNY();
          int newz=(z+offsetZ)%topoMgr->getDimNZ();
          destpe =  topoMgr->coordinatesToRank(newx, newy, newz, t);
          maptable->set(state, plane, destpe);
        }
    }
    else
    {
      CkPrintf("WARNING: using co-mapping for instances because I'm too lazy to partition the processors in the non topo case\n");
      for(int state=0; state<nstates; state++)
        for(int plane=0; plane<sizeZ; plane++) {
          maptable->set(state,plane,_frommap->get(state, plane));
        }
    }
  }
}


RPPMapTable::RPPMapTable(MapType2  *_map, 
    PeList *_availprocs, PeList *exclusion,
    int _nstates, int _sizeZNL, int _Rstates_per_pe,
    int boxSize, bool usePPmap, int nchareG,
    MapType2 *pp_map) :
  nstates(_nstates), sizeZNL(_sizeZNL),
  Rstates_per_pe(_Rstates_per_pe)
{
  int states_per_pe=Rstates_per_pe;
  int totalChares=nstates*sizeZNL;
  reverseMap=NULL;
  maptable=_map;
  availprocs=_availprocs;
  //        if(useCuboidMap)
  //	  states_per_pe=nstates/boxSize;		// no of states in one chunk
  //	else
  //	  states_per_pe=nstates/boxSize;		// no of states in one chunk
  //	else
  bool useExclusion=true;
  PeList *RPPlist=availprocs;
  int chares_per_pe=totalChares/config.numPesPerInstance;
  //  pl = nstates / states_per_pe;
  CkPrintf("CharesPerPe %d states_per_pe %d\n",chares_per_pe, states_per_pe);
  if(exclusion==NULL || exclusion->count()==0 || config.numPesPerInstance <=exclusion->count() )
    useExclusion=false;
  int afterExclusion=availprocs->count();
  if(useExclusion)
    afterExclusion=availprocs->count() - exclusion->count();
  //  CkPrintf("RPP excluded list\n");
  //  exclusion->dump();
  if(useExclusion && afterExclusion > chares_per_pe*sizeZNL)
  { // we can fit the exclusion without blinking
    CkPrintf("RPP using density exclusion to avoid %d processors\n",exclusion->count());
    *RPPlist-*exclusion;
    RPPlist->reindex();
  }
  else
  {// so an rstates_per_pe chosen for realstate might be too big
    if(useExclusion && afterExclusion > (sizeZNL+nstates)*2)
    { // set states_per_pe to fit the exclusion

      states_per_pe=(int) ((float) (totalChares/afterExclusion))*0.75;
      CkPrintf("RPP adjusting states per pe from %d to %d to use density exclusion to stay within %d processors\n",Rstates_per_pe, states_per_pe, exclusion->count());
      *RPPlist-*exclusion;
      RPPlist->reindex();
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
  // CkPrintf("nstates %d sizeZNL %d Pes %d\n", nstates, sizeZNL, RPPlist->count());	
  // RPPlist->dump();
  if(usePPmap)
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
      CkAssert(state_map->count()>0);
      if(useExclusion)
      {
        *state_map-*exclusion;
        state_map->reindex();
        if(state_map->count()==0 || sizeZNL/state_map->count() > states_per_pe)
        { //not enough for exclusion
          delete state_map;
          state_map=subListState( state, nchareG, pp_map);
          useExclude[state]=false;
        }
        else
        {
          useExclude[state]=true;
        }
      }
      else
      {
        useExclude[state]=false;
      }
      int thischaresperpe=sizeZNL/state_map->count() + 1;
      //	  CkPrintf("state %d has %d charesperpe in %d pemap \n",state,thischaresperpe,state_map->count());


      maxcharesperpe=(thischaresperpe>maxcharesperpe) ? thischaresperpe : maxcharesperpe;

      // if(states_per_pe>maxcharesperpe)
      //	    maxcharesperpe=states_per_pe;
      delete state_map;
    }
    int totcharesperpe=sizeZNL*nstates/config.numPesPerInstance+1;
    maxcharesperpe=(maxcharesperpe>totcharesperpe) ? maxcharesperpe : totcharesperpe;
    //      maxcharesperpe++;
    PeList *usedbyRPP=NULL;
    CkPrintf("RPP maxcharesperpe is %d\n",maxcharesperpe);
    int origmaxcharesperpe=maxcharesperpe;
    PeList *excludedBigmap=NULL;
    for(int state=0; state < nstates ; state++)
    {
      PeList *state_map=subListState( state, nchareG, pp_map);
      if(useExclude[state])
      {
        *state_map-*exclusion;
        state_map->reindex();
      }
      while(state_map->count()<=0 && maxcharesperpe<=sizeZNL*nstates)
      {
        if(config.useStrictCuboid)
        {
          CkPrintf("State %d  Ran out of procs in RPP centroid map increasing rpp objects per proc to %d\n",state,maxcharesperpe);
          maxcharesperpe++;
          if(usedbyRPP!=NULL)
            delete usedbyRPP;
          usedbyRPP=rebuildExclusion(usedPes, maxcharesperpe);
          delete state_map;
          state_map = subListState(state, nchareG, pp_map);
          if(usedbyRPP!=NULL){
            *state_map - *usedbyRPP;
            state_map->reindex();
            state_map->reset();
          }
        }
        else
        {
          // ditch topo scheme for overflow
          // use the RPPlist
          CkPrintf("State %d  Ran out of procs in RPP centroid using full RPPlist\n",state);
          delete state_map;
          if(!neednewexc && excludedBigmap!=NULL)
          {
            state_map= new  PeList(*excludedBigmap);
          }
          else
          {
            state_map = new PeList(*RPPlist);
            if(usedbyRPP!=NULL){
              *state_map - *usedbyRPP;
              state_map->reindex();
              state_map->reset();
              if(excludedBigmap!=NULL)
                delete excludedBigmap;
              excludedBigmap=new PeList(*state_map);
            }
          }
          if(state_map->count()<=0)
          {
            maxcharesperpe++;
            CkPrintf("plane %d  Ran out of procs in RPP centroid using full RPPlist and bumping maxcharesperpe to %d\n",state, maxcharesperpe);
            if(usedbyRPP!=NULL)
              delete usedbyRPP;
            usedbyRPP=rebuildExclusion(usedPes, maxcharesperpe);
            delete state_map;
            state_map = new PeList(*RPPlist);
            if(usedbyRPP!=NULL && usedbyRPP->count()>0){
              *state_map - *usedbyRPP;
              state_map->reindex();
              state_map->reset();
            }
            if(excludedBigmap!=NULL)
              delete excludedBigmap;
            excludedBigmap= new PeList(*state_map);
            neednewexc=false;
          }
          // ADD SORT here
        }
        state_map->sortSource(pp_map->get(state,0));
      }
      for(int plane=0; plane < sizeZNL; plane++)
      {
        if(state_map->count()==0)
          state_map->reset();
        while(state_map->count()<=0 && maxcharesperpe<=sizeZNL*nstates)
        {
          if(config.useStrictCuboid)
          {
            while(state_map->count()<=0 && maxcharesperpe<=sizeZNL*nstates)
            {
              CkPrintf("State %d  Ran out of procs in RPP centroid map increasing rpp objects per proc to %d\n",state,maxcharesperpe);
              maxcharesperpe++;
              if(usedbyRPP!=NULL)
                delete usedbyRPP;
              usedbyRPP=rebuildExclusion(usedPes, maxcharesperpe);
              delete state_map;
              state_map = subListState(state, nchareG, pp_map);
              if(usedbyRPP!=NULL){
                *state_map - *usedbyRPP;
                state_map->reindex();
                state_map->reset();
              }

            }
          }
          else
          {

            CkPrintf("State %d  Ran out of procs in RPP centroid using full RPPlist\n",state);
            delete state_map;
            if(!neednewexc && excludedBigmap!=NULL&& excludedBigmap->count()>0)
            {
              state_map= new PeList(*excludedBigmap);
            }
            else
            {
              state_map = new PeList(*RPPlist);
              if(usedbyRPP!=NULL){
                *state_map - *usedbyRPP;
                state_map->reindex();
                state_map->reset();
                if(excludedBigmap!=NULL)
                  delete excludedBigmap;
                if(state_map->count()==0)
                { // man we're totally dry here.
                  // this should be handled in the next block
                  excludedBigmap=NULL;
                }
                else
                {
                  excludedBigmap= new PeList(*state_map);
                  neednewexc=false;
                }
              }
            }
            if(state_map->count()<=0)
            {
              maxcharesperpe++;
              CkPrintf("plane %d  Ran out of procs in RPP centroid using full RPPlist and bumping maxcharesperpe to %d\n",state, maxcharesperpe);

              if(usedbyRPP!=NULL)
                delete usedbyRPP;
              usedbyRPP=rebuildExclusion(usedPes, maxcharesperpe);
              delete state_map;
              state_map = new PeList(*RPPlist);
              if(usedbyRPP!=NULL)
              {
                *state_map - *usedbyRPP;
                state_map->reindex();
                state_map->reset();
              }
              if(state_map->count()==0)
              { // man we're totally dry here.
                // reboot
                CkPrintf("plane %d  Ran out of procs in RPP centroid using full RPPlist after bumping maxcharesperpe to %d, clearing used list, resetting maxcharesperpe to  %d\n",state, maxcharesperpe, origmaxcharesperpe);


                maxcharesperpe=origmaxcharesperpe;
                bzero(usedPes, config.numPes * sizeof(int));
                if(usedbyRPP!=NULL)
                  delete usedbyRPP;
                usedbyRPP=rebuildExclusion(usedPes, maxcharesperpe);
                delete state_map;
                state_map = new PeList(*RPPlist);
                if(usedbyRPP!=NULL)
                {
                  *state_map - *usedbyRPP;
                  state_map->reindex();
                  state_map->reset();
                }
              }

              if(excludedBigmap!=NULL)
                delete excludedBigmap;
              excludedBigmap= new PeList(*state_map);
              neednewexc=false;
            }
          }
          state_map->sortSource(pp_map->get(state,0));
          state_map->reset();
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
          if(usedbyRPP==NULL)
          { 
            usedbyRPP= new PeList(1);
            usedbyRPP->TheList[0]=destpe;
          }
          else
            usedbyRPP->mergeOne(destpe);
          usedbyRPP->reindex();
          exclusion->mergeOne(destpe);
          exclusion->reindex();
          PeList thisOne(1);
          thisOne.TheList[0]=destpe;
          if(excludedBigmap!=NULL)
          {
            neednewexc=false;
            *excludedBigmap-thisOne;
          }
          *state_map - thisOne;
          state_map->reindex();
          state_map->reset();
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
        //	      CkPrintf("RPP states_per_pe %d setting %d %d to pe %d\n", states_per_pe, state,plane,destpe);		    
#ifdef USE_INT_MAP
        maptable->set(state, plane, destpe);
#else
        maptable->put(intdual(state, plane))=destpe;
#endif
        if(++charesOnThisPe==chares_per_pe)
        {
          destpe=RPPlist->findNext();
          charesOnThisPe=0;
          if(RPPlist->count()==0)
            RPPlist->reset();
        }
      }
      if(rem && RPPlist->count()<2*nstates)
      {
        chares_per_pe=starting_cpp+1;
        //	      CkPrintf("At state %d count %d bumping chares_per_pe to %d\n",state, RPPlist->count(),chares_per_pe);
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
  oobjs_per_pe = northo/(config.numPesPerInstance);
  int *Pecount= new int [config.numPesPerInstance];
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
        *thisStateBox - *exclusionList;
        thisStateBox->reindex();
        thisStateBox->reset();
      }
      else
      {
        thisStateBox=availprocs;
      }
      if(thisStateBox->count() == 0)
      {  // the sublist scheme failed 
        CkPrintf("Ortho %d %d ignoring SubList\n", state1, state2);
        if(thisStateBox!=availprocs)
          delete thisStateBox;
        thisStateBox = availprocs;
        //*thisStateBox - *exclusionList;
        //thisStateBox->reindex();
        //useSublist=false;
      }

      if(thisStateBox->count() == 0)
      {
        CkPrintf("Ortho %d %d ignoring exclusion\n", state1, state2);
        if(thisStateBox!=availprocs)
          delete thisStateBox;
        thisStateBox = subListState2(s1, s2, nplanes, numChunks, scalcmap);
        useExclude = false;
      }

      /*sortByCentroid(thisStateBox, nplanes, s1, s2, numChunks, scalcmap);
        destpe=thisStateBox->findNext();
        if(thisStateBox->count()==0)
        thisStateBox->reset();
       */

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
        exclusionList->mergeOne(destpe);
        if(useExclude)
        {
          PeList one(1);
          one.TheList[0]=destpe;
          *availprocs - one;
          availprocs->reindex();
          availprocs->reset();
        }
      }
      if(thisStateBox!=availprocs)
        delete thisStateBox;
    }
  delete [] Pecount;
}

OrthoHelperMapTable::OrthoHelperMapTable(MapType2 *_map, int _nstates, int _orthograinsize, MapType2 *omap, PeList *_avail, PeList *exclude): nstates(_nstates), orthoGrainSize(_orthograinsize) 
{
  maptable = _map;
  int destpe = 0;
  // map orthohelper near but not on ortho by removing all ortho procs
  availprocs= new PeList(*_avail);
  *availprocs-*exclude;
  //  CkPrintf("exlude for helpers is \n");
  //  exclude->dump();
  bool useExclude=true;
  if(availprocs->count() < ((nstates/orthoGrainSize)*(nstates/orthoGrainSize)))
  {
    CkPrintf("There aren't enough processors to effectively parallelize orthoHelpers, you should disable orthoHelpers!!!\n");
    delete availprocs;
    availprocs= new PeList(*_avail);
    availprocs->reset();
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

        PeList one(1);
        one.TheList[0]=destpe;
        *availprocs - one;
        availprocs->reindex();
      }

    }
  delete availprocs;
}

RhoRSMapTable::RhoRSMapTable(MapType2  *_map, PeList *_availprocs, int _nchareRhoR, int _rhoRsubplanes, int max_states, bool useCentroid, MapType2 *rsmap, PeList *exclude): nchareRhoR(_nchareRhoR), rhoRsubplanes(_rhoRsubplanes)
{
  reverseMap=NULL;
  maptable=_map;
  availprocs=_availprocs;
  int rrsobjs_per_pe, rem;
  int srcpe=0;
  int numChares=nchareRhoR*rhoRsubplanes;
  int pesused=0;
  if(availprocs->count()==1)
  {
    rrsobjs_per_pe= numChares;
    rem=0;
  }
  else
  {
    rrsobjs_per_pe= numChares/(availprocs->count());
    rem = numChares % (availprocs->count());
    if(numChares<availprocs->count())
    {
      rem=0;
      rrsobjs_per_pe=1;
    }
    if(rem!=0)
      rrsobjs_per_pe += 1;
  }

  if(availprocs->count()==0)
    availprocs->reset();
  int destpe;
  int *Pecount = new int [config.numPes];
  bzero(Pecount, config.numPes*sizeof(int)); 

  //if(CkMyPe()==0) CkPrintf("nchareRhoR %d rrsobjs_per_pe %d rem %d\n", nchareRhoR, rrsobjs_per_pe, rem);   
  if(useCentroid) 
  {
    PeList *exclusionList = NULL;
    //CkAssert(numChares<availprocs->count());

    for(int chunk=0; chunk<nchareRhoR; chunk++)
    {
      PeList *thisPlaneBox = subListPlane(chunk, max_states, rsmap);
      bool useExclude=true;
      if(exclusionList!=NULL) {
        *thisPlaneBox - *exclusionList;
        thisPlaneBox->reindex();
      }

      if(thisPlaneBox->count()==0)
      {	
        CkPrintf("Rho RS %d ignoring plane sublist\n",chunk);
        thisPlaneBox = new PeList(*availprocs);
        if(exclusionList!=NULL) {
          *thisPlaneBox - *exclusionList;
          thisPlaneBox->reindex();
        }
        if(thisPlaneBox->count()==0)
        {
          CkPrintf("Rho RS %d ignoring plane sublist and exclusion\n",chunk);
          thisPlaneBox = new PeList(*availprocs);
        }
      }
      sortByCentroid(thisPlaneBox, chunk, max_states, rsmap);
      // CkPrintf("RhoR %d has %d procs from RS plane\n",chunk,thisPlaneBox->count());

      destpe=thisPlaneBox->findNext();
      if(thisPlaneBox->count()==0)
        thisPlaneBox->reset();
      for(int subplane=0 ; subplane<rhoRsubplanes ; subplane++)
      {
#ifdef USE_INT_MAP
        maptable->set(chunk, subplane, destpe);
#else
        maptable->put(intdual(chunk, subplane))=destpe;
#endif
        Pecount[destpe]++;	
        if(Pecount[destpe]>=rrsobjs_per_pe)
        {
          if(exclusionList==NULL)
          {
            exclusionList=new PeList(1);
            exclusionList->TheList[0]=destpe;
          }
          else
            exclusionList->mergeOne(destpe);
          if(useExclude && thisPlaneBox->size>1)
          {
            *thisPlaneBox - *exclusionList;
            thisPlaneBox->reindex();
          }
          sortByCentroid(thisPlaneBox, chunk, max_states, rsmap);
        }
        destpe=thisPlaneBox->findNext();
        if(thisPlaneBox->count()==0)
          thisPlaneBox->reset();
      }
      delete thisPlaneBox;
    }
    // now include the partially filled processors
    for(int i=0; i<config.numPesPerInstance;i++)
    {
      if(Pecount[i]>0)
      { 
        pesused++;
        if(Pecount[i] <rrsobjs_per_pe)
        {
          if(exclusionList==NULL)
          {
            exclusionList=new PeList(1);
            exclusionList->TheList[0]=i;
          }
          else
            exclusionList->mergeOne(i);
        }
      }
      exclusionList->reindex();
    }
    if(exclusionList!=NULL)
    {
      exclude->append(*exclusionList);
      delete exclusionList;
    }
  }
  else
  {
    int nprocs=0, objs=0;
    destpe=availprocs->findNext();
    if(availprocs->count()==0)
      availprocs->reset();
    for(int chunk=0; chunk<nchareRhoR; chunk++)
    {
      for(int subplane=0; subplane<rhoRsubplanes; subplane++)
      {
        if(rem!=0)
          if(nprocs==rem)
            rrsobjs_per_pe -= 1;
#ifdef USE_INT_MAP
        maptable->set(chunk, subplane, destpe);
#else
        maptable->put(intdual(chunk, subplane))=destpe;
#endif
        objs++;
        if(objs==rrsobjs_per_pe)
        {
          destpe=availprocs->findNext();
          if(availprocs->count()==0)
            availprocs->reset();
          objs=0;
          nprocs++;
        }
      }
    }
  }
  delete [] Pecount;
  CkPrintf("Built RhoRS Map [%d, %d] on %d processors\n",nchareRhoR,rhoRsubplanes, pesused ); 
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

RhoGSMapTable::RhoGSMapTable(MapType2  *_map, PeList *_availprocs, int _nchareRhoG,  bool useCentroid, MapType2 *rhorsmap, PeList *exclude): nchareRhoG(_nchareRhoG)
{
  reverseMap=NULL;
  maptable=_map;
  availprocs=_availprocs;
  int rgsobjs_per_pe, rem;

  if(availprocs->count()==0)
    availprocs->reset();

  //if(CkMyPe()==0) CkPrintf("nchareRhoG %d rgsobjs_per_pe %d rem
  //%d\n", nchareRhoG, rgsobjs_per_pe, rem); 
  if(useCentroid)
  {
    // get centroid of rsmap  use it to sort the avail list
    availprocs->sortSource(rhorsmap->getCentroid(config.torusMap));
    availprocs->reset();
  }
  PeList *avail= new PeList(*availprocs);
  *avail-*exclude;
  if(avail->count()>nchareRhoG)
  {
    // try an exclusion
    CkPrintf("RhoG excluding %d from avail %d\n",exclude->count(), availprocs->count());
    *availprocs-*exclude;
    availprocs->reindex();
    // CkPrintf("avail now %d\n", availprocs->count());
  }
  else
  {
    CkPrintf("cannot use exclusion in rhog\n");
    availprocs->reset();
  }
  delete avail;

  if(availprocs->count()==1)
  {
    rgsobjs_per_pe= nchareRhoG;
    rem=0;
  }
  else
  {
    rgsobjs_per_pe= nchareRhoG/(availprocs->count());
    rem = nchareRhoG % (availprocs->count());
    if(nchareRhoG<availprocs->count())
    {
      rem=0;
      rgsobjs_per_pe=1;
    }
    if(rem!=0)
      rgsobjs_per_pe += 1;
  }

  int destpe=availprocs->findNext();
  if(availprocs->count()==0)
    availprocs->reset();

  for(int chunk=0; chunk<nchareRhoG; chunk+=rgsobjs_per_pe)
  {
    if(rem>1)
      if(chunk==rem*rgsobjs_per_pe)
        rgsobjs_per_pe -= 1;  
    for(int i=chunk;((i<chunk+rgsobjs_per_pe)&&(i<nchareRhoG));i++)
    {
#ifdef USE_INT_MAP
      maptable->set(i, 0, destpe);
#else
      maptable->put(intdual(i, 0))=destpe;
#endif
    } 
    exclude->mergeOne(destpe);
    if(chunk+1<nchareRhoG)
      destpe=availprocs->findNext();
    if(availprocs->count()==0)
      availprocs->reset();
  }
#ifdef _MAP_DEBUG_
  CkPrintf("RhoGSMap created on processor %d\n", CkMyPe());
  dump();
#endif
}


RhoRHartMapTable::RhoRHartMapTable(MapType3  *_map, PeList *_availprocs, int _nchareRhoRHart, int rhoRsubplanes, int nchareHartAtmT,PeList *exclude ): nchareRhoRHart(_nchareRhoRHart)
{
  reverseMap=NULL;
  maptable=_map;
  availprocs=_availprocs;
  int rrsobjs_per_pe, rem;
  int srcpe=0;
  int numChares=nchareRhoRHart*rhoRsubplanes*nchareHartAtmT;
  if(availprocs->count()==0)
    availprocs->reset();
  PeList *avail= new PeList(*availprocs);
  *avail-*exclude;
  if(avail->count()>numChares)
  {
    // try an exclusion
    CkPrintf("RhoRHart excluding %d from avail %d\n",exclude->count(), availprocs->count());
    *availprocs-*exclude;
    availprocs->reindex();
    // CkPrintf("avail now %d\n", availprocs->count());
  }
  else
  {
    // CkPrintf("cannot use exclusion in rhoRhart\n");
    availprocs->reset();
  }
  delete avail;

  if(availprocs->count()==1)
  {
    rrsobjs_per_pe= numChares;
    rem=0;
  }
  else
  {
    rrsobjs_per_pe= numChares/(availprocs->count());
    rem = numChares % (availprocs->count());
    if(numChares<availprocs->count())
    {
      rrsobjs_per_pe=1;
      rem=0;
    }
    if(rem!=0)
      rrsobjs_per_pe += 1;
  }

  int destpe=availprocs->findNext(); 
  srcpe=destpe;
  if(availprocs->count()==0)
    availprocs->reset();

  //if(CkMyPe()==0) CkPrintf("nchareRhoR %d rrsobjs_per_pe %d rem %d\n", nchareRhoRHart, rrsobjs_per_pe, rem);   

  int nprocs=0, objs=0;
  destpe=availprocs->findNext();
  if(availprocs->count()==0)
    availprocs->reset();

  for(int atmtype=0; atmtype< nchareHartAtmT;atmtype++)
  {
    for(int chunk=0; chunk<nchareRhoRHart; chunk++)
    {
      for(int subplane=0; subplane<rhoRsubplanes; subplane++)
      {
        if(rem!=0)
          if(nprocs==rem)
            rrsobjs_per_pe -= 1;
#ifdef USE_INT_MAP
        maptable->set(chunk, subplane, atmtype,destpe);
#else
        maptable->put(inttriple(chunk, subplane, atmtype))=destpe;
#endif
        objs++;
        exclude->mergeOne(destpe);
        if(objs>=rrsobjs_per_pe)
        {

          destpe=availprocs->findNext();
          if(availprocs->count()==0)
            availprocs->reset();
          objs=0;
          nprocs++;
        }
      }
    }
  }
#ifdef _MAP_DEBUG_
  CkPrintf("RhoRHartMap created on processor %d\n", CkMyPe());
  dump();
#endif

}

RhoGHartMapTable::RhoGHartMapTable(MapType2  *_map, PeList *_availprocs, int _nchareRhoGHart, int nchareHartAtmT, int useCentroid, MapType3 *rhartmap, PeList *exclude): nchareRhoGHart(_nchareRhoGHart)
{
  int npes;
  reverseMap=NULL;
  maptable=_map;
  PeList *availprocs=_availprocs;
  PeList *avail= new PeList(*availprocs);
  avail->reset();
  exclude->reset();
  *avail-*exclude;
  int numchares=nchareRhoGHart*nchareHartAtmT;
  bool excluded=true;
  if(avail->count()>numchares)
  {
    // try an exclusion
    CkPrintf("RhoGHart excluding %d from avail %d\n",exclude->count(), availprocs->count());
    *availprocs-*exclude;
    availprocs->reindex();
  }
  else
  {
    excluded=false;
    CkPrintf("cannot use exclusion in rhoghart\n");
    availprocs->reset();      
  }
  delete avail;
  if(availprocs->count()==0)
    availprocs->reset();
  npes=availprocs->count();
  int rghobjs_per_pe, rem;

  if(availprocs->count()==1)
  {
    rghobjs_per_pe= numchares;
    rem=0;
  }
  else
  {
    rghobjs_per_pe= numchares/npes;
    rem = numchares % npes;
    if(numchares < npes)
    {
      rem=0;
      rghobjs_per_pe=1;
    }
    if(rem!=0)
      rghobjs_per_pe += 1;
  }
  if(useCentroid && excluded && rhartmap!=NULL)
  {
    // get centroid of rhartmap  use it to sort the avail list
    availprocs->reset();
    availprocs->sortSource(rhartmap->getCentroid(config.torusMap));
    availprocs->reset();
  }
  int destpe=availprocs->findNext();
  if(availprocs->count()==0)
    availprocs->reset();
  for(int atmtype=0; atmtype< nchareHartAtmT;atmtype++)
  {
    for(int chunk=0; chunk<nchareRhoGHart; chunk+=rghobjs_per_pe)
    {
      //      CkAssert(exclude->exists(destpe)<0);
      if(rem>1 && rghobjs_per_pe>1)
        if(chunk==rem*rghobjs_per_pe)
          rghobjs_per_pe -= 1; 
      for(int i=chunk;((i<chunk+rghobjs_per_pe)&&(i<nchareRhoGHart));i++)
      {
#ifdef USE_INT_MAP
        maptable->set(i, atmtype, destpe);
#else
        maptable->put(intdual(i, atmtype))=destpe;
#endif
      } 
      exclude->mergeOne(destpe);
      //	  availprocs->sortSource(srcpe);
      if(chunk+1<nchareRhoGHart)
        destpe=availprocs->findNext();
      else
        chunk=nchareRhoGHart;  //get us out of here
      if(availprocs->count()==0)
        availprocs->reset();

    }
  }
#ifdef _MAP_DEBUG_
  CkPrintf("RhoGHartMap created on processor %d\n", CkMyPe());
  dump();
#endif

}

/* We receive a copy of the RhoRS and then act independently with the
   bulk of our communication with VdWGS. Exclusion mapping will try to
   put us on different processors from RhoRS. Our connection to RS is
   more tenuous, so centroid mapping to RS is not as compelling. */

VdWRSMapTable::VdWRSMapTable(MapType3  *_map, PeList *_availprocs, int _nchareRhoR, int _rhoRsubplanes, int _nchareVdW, int max_states, PeList *exclude): nchareRhoR(_nchareRhoR), rhoRsubplanes(_rhoRsubplanes)
                                                                                                                                                          ,nchareVdW(_nchareVdW) {
                                                                                                                                                            reverseMap=NULL;
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

                                                                                                                                                            if(availprocs->count()==0)
                                                                                                                                                              availprocs->reset();
                                                                                                                                                            int destpe;
                                                                                                                                                            int nprocs=0, objs=0;
                                                                                                                                                            destpe=availprocs->findNext();
                                                                                                                                                            if(availprocs->count()==0)
                                                                                                                                                              availprocs->reset();
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
                                                                                                                                                                    if(availprocs->count()==0)
                                                                                                                                                                      availprocs->reset();
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
  reverseMap=NULL;
  maptable=_map;
  availprocs=_availprocs;
  int rgsobjs_per_pe, rem;

  if(availprocs->count()==0)
    availprocs->reset();
  PeList *avail= new PeList(*availprocs);
  *avail-*exclude;
  int numVdWObjs=nchareRhoG*nchareVdW;
  if(avail->count()>numVdWObjs)
  {
    // try an exclusion
    CkPrintf("RhoG excluding %d from avail %d\n",exclude->count(), availprocs->count());
    *availprocs-*exclude;
    availprocs->reindex();
    // CkPrintf("avail now %d\n", availprocs->count());
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
  if(availprocs->count()==0)
    availprocs->reset();

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
        exclude->mergeOne(destpe);
        destpe=availprocs->findNext();
        if(availprocs->count()==0)
          availprocs->reset();
      }
    }
  }
#ifdef _MAP_DEBUG_
  CkPrintf("VdWGSMap created on processor %d\n", CkMyPe());
  dump();
#endif
}



void MapTable2::makeReverseMap()
{
#ifndef USE_INT_MAP
  CkHashtableIterator *it=maptable->iterator();
  it->seekStart();
  intdual *key;
  reverseMap= new CkVec <intdual> [config.numPesPerInstance];
  while(it->hasNext())
  {
    it->next((void **) &key);
    int proc =maptable->get(key[0]);
    reverseMap[proc].push_back(key[0]);
  }
  delete it;
#endif
}



PeList *subListPlane(int plane, int nstates, MapType2 *smap)
{

  PeList *thisPlane = new PeList(1);
  thisPlane->TheList[0]=smap->get(0,plane);
  for(int state=1; state<nstates; state++)
  {
    thisPlane->mergeOne(smap->get(state,plane));
  }
  thisPlane->reset();
  thisPlane->reindex();
  return(thisPlane);
}

PeList *subListState(int state, int nplanes, MapType2 *smap)
{
  //      CkPrintf("in sublist state %d\n",state);
  PeList *thisState= new PeList(1);
  thisState->TheList[0]=smap->get(state,0);
  thisState->sortIdx[0]=0;
  for(int plane=1; plane < nplanes; plane++)
  {
    int pe=smap->get(state, plane);
    thisState->mergeOne(pe);
  }
  thisState->reindex();
  thisState->reset();
  return(thisState);
}

PeList *subListState2(int state1, int state2, int nplanes, int numChunks, MapType4 *smap)
{
  PeList *thisState = new PeList(1);
  thisState->TheList[0]=smap->get(0, state1, state2, 0);
  thisState->sortIdx[0]=0;      
  for(int plane=0; plane<nplanes; plane++)
    for(int chunk=0; chunk<numChunks; chunk++)
    {
      int pe = smap->get(plane, state1, state2, chunk);
      thisState->mergeOne(pe);
    }
  thisState->reindex();
  thisState->reset();
  return(thisState);
}



void RhoRSMapTable::sortByCentroid(PeList *avail, int plane, int nstates, MapType2 *rsmap)
{
  int points=0, bestPe;
  if(config.torusMap==1) {
    int sumX=0, sumY=0, sumZ=0;
    for(int state=0;state<nstates;state++)
    {
      int X, Y, Z, T;
      topoMgr->rankToCoordinates(rsmap->get(state, plane), X, Y, Z, T);
      sumX+=X;
      sumY+=Y;
      sumZ+=Z;
      points++;
    }
    int avgX=sumX/points;
    int avgY=sumY/points;
    int avgZ=sumZ/points;
    bestPe=topoMgr->coordinatesToRank(avgX, avgY, avgZ, 0);
    avail->sortSource(bestPe);
    avail->reset();
  }
  else {
    int sumPe=0;
    for(int state=0;state<nstates;state++)
    {
      sumPe+=rsmap->get(state,plane);
      points++;
    }
    bestPe=sumPe/points;
    avail->sortSource(bestPe);
    avail->reset();
  }
}


void SCalcMapTable::sortByCentroid(PeList *avail, int plane, int stateX, int stateY, int grainsize, MapType2 *gsmap)
{
  int points=0, bestPe;

  if(config.torusMap==1) {
    int sumX=0, sumY=0, sumZ=0;
    for(int state=stateX;(state<stateX+grainsize)&&(state<config.nstates);state++)
    {
      int X, Y, Z, T;
      topoMgr->rankToCoordinates(gsmap->get(state, plane), X, Y, Z, T);
      sumX+=X;
      sumY+=Y;
      sumZ+=Z;
      points++;
    }
    for(int state=stateY;(state<stateY+grainsize)&&(state<config.nstates);state++)
    {
      int X, Y, Z, T;
      topoMgr->rankToCoordinates(gsmap->get(state, plane), X, Y, Z, T);
      sumX+=X;
      sumY+=Y;
      sumZ+=Z;
      points++;
    }
    int avgX=sumX/points;
    int avgY=sumY/points;
    int avgZ=sumZ/points;
    bestPe=topoMgr->coordinatesToRank(avgX, avgY, avgZ, 0);
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
  }
  avail->sortSource(bestPe);
  avail->reset();

}


void OrthoMapTable::sortByCentroid(PeList *avail, int nplanes, int state1, int state2, int numChunks, MapType4 *smap)
{
  int points=0, bestPe;
  if(config.torusMap==1) {
    int sumX = 0, sumY = 0, sumZ = 0;

    for(int plane=0; plane<nplanes; plane++)
      for(int chunk=0; chunk<numChunks; chunk++)
      {  
        int X, Y, Z, T;
        topoMgr->rankToCoordinates(smap->get(plane, state1, state2, chunk), X, Y, Z, T);
        sumX += X;
        sumY += Y;
        sumZ += Z;
        points++;
      }
    int avgX = sumX/points;
    int avgY = sumY/points;
    int avgZ = sumZ/points;
    bestPe = topoMgr->coordinatesToRank(avgX, avgY, avgZ, 0);
    avail->sortSource(bestPe);
    avail->reset();
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
    avail->sortSource(bestPe);
    avail->reset();
  }
}

int OrthoMapTable::minDistCentroid(PeList *avail, int nplanes, int state1, int state2, int numChunks, MapType4 *smap)
{
  int points=0, bestPe;
  if(config.torusMap==1) {
    int sumX = 0, sumY = 0, sumZ = 0;

    for(int plane=0; plane<nplanes; plane++)
      for(int chunk=0; chunk<numChunks; chunk++)
      {    
        int X, Y, Z, T;
        topoMgr->rankToCoordinates(smap->get(plane, state1, state2, chunk), X, Y, Z, T);
        sumX += X;
        sumY += Y;
        sumZ += Z;
        points++;
      }
    int avgX = sumX/points;
    int avgY = sumY/points;
    int avgZ = sumZ/points;
    bestPe = topoMgr->coordinatesToRank(avgX, avgY, avgZ, 0);
    return(avail->minDist(bestPe));
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

