#include "charm++.h"
#include "PeList.h"
#include "MapFile.h"

#define USE_INT_MAP
#ifdef USE_INT_MAP
#include "IntMap.h"
//#define USE_INT_2on1
#ifdef USE_INT_2on1
typedef IntMap2on1 MapType2;
#else
//class MapType2;
#endif
typedef IntMap4 MapType4;
typedef IntMap3 MapType3;
#endif
#include "MapTable.h"
#ifdef CMK_VERSION_BLUEGENE
    extern 	BGLTorusManager *bgltm;
#endif      

int MapType2::getCentroid(){
#ifdef CMK_VERSION_BLUEGENE
  int points=0;
  int  sumX=0;
  int sumY=0;
  int sumZ=0;
  int X=0,Y=0,Z=0;
    for(int i=0;i<getXmax();i++)
      for(int j=0;j<getYmax();j++){
	  bgltm->getCoordinatesByRank(get(i,j),X, Y, Z);
	  sumX+=X;
	  sumY+=Y;
	  sumZ+=Z;
	  points++;
	}
  int avgX=sumX/points;
  int avgY=sumY/points;
  int avgZ=sumZ/points;
  int bestPe=bgltm->coords2rank(avgX, avgY, avgZ);
#else
  int points=0, sum=0;
  for(int i=0;i<getXmax();i++)
    for(int j=0;j<getYmax();j++){
	sum+=get(i,j);
	points++;
      }
  int bestPe=sum/points;
#endif      
  return(bestPe);
  }

int MapType3::getCentroid(){
#ifdef CMK_VERSION_BLUEGENE
  int points=0;
  int  sumX=0;
  int sumY=0;
  int sumZ=0;
  int X=0,Y=0,Z=0;
    for(int i=0;i<getXmax();i++)
      for(int j=0;j<getYmax();j++){
	for(int k=0;k<getZmax();k++){
	  bgltm->getCoordinatesByRank(get(i,j,k),X, Y, Z);
	  sumX+=X;
	  sumY+=Y;
	  sumZ+=Z;
	  points++;
	}
    }
  int avgX=sumX/points;
  int avgY=sumY/points;
  int avgZ=sumZ/points;
  int bestPe=bgltm->coords2rank(avgX, avgY, avgZ);
#else
  int points=0, sum=0;
  for(int i=0;i<getXmax();i++)
    for(int j=0;j<getYmax();j++){
      for(int k=0;k<getZmax();k++){
	sum+=get(i,j,k);
	points++;
      }
    }
  int bestPe=sum/points;
#endif      
  return(bestPe);
  }

GSMapTable::GSMapTable(MapType2  *_map, PeList *_availprocs, 
		       int _nchareG, int _nstates,  
		       int _Gstates_per_pe, bool useCuboidMap)    : 
      nchareG(_nchareG), 
    nstates(_nstates), Gstates_per_pe(_Gstates_per_pe)
{ 
  reverseMap=NULL;
  maptable=_map;
  availprocs=_availprocs;
  state_load = 0.0;
  int l, m, pl, pm, srem, rem, i=0;

  l=Gstates_per_pe;		// no of states in one chunk
  pl = nstates / l;		// no of procs on y axis
  if(nstates % l == 0)
    srem = 0;
  else
  {
    while(pow(2.0, (double)i) < pl)
      i++;
    pl = (int) pow(2.0, (double)(i-1));             // make it same as the nearest smaller power of 2
    srem = nstates % pl;
  }
  pm = availprocs->count() / pl;		// no of procs on x axis
        
  if(pm==0)
  CkAbort("Choose a larger Gstates_per_pe\n");
        
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
  int destpe=availprocs->findNext();

  if(useCuboidMap)
  {
    for(int plane=0;plane<nchareG;plane++)
      for(int state=0;state<nstates;state+=Gstates_per_pe)
      {
	for(int stateperpe=0;stateperpe<Gstates_per_pe;stateperpe++)
	{
#ifdef USE_INT_MAP
	  maptable->set(state+stateperpe, plane,destpe);
#else
	  maptable->put(intdual(state+stateperpe, plane))=destpe;
#endif
	}
	if(availprocs->count()==0)
	  availprocs->reset();
	destpe=availprocs->findNext();
	/*if(availprocs->count()==0)
	{
	  CkPrintf("GSMap created on processor %d\n", CkMyPe());
	  dump();
	  CkPrintf("GSMap ran out of nodes on plane %d state %d\n", plane, state);
	  CkExit();
	}
	*/
      }
  }
  else
  {
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
  MapFile *mf = new MapFile("GSMap", 2, size, CkNumPes(), "TXYZ", 2, 1, 1, 1);
  mf->dumpMap(maptable);
#endif
}

SCalcMapTable::SCalcMapTable(MapType4  *_map, PeList *_availprocs, 
			     int _nstates, int _nchareG,  int _grainsize, 
			     CmiBool _flag, int _scalc_per_plane,  
			     int _planes_per_pe, 
			     int _numChunksA, 
			     int _numChunksS, 
			     MapType2  *gsmap, bool useCuboidMap, bool useCentroid, int boxSize) : 
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

  if(planes_per_pe==0)
    CkAbort("Choose a smaller Gstates_per_pe\n");
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
		  //		  CkPrintf("plane %d making pelist from plane's pes\n",plane);
		  thisPlaneBox= subListPlane(plane, max_states, gsmap);
		}
	      else
		{
		  //		  CkPrintf("plane %d making pelist from boxSize %d\n",plane,boxSize);
		  thisPlaneBox= new PeList(availprocs, plane*boxSize, boxSize);
		}

	      //	      PeList *thisPlaneBox= availprocs;
	      //	      thisPlaneBox->dump();
	      if(!useCentroid)
		destpe=thisPlaneBox->findNext();
	      for(int xchunk=0; xchunk<max_states; xchunk=xchunk+grainsize)
		for(int ychunk=xchunk; ychunk<max_states; ychunk=ychunk+grainsize)
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
			    if(rem!=0)
			      if(procno==rem)
				scobjs_per_pe-=1;
			    count=0;

			  }
			//			CkPrintf("%d %d %d %d mapped to %d\n",plane,xchunk,ychunk,newdim, destpe);
#ifdef USE_INT_MAP
			maptable->set(plane, xchunk, ychunk, newdim,destpe);
#else
			CkArrayIndex4D idx4d(plane, xchunk, ychunk, newdim);
			memcpy(intidx,idx4d.index,2*sizeof(int));
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
	      for(int xchunk=0; xchunk<max_states; xchunk=xchunk+grainsize)
		for(int ychunk=xchunk; ychunk<max_states; ychunk=ychunk+grainsize)
		  //for(int newdim=0; newdim<numChunksSym; newdim++)
		  for(int plane=pchunk; plane<pchunk+planes_per_pe && plane<nchareG; plane++)
		    {
		      CkArrayIndex4D idx4d(plane, xchunk, ychunk, newdim);
		      memcpy(intidx,idx4d.index,2*sizeof(int));  // our 4 shorts are now 2 ints

		      if(count<scobjs_per_pe)
			{
			  //if(CkMyPe()==0) CkPrintf("plane %d x %d y %d newdim %d = proc %d\n", plane, xchunk, ychunk, newdim, fp.next[0]*x*y+fp.next[1]*x+fp.next[2]);
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
		  //		  CkPrintf("plane %d making pelist from plane's pes\n",plane);
		  thisPlaneBox= subListPlane(plane, max_states, gsmap);
		}
	      else
		{
		  //		  CkPrintf("plane %d making pelist from boxSize %d\n",plane,boxSize);
		  thisPlaneBox= new PeList(availprocs, plane*boxSize, boxSize);
		}

	      //	      PeList *thisPlaneBox= availprocs;
	      //	      thisPlaneBox->dump();
	      if(!useCentroid)
		destpe=thisPlaneBox->findNext();
	      for(int xchunk=0; xchunk<max_states; xchunk=xchunk+grainsize)
		for(int ychunk=0; ychunk<max_states; ychunk=ychunk+grainsize)
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
		    for(int newdim=0; newdim<numChunksAsym; newdim++)
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
			    if(rem!=0)
			      if(procno==rem)
				scobjs_per_pe-=1;
			    count=0;

			  }
			//			CkPrintf("%d %d %d %d mapped to %d\n",plane,xchunk,ychunk,newdim, destpe);
#ifdef USE_INT_MAP
			maptable->set(plane, xchunk, ychunk, newdim,destpe);
#else
			CkArrayIndex4D idx4d(plane, xchunk, ychunk, newdim);
			memcpy(intidx,idx4d.index,2*sizeof(int));
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
	      for(int xchunk=0; xchunk<max_states; xchunk=xchunk+grainsize)
		for(int ychunk=0; ychunk<max_states; ychunk=ychunk+grainsize)
		  for(int plane=pchunk; plane<pchunk+planes_per_pe && plane<nchareG; plane++)
		    {
		      CkArrayIndex4D idx4d(plane, xchunk, ychunk, newdim);
		      memcpy(intidx,idx4d.index,2*sizeof(int));  // our 4 shorts are now 2 ints

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

RSMapTable::RSMapTable(MapType2  *_map, PeList *_availprocs,
	int _nstates, int _sizeZ, int _Rstates_per_pe, bool useCuboidMap,
	MapType2 *gsmap, int nchareG) :
   nstates(_nstates), sizeZ(_sizeZ),
  Rstates_per_pe(_Rstates_per_pe)
{
  int l, m, pl, pm, srem, rem, i=0, rsobjs_per_pe;
  reverseMap=NULL;
  maptable=_map;
  availprocs=_availprocs;
  int *Pecount= new int [CkNumPes()];

  bzero(Pecount, CkNumPes() *sizeof(int));

  rsobjs_per_pe = nstates*sizeZ/availprocs->count();
  l=Rstates_per_pe;		// no of states in one chunk
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
        
  if(pm==0)
    CkAbort("Choose a larger Rstates_per_pe\n");

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
      // place all the states box by box
      PeList *exclusionList = NULL;
      for(int state=0; state < nstates; state++)
	{
	  PeList *thisStateBox, *single;
	  thisStateBox = subListState(state, nchareG, gsmap);
	  bool useExclude=true;
	  if(exclusionList!=NULL)
	    {
	      *thisStateBox - *exclusionList;
	      thisStateBox->reindex();
	    }
	  if(thisStateBox->count()==0)
	    {
	      useExclude=false;
	      CkPrintf("RS state %d excluded to nil, ignoring exclusions\n",state);
	      delete thisStateBox;
	      thisStateBox = subListState(state, nchareG, gsmap);
	    }
	  
	  // sort by centroid (not necessary)
	  if(thisStateBox->count()==0)
	    thisStateBox->reset();
	  destpe = thisStateBox->findNext();

	  for(int plane=0; plane < sizeZ; plane++)
	    {
	 	
#ifdef USE_INT_MAP
	      maptable->set(state, plane, destpe);
#else						
	      maptable->put(intdual(state, plane))= destpe;
#endif
	      //		if(CkMyPe()==0) CkPrintf("%d %d [%d]\n", state, plane, destpe);
	      Pecount[destpe]++;
	      if(Pecount[destpe]>=rsobjs_per_pe+1)
		{
		  if(exclusionList==NULL)
		    {
		      exclusionList=new PeList(1);
		      exclusionList->TheList[0]=destpe;
		    }
		  else
		    exclusionList->mergeOne(destpe);
		  if(useExclude && thisStateBox->size>1)
		    {
		      *thisStateBox - *exclusionList;
		      thisStateBox->reindex();
		    }
		}
	      if(useExclude && Pecount[destpe]>rsobjs_per_pe+1)
		{
		  CkPrintf("[%d] %d\n", destpe, Pecount[destpe]);
		  exclusionList->dump();
		  thisStateBox->dump();
		  CkAbort("RS has too many chares\n");
		}
	      if(thisStateBox->count()==0)
		thisStateBox->reset();
	      destpe = thisStateBox->findNext();
	    }
	  delete thisStateBox;
	}
      if(exclusionList!=NULL)
	delete exclusionList;

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
}


RPPMapTable::RPPMapTable(MapType2  *_map, 
			 PeList *_availprocs, PeList *exclusion,
			 int _nstates, int _sizeZNL, int _Rstates_per_pe,
			 int boxSize, bool usePPmap, int nchareG,
			 MapType2 *pp_map) :
  nstates(_nstates), sizeZNL(_sizeZNL),
  Rstates_per_pe(_Rstates_per_pe)
{
  int states_per_pe;
  int totalChares=nstates*sizeZNL;
  int m, pl, pm, srem, rem, i=0;
  reverseMap=NULL;
  maptable=_map;
  availprocs=_availprocs;
  //        if(useCuboidMap)
  //	  states_per_pe=nstates/boxSize;		// no of states in one chunk
  //	else
  //	  states_per_pe=nstates/boxSize;		// no of states in one chunk
  //	else
  PeList *RPPlist=availprocs;
  states_per_pe=Rstates_per_pe;		// no of states in one chunk
  pl = nstates / states_per_pe;
  int afterExclusion=availprocs->count() - exclusion->count();
  //  CkPrintf("RPP excluded list\n");
  //  exclusion->dump();
  if(afterExclusion > pl*sizeZNL)
    { // we can fit the exclusion without blinking
      CkPrintf("RPP using density exclusion to avoid %d processors\n",exclusion->count());
      RPPlist-exclusion;
    }
  else
    {// so an rstates_per_pe chosen for realstate might be too big
      if(afterExclusion > (sizeZNL+nstates)*2)
	{ // set states_per_pe to fit the exclusion

	  states_per_pe=totalChares/afterExclusion/2;
	  CkPrintf("RPP adjusting states per pe from %d to %d to use density exclusion to stay within %d processors\n",Rstates_per_pe, states_per_pe, exclusion->count());
	  RPPlist-exclusion;
	}
      else
	{
	  CkPrintf("RPP with %d chares ignoring density exclusion which left %d out of %d processors\n",totalChares, exclusion->count(), availprocs->count());
	}
    }
  if(nstates % states_per_pe == 0)
    srem = 0;
  else
    {
      while(pow(2.0, (double)i) < pl)
	i++;
      pl = (int) pow(2.0, (double)(i-1));		// make it same as the nearest smaller power of 2
      srem = nstates % pl;
    }
  pm = RPPlist->count() / pl;
        
  if(pm==0)
    CkAbort("Choose a larger Rstates_per_pe\n");

  m = sizeZNL / pm;
  rem = sizeZNL % pm;

  // CkPrintf("nstates %d sizeZNL %d Pes %d\n", nstates, sizeZNL, RPPlist->count());	
  // CkPrintf("l %d, m %d pl %d pm %d srem %d rem %d\n", l, m,
  // pl, pm, srem, rem);
  int srcpe=0;
  // CkPrintf("nstates %d sizeZNL %d Pes %d\n", nstates, sizeZNL, RPPlist->count());	
  // RPPlist->dump();
  if(usePPmap)
    {
      /*  
       * RPP(s,*) <-> PP(s,*)  
       */

      // make all the maps
      // subtract the exclusion from each 
      // keep the smallest count as the max charesperpe
      PeList **maps= new PeList* [nstates];
      int maxcharesperpe=0;
      int *usedPes= new int[CkNumPes()];
      bzero(usedPes, CkNumPes() * sizeof(int));
      for(int state=0; state < nstates ; state++)
	{
	  // have variable number of exclusions per list
	  // need to figure out what the max should be
	  // for the border cases
	  maps[state]= subListState( state, nchareG, pp_map);
	  //	  *maps[state]-*exclusion;
	  //	  maps[state]->reindex();
	  int thischaresperpe=sizeZNL/maps[state]->count() + 1;
	  //	  CkPrintf("state %d has %d charesperpe in %d pemap \n",state,thischaresperpe,maps[state]->count());
	  //	  maxcharesperpe=(thischaresperpe>maxcharesperpe) ? thischaresperpe : maxcharesperpe;
	  maxcharesperpe=sizeZNL*nstates/CkNumPes();
	}
      PeList *usedbyRPP=NULL;
      //      CkPrintf("maxcharesperpe is %d\n",maxcharesperpe);
      for(int state=0; state < nstates ; state++)
	{
	  if(usedbyRPP!=NULL)
	    {
	      *(maps[state])-*usedbyRPP;
	      maps[state]->reindex();
	    }
	  for(int plane=0; plane < sizeZNL; plane++)
	    {
	      if(maps[state]->count()==0)
		maps[state]->reset();
	      int destpe=maps[state]->findNext(); 
#ifdef USE_INT_MAP
	      maptable->set(state, plane, destpe);
#else
	      maptable->put(intdual(state, plane))=destpe;
#endif
	      usedPes[destpe]++;
	      if(usedPes[destpe]>maxcharesperpe)
		{
		  if(usedbyRPP==NULL)
		    { 
		      usedbyRPP= new PeList(1);
		      usedbyRPP->TheList[0]=destpe;
		    }
		  else
		    usedbyRPP->mergeOne(destpe);
		  exclusion->mergeOne(destpe);
		  exclusion->reindex();
		  //		  maps[state]->dump();
		  *(maps[state]) - *usedbyRPP;
		  //		  maps[state]->dump();
		  maps[state]->reindex();
		  maps[state]->reset();
		  //		  maps[state]->dump();
		  //		  CkPrintf("removed %d\n",destpe);
		  
		}
	    }
	}
      delete usedbyRPP;
      delete [] usedPes;
      for(int state=0;state<nstates;state++)
	{
	  delete maps[state];
	  maps[state]=NULL;
	}
      delete [] maps;
    }
  else
    {
      int destpe=RPPlist->findNext();
      for(int ychunk=0; ychunk<sizeZNL; ychunk=ychunk+m)
	{
	  if(ychunk==(pm-rem)*m)
	    m=m+1;
	  for(int xchunk=0; xchunk<nstates; xchunk=xchunk+states_per_pe)
	    {
	      if(xchunk==(pl-srem)*states_per_pe)
		states_per_pe=states_per_pe+1;
	      for(int state=xchunk; state<xchunk+states_per_pe && state<nstates; state++)
		{
		  for(int plane=ychunk; plane<ychunk+m && plane<sizeZNL; plane++)
		    {
		      //		  CkPrintf("RPP setting %d %d to pe %d\n",state,plane,destpe);		    
#ifdef USE_INT_MAP
		      maptable->set(state, plane, destpe);
#else
		      maptable->put(intdual(state, plane))=destpe;
#endif

		    }
		}
	      if(RPPlist->count()==0)
		RPPlist->reset();
	      destpe=RPPlist->findNext();
	    }
	}
    }
#ifdef _MAP_DEBUG_
  CkPrintf("RPPMap created on processor %d\n", CkMyPe());
  dump();
#endif
}

OrthoMapTable::OrthoMapTable(MapType2 *_map, PeList *_availprocs, int _nstates, int _orthograinsize, MapType4 *scalcmap, int nplanes, int numChunks, int sGrainSize): nstates(_nstates), orthoGrainSize(_orthograinsize)
{
  maptable = _map;
  availprocs = _availprocs;
  int oobjs_per_pe;
  int srcpe = 0, destpe;
  oobjs_per_pe = (nstates/orthoGrainSize)*(nstates/orthoGrainSize)/(availprocs->count()) + 1;
  int *Pecount= new int [CkNumPes()];
  bzero(Pecount, CkNumPes()*sizeof(int)); 
  int s1 = 0, s2 = 0;
  PeList *exclusionList = NULL;
  for(int state1 = 0; state1 < nstates; state1 += orthoGrainSize)
    for(int state2 = 0; state2 < nstates; state2 += orthoGrainSize)
    {
      s1 = (state1/sGrainSize);
      s1 = s1 * sGrainSize;
      s2 = (state2/sGrainSize);
      s2 = s2 * sGrainSize;
      PeList *thisStateBox = subListState2(s1, s2, nplanes, numChunks, scalcmap);
      bool useExclude = true;
      if(exclusionList != NULL) 
      {
	*thisStateBox - *exclusionList;
	thisStateBox->reindex();
      }

      if(thisStateBox->count() == 0)
      {
	CkPrintf("Ortho %d %d ignoring exclusion\n", state1, state2);
	delete thisStateBox;
	thisStateBox = subListState2(s1, s2, nplanes, numChunks, scalcmap);
	useExclude = false;
      }
	
      sortByCentroid(thisStateBox, nplanes, s1, s2, numChunks, scalcmap);
	
      destpe=thisStateBox->findNext();
      if(thisStateBox->count()==0)
	thisStateBox->reset();

#ifdef USE_INT_MAP
      maptable->set(state1/orthoGrainSize, state2/orthoGrainSize, destpe);
#else
      maptable->put(intdual(state1/orthoGrainSize, state2/orthoGrainSize))=destpe;
#endif
      Pecount[destpe]++;	
      if(Pecount[destpe]>=oobjs_per_pe)
      {
	if(exclusionList==NULL)
	{
	  exclusionList=new PeList(1);
	  exclusionList->TheList[0]=destpe;
        }
        else
	  exclusionList->mergeOne(destpe);
	if(useExclude && thisStateBox->size>1)
	{
	  *thisStateBox - *exclusionList;
	  thisStateBox->reindex();
	}
	sortByCentroid(thisStateBox, nplanes, s1, s2, numChunks, scalcmap);
      }
      destpe=thisStateBox->findNext();
      if(thisStateBox->count()==0)
         thisStateBox->reset();
      delete thisStateBox;
    }
  delete [] Pecount;
}

OrthoHelperMapTable::OrthoHelperMapTable(MapType2 *_map, int _nstates, int _orthograinsize, MapType2 *omap): nstates(_nstates), orthoGrainSize(_orthograinsize)
{
  maptable = _map;
  int destpe = 0;

  for(int state1 = 0; state1 < nstates/orthoGrainSize; state1++)
    for(int state2 = 0; state2 < nstates/orthoGrainSize; state2++)
    {
#ifdef USE_INT_MAP
      destpe = omap->get(state1, state2);
      if(destpe	== 0)
        maptable->set(state1, state2, (destpe+1)%CkNumPes());
      else
        maptable->set(state1, state2, destpe-1);
#else
      destpe = omap->get(intdual(state1, state2));
      if(destpe	== 0)
        maptable->put(intdual(state1, state2))=(destpe+1)%CkNumPes();
      else
        maptable->put(intdual(state1, state2))=destpe-1;
#endif
    }
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
      if(rem!=0)
	rrsobjs_per_pe += 1;
    }

  if(availprocs->count()==0)
    availprocs->reset();
  int destpe;
  int *Pecount= new int [CkNumPes()];
  bzero(Pecount, CkNumPes()*sizeof(int)); 

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
	    CkPrintf("Rho RS %d ignoring exclusion\n",chunk);
	    delete thisPlaneBox;
	    thisPlaneBox = subListPlane(chunk, max_states, rsmap);
	    useExclude=false;
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
      for(int i=0; i<CkNumPes();i++)
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

  if(availprocs->count()==1)
    {
      rgsobjs_per_pe= nchareRhoG;
      rem=0;
    }
  else
    {
      rgsobjs_per_pe= nchareRhoG/(availprocs->count());
      rem = nchareRhoG % (availprocs->count());
      if(rem!=0)
	rgsobjs_per_pe += 1;
    }

  if(availprocs->count()==0)
    availprocs->reset();

  //if(CkMyPe()==0) CkPrintf("nchareRhoG %d rgsobjs_per_pe %d rem
  //%d\n", nchareRhoG, rgsobjs_per_pe, rem); 
  if(useCentroid)
    {
      // get centroid of rsmap  use it to sort the avail list
      availprocs->sortSource(rhorsmap->getCentroid());
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
      CkPrintf("avail now %d\n", availprocs->count());
    }
  else
    {
      CkPrintf("cannot use exclusion in rhog\n");
    }
  delete avail;

  int destpe=availprocs->findNext();
  if(availprocs->count()==0)
    availprocs->reset();

  for(int chunk=0; chunk<nchareRhoG; chunk+=rgsobjs_per_pe)
    {
      if(rem!=0)
	if(chunk==rem*rgsobjs_per_pe)
	  rgsobjs_per_pe -= 1;  
      for(int i=chunk;i<chunk+rgsobjs_per_pe;i++)
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

  if(availprocs->count()==1)
    {
      rrsobjs_per_pe= numChares;
      rem=0;
    }
  else
    {
      rrsobjs_per_pe= numChares/(availprocs->count());
      rem = numChares % (availprocs->count());
      if(rem!=0)
	rrsobjs_per_pe += 1;
    }
  PeList *avail= new PeList(*availprocs);
  *avail-*exclude;
  if(avail->count()>numChares)
    {
      // try an exclusion
      CkPrintf("RhoRHart excluding %d from avail %d\n",exclude->count(), availprocs->count());
      *availprocs-*exclude;
      availprocs->reindex();
      CkPrintf("avail now %d\n", availprocs->count());
    }
  else
    {
      CkPrintf("cannot use exclusion in rhoRhart\n");
    }
  delete avail;

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
  availprocs=_availprocs;
  npes=availprocs->count();

  int rghobjs_per_pe, rem;
  int numchares=nchareRhoGHart*nchareHartAtmT;
  if(availprocs->count()==1)
    {
      rghobjs_per_pe= numchares;
      rem=0;
    }
  else
    {
      rghobjs_per_pe= numchares/npes;
      rem = numchares % npes;
      if(rem!=0)
	rghobjs_per_pe += 1;
    }
  PeList *avail= new PeList(*availprocs);
  *avail-*exclude;
  if(avail->count()>numchares)
    {
      // try an exclusion
      CkPrintf("RhoGHart excluding %d from avail %d\n",exclude->count(), availprocs->count());
      *availprocs-*exclude;
      availprocs->reindex();
      CkPrintf("avail now %d\n", availprocs->count());
    }
  else
    {
      CkPrintf("cannot use exclusion in rhoghart\n");
    }
  delete avail;
    

  if(availprocs->count()==0)
    availprocs->reset();
  if(useCentroid)
    {
      // get centroid of rsmap  use it to sort the avail list
      availprocs->sortSource(rhartmap->getCentroid());
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
	  if(rem!=0)
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
	  if(chunk+1<numchares)
	    destpe=availprocs->findNext();
	  if(availprocs->count()==0)
	    availprocs->reset();

	}
    }
#ifdef _MAP_DEBUG_
  CkPrintf("RhoGHartMap created on processor %d\n", CkMyPe());
  dump();
#endif

}


void MapTable2::makeReverseMap()
{
#ifndef USE_INT_MAP
    CkHashtableIterator *it=maptable->iterator();
    it->seekStart();
    intdual *key;
    reverseMap= new CkVec <intdual> [CkNumPes()];
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

      PeList *thisPlane = new PeList(nstates);
      int count = 0;
      for(int state=0; state<nstates; state++)
      {
	bool newPe = true;
	int pe = smap->get(state,plane);
	for(int i=0; i<count; i++)
	{
	    if(thisPlane->TheList[i] == pe)
	    newPe = false;
	}
	if(newPe)
	{
	    thisPlane->sortIdx[count] = count;
	    thisPlane->TheList[count] = pe;
	    count++;
	}
      }
      thisPlane->size = count;
      thisPlane->current = 0;
      return(thisPlane);
}

PeList *subListState(int state, int nplanes, MapType2 *smap)
{
  //      CkPrintf("in sublist state %d\n",state);
      PeList *thisState= new PeList(nplanes);
      int count=0;
      for(int plane=0; plane < nplanes; plane++)
	{
	  bool newPe=true;
	  int pe=smap->get(state, plane);
	  for(int i=0; i < count; i++)
	    {
	      if(thisState->TheList[i] == pe)
		newPe=false;
	    }
	  if(newPe)
	    {
	      thisState->sortIdx[count]=count;
	      thisState->TheList[count]=pe;
	      count++;
	    }
	}
      thisState->size=count;
      thisState->current=0;
      return(thisState);
}

PeList *subListState2(int state1, int state2, int nplanes, int numChunks, MapType4 *smap)
{
      PeList *thisState = new PeList(nplanes*numChunks);
      int count = 0;
      for(int plane=0; plane<nplanes; plane++)
        for(int chunk=0; chunk<numChunks; chunk++)
	{
	  bool newPe = true;
	  int pe = smap->get(plane, state1, state2, chunk);
	  for(int i=0; i<count; i++)
	    {
	      if(thisState->TheList[i] == pe)
		newPe = false;
	    }
	  if(newPe)
	    {
	      thisState->sortIdx[count] = count;
	      thisState->TheList[count] = pe;
	      count++;
	    }
	}
      thisState->size = count;
      thisState->current = 0;
      return(thisState);
}

#ifdef CMK_VERSION_BLUEGENE


void RhoRSMapTable::sortByCentroid(PeList *avail, int plane, int nstates, MapType2 *rsmap)
{
  int sumX=0, sumY=0, sumZ=0;
  int points=0;
      
  for(int state=0;state<nstates;state++)
    {
      int X, Y, Z;
      bgltm->getCoordinatesByRank(rsmap->get(state,plane),X, Y, Z);
      sumX+=X;
      sumY+=Y;
      sumZ+=Z;
      points++;
    }
  int avgX=sumX/points;
  int avgY=sumY/points;
  int avgZ=sumZ/points;
  int bestPe=bgltm->coords2rank(avgX, avgY, avgZ);
  avail->sortSource(bestPe);
  avail->reset();
}

void SCalcMapTable::sortByCentroid(PeList *avail, int plane, int stateX, int stateY, int grainsize, MapType2 *gsmap)
{
  int sumX=0, sumY=0, sumZ=0;
  int points=0;
      
  for(int state=stateX;state<stateX+grainsize;state++)
    {
      int X, Y, Z;
      bgltm->getCoordinatesByRank(gsmap->get(state,plane), X, Y, Z);
      sumX+=X;
      sumY+=Y;
      sumZ+=Z;
      points++;
    }
  for(int state=stateY;state<stateY+grainsize;state++)
    {
      int X, Y, Z;
      bgltm->getCoordinatesByRank(gsmap->get(state,plane), X, Y, Z);
      sumX+=X;
      sumY+=Y;
      sumZ+=Z;
      points++;
    }
  int avgX=sumX/points;
  int avgY=sumY/points;
  int avgZ=sumZ/points;
  int bestPe=bgltm->coords2rank(avgX, avgY, avgZ);
  avail->sortSource(bestPe);
  avail->reset();
}

void OrthoMapTable::sortByCentroid(PeList *avail, int nplanes, int state1, int state2, int numChunks, MapType4 *smap)
{
  int sumX = 0, sumY = 0, sumZ = 0;
  int points = 0;

  for(int plane=0; plane<nplanes; plane++)
    for(int chunk=0; chunk<numChunks; chunk++)
    {    
      int X, Y, Z;
      bgltm->getCoordinatesByRank(smap->get(plane, state1, state2, chunk), X, Y, Z);
      sumX += X;
      sumY += Y;
      sumZ += Z;
      points++;
    }
  int avgX = sumX/points;
  int avgY = sumY/points;
  int avgZ = sumZ/points;
  int bestPe = bgltm->coords2rank(avgX, avgY, avgZ);
  avail->sortSource(bestPe);
  avail->reset();
}

#else
// arguably meaningless in non torus case but it was easy to implement
void SCalcMapTable::sortByCentroid(PeList *avail, int plane, int stateX, int stateY, int grainsize, MapType2 *gsmap)
{
  int sumPe=0;
  int points=0;
  for(int state=stateX;state<stateX+grainsize;state++)
    {
      sumPe+=gsmap->get(state,plane);
      points++;
    }
  for(int state=stateY;state<stateY+grainsize;state++)
    {
      sumPe+=gsmap->get(state,plane);
      points++;
    }
  int bestPe=sumPe/points;
  avail->sortSource(bestPe);
  avail->reset();
}


void RhoRSMapTable::sortByCentroid(PeList *avail, int plane, int nstates, MapType2 *rsmap)
{
  int sumPe=0;
  int points=0;
  for(int state=0;state<nstates;state++)
    {
      sumPe+=rsmap->get(state,plane);
      points++;
    }
  int bestPe=sumPe/points;
  avail->sortSource(bestPe);
  avail->reset();
}

void OrthoMapTable::sortByCentroid(PeList *avail, int nplanes, int state1, int state2, int numChunks, MapType4 *smap)
{
  int sumPe = 0;
  int points = 0;

  for(int plane=0; plane<nplanes; plane++)
    for(int chunk=0; chunk<numChunks; chunk++)
    {    
      sumPe += smap->get(plane, state1, state2, chunk);
      points++;
    }
  int bestPe = sumPe/points;
  avail->sortSource(bestPe);
  avail->reset();
}

#endif
