#include "charm++.h"
#include "../../include/debug_flags.h"
#include "PeList.h"
#define USE_INT_MAP
#ifdef USE_INT_MAP
#include "IntMap.h"
//#define USE_INT_2on1
#ifdef USE_INT_2on1
typedef IntMap2on1 MapType2;
#else
typedef IntMap2on2 MapType2;
#endif
typedef IntMap4 MapType4;
#endif
#include "MapTable.h"

#define MAP_DEBUG


GSMapTable::GSMapTable(MapType2  *_map, PeList *_availprocs, 
		       int _nchareG, double *_lines_per_chareG, double *_pts_per_chareG, 
		       int _nstates,  int _Gstates_per_pe, bool useCuboidMap)    : 
      nchareG(_nchareG), 
    lines_per_chareG(_lines_per_chareG),  pts_per_chareG(_pts_per_chareG), 
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
		  destpe=availprocs->findNext();
		  if(availprocs->count()==0)
		    availprocs->reset();
		  /*		  if(availprocs->count()==0)
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
	    //}
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
		    destpe=availprocs->findNext();
		  }
	      }
	  }
#ifdef MAP_DEBUG
	CkPrintf("GSMap created on processor %d\n", CkMyPe());
	dump();
#endif
      }

SCalcMapTable::SCalcMapTable(MapType4  *_map, PeList *_availprocs, 
			     int _nstates, int _nchareG,  int _grainsize, CmiBool _flag, 
			     int _nplanes,  double *_lines_per_chareG, 
			     double *_pts_per_chareG, int _scalc_per_plane,  
			     int _planes_per_pe, 
			     int _numChunksA, 
			     int _numChunksS, 
			     MapType2  *gsmap, bool useCuboidMap, bool useCentroid, int boxSize) : 
  max_states(_nstates), nchareG(_nchareG),  
  grainsize(_grainsize), symmetric(_flag), max_planes(_nplanes), 
  lines_per_chareG(_lines_per_chareG), pts_per_chareG(_pts_per_chareG),
  scalc_per_plane(_scalc_per_plane), planes_per_pe(_planes_per_pe), numChunksAsym(_numChunksA), numChunksSym(_numChunksS)
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
#ifdef MAP_DEBUG
      CkPrintf(" lesser_scalc %d *nchareG %d *numChunksSym %d %% availprocs->count() %d = rem %d and scobjs_per_pe is %d\n", lesser_scalc,scalc_per_plane,nchareG,numChunksAsym , availprocs->count(),rem, scobjs_per_pe);
#endif

      int srcpe=0,destpe=availprocs->findNext();
      if(availprocs->count()==0)
	availprocs->reset();

      //if(CkMyPe()==0) CkPrintf("scobjs_per_pe %d grainsize %d nchareG %d scalc_per_plane %d planes_per_pe %d numChunks %d rem %d\n", scobjs_per_pe, grainsize, nchareG, scalc_per_plane, planes_per_pe, numChunksSym, rem);
			
      for(int pchunk=0; pchunk<nchareG; pchunk=pchunk+planes_per_pe)
	for(int newdim=0; newdim<numChunksSym; newdim++)
	  for(int xchunk=0; xchunk<max_states; xchunk=xchunk+grainsize)
	    for(int ychunk=xchunk; ychunk<max_states; ychunk=ychunk+grainsize)
	      //for(int newdim=0; newdim<numChunksSym; newdim++)
	      for(int plane=pchunk; plane<pchunk+planes_per_pe && plane<nchareG; plane++)
		{
		  CkArrayIndex4D idx4d(plane, xchunk, ychunk, newdim);
		  CmiMemcpy(intidx,idx4d.index,2*sizeof(int));  // our 4 shorts are now 2 ints

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
		      if(availprocs->noPes())
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
#ifdef MAP_DEBUG
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
#ifdef MAP_DEBUG
      CkPrintf(" scalc_per_plane %d *nchareG %d *numChunksAsym %d %% availprocs->count() %d = rem %d and scobjs_per_pe is %d\n", scalc_per_plane,nchareG,numChunksAsym , availprocs->count(),rem, scobjs_per_pe);
#endif
      //if(CkMyPe()==0) CkPrintf("scobjs_per_pe %d grainsize %d nchareG %d scalc_per_plane %d planes_per_pe %d numChunksAsym %d rem %d\n", scobjs_per_pe, grainsize, nchareG, scalc_per_plane, planes_per_pe, numChunksAsym, rem);
      int srcpe=0,destpe=0;
      if(!useCentroid && !useCuboidMap)
	destpe=availprocs->findNext();
      if(availprocs->count()==0)
	availprocs->reset();
      if(useCuboidMap)
	{ // in the cuboid map case we place all planes box by box
	  for(int plane=0; plane<nchareG; plane++)
	    { // could restrict list to the gs box here
	      PeList *thisPlaneBox= new PeList(availprocs, plane*boxSize, boxSize);
	      //	      PeList *thisPlaneBox= availprocs;
	      if(!useCentroid)
		destpe=thisPlaneBox->findNext();
	      for(int xchunk=0; xchunk<max_states; xchunk=xchunk+grainsize)
		for(int ychunk=0; ychunk<max_states; ychunk=ychunk+grainsize)
		  { // could find centroid here
		    if(useCentroid)
		      {
			sortByCentroid(thisPlaneBox, plane, xchunk, ychunk, grainsize, gsmap);
			destpe=thisPlaneBox->findNext();
			if(thisPlaneBox->noPes())
			  thisPlaneBox->reset();
			count=0;
		      }
		    for(int newdim=0; newdim<numChunksAsym; newdim++)
		      {
			CkArrayIndex4D idx4d(plane, xchunk, ychunk, newdim);
			CmiMemcpy(intidx,idx4d.index,2*sizeof(int));  // our 4 shorts are now 2 ints
			if(count<scobjs_per_pe)
			  {
			    //if(CkMyPe()==0) CkPrintf("plane %d x %d y %d newdim %d= proc %d\n", plane, xchunk, ychunk, newdim, assign[0]*x*y+assign[1]*x+assign[2]);
#ifdef USE_INT_MAP		      
			    maptable->set(plane, xchunk, ychunk, newdim,destpe);
#else
			    maptable->put(intdual(intidx[0], intidx[1]))=destpe;
#endif
			    count++;
			  }
			else
			  {  // new partition
			    procno++;
			    srcpe=destpe;
			    if(thisPlaneBox->noPes())
			      thisPlaneBox->reset();

			    destpe=thisPlaneBox->findNext();
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
	      //	      delete thisPlaneBox;
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
		      CmiMemcpy(intidx,idx4d.index,2*sizeof(int));  // our 4 shorts are now 2 ints

		      if(count<scobjs_per_pe)
			{
			  //if(CkMyPe()==0) CkPrintf("plane %d x %d y %d newdim %d= proc %d\n", plane, xchunk, ychunk, newdim, assign[0]*x*y+assign[1]*x+assign[2]);
#ifdef USE_INT_MAP		      
			  maptable->set(plane, xchunk, ychunk, newdim,destpe);
#else
			  maptable->put(intdual(intidx[0], intidx[1]))=destpe;
#endif
			  count++;
			}
		      else
			{  // new partition
			  procno++;
			  srcpe=destpe;
			  if(availprocs->noPes())
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
#ifdef MAP_DEBUG
      CkPrintf("Asymmetric SCalcMap created on processor %d\n", CkMyPe());
      dump();
#endif
    }

}

RSMapTable::RSMapTable(MapType2  *_map, PeList *_availprocs,
	int _nstates, int _sizeZ, int _Rstates_per_pe) :
   nstates(_nstates), sizeZ(_sizeZ),
  Rstates_per_pe(_Rstates_per_pe)
{
        int c = 0;
	
	int l, m, pl, pm, srem, rem, i=0;
	reverseMap=NULL;
	maptable=_map;
	availprocs=_availprocs;
        
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

        //CkPrintf("nstates %d sizeZ %d Pes %d\n", nstates, sizeZ, availprocs->count());	
	//CkPrintf("l %d, m %d pl %d pm %d srem %d rem %d\n", l, m,
	//pl, pm, srem, rem);
	int srcpe=0;
	int destpe=availprocs->findNext();

	if(availprocs->count()==0)
	  availprocs->reset();	

        for(int ychunk=0; ychunk<sizeZ; ychunk=ychunk+m)
        {
                if(ychunk==(pm-rem)*m)
              		m=m+1;
        	for(int xchunk=0; xchunk<nstates; xchunk=xchunk+l)
		{
                	if(xchunk==(pl-srem)*l)
				l=l+1;
			if(xchunk==0 && ychunk==0) {}
			else
			{
			    srcpe=destpe;
			    //			    availprocs->sortSource(srcpe);
			    destpe=availprocs->findNext();
			}
                        c=0;
			for(int state=xchunk; state<xchunk+l && state<nstates; state++)
			{
				for(int plane=ychunk; plane<ychunk+m && plane<sizeZ; plane++)
				{
					if(xchunk==0 && ychunk==0)
					{
                                                c++;
#ifdef USE_INT_MAP
						maptable->set(state, plane,0);
#else
						maptable->put(intdual(state, plane))=0;
#endif
						//CkPrintf("%d %d on 0\n", state, plane);
					}
					else
					{
                                                c++;
#ifdef USE_INT_MAP
						maptable->set(state, plane, destpe);

#else						
						maptable->put(intdual(state, plane))=destpe;
#endif
					}
				}
			}
                }
        }
#ifdef MAP_DEBUG
	CkPrintf("RSMap created on processor %d\n", CkMyPe());
	dump();
#endif
}


RSPMapTable::RSPMapTable(MapType2  *_map, 
			 PeList *_availprocs, PeList *exclusion,
			 int _nstates, int _sizeZNL, int _Rstates_per_pe,
			 int boxSize, bool useCuboidMap) :
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
  PeList *RPPlist=availprocs;
  states_per_pe=Rstates_per_pe;		// no of states in one chunk
  pl = nstates / states_per_pe;
  int afterExclusion=exclusion->count();
  
  if(afterExclusion > pl*sizeZNL)
    { // we can fit the exclusion without blinking
      CkPrintf("RPP using density exclusion to avoid %d processors\n",exclusion->count());
      RPPlist=exclusion;
    }
  else
    {// so an rstates_per_pe chosen for realstate might be too big
      if(afterExclusion > (sizeZNL+nstates)*2)
	{ // set states_per_pe to fit the exclusion

	  states_per_pe=totalChares/afterExclusion/2;
	  CkPrintf("RPP adjusting states per pe from %d to %d to use density exclusion to stay within %d processors\n",Rstates_per_pe, states_per_pe, exclusion->count());
	  RPPlist=exclusion;
	}
      else
	{
	  CkPrintf("RPP with %d chares ignoring density exclusion of %d out of %d processors\n",totalChares, exclusion->count(), availprocs->count());
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

  //CkPrintf("nstates %d sizeZNL %d Pes %d\n", nstates, sizeZNL, RPPlist->count());	
  //CkPrintf("l %d, m %d pl %d pm %d srem %d rem %d\n", l, m,
  //pl, pm, srem, rem);
  int srcpe=0;

  /*  RPP isn't planewise with GPP.  Need different cuboid code
      if(useCuboidMap)
      {
      CkPrintf("nstates %d sizeZNL %d Pes %d states_per_pe %d boxSize %d\n", nstates, sizeZNL, RPPlist->count(), states_per_pe, boxSize);	
      for(int plane=0;plane<sizeZNL;plane++)
      {
      // get the cube for this plane
      PeList *thisPlaneBox= new PeList(RPPlist, plane*boxSize, boxSize);
      // try a subtraction from the exclusion
      if(exclusion->count()>0)
      *thisPlaneBox-*exclusion;
      if(thisPlaneBox->count()<nstates/states_per_pe)
      {
      CkPrintf("RSP ignoring exclusion map for plane %d count %d\n",plane,thisPlaneBox->count());
      // if not enough ignore the exclusion map
      delete thisPlaneBox;
      thisPlaneBox= new PeList(RPPlist,plane*bSize,bSize);
      }
      else
      { // repair the index
      thisPlaneBox->reindex();
      }
      destpe=thisPlaneBox->findNext();
      for(int state=0;state<nstates;state+=states_per_pe)
      {
      for(int stateperpe=0;stateperpe<states_per_pe;stateperpe++)
      {
      maptable->put(intdual(state+stateperpe, plane))=destpe;
      }
      destpe=thisPlaneBox->findNext();
      if(thisPlaneBox->count()==0)
      thisPlaneBox->reset();
      }
      delete thisPlaneBox;
      }
      }
      else
  */
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
	    srcpe=destpe;
		//			    RPPlist->sortSource(srcpe);
	    for(int state=xchunk; state<xchunk+states_per_pe && state<nstates; state++)
	      {
		for(int plane=ychunk; plane<ychunk+m && plane<sizeZNL; plane++)
		  {
#ifdef USE_INT_MAP
			maptable->set(state, plane,destpe);
#else
			maptable->put(intdual(state, plane))=destpe;
#endif

		  }
	      }
	    destpe=RPPlist->findNext();
	    if(RPPlist->count()==0)
	      RPPlist->reset();
	  }
      }
  }
#ifdef MAP_DEBUG
  CkPrintf("RSPMap created on processor %d\n", CkMyPe());
  dump();
#endif
}


RhoRSMapTable::RhoRSMapTable(MapType2  *_map, PeList *_availprocs, int _nchareRhoR): nchareRhoR(_nchareRhoR)
{
  reverseMap=NULL;
  maptable=_map;
  availprocs=_availprocs;
  int rrsobjs_per_pe, rem;
  int srcpe=0;

  if(availprocs->count()==1)
    {
      rrsobjs_per_pe= nchareRhoR;
      rem=0;
    }
  else
    {
      rrsobjs_per_pe= nchareRhoR/(availprocs->count());
      rem = nchareRhoR % (availprocs->count());
      if(rem!=0)
	rrsobjs_per_pe += 1;
    }
  int destpe=availprocs->findNext(); 
  if(availprocs->count()==0)
    availprocs->reset();

  //if(CkMyPe()==0) CkPrintf("nchareRhoR %d rrsobjs_per_pe %d rem %d\n", nchareRhoR, rrsobjs_per_pe, rem);   
  for(int chunk=0; chunk<nchareRhoR; chunk+=rrsobjs_per_pe)
    {
      if(rem!=0)
	if(chunk==rem*rrsobjs_per_pe)
	  rrsobjs_per_pe -= 1;
      if(chunk==0) {}
      else
	{
	  srcpe=destpe;
	  //	  availprocs->sortSource(srcpe);
	  destpe=availprocs->findNext();
	}
      for(int i=chunk;i<chunk+rrsobjs_per_pe;i++)
	{
	  if(chunk==0)
            {
#ifdef USE_INT_MAP
              maptable->set(i, 0,0);
#else
              maptable->put(intdual(i, 0))=0;
#endif
              //CkPrintf("%d on %d\n", i, 0);
            }
	  else
            {
#ifdef USE_INT_MAP
	      maptable->set(i, 0, destpe);
#else
	      maptable->put(intdual(i, 0))=destpe;
#endif
	    }
	}
    }
#ifdef MAP_DEBUG
	CkPrintf("RhoRSMap created on processor %d\n", CkMyPe());
	dump();
#endif

}

RhoGSMapTable::RhoGSMapTable(MapType2  *_map, PeList *_availprocs, int _nchareRhoG): nchareRhoG(_nchareRhoG)
{
  reverseMap=NULL;
  maptable=_map;
  availprocs=_availprocs;
  int rgsobjs_per_pe, rem;
  int srcpe=0;

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
  int destpe=availprocs->findNext();

  if(availprocs->count()==0)
    availprocs->reset();

  //if(CkMyPe()==0) CkPrintf("nchareRhoG %d rgsobjs_per_pe %d rem %d\n", nchareRhoG, rgsobjs_per_pe, rem);   
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
      srcpe=destpe;
      if(chunk+1<nchareRhoG)
	destpe=availprocs->findNext();
    }
#ifdef MAP_DEBUG
  CkPrintf("RhoGSMap created on processor %d\n", CkMyPe());
  dump();
#endif
}


RhoRHartMapTable::RhoRHartMapTable(MapType2  *_map, PeList *_availprocs, int _nchareRhoRHart): nchareRhoRHart(_nchareRhoRHart)
{
  reverseMap=NULL;
  maptable=_map;
  availprocs=_availprocs;
  int rrsobjs_per_pe, rem;
  int srcpe=0;
  if(availprocs->count()==0)
    availprocs->reset();

  if(availprocs->count()==1)
    {
      rrsobjs_per_pe= nchareRhoRHart;
      rem=0;
    }
  else
    {
      rrsobjs_per_pe= nchareRhoRHart/(availprocs->count());
      rem = nchareRhoRHart % (availprocs->count());
      if(rem!=0)
	rrsobjs_per_pe += 1;
    }
  int destpe=availprocs->findNext(); 
  srcpe=destpe;
  if(availprocs->count()==0)
    availprocs->reset();

  //if(CkMyPe()==0) CkPrintf("nchareRhoR %d rrsobjs_per_pe %d rem %d\n", nchareRhoRHart, rrsobjs_per_pe, rem);   
  for(int chunk=0; chunk<nchareRhoRHart; chunk+=rrsobjs_per_pe)
    {
      if(rem!=0)
	if(chunk==rem*rrsobjs_per_pe)
	  rrsobjs_per_pe -= 1;
      for(int i=chunk;i<chunk+rrsobjs_per_pe;i++)
	{
#ifdef USE_INT_MAP
	  maptable->set(i, 0, destpe);
#else
	  maptable->put(intdual(i, 0))=destpe;
#endif
	}
      srcpe=destpe;
      //	  availprocs->sortSource(srcpe);
      if(chunk+1<nchareRhoRHart)
	destpe=availprocs->findNext();

    }
#ifdef MAP_DEBUG
	CkPrintf("RhoRHartMap created on processor %d\n", CkMyPe());
	dump();
#endif

}

RhoGHartMapTable::RhoGHartMapTable(MapType2  *_map, PeList *_availprocs, int _nchareRhoGHart): nchareRhoGHart(_nchareRhoGHart)
{
  int npes;
  int srcpe=0;
  reverseMap=NULL;
  maptable=_map;
  availprocs=_availprocs;
  npes=availprocs->count();

  int rghobjs_per_pe, rem;

  if(availprocs->count()==1)
    {
      rghobjs_per_pe= nchareRhoGHart;
      rem=0;
    }
  else
    {
      rghobjs_per_pe= nchareRhoGHart/npes;
      rem = nchareRhoGHart % npes;
      if(rem!=0)
	rghobjs_per_pe += 1;
    }
  int destpe=availprocs->findNext();
  if(availprocs->count()==0)
    availprocs->reset();
  for(int chunk=0; chunk<nchareRhoGHart; chunk+=rghobjs_per_pe)
    {
      if(rem!=0)
	if(chunk==rem*rghobjs_per_pe)
	  rghobjs_per_pe -= 1; 
      for(int i=chunk;i<chunk+rghobjs_per_pe;i++)
	{
#ifdef USE_INT_MAP
	  maptable->set(i, 0, destpe);
#else
	  maptable->put(intdual(i, 0))=destpe;
#endif
	} 
      srcpe=destpe;
      //	  availprocs->sortSource(srcpe);
      if(chunk+1<nchareRhoGHart)
	destpe=availprocs->findNext();
    }
#ifdef MAP_DEBUG
  CkPrintf("RhoGHartMap created on processor %d\n", CkMyPe());
  dump();
#endif

}


void MapTable::makeReverseMap()
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

#ifdef CMK_VERSION_BLUEGENE
extern 	BGLTorusManager *bgltm;

void SCalcMapTable::sortByCentroid(PeList *avail, int plane, int stateX, int stateY, int grainsize, MapType2 *gsmap)
{
      int sumX=0, sumY=0, sumZ=0;
      int points=0;
      
      for(int state=stateX;state<stateX+grainsize;state++)
	{
	  int X, Y, Z;
	  bgltm->getCoordinatesByRank(gsmap->get(state,plane),X, Y, Z);
	  sumX+=X;
	  sumY+=Y;
	  sumZ+=Z;
	  points++;
	}
      for(int state=stateY;state<stateY+grainsize;state++)
	{
	  int X, Y, Z;
	  bgltm->getCoordinatesByRank(gsmap->get(state,plane),X, Y, Z);
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
    }
#else
// arguably meaningless in non torus case but it was easy to implement
void SCalcMapTable::sortByCentroid(PeList *avail, int plane, int stateX, int stateY, int grainsize, MapType2 *gsmap)
    {
      int sumPe=0;
      int points=0;
      for(int state=stateX;state<stateX+grainsize;state++)
	{
	  sumPe=gsmap->get(state,plane);
	  points++;
	}
      for(int state=stateY;state<stateY+grainsize;state++)
	{
	  sumPe=gsmap->get(state,plane);
	  points++;
	}
      int bestPe=sumPe/points;
      avail->sortSource(bestPe);
    }
#endif
