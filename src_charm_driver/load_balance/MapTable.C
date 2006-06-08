#include "charm++.h"
#include "PeList.h"
#include "MapTable.h"


GSMapTable::GSMapTable(CkHashtableT <intdual, int > *_map, PeList *_availprocs, 
		       int _nchareG, double *_lines_per_chareG, double *_pts_per_chareG, 
		       int _nstates,  int _Gstates_per_pe)    : 
      nchareG(_nchareG), 
    lines_per_chareG(_lines_per_chareG),  pts_per_chareG(_pts_per_chareG), 
    nstates(_nstates), Gstates_per_pe(_Gstates_per_pe)
      { 
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
        for(int ychunk=0; ychunk<nchareG; ychunk=ychunk+m)
        {
                if(ychunk==(pm-rem)*m)
		    m=m+1;
                for(int xchunk=0; xchunk<nstates; xchunk=xchunk+l)
		{
			if(xchunk==(pl-srem)*l)
				l=l+1;
			if(xchunk==0 && ychunk==0) {}
			else
			  { // shift to next proc
			    srcpe=destpe;
			    destpe=availprocs->findNext();

			  }
			for(int state=xchunk; state<xchunk+l && state<nstates; state++)
			{
				for(int plane=ychunk; plane<ychunk+m && plane<nchareG; plane++)
				{
				    maptable->put(intdual(state, plane))=destpe;
				}
			}
                }
        }

#ifdef MAP_DEBUG
	CkPrintf("GSMap created on processor %d\n", CkMyPe());
	dump();
#endif
      }

SCalcMapTable::SCalcMapTable(CkHashtableT <intdual, int > *_map, PeList *_availprocs, 
			     int _nstates, int _nchareG,  int _grainsize, CmiBool _flag, 
			     int _nplanes,  double *_lines_per_chareG, 
			     double *_pts_per_chareG, int _scalc_per_plane,  
			     int _planes_per_pe, int _numChunks) : 
  max_states(_nstates), nchareG(_nchareG),  
  grainsize(_grainsize), symmetric(_flag), max_planes(_nplanes), 
  lines_per_chareG(_lines_per_chareG), pts_per_chareG(_pts_per_chareG),
  scalc_per_plane(_scalc_per_plane), planes_per_pe(_planes_per_pe), numChunks(_numChunks)
{ 
  
  int scobjs_per_pe, rem;
  int count=0, procno=0;
  int intidx[2];
  int lesser_scalc = 0;
  maptable=_map;
  availprocs=_availprocs;
  int srcpe=0,destpe=availprocs->findNext();

  if(planes_per_pe==0)
    CkAbort("Choose a smaller Gstates_per_pe\n");
  if(symmetric)
    {
      for(int i=1; i<=max_states/grainsize; i++)
	lesser_scalc += i;
      scobjs_per_pe = lesser_scalc*nchareG*numChunks/availprocs->count();
      rem = lesser_scalc*nchareG*numChunks % availprocs->count();
      if(rem!=0)
	scobjs_per_pe+=1;
      //if(CkMyPe()==0) CkPrintf("scobjs_per_pe %d grainsize %d nchareG %d scalc_per_plane %d planes_per_pe %d numChunks %d rem %d\n", scobjs_per_pe, grainsize, nchareG, scalc_per_plane, planes_per_pe, numChunks, rem);
			
      for(int pchunk=0; pchunk<nchareG; pchunk=pchunk+planes_per_pe)
	for(int newdim=0; newdim<numChunks; newdim++)
	  for(int xchunk=0; xchunk<max_states; xchunk=xchunk+grainsize)
	    for(int ychunk=xchunk; ychunk<max_states; ychunk=ychunk+grainsize)
	      //for(int newdim=0; newdim<numChunks; newdim++)
	      for(int plane=pchunk; plane<pchunk+planes_per_pe && plane<nchareG; plane++)
		{
		  CkArrayIndex4D idx4d(plane, xchunk, ychunk, newdim);
		  CmiMemcpy(intidx,idx4d.index,2*sizeof(int));  // our 4 shorts are now 2 ints
		  if(xchunk==0 && ychunk==0 && newdim==0 && plane==0)
		    {

		      maptable->put(intdual(intidx[0], intidx[1]))=destpe;
		      count++;
		    }
		  else
		    {
		      if(count<scobjs_per_pe)
			{
			  //if(CkMyPe()==0) CkPrintf("plane %d x %d y %d newdim %d = proc %d\n", plane, xchunk, ychunk, newdim, fp.next[0]*x*y+fp.next[1]*x+fp.next[2]);
			  maptable->put(intdual(intidx[0], intidx[1]))=destpe;
			  count++;
			}
		      else
			{
			  // new partition
			  procno++;
			  srcpe=destpe;
			  if(availprocs->noPes())
			    availprocs->rebuild();
			  //			  availprocs->sortSource(srcpe);
			  destpe=availprocs->findNext();
			  if(rem!=0)
			    if(procno==rem)
			      scobjs_per_pe-=1;
			  maptable->put(intdual(intidx[0], intidx[1]))=destpe;
			  count=0;
			  count++;
			}
		    }
		}
#ifdef MAP_DEBUG
      CkPrintf("Symmetric SCalcMap created on processor %d\n", CkMyPe());
      dump();
#endif
    }
  else
    {
      scobjs_per_pe = scalc_per_plane*nchareG*numChunks/availprocs->count();
      rem = scalc_per_plane*nchareG*numChunks % availprocs->count();
      if(rem!=0)
	scobjs_per_pe+=1;

      //if(CkMyPe()==0) CkPrintf("scobjs_per_pe %d grainsize %d nchareG %d scalc_per_plane %d planes_per_pe %d numChunks %d rem %d\n", scobjs_per_pe, grainsize, nchareG, scalc_per_plane, planes_per_pe, numChunks, rem);
			
      for(int pchunk=0; pchunk<nchareG; pchunk=pchunk+planes_per_pe)
	for(int newdim=0; newdim<numChunks; newdim++)
	  for(int xchunk=0; xchunk<max_states; xchunk=xchunk+grainsize)
	    for(int ychunk=0; ychunk<max_states; ychunk=ychunk+grainsize)
	      for(int plane=pchunk; plane<pchunk+planes_per_pe && plane<nchareG; plane++)
		{
		  CkArrayIndex4D idx4d(plane, xchunk, ychunk, newdim);
		  CmiMemcpy(intidx,idx4d.index,2*sizeof(int));  // our 4 shorts are now 2 ints
		  if(xchunk==0 && ychunk==0 && newdim==0 && plane==0)
		    {
		      //if(CkMyPe()==0) CkPrintf("plane %d x %d y %d newdim %d= proc 0\n", plane, xchunk, ychunk, newdim); 
		      maptable->put(intdual(intidx[0], intidx[1]))=0;
		      count++;
		    }
		  else
		    {
		      if(count<scobjs_per_pe)
			{
			  //if(CkMyPe()==0) CkPrintf("plane %d x %d y %d newdim %d= proc %d\n", plane, xchunk, ychunk, newdim, assign[0]*x*y+assign[1]*x+assign[2]);
			  maptable->put(intdual(intidx[0], intidx[1]))=destpe;
			  count++;
			}
		      else
			{  // new partition
			  procno++;
			  srcpe=destpe;
			  if(availprocs->noPes())
			      availprocs->rebuild();
			  //			  availprocs->sortSource(srcpe);
			  destpe=availprocs->findNext();
			  if(rem!=0)
			    if(procno==rem)
			      scobjs_per_pe-=1;
			  maptable->put(intdual(intidx[0], intidx[1]))=destpe;
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

RSMapTable::RSMapTable(CkHashtableT <intdual, int > *_map, PeList *_availprocs,
	int _nstates, int _sizeY, int _Rstates_per_pe) :
   nstates(_nstates), sizeY(_sizeY),
  Rstates_per_pe(_Rstates_per_pe)
{
        int c = 0;
	
	int l, m, pl, pm, srem, rem, i=0;
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

        m = sizeY / pm;
        rem = sizeY % pm;

        //CkPrintf("nstates %d sizeY %d Pes %d\n", nstates, sizeY, availprocs->count());	
	//CkPrintf("l %d, m %d pl %d pm %d srem %d rem %d\n", l, m,
	//pl, pm, srem, rem);
	int srcpe=0;
	int destpe=availprocs->findNext();
	
        for(int ychunk=0; ychunk<sizeY; ychunk=ychunk+m)
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
				for(int plane=ychunk; plane<ychunk+m && plane<sizeY; plane++)
				{
					if(xchunk==0 && ychunk==0)
					{
                                                c++;
						maptable->put(intdual(state, plane))=0;
						//CkPrintf("%d %d on 0\n", state, plane);
					}
					else
					{
                                                c++;
						maptable->put(intdual(state, plane))=destpe;
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

RhoRSMapTable::RhoRSMapTable(CkHashtableT <intdual, int > *_map, PeList *_availprocs, int _nchareRhoR): nchareRhoR(_nchareRhoR)
{
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
      rrsobjs_per_pe= nchareRhoR/(availprocs->count()/2);
      rem = nchareRhoR % (availprocs->count()/2);
      if(rem!=0)
	rrsobjs_per_pe += 1;
    }
  int destpe=availprocs->findNext(); 
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
              maptable->put(intdual(i, 0))=0;
              //CkPrintf("%d on %d\n", i, 0);
            }
	  else
            {

	      maptable->put(intdual(i, 0))=destpe;
	    }
	}
    }
#ifdef MAP_DEBUG
	CkPrintf("RhoRSMap created on processor %d\n", CkMyPe());
	dump();
#endif

}

RhoGSMapTable::RhoGSMapTable(CkHashtableT <intdual, int > *_map, PeList *_availprocs, int _nchareRhoG): nchareRhoG(_nchareRhoG)
{
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
      rgsobjs_per_pe= nchareRhoG/(availprocs->count()/2);
      rem = nchareRhoG % (availprocs->count()/2);
      if(rem!=0)
	rgsobjs_per_pe += 1;
    }
  int destpe=availprocs->findNext();
  //if(CkMyPe()==0) CkPrintf("nchareRhoG %d rgsobjs_per_pe %d rem %d\n", nchareRhoG, rgsobjs_per_pe, rem);   
  for(int chunk=0; chunk<nchareRhoG; chunk+=rgsobjs_per_pe)
    {
      if(rem!=0)
	if(chunk==rem*rgsobjs_per_pe)
	  rgsobjs_per_pe -= 1;  
      if(chunk==0) {}
      else
	{
	  srcpe=destpe;
	  //	  availprocs->sortSource(srcpe);
	  destpe=availprocs->findNext();
	}
      for(int i=chunk;i<chunk+rgsobjs_per_pe;i++)
	{
	  maptable->put(intdual(i, 0))=destpe;
	} 
    }
#ifdef MAP_DEBUG
  CkPrintf("RhoGSMap created on processor %d\n", CkMyPe());
  dump();
#endif
}


RhoGHartMapTable::RhoGHartMapTable(CkHashtableT <intdual, int > *_map, PeList *_availprocs, int _nchareRhoGHart): nchareRhoGHart(_nchareRhoGHart)
{
  int npes, procno=2, normal=0;
  int srcpe=0;

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
  for(int chunk=0; chunk<nchareRhoGHart; chunk+=rghobjs_per_pe)
    {
      if(rem!=0)
	if(chunk==rem*rghobjs_per_pe)
	  rghobjs_per_pe -= 1; 
      if(chunk==0) {}
      else
	{
	  srcpe=destpe;
	  //	  availprocs->sortSource(srcpe);
	  destpe=availprocs->findNext();
	}
      for(int i=chunk;i<chunk+rghobjs_per_pe;i++)
	{
	  maptable->put(intdual(i, 0))=destpe;
	} 
    }
#ifdef MAP_DEBUG
  CkPrintf("RhoGHartMap created on processor %d\n", CkMyPe());
  dump();
#endif

}
