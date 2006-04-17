/** \file TopoMap.C
 * The functions in this file generate a topology sensitive
 * mapping scheme for the GSpace, PairCalc and RealSpace objects.
 * 
 */
#include "charm++.h"
#include "cpaimd.h"
#include "util.h"

#ifdef USE_TOPOMAP
#include "FindProcessor.h"
#ifdef CMK_VERSION_BLUEGENE
#include "bgltorus.h"
#endif

/** 
 * Function for GSpace objects
 * Makes the map for the first time and places
 * all the objects.
 */
void GSMap::makemap()
{
	FindProcessor fp;
	int x = CkNumPes();
	int y = 1;
	int z = 1;
	int vn = 0;
        int c = 0;
        
#ifdef CMK_VERSION_BLUEGENE
	BGLTorusManager *bgltm = BGLTorusManager::getObject();
	x = bgltm->getXSize();
	y = bgltm->getYSize();
	z = bgltm->getZSize();
	vn = bgltm->isVnodeMode();
#endif

#ifdef MAP_DEBUG
        /*char name[30];
        sprintf(name, "proc%d", CkMyPe());
        f=fopen(name, "w");*/
#endif

	fp.count=1;
	fp.nopX=x;
	fp.nopY=y;
	fp.nopZ=z;
	
	int assign[3]={0, 0, 0};
	int w = 0;
	for(int i=0;i<3;i++)
		fp.start[i]=fp.next[i]=0;
	/*int gsobjs_per_pe;
	
	if((nstates*nchareG) % CkNumPes() == 0)
	    gsobjs_per_pe = (nstates*nchareG)/CkNumPes();
	else
	    gsobjs_per_pe = (nstates*nchareG)/CkNumPes()+1;*/
	int l=Gstates_per_pe;
	int m, pl, pm, rem;
        
        pl = nstates / l;
        pm = CkNumPes() / pl;
        
	if(pm==0)
	  CkAbort("Choose a larger Gstates_per_pe\n");
        
	m = nchareG / pm;
        rem = nchareG % pm;
        
        planes_per_pe=m;

        //CkPrintf("nstates %d nchareG %d Pes %d, gsobjs_per_pe %d\n", nstates, nchareG, CkNumPes(), gsobjs_per_pe);	
	//CkPrintf("l %d, m %d pl %d pm %d rem %d\n", l, m, pl, pm, rem);
	
        for(int ychunk=0; ychunk<nchareG; ychunk=ychunk+m)
        {
                if(ychunk==(pm-rem)*m)
                  m=m+1;
                for(int xchunk=0; xchunk<nstates; xchunk=xchunk+l)
		{
			if(xchunk==0 && ychunk==0) {}
			else
			{
				for(int i=0;i<3;i++)
					fp.start[i]=fp.next[i];
				if(fp.start[2]>x/2)
					assign[2]=fp.start[2]-x;
				else
					assign[2]=fp.start[2];
				if(fp.start[1]>y/2)
					assign[1]=fp.start[1]-y;
				else
					assign[1]=fp.start[1];
				if(fp.start[0]>z/2)
					assign[0]=fp.start[0]-z;
				else
					assign[0]=fp.start[0];
				if(vn==0)
					fp.findNextInTorus(assign);
                                else
                                {
                                	fp.findNextInTorusV(w, assign);
                                        w = fp.w;
                                }	
			}
                        c=0;
			for(int state=xchunk; state<xchunk+l && state<nstates; state++)
			{
				for(int plane=ychunk; plane<ychunk+m && plane<nchareG; plane++)
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
						if(vn==0)
						  maptable->put(intdual(state, plane))=fp.next[0]*x*y+fp.next[1]*x+fp.next[2];
                                                else
                                                  maptable->put(intdual(state, plane))=(fp.next[0]*x*y+fp.next[1]*x+fp.next[2])*2+fp.w;
						//CkPrintf("%d %d on %d\n", state, plane, fp.next[0]*x*y+fp.next[1]*x+fp.next[2]);
						
					}
				}
			}
#ifdef MAP_DEBUG
                        /*if(xchunk==0 && ychunk==0)
                          fprintf(f, "objs %d pe %d\n", c, 0);
                        else
                          fprintf(f, "objs %d pe %d\n", c, fp.next[0]*x*y+fp.next[1]*x+fp.next[2]);
*/
#endif
                }
        }

#ifdef MAP_DEBUG
        //fclose(f);
	CkPrintf("GSMap created on processor %d\n", CkMyPe());
#endif
}

/**
 * Function which returns the processor no. for a GSPace object.
 * Uses the hash table built by makemap
 */
int GSMap::procNum(int handle, const CkArrayIndex &index)
{
	CkArrayIndex2D idx2d = *(CkArrayIndex2D *) &index;
	if(maptable==NULL)
	{
	        maptable= new CkHashtableT<intdual, int> (nstates*nchareG);
		makemap();
	}
	return maptable->get(intdual(idx2d.index[0], idx2d.index[1]));
}
   
/** Function for PairCalc objects
 *
 */
void SCalcMap::makemap()
{
	FindProcessor fp;
	int n = scalc_per_plane;
	int scobjs_per_pe;
	int grainsize = gs;
	int count=0;
	int intidx[2];
	int lesser_scalc = 0;
	int x = CkNumPes();
	int y = 1;
	int z = 1;
	int vn = 0;
	
#if CMK_VERSION_BLUEGENE
	BGLTorusManager *bgltm = BGLTorusManager::getObject();
	x = bgltm->getXSize();
	y = bgltm->getYSize();
	z = bgltm->getZSize();
	vn = bgltm->isVnodeMode();
#endif
	fp.count=1;
	fp.nopX=x;
	fp.nopY=y;
	fp.nopZ=z;
	fp.w = 0;

	if(planes_per_pe==0)
		CkAbort("Choose a smaller Gstates_per_pe\n");
	int assign[3]={0, 0, 0};
	int w = 0;
        for(int i=0;i<3;i++)
        	fp.start[i]=fp.next[i]=0;
		
	if(symmetric)
	{
		for(int i=1; i<=max_states/grainsize; i++)
			lesser_scalc += i;
		if(lesser_scalc*nchareG*numChunks % CkNumPes() == 0)
			scobjs_per_pe = lesser_scalc*nchareG*numChunks/CkNumPes();
		else
			scobjs_per_pe = lesser_scalc*nchareG*numChunks/CkNumPes() + 1;
		//CkPrintf("scobjs_per_pe %d grainsize %d nchareG %d scalc_per_plane %d planes_per_pe %d numChunks %d\n", scobjs_per_pe, grainsize, nchareG, scalc_per_plane, planes_per_pe, numChunks);
			
		for(int pchunk=0; pchunk<nchareG; pchunk=pchunk+planes_per_pe)
			for(int xchunk=0; xchunk<max_states; xchunk=xchunk+grainsize)
				for(int ychunk=xchunk; ychunk<max_states; ychunk=ychunk+grainsize)
				    for(int newdim=0; newdim<numChunks; newdim++)
					for(int plane=pchunk; plane<pchunk+planes_per_pe && plane<nchareG; plane++)
					{
						CkArrayIndex4D idx4d(plane, xchunk, ychunk, newdim);
						CmiMemcpy(intidx,idx4d.index,2*sizeof(int));  // our 4 shorts are now 2 ints
						if(xchunk==0 && ychunk==0 && newdim==0 && plane==0)
						{
							//CkPrintf("plane %d x %d y %d newdim %d = proc 0\n", plane, xchunk, ychunk, newdim); 
							maptable->put(intdual(intidx[0], intidx[1]))=0;
							count++;
						}
						else
						{
							if(count<scobjs_per_pe)
							{
								//CkPrintf("plane %d x %d y %d newdim %d = proc %d\n", plane, xchunk, ychunk, newdim, fp.next[0]*x*y+fp.next[1]*x+fp.next[2]);
								if(vn==0)
								  maptable->put(intdual(intidx[0], intidx[1]))=fp.next[0]*x*y+fp.next[1]*x+fp.next[2];
                                                		else
                                                		  maptable->put(intdual(intidx[0], intidx[1]))=(fp.next[0]*x*y+fp.next[1]*x+fp.next[2])*2+fp.w;
								count++;
							}
							else
							{
								count=0;
								for(int i=0;i<3;i++)
									fp.start[i]=fp.next[i];
								if(fp.start[2]>x/2)
									assign[2]=fp.start[2]-x;
								else
									assign[2]=fp.start[2];
								if(fp.start[1]>y/2)
									assign[1]=fp.start[1]-y;
								else
									assign[1]=fp.start[1];
								if(fp.start[0]>z/2)
									assign[0]=fp.start[0]-z;
								else
									assign[0]=fp.start[0];
								if(vn==0)
									fp.findNextInTorus(assign);
                                				else
                                                                {
                                					fp.findNextInTorusV(w, assign);
                                                                        w = fp.w;
                                                                }
								for(int i=0;i<3;i++)
									assign[i]=fp.next[i];
								
								//CkPrintf("plane %d x %d y %d newdim %d = proc %d\n", plane, xchunk, ychunk, newdim, fp.next[0]*x*y+fp.next[1]*x+fp.next[2]);
								if(vn==0)
								  maptable->put(intdual(intidx[0], intidx[1]))=fp.next[0]*x*y+fp.next[1]*x+fp.next[2];
                                                		else
                                                		  maptable->put(intdual(intidx[0], intidx[1]))=(fp.next[0]*x*y+fp.next[1]*x+fp.next[2])*2+fp.w;
								count++;
							}
						}
					}
#ifdef MAP_DEBUG
	CkPrintf("Symmetric SCalcMap created on processor %d\n", CkMyPe());
#endif
	}
	else
	{
		if(n*nchareG*numChunks % CkNumPes() == 0)
		  scobjs_per_pe = n*nchareG*numChunks/CkNumPes();
		else
		  scobjs_per_pe = n*nchareG*numChunks/CkNumPes() + 1;
		//CkPrintf("scobjs_per_pe %d grainsize %d nchareG %d scalc_per_plane %d planes_per_pe %d numChunks %d\n", scobjs_per_pe, grainsize, nchareG, scalc_per_plane, planes_per_pe, numChunks);
			
		for(int pchunk=0; pchunk<nchareG; pchunk=pchunk+planes_per_pe)
			for(int xchunk=0; xchunk<max_states; xchunk=xchunk+grainsize)
				for(int ychunk=0; ychunk<max_states; ychunk=ychunk+grainsize)
				    for(int newdim=0; newdim<numChunks; newdim++)
					for(int plane=pchunk; plane<pchunk+planes_per_pe && plane<nchareG; plane++)
					{
						CkArrayIndex4D idx4d(plane, xchunk, ychunk, newdim);
						CmiMemcpy(intidx,idx4d.index,2*sizeof(int));  // our 4 shorts are now 2 ints
						if(xchunk==0 && ychunk==0 && newdim==0 && plane==0)
						{
							//CkPrintf("plane %d x %d y %d newdim %d= proc 0\n", plane, xchunk, ychunk, newdim); 
							maptable->put(intdual(intidx[0], intidx[1]))=0;
							count++;
						}
						else
						{
							if(count<scobjs_per_pe)
							{
								//CkPrintf("plane %d x %d y %d newdim %d= proc %d\n", plane, xchunk, ychunk, newdim, assign[0]*x*y+assign[1]*x+assign[2]);
								if(vn==0)
								  maptable->put(intdual(intidx[0], intidx[1]))=fp.next[0]*x*y+fp.next[1]*x+fp.next[2];
                                                		else
                                                		  maptable->put(intdual(intidx[0], intidx[1]))=(fp.next[0]*x*y+fp.next[1]*x+fp.next[2])*2+fp.w;
								count++;
							}
							else
							{
								count=0;
								for(int i=0;i<3;i++)
									fp.start[i]=fp.next[i];
								if(fp.start[2]>x/2)
									assign[2]=fp.start[2]-x;
								else
									assign[2]=fp.start[2];
								if(fp.start[1]>y/2)
									assign[1]=fp.start[1]-y;
								else
									assign[1]=fp.start[1];
								if(fp.start[0]>z/2)
									assign[0]=fp.start[0]-z;
								else
									assign[0]=fp.start[0];
								if(vn==0)
									fp.findNextInTorus(assign);
                                				else
                                                                {
                                					fp.findNextInTorusV(w, assign);
                                                                        w = fp.w;
                                                                }
								for(int i=0;i<3;i++)
									assign[i]=fp.next[i];
								//CkPrintf("plane %d x %d y %d newdim %d= proc %d\n", plane, xchunk, ychunk, newdim, assign[0]*x*y+assign[1]*x+assign[2]);
								if(vn==0)
								  maptable->put(intdual(intidx[0], intidx[1]))=fp.next[0]*x*y+fp.next[1]*x+fp.next[2];
                                                		else
                                                		  maptable->put(intdual(intidx[0], intidx[1]))=(fp.next[0]*x*y+fp.next[1]*x+fp.next[2])*2+fp.w;
								count++;
							}
						}
					}
#ifdef MAP_DEBUG
	CkPrintf("Asymmetric SCalcMap created on processor %d\n", CkMyPe());
#endif
	}
}

/**
 *
 */
int SCalcMap::procNum(int handle, const CkArrayIndex &index)
{
	CkArrayIndex4D &idx4d = *(CkArrayIndex4D *)&index;
	int intidx[2];
	CmiMemcpy(intidx,idx4d.index,2*sizeof(int));  // our 4 shorts are now 2 ints

	if(maptable==NULL)
	{
		maptable= new CkHashtableT<intdual, int> (scalc_per_plane*nchareG); 
		makemap();
	}
	return maptable->get(intdual(intidx[0], intidx[1]));
}

/** Function for RealSpace objects
 *
 */
void RSMap::makemap()
{
	//CkPrintf("%d start rsmap\n", CkMyPe());
	FindProcessor fp;
	int x = CkNumPes();
	int y = 1;
	int z = 1;
	int vn = 0;
        int c = 0;
	
#ifdef CMK_VERSION_BLUEGENE
	BGLTorusManager *bgltm = BGLTorusManager::getObject();
	x = bgltm->getXSize();
	y = bgltm->getYSize();
	z = bgltm->getZSize();
	vn = bgltm->isVnodeMode();
#endif
	
	fp.count=1;
	fp.nopX=x;
	fp.nopY=y;
	fp.nopZ=z;
	
	int assign[3]={0, 0, 0};
	int w = 0;
	for(int i=0;i<3;i++)
		fp.start[i]=fp.next[i]=0;
	
	/*int rsobjs_per_pe;
	if(nstates*sizeY % CkNumPes() == 0)
	    rsobjs_per_pe = nstates*sizeY/CkNumPes();
	else
	    rsobjs_per_pe = nstates*sizeY/CkNumPes()+1;*/
	int l=Rstates_per_pe;
	int m, pl, pm, rem;
        
        pl = nstates / l;
        pm = CkNumPes() / pl;
        
	if(pm==0)
	  CkAbort("Choose a larger Rstates_per_pe\n");

        m = sizeY / pm;
        rem = sizeY % pm;

        //CkPrintf("nstates %d sizeY %d Pes %d, rsobjs_per_pe %d\n", nstates, sizeY, CkNumPes(), rsobjs_per_pe);	
	//CkPrintf("l %d, m %d pl %d pm %d rem %d\n", l, m, pl, pm, rem);
	
        for(int ychunk=0; ychunk<sizeY; ychunk=ychunk+m)
        {
                if(ychunk==(pm-rem)*m)
                  m=m+1;
                for(int xchunk=0; xchunk<nstates; xchunk=xchunk+l)
		{
			if(xchunk==0 && ychunk==0) {}
			else
			{
				for(int i=0;i<3;i++)
					fp.start[i]=fp.next[i];
				if(fp.start[2]>x/2)
					assign[2]=fp.start[2]-x;
				else
					assign[2]=fp.start[2];
				if(fp.start[1]>y/2)
					assign[1]=fp.start[1]-y;
				else
					assign[1]=fp.start[1];
				if(fp.start[0]>z/2)
					assign[0]=fp.start[0]-z;
				else
					assign[0]=fp.start[0];
				if(vn==0)
                                  fp.findNextInTorus(assign);
                                else
                                {
                                        fp.findNextInTorusV(w, assign);
                                        w = fp.w;
                                }	
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
						if(vn==0)
						  maptable->put(intdual(state, plane))=fp.next[0]*x*y+fp.next[1]*x+fp.next[2];
                                                else
                                                {
                                                  maptable->put(intdual(state, plane))=(fp.next[0]*x*y+fp.next[1]*x+fp.next[2])*2+fp.w;
                                                  //CkPrintf("%d %d on %d\n", state, plane, (fp.next[0]*x*y+fp.next[1]*x+fp.next[2])*2+fp.w);
                                                }
						//CkPrintf("%d %d on %d\n", state, plane, fp.next[0]*x*y+fp.next[1]*x+fp.next[2]);
						
					}
				}
			}
                }
        }
#ifdef MAP_DEBUG
	CkPrintf("RSMap created on processor %d\n", CkMyPe());
#endif
}

/**
 *
 */
int RSMap::procNum(int handle, const CkArrayIndex &index)
{
	CkArrayIndex2D idx2d = *(CkArrayIndex2D *) &index;
	if(maptable==NULL)
	{
		maptable= new CkHashtableT<intdual, int> (nstates*sizeY); 
		makemap();
	}
	return maptable->get(intdual(idx2d.index[0], idx2d.index[1]));
}

/** 
 * New functions beind added for the topology mapping of the
 * density objects - RhoR, RhoG, RhoGHartExt
 */

void RhoRSMap::makemap()
{
	FindProcessor fp;
	int x = CkNumPes();
	int y = 1;
	int z = 1;
	int vn = 0;
	
#ifdef CMK_VERSION_BLUEGENE
	BGLTorusManager *bgltm = BGLTorusManager::getObject();
	x = bgltm->getXSize();
	y = bgltm->getYSize();
	z = bgltm->getZSize();
	vn = bgltm->isVnodeMode();
#endif
	
	fp.count=1;
	fp.nopX=x;
	fp.nopY=y;
	fp.nopZ=z;
	
	int assign[3]={0, 0, 0};
	int w = 0;
	for(int i=0;i<3;i++)
		fp.start[i]=fp.next[i]=0;
                
        int rrsobjs_per_pe= nchareRhoR/(CkNumPes()/2);
        int rem;
        
        rem = nchareRhoR % (CkNumPes()/2);
        if(rem!=0)
          rrsobjs_per_pe += 1;
          
        for(int chunk=0; chunk<nchareRhoR; chunk+=rrsobjs_per_pe)
        {
          if(rem!=0)
            if(chunk==rem*rrsobjs_per_pe)
              rrsobjs_per_pe -= 1;
          if(chunk==0) {}
          else
          {
            for(int i=0;i<3;i++)
                  fp.start[i]=fp.next[i];
            if(fp.start[2]>x/2)
                  assign[2]=fp.start[2]-x;
            else
                  assign[2]=fp.start[2];
            if(fp.start[1]>y/2)
                  assign[1]=fp.start[1]-y;
            else
                  assign[1]=fp.start[1];
            if(fp.start[0]>z/2)
                  assign[0]=fp.start[0]-z;
            else
                  assign[0]=fp.start[0];
            if(vn==0)
            {
              fp.findNextInTorus(assign);
              for(int i=0;i<3;i++)
                  fp.start[i]=fp.next[i];
              if(fp.start[2]>x/2)
                  assign[2]=fp.start[2]-x;
              else
                  assign[2]=fp.start[2];
              if(fp.start[1]>y/2)
                  assign[1]=fp.start[1]-y;
              else
                  assign[1]=fp.start[1];
              if(fp.start[0]>z/2)
                  assign[0]=fp.start[0]-z;
              else
                  assign[0]=fp.start[0];
              fp.findNextInTorus(assign);
            }
            else
            {
              fp.findNextInTorusV(w, assign);
              w = fp.w;
              for(int i=0;i<3;i++)
                  fp.start[i]=fp.next[i];
              if(fp.start[2]>x/2)
                  assign[2]=fp.start[2]-x;
              else
                  assign[2]=fp.start[2];
              if(fp.start[1]>y/2)
                  assign[1]=fp.start[1]-y;
              else
                  assign[1]=fp.start[1];
              if(fp.start[0]>z/2)
                  assign[0]=fp.start[0]-z;
              else
                  assign[0]=fp.start[0];
              fp.findNextInTorusV(w, assign);
            }
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
              if(vn==0)
		maptable->put(intdual(i, 0))=fp.next[0]*x*y+fp.next[1]*x+fp.next[2];
              else
                maptable->put(intdual(i, 0))=(fp.next[0]*x*y+fp.next[1]*x+fp.next[2])*2+fp.w;
              //CkPrintf("%d on %d\n", i, fp.next[0]*x*y+fp.next[1]*x+fp.next[2]);
            }
          }   
        }
#ifdef MAP_DEBUG
	CkPrintf("RhoRSMap created on processor %d\n", CkMyPe());
#endif
}

int RhoRSMap::procNum(int arrayHdl, const CkArrayIndex &idx)
{
      CkArrayIndex2D idx2d = *(CkArrayIndex2D *) &idx;
      if(maptable==NULL)
      {
        maptable= new CkHashtableT<intdual, int> (nchareRhoR*1); 
	makemap();
      }  
      return maptable->get(intdual(idx2d.index[0], idx2d.index[1]));
}

void RhoGSMap::makemap()
{
	FindProcessor fp;
	int x = CkNumPes();
	int y = 1;
	int z = 1;
	int vn = 0;
	
#ifdef CMK_VERSION_BLUEGENE
	BGLTorusManager *bgltm = BGLTorusManager::getObject();
	x = bgltm->getXSize();
	y = bgltm->getYSize();
	z = bgltm->getZSize();
	vn = bgltm->isVnodeMode();
#endif
	
	fp.count=1;
	fp.nopX=x;
	fp.nopY=y;
	fp.nopZ=z;
	
        int assign[3]={0, 0, 0};
        int w;
        
        if(vn==0)
          assign[2]=1;
        else
          w=1;
        
	for(int i=0;i<3;i++)
		fp.start[i]=fp.next[i]=assign[i];
        
        int rgsobjs_per_pe= nchareRhoG/(CkNumPes()/2);
        int rem;
        
        rem = nchareRhoG % (CkNumPes()/2);
        if(rem!=0)
          rgsobjs_per_pe += 1;
        
        //CkPrintf("rem %d rgsobjs_per_pe %d\n", rem, rgsobjs_per_pe);
        for(int chunk=0; chunk<nchareRhoG; chunk+=rgsobjs_per_pe)
        {
          if(rem!=0)
            if(chunk==rem*rgsobjs_per_pe)
              rgsobjs_per_pe -= 1;  
          if(chunk==0) {}
          else
          {
            for(int i=0;i<3;i++)
                  fp.start[i]=fp.next[i];
            if(fp.start[2]>x/2)
                  assign[2]=fp.start[2]-x;
            else
                  assign[2]=fp.start[2];
            if(fp.start[1]>y/2)
                  assign[1]=fp.start[1]-y;
            else
                  assign[1]=fp.start[1];
            if(fp.start[0]>z/2)
                  assign[0]=fp.start[0]-z;
            else
                  assign[0]=fp.start[0];
            if(vn==0)
            {
              fp.findNextInTorus(assign);
              for(int i=0;i<3;i++)
                  fp.start[i]=fp.next[i];
              if(fp.start[2]>x/2)
                  assign[2]=fp.start[2]-x;
              else
                  assign[2]=fp.start[2];
              if(fp.start[1]>y/2)
                  assign[1]=fp.start[1]-y;
              else
                  assign[1]=fp.start[1];
              if(fp.start[0]>z/2)
                  assign[0]=fp.start[0]-z;
              else
                  assign[0]=fp.start[0];
              fp.findNextInTorus(assign);
            }
            else
            {
              fp.findNextInTorusV(w, assign);
              w = fp.w;
              for(int i=0;i<3;i++)
                  fp.start[i]=fp.next[i];
              if(fp.start[2]>x/2)
                  assign[2]=fp.start[2]-x;
              else
                  assign[2]=fp.start[2];
              if(fp.start[1]>y/2)
                  assign[1]=fp.start[1]-y;
              else
                  assign[1]=fp.start[1];
              if(fp.start[0]>z/2)
                  assign[0]=fp.start[0]-z;
              else
                  assign[0]=fp.start[0];
              fp.findNextInTorusV(w, assign);
            }
          }
          for(int i=chunk;i<chunk+rgsobjs_per_pe;i++)
          {
            if(chunk==0)
            {
              if(vn==0)
                maptable->put(intdual(i, 0))=assign[0]*x*y+assign[1]*x+assign[2];
              else
              maptable->put(intdual(i, 0))=(assign[0]*x*y+assign[1]*x+assign[2])*2+w;
              //CkPrintf("%d on %d\n", i, assign[0]*x*y+assign[1]*x+assign[2]);
            }
            else
            {
              if(vn==0)
		maptable->put(intdual(i, 0))=fp.next[0]*x*y+fp.next[1]*x+fp.next[2];
              else
                maptable->put(intdual(i, 0))=(fp.next[0]*x*y+fp.next[1]*x+fp.next[2])*2+fp.w;
              //CkPrintf("%d on %d\n", i, fp.next[0]*x*y+fp.next[1]*x+fp.next[2]);
            }
          } 
        }
#ifdef MAP_DEBUG
	CkPrintf("RhoGSMap created on processor %d\n", CkMyPe());
#endif
}
    
int RhoGSMap::procNum(int arrayHdl, const CkArrayIndex &idx)
{
      CkArrayIndex2D idx2d = *(CkArrayIndex2D *) &idx;
      if(maptable==NULL)
      {
        maptable= new CkHashtableT<intdual, int> (nchareRhoG*1); 
	makemap();
      }  
      return maptable->get(intdual(idx2d.index[0], idx2d.index[1]));
}

void RhoGHartMap::makemap()
{
	FindProcessor fp;
	int x = CkNumPes();
	int y = 1;
	int z = 1;
	int vn = 0;
	
#ifdef CMK_VERSION_BLUEGENE
	BGLTorusManager *bgltm = BGLTorusManager::getObject();
	x = bgltm->getXSize();
	y = bgltm->getYSize();
	z = bgltm->getZSize();
	vn = bgltm->isVnodeMode();
#endif
	
	fp.count=1;
	fp.nopX=x;
	fp.nopY=y;
	fp.nopZ=z;
	
        int assign[3]={0, 0, 0};
        int w;
        
        if(vn==0)
          assign[2]=1;
        else
          w=1;
          
	for(int i=0;i<3;i++)
		fp.start[i]=fp.next[i]=assign[i];
        
        int rghobjs_per_pe= nchareRhoGHart/(CkNumPes()/2);
        int rem;
        
        rem = nchareRhoGHart % (CkNumPes()/2);
        if(rem!=0)
          rghobjs_per_pe += 1;
          
        for(int chunk=0; chunk<nchareRhoGHart; chunk+=rghobjs_per_pe)
        {
          if(rem!=0)
            if(chunk==rem*rghobjs_per_pe)
              rghobjs_per_pe -= 1; 
          if(chunk==0) {}
          else
          {
            for(int i=0;i<3;i++)
                  fp.start[i]=fp.next[i];
            if(fp.start[2]>x/2)
                  assign[2]=fp.start[2]-x;
            else
                  assign[2]=fp.start[2];
            if(fp.start[1]>y/2)
                  assign[1]=fp.start[1]-y;
            else
                  assign[1]=fp.start[1];
            if(fp.start[0]>z/2)
                  assign[0]=fp.start[0]-z;
            else
                  assign[0]=fp.start[0];
            if(vn==0)
            {
              fp.findNextInTorus(assign);
              for(int i=0;i<3;i++)
                  fp.start[i]=fp.next[i];
              if(fp.start[2]>x/2)
                  assign[2]=fp.start[2]-x;
              else
                  assign[2]=fp.start[2];
              if(fp.start[1]>y/2)
                  assign[1]=fp.start[1]-y;
              else
                  assign[1]=fp.start[1];
              if(fp.start[0]>z/2)
                  assign[0]=fp.start[0]-z;
              else
                  assign[0]=fp.start[0];
              fp.findNextInTorus(assign);
            }
            else
            {
              fp.findNextInTorusV(w, assign);
              w = fp.w;
              for(int i=0;i<3;i++)
                  fp.start[i]=fp.next[i];
              if(fp.start[2]>x/2)
                  assign[2]=fp.start[2]-x;
              else
                  assign[2]=fp.start[2];
              if(fp.start[1]>y/2)
                  assign[1]=fp.start[1]-y;
              else
                  assign[1]=fp.start[1];
              if(fp.start[0]>z/2)
                  assign[0]=fp.start[0]-z;
              else
                  assign[0]=fp.start[0];
              fp.findNextInTorusV(w, assign);
            }
          }
          for(int i=chunk;i<chunk+rghobjs_per_pe;i++)
          {
            if(chunk==0)
            {
              if(vn==0)
                maptable->put(intdual(i, 0))=assign[0]*x*y+assign[1]*x+assign[2];
              else
                maptable->put(intdual(i, 0))=(assign[0]*x*y+assign[1]*x+assign[2])*2+w;
              //CkPrintf("%d on %d\n", i, assign[0]*x*y+assign[1]*x+assign[2]);
            }
            else
            {
              if(vn==0)
		maptable->put(intdual(i, 0))=fp.next[0]*x*y+fp.next[1]*x+fp.next[2];
              else
                maptable->put(intdual(i, 0))=(fp.next[0]*x*y+fp.next[1]*x+fp.next[2])*2+fp.w;
              //CkPrintf("%d on %d\n", i, fp.next[0]*x*y+fp.next[1]*x+fp.next[2]);
            }
          }  
        }
#ifdef MAP_DEBUG
	CkPrintf("RhoGHartMap created on processor %d\n", CkMyPe());
#endif
}

int RhoGHartMap::procNum(int arrayHdl, const CkArrayIndex &idx)
{
      CkArrayIndex2D idx2d = *(CkArrayIndex2D *) &idx;
      if(maptable==NULL)
      {
        maptable= new CkHashtableT<intdual, int> (nchareRhoGHart*1); 
	makemap();
      }  
      return maptable->get(intdual(idx2d.index[0], idx2d.index[1]));
}    
/**
 * End of the topology sensitive map functions.
 * Below is the alternative mapping scheme which was being used earlier.
 */
 
#else


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** 
 * The maps are intended to spread load around evenly and keep plane
 * and state objects which communicate across phases close together.

 * It basically uses this key idea pe = cumload /  average_load so
 * cumload depends on the object list traversal.  For  plane wise we
 * traverse it plane wise and for state wise do it states in the
 * first loop and planes in the second loop.  For hybrid map, it
 * partitions states and processors in to K partitions  
 * k = numpes/pesperstate so apply the cumulative load logic in each 
 * partition, and add the start pe of the partition to it

 * cumload is the cumulative load till that object average load
 * is the total_load/numpes

 * Assuming object s1,p1 is being inserted s1 is a state and p1 is a
 * plane the map will use a certain  traversal through object space
 * for state map it will traverse all  objects state by state hence
 * cumload = load_sum { s1 * numplanesperstate + p1} for gspace it
 * will be a loop as each plane has a different load.
 */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================

#include <math.h>
#include "ckarray.h"
#include "groups.h"
#include "sim_subroutines.h"
#include "CP_State_Plane.h"

//============================================================================

extern int sizeX;
extern Config config;
extern CProxy_CPcharmParaInfoGrp scProxy;

//void hackGSpacePlaneLoad(CPcharmParaInfo *sim,int idx, double *line_load, 
//                         double *pt_load);
void GSpacePlaneLoad(int , double *, double *);

//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void GSMap::GSpacePlaneLoad(int idx, double *line_load, double *pt_load){
//============================================================================

  if( (idx < 0) || (idx >= nchareG) ){
   CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
   CkPrintf("GSpace plane index %d out of range: 0 < idx > %d\n",idx,nchareG);
   CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
   CkExit();
  }//endif

  line_load[0] = lines_per_chareG[idx];
  pt_load[0]   =   pts_per_chareG[idx];

//============================================================================
  }//end routine
//============================================================================

void GSMap::makemap() {}

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void SCalcMap::GSpacePlaneLoad(int idx, double *line_load, double *pt_load){
//============================================================================


  if( (idx < 0) || (idx >= nchareG) ){
   CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
   CkPrintf("GSpace plane index %d out of range: 0 < idx > %d\n",idx,nchareG);
   CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
   CkExit();
  }//endif

  line_load[0] = lines_per_chareG[idx];
  pt_load[0]   =   pts_per_chareG[idx];

//============================================================================
  }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** 
 * Void allocate each state to more than one processor. Hybrid between 
 * state allocation and plane allocation
 */
//============================================================================
int basicMap(CkArrayIndex2D &idx2d, int numPlanes) {
//============================================================================
  int pes_per_state = config.RpesPerState;
  int np = CkNumPes()/pes_per_state;
 
  if(np < 1)
    np = 1;

  int partition_nstates = config.nstates / np;
  if(config.nstates % np != 0)
    partition_nstates ++;
  
  double total_planes = config.nstates * numPlanes;
  double numPerProc = total_planes / CkNumPes();
    
  //My state in the current partition
  int s = idx2d.index[0]%partition_nstates; 

  int pe = (int) floor((idx2d.index[1]*partition_nstates + s)/numPerProc); 
    
  pe = pe % pes_per_state;
  pe +=  pes_per_state *  (idx2d.index[0]/partition_nstates);

  if(pe >= CkNumPes())
    pe = pe % CkNumPes();

  return pe; 
//============================================================================
   }// end routine
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//int GSMap::slowprocNum(int arrayHdl, const CkArrayIndex2D &idx2d)
//============================================================================
int GSMap::procNum(int arrayHdl, const CkArrayIndex &idx){
//============================================================================

  CkArrayIndex2D idx2d = *(CkArrayIndex2D *) &idx;

  int pe        = 0;
  int numChareG = 0;
  
  if(config.doublePack) {
     numChareG = nchareG;
  }else{
     CkAbort("not doublepack broken\n");
  }//endif
  
  if(state_load <= 0.0 ) {
    for(int x = 0; x  < numChareG; x++) {
      double curload = 0.0;
      double sload = 0.0;
      
      GSpacePlaneLoad(x, &curload, &sload);
      
      state_load += curload;
    }//endfor
  }//endif

  int pes_per_state = config.GpesPerState;
  int np = CkNumPes()/pes_per_state;
  
  if(np < 1){np = 1;}
  
  int partition_nstates = config.nstates / np;
  if(config.nstates % np != 0){partition_nstates ++;}
  
  int start_pe = (idx2d.index[0]/partition_nstates) * pes_per_state;
  int start_state = (idx2d.index[0]/partition_nstates) * partition_nstates;
  
  double cum_load = 0.0;
  double average_load = state_load * config.nstates / CkNumPes();

  for(int x = 0; x < numChareG; x++) {
    for(int s = 0; s < partition_nstates; s++) {
      double curload = 0.0;
      double sload = 0.0;
      
      GSpacePlaneLoad(x, &curload, &sload);

      cum_load += curload;

      if((idx2d.index[0] == s + start_state) && idx2d.index[1] == x) {
	double dpe = 0.0;
	dpe = cum_load / average_load;

	pe = (int)dpe;
	pe = pe % pes_per_state;
	pe += start_pe;

	return (pe % CkNumPes());
      }//endif
    }//endfor
  }//endfor
  return (idx2d.index[0]*1037+idx2d.index[1])%CkNumPes();

//============================================================================
   }//end routine
//============================================================================

void RSMap::makemap() {}

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
int RSMap::procNum(int arrayHdl, const CkArrayIndex &idx){
    CkArrayIndex2D idx2d = *(CkArrayIndex2D *) &idx;

    int numPlanes = 0;
    numPlanes = sizeX;
    
    int pe  = 0;
    
    double total_planes = config.nstates * numPlanes;
    double planes_per_proc = total_planes / CkNumPes();
    
    pe = basicMap(idx2d, numPlanes);
    return pe;
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
inline double scalc_load(int x, int numplanes, double *load) {
  *load = numplanes * numplanes - x*x;
  return *load;
} 
//============================================================================


//============================================================================
/**
 * this one uses a lookup table built by calling the slow version
 */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
int SCalcMap::procNum(int hdl, const CkArrayIndex &idx){
    CkArrayIndex4D &idx4d = *(CkArrayIndex4D *)&idx;
    int intidx[2];
    CmiMemcpy(intidx,idx4d.index,2*sizeof(int));  // our 4 shorts are now 2 ints
    if(maptable==NULL){
	int calcchunks=max_states/gs;
	maptable= new CkHashtableT<intdual,int> (nchareG*calcchunks*calcchunks); 
                  // times blkSize, but its always 1
	makemap();
    }//endif
    return maptable->get(intdual(intidx[0], intidx[1]));
}
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void SCalcMap::makemap(){

  if(symmetric){
    for(int z= 0; z<numChunks; z++)
    for(int numX = 0; numX < nchareG; numX++){
      for (int s1 = 0; s1 < max_states; s1 += gs) {
	for (int s2 = s1; s2 < max_states; s2 += gs) {
	    CkArrayIndex4D idx4d(numX,s1,s2,z);
	    int intidx[2];
	    CmiMemcpy(intidx,idx4d.index,2*sizeof(int));  // our 4 shorts are now 2 ints
	    maptable->put(intdual(intidx[0],intidx[1]))=slowprocNum(0,idx4d);
#ifdef MAP_DEBUG
	    //CkPrintf("SYM: plane: %d x: %d y: %d pe %d\n", numX, s1, s2, slowprocNum(0,idx4d));
#endif
	}//endfor
      }//endfor
    }//endfor

  }else{

    for(int z= 0; z<numChunks; z++)
      for(int numX = 0; numX < nchareG; numX++){
	  for (int s1 = 0; s1 < max_states; s1 += gs) {
	      for (int s2 = 0; s2 < max_states; s2 += gs) {
		  CkArrayIndex4D idx4d(numX,s1,s2,z);
		  int intidx[2];
		  CmiMemcpy(intidx,idx4d.index,2*sizeof(int));  // our 4 shorts now 2 ints
		  maptable->put(intdual(intidx[0],intidx[1]))=slowprocNum(0,idx4d);
#ifdef MAP_DEBUG
		  //CkPrintf("ASYM: plane: %d x: %d y: %d pe %d\n", numX, s1, s2, slowprocNum(0,idx4d));
#endif
	      }//endfor
	  }//endfor
      }//endfor

  }//endif

//============================================================================
   }//end routine
//============================================================================

int SCalcMap::slowprocNum2(int hdl, const CkArrayIndex4D &idx4d){
  // Just use gspace as our guide for (w,x,y,z) use gsp(y+x/grainsize,w);
  return cheesyhackgsprocNum(scProxy.ckLocalBranch()->cpcharmParaInfo, idx4d.index[2]+idx4d.index[3]/gs,idx4d.index[0]);
}

/*
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
int SCalcMap::slowprocNum(int hdl, const CkArrayIndex4D &idx4d){
//============================================================================
//    CkArrayIndex4D &idx4d = *(CkArrayIndex4D *)&idx;
//    CkPrintf("scalc map call for [%d %d %d %d] on pe %d\n", 
//       idx4d.index[0],idx4d.index[1],idx4d.index[2],idx4d.index[3],CkMyPe());
//============================================================================
#ifndef CMK_OPTIMIZE
      double StartTime=CmiWallTimer();
#endif
    //Here maxY is the max number of planes;
    int planeid = idx4d.index[0];
    int numChareG = 0;

    if(config.doublePack){
        numChareG = config.nchareG;
    }else{
        numChareG = max_planes/2;
    }//endif

    //for asymetric
    double *load = new double[CkNumPes()];
    memset(load, 0, CkNumPes() * sizeof(double));

    int w=0, x=0, y = 0, z =0; 

    if(totalload <= 0.0) { 
      for(w = 0; w < numChareG; w++) {
      for(z = 0; z< numChunks; z++ )
	for(x = 0; x < max_states; x += gs) {
	  if (symmetric){
	    y = x;
	  }else{
	    y = 0;
	  }//endif
	  
	  for(; y < max_states; y += gs) {
	    
	    double curload = 0.0;
	    double gload = 0.0;
	    
	    //scalc_load(w, numChareG, &curload);
	    GSpacePlaneLoad(w, &gload, &curload);
	    
	    totalload += curload;
	  }//endfor
	}//endfor
    }//endif
    
    int pe = 0;
    
    for(w = 0; w < numChareG; w ++) {
    for(z= 0; z<numChunks; z++)
      for(x = 0; x < max_states; x += gs){
	if (symmetric){
	  y = x;
	}else{
	  y = 0;
	}//endif

	for(; y < max_states; y += gs) {
	  double curload = 0.0;
	  double gload = 0.0;
	  
	  //scalc_load(w, numChareG, &curload);
	  GSpacePlaneLoad(w, &gload, &curload);
	  
	  curload /=  totalload;
          
          if(load[pe] + curload > 0.3/CkNumPes()){pe ++;}
	  
	  if(pe >= CkNumPes()){pe = 0;}
          
	  load[pe] += curload;
	  
	  if((w == idx4d.index[0]) && (x == idx4d.index[1]) &&
	     (y == idx4d.index[2]) && (z == idx4d.index[3])) {
#ifndef CMK_OPTIMIZE
             traceUserBracketEvent(Scalcmap_, StartTime, CmiWallTimer());    
#endif
             delete [] load;
             return pe;
	  }//endif
        }//endfor
      }//endfor
    }//endfor : w

    delete [] load;
#ifndef CMK_OPTIMIZE
    traceUserBracketEvent(Scalcmap_, StartTime, CmiWallTimer());    
#endif

//============================================================================
    return (idx4d.index[0]*197+idx4d.index[1]*23+idx4d.index[2]*7
            +idx4d.index[3])%CkNumPes();    
//============================================================================
   }//end routine
//============================================================================
*/


int SCalcMap::slowprocNum(int hdl, const CkArrayIndex4D &idx4d)
{
  
    //Here maxY is the max number of planes;
    int planeid = idx4d.index[0];
    int numChareG = 0;


    if(config.doublePack){
        numChareG = config.nchareG;
    }else{
        numChareG = max_planes/2;
    }//endif

    //for asymetric
    double *load = new double[CkNumPes()];
    memset(load, 0, CkNumPes() * sizeof(double));

    int w=0, x=0, y = 0; 

    if(totalload <= 0.0) { 
      for(w = 0; w < numChareG; w ++) 
	for(x = 0; x < max_states; x += gs) {
	  if (symmetric)
	    y = x;
	  else
	    y = 0;
	  
	  for(; y < max_states; y += gs) {
	    
	    double curload = 0.0;
	    double gload = 0.0;
	    
	    //scalc_load(w, numChareG, &curload);
	    GSpacePlaneLoad(w, &gload, &curload);
	    
	    totalload += curload;
	  }
	}
    }
    
    int pe = 0;
    
    for(w = 0; w < numChareG; w ++)  
      for(x = 0; x < max_states; x += gs) {
	if (symmetric)
	  y = x;
	else
	  y = 0;
	
	for(; y < max_states; y += gs) {
	  double curload = 0.0;
	  double gload = 0.0;
	  
	  //scalc_load(w, numChareG, &curload);
	  GSpacePlaneLoad(w, &gload, &curload);
	  
	  curload /=  totalload;
          
          if(load[pe] + curload > 0.3/CkNumPes()) 
              pe ++;
	  
	  if(pe >= CkNumPes())
              pe = 0;
          
	  load[pe] += curload;
	  
	  if((w == idx4d.index[0]) && (x == idx4d.index[1]) &&
	     (y == idx4d.index[2])) {

              //if(CkMyPe() == 0)
              //  CkPrintf ("scalc %d %d %d %d assigned to pe %d and curload = %f, load = %f\n", w, x ,y, symmetric, pe, curload, load[pe]);
              
              delete [] load;
	    return pe;
	  }
        }
      }
    
    delete [] load;
    return (idx4d.index[0]*197+idx4d.index[1]*23+idx4d.index[2]*7+idx4d.index[3])%CkNumPes();    

}

int RhoRSMap::procNum(int arrayHdl, const CkArrayIndex &idx){
      CkArrayIndex2D idx2d = *(CkArrayIndex2D *) &idx;
      return (((N * idx2d.index[0]) + idx2d.index[1] + off) % CkNumPes());
    }

int RhoGSMap::procNum(int arrayHdl, const CkArrayIndex &idx){
      CkArrayIndex2D idx2d = *(CkArrayIndex2D *) &idx;
//      return (((N * idx2d.index[0]) + idx2d.index[1] + off) % CkNumPes());
      int pe=(((N * idx2d.index[0])  + off) % CkNumPes());
      // avoid PEs favored by the array characterized by the avoid and avoid_off parms
      if(avoid>1 && (pe-avoid_off)%avoid==0)
      {
	  pe=(pe+1) %CkNumPes();
	  
      }
      return pe;
    }

int RhoGHartMap::procNum(int arrayHdl, const CkArrayIndex &idx){
      CkArrayIndex2D idx2d = *(CkArrayIndex2D *) &idx;
//      return (((N * idx2d.index[0]) + idx2d.index[1] + off) % CkNumPes());
      int pe=(((N * idx2d.index[0])  + off) % CkNumPes());
      // avoid PEs favored by the array characterized by the avoid and avoid_off parms
      if(avoid>1 && (pe-avoid_off)%avoid==0)
      {
	  pe=(pe+1) %CkNumPes();
	  
      }
      return pe;
    }    
#endif

