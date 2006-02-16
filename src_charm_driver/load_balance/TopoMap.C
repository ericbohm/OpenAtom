/** \file TopoMap.C
 * The functions in this file generate a topology sensitive
 * mapping scheme for the GSpace, PairCalc and RealSpace objects.
 * 
 */
#include "charm++.h"
#include "FindProcessor.h"
#include "cpaimd.h"

#ifdef CMK_VERSION_BLUEGENE
#include "bgltorus.h"

#endif

#include "util.h" 
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
	
#if CMK_VERSION_BLUEGENE
	BGLTorusManager *bgltm = BGLTorusManager::getObject();
	x = bgltm->getXSize();
	y = bgltm->getYSize();
	z = bgltm->getZSize();
#endif
	
	int procs[x][y][z];      // get from BGLTorusManager
	for(int i=0; i<x; i++)
	    for(int j=0; j<y; j++)
		for(int k=0; k<z; k++)
		    procs[i][j][k]=0;
		    
	//CkPrintf("x %d, y %d, z %d no of pe's %d\n", x, y, z, CkNumPes());
	//char fname[100];
	//sprintf(fname, "proc%d", CkMyPe()); 
	//FILE *f = fopen(fname, "w");
		
	fp.count=1;
	fp.nopX=x;
	fp.nopY=y;
	fp.nopZ=z;
	//CkPrintf("fp -> %d %d %d\n", fp.nopX, fp.nopY, fp.nopZ);
	int assign[3]={0, 0, 0};
	for(int i=0;i<3;i++)
		fp.start[i]=fp.next[i]=0;
	int gsobjs_per_pe;
	
	//CkPrintf("[%d] nstates %d, nchareG %d\n", CkMyPe(), nstates, nchareG);
	
	if((nstates*nchareG) % CkNumPes() == 0)
	    gsobjs_per_pe = (nstates*nchareG)/CkNumPes();
	else
	    gsobjs_per_pe = (nstates*nchareG)/CkNumPes()+1;
	int l=states_per_pe;
	int m;

	/*while(gsobjs_per_pe%l!=0)
		l++;
	// l--;                     l now divides gobjs_per_pe exactly 
	m = gsobjs_per_pe/l;  // each chunk will be of size l states by m planes*/
	
	if(gsobjs_per_pe%l==0)
		m = gsobjs_per_pe/l;
	else
		m = gsobjs_per_pe/l + 1;
	planes_per_pe=m;
	
	//CkPrintf("gsobjs_per_pe %d, l %d, m %d\n", gsobjs_per_pe, l, m);
	
	for(int ychunk=0; ychunk<nchareG; ychunk=ychunk+m)
		for(int xchunk=0; xchunk<nstates; xchunk=xchunk+l)
		{
			if(xchunk==0 && ychunk==0) {}
			else
			{
				procs[fp.next[2]][fp.next[1]][fp.next[0]]=1;
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
				
				while(procs[fp.next[2]][fp.next[1]][fp.next[0]]==1)
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
					fp.findNextInTorus(assign);
				}
				
			}
			for(int state=xchunk; state<xchunk+l && state<nstates; state++)
			{
				for(int plane=ychunk; plane<ychunk+m && plane<nchareG; plane++)
				{
					if(xchunk==0 && ychunk==0)
					{
						maptable->put(intdual(state, plane))=0;
						//CkPrintf("%d %d on 0\n", state, plane);
					}
					else
					{
						maptable->put(intdual(state, plane))=fp.next[0]*x*y+fp.next[1]*x+fp.next[2];
						//CkPrintf("%d %d on %d\n", state, plane, fp.next[0]*x*y+fp.next[1]*x+fp.next[2]);
						
					}
				}
			}
		}
	CkPrintf("GSMap created on processor %d\n", CkMyPe());
	//fclose(f);
	/*for(int i=0; i<nstates; i++)
		for(int j=0; j<nchareG; j++)
			maptable->put(intdual(i, j))=0;*/
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
	
#if CMK_VERSION_BLUEGENE
	BGLTorusManager *bgltm = BGLTorusManager::getObject();
	x = bgltm->getXSize();
	y = bgltm->getYSize();
	z = bgltm->getZSize();
#endif

	int procs[x][y][z];      // get from BGLTorusManager
	for(int i=0; i<x; i++)
		for(int j=0; j<y; j++)
			for(int k=0; k<z; k++)
				procs[i][j][k]=0;
			
	fp.count=1;
	fp.nopX=x;
	fp.nopY=y;
	fp.nopZ=z;
		
	if(symmetric)
	{
		for(int i=1; i<=max_states/grainsize; i++)
			lesser_scalc += i;
		if(lesser_scalc*nchareG % CkNumPes() == 0)
			scobjs_per_pe = lesser_scalc*nchareG/CkNumPes();
		else
			scobjs_per_pe = lesser_scalc*nchareG/CkNumPes() + 1;
		//CkPrintf("scobjs_per_pe %d grainsize %d nchareG %d scalc_per_plane %d planes_per_pe %d\n", scobjs_per_pe, grainsize, nchareG, scalc_per_plane, planes_per_pe);
				
		int assign[3]={0, 0, 0};
		for(int i=0;i<3;i++)
			fp.start[i]=fp.next[i]=0;
			
		for(int pchunk=0; pchunk<nchareG; pchunk=pchunk+planes_per_pe)
			for(int xchunk=0; xchunk<max_states; xchunk=xchunk+grainsize)
				for(int ychunk=xchunk; ychunk<max_states; ychunk=ychunk+grainsize)
					for(int plane=pchunk; plane<pchunk+planes_per_pe && plane<nchareG; plane++)
					{
						CkArrayIndex4D idx4d(plane, xchunk, ychunk, 0);
						CmiMemcpy(intidx,idx4d.index,2*sizeof(int));  // our 4 shorts are now 2 ints
						if(xchunk==0 && ychunk==0 && plane==0)
						{
							//CkPrintf("plane %d x %d y %d = proc 0\n", plane, xchunk, ychunk); 
							maptable->put(intdual(intidx[0], intidx[1]))=0;
							count++;
						}
						else
						{
							if(count<scobjs_per_pe)
							{
								//CkPrintf("plane %d x %d y %d = proc %d\n", plane, xchunk, ychunk, assign[0]*x*y+assign[1]*x+assign[2]);
								maptable->put(intdual(intidx[0], intidx[1]))=assign[0]*x*y+assign[1]*x+assign[2];
								count++;
							}
							else
							{
								count=0;
								procs[assign[2]][assign[1]][assign[0]]=1;
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
								while(procs[fp.next[2]][fp.next[1]][fp.next[0]]==1)
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
									fp.findNextInTorus(assign);
								}
								for(int i=0;i<3;i++)
									assign[i]=fp.next[i];
								//CkPrintf("plane %d x %d y %d = proc %d\n", plane, xchunk, ychunk, assign[0]*x*y+assign[1]*x+assign[2]);
								maptable->put(intdual(intidx[0], intidx[1]))=assign[0]*x*y+assign[1]*x+assign[2];
								count++;
							}
						}
					}
	CkPrintf("Symmetric SCalcMap created on processor %d\n", CkMyPe());
	}
	else
	{
		if(n*nchareG % CkNumPes() == 0)
		scobjs_per_pe = n*nchareG/CkNumPes();
		else
		scobjs_per_pe = n*nchareG/CkNumPes() + 1;
		//CkPrintf("scobjs_per_pe %d grainsize %d nchareG %d scalc_per_plane %d planes_per_pe %d\n", scobjs_per_pe, grainsize, nchareG, scalc_per_plane, planes_per_pe);
			
		int assign[3]={0, 0, 0};
		for(int i=0;i<3;i++)
			fp.start[i]=fp.next[i]=0;
			
		for(int pchunk=0; pchunk<nchareG; pchunk=pchunk+planes_per_pe)
			for(int xchunk=0; xchunk<max_states; xchunk=xchunk+grainsize)
				for(int ychunk=0; ychunk<max_states; ychunk=ychunk+grainsize)
					for(int plane=pchunk; plane<pchunk+planes_per_pe && plane<nchareG; plane++)
					{
						CkArrayIndex4D idx4d(plane, xchunk, ychunk, 0);
						CmiMemcpy(intidx,idx4d.index,2*sizeof(int));  // our 4 shorts are now 2 ints
						if(xchunk==0 && ychunk==0 && plane==0)
						{
							//CkPrintf("plane %d x %d y %d = proc 0\n", plane, xchunk, ychunk); 
							maptable->put(intdual(intidx[0], intidx[1]))=0;
							count++;
						}
						else
						{
							if(count<scobjs_per_pe)
							{
								//CkPrintf("plane %d x %d y %d = proc %d\n", plane, xchunk, ychunk, assign[0]*x*y+assign[1]*x+assign[2]);
								maptable->put(intdual(intidx[0], intidx[1]))=assign[0]*x*y+assign[1]*x+assign[2];
								count++;
							}
							else
							{
								count=0;
								procs[assign[2]][assign[1]][assign[0]]=1;
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
								while(procs[fp.next[2]][fp.next[1]][fp.next[0]]==1)
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
									fp.findNextInTorus(assign);
								}
								for(int i=0;i<3;i++)
									assign[i]=fp.next[i];
								//CkPrintf("plane %d x %d y %d = proc %d\n", plane, xchunk, ychunk, assign[0]*x*y+assign[1]*x+assign[2]);
								maptable->put(intdual(intidx[0], intidx[1]))=assign[0]*x*y+assign[1]*x+assign[2];
								count++;
							}
						}
					}
	CkPrintf("Asymmetric SCalcMap created on processor %d\n", CkMyPe());
	}
	
	/*if(symmetric)
	{
		for(int i=0; i<nchareG; i++)
			for(int j=0; j<max_states; j=j+grainsize)
				for(int k=j; k<max_states; k=k+grainsize)
				{
					CkArrayIndex4D idx4d(i, j, k, 0);
					CmiMemcpy(intidx,idx4d.index,2*sizeof(int));
					maptable->put(intdual(intidx[0], intidx[1]))=0;
				}
	}
	else
	{
		for(int i=0; i<nchareG; i++)
			for(int j=0; j<max_states; j=j+grainsize)
				for(int k=0; k<max_states; k=k+grainsize)
				{
					CkArrayIndex4D idx4d(i, j, k, 0);
					CmiMemcpy(intidx,idx4d.index,2*sizeof(int));
					maptable->put(intdual(intidx[0], intidx[1]))=0;
				}
	}*/
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
	
#ifdef CMK_VERSION_BLUEGENE
	BGLTorusManager *bgltm = BGLTorusManager::getObject();
	x = bgltm->getXSize();
	y = bgltm->getYSize();
	z = bgltm->getZSize();
#endif

	int procs[x][y][z];
	for(int i=0; i<x; i++)
	    for(int j=0; j<y; j++)
		for(int k=0; k<z; k++)
		    procs[i][j][k]=0;
	
	fp.count=1;
	fp.nopX=x;
	fp.nopY=y;
	fp.nopZ=z;
	int assign[3]={0, 0, 0};
	for(int i=0;i<3;i++)
		fp.start[i]=fp.next[i]=0;
	
	int rsobjs_per_pe;
	if(nstates*nchareG % CkNumPes() == 0)
	    rsobjs_per_pe = nstates*nchareG/CkNumPes();
	else
	    rsobjs_per_pe = nstates*nchareG/CkNumPes()+1;
	int count=0;
	
	//CkPrintf("[%d] nstates %d, nchareG %d\n", CkMyPe(), nstates, nchareG);
	//CkPrintf("rsobjs_per_pe %d\n", rsobjs_per_pe);
	
	for(int state=0; state<nstates; state++)
		for(int plane=0; plane<nchareG; plane++)
		{
			if(plane==0 && state==0)
			{
			        //CkPrintf("state %d plane %d = proc 0\n", state, plane);
				maptable->put(intdual(state, plane))=0;
				count++;
			}
			else
			{
				if(count<rsobjs_per_pe)
				{
				        //CkPrintf("state %d plane %d = proc %d\n", state, plane, fp.next[0]*x*y+fp.next[1]*x+fp.next[2]);
					maptable->put(intdual(state, plane))=assign[0]*x*y+assign[1]*x+assign[2];
					count++;
				}
				else
				{
					count=0;
					procs[assign[2]][assign[1]][assign[0]]=1;

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
				
					while(procs[fp.next[2]][fp.next[1]][fp.next[0]]==1)
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
						fp.findNextInTorus(assign);
					}
					for(int i=0;i<3;i++)
						assign[i]=fp.next[i];
					//CkPrintf("state %d plane %d = proc %d\n", state, plane, fp.next[0]*x*y+fp.next[1]*x+fp.next[2]);
					maptable->put(intdual(state, plane))=assign[0]*x*y+assign[1]*x+assign[2];
					count++;
				}
			}
		}
	CkPrintf("RSMap created on processor %d\n", CkMyPe());
	/*for(int i=0; i<nstates; i++)
		for(int j=0; j<nchareG; j++)
			maptable->put(intdual(i, j))=0;*/
}

/**
 *
 */
int RSMap::procNum(int handle, const CkArrayIndex &index)
{
	CkArrayIndex2D idx2d = *(CkArrayIndex2D *) &index;
	if(maptable==NULL)
	{
		maptable= new CkHashtableT<intdual, int> (nstates*nchareG); 
		makemap();
	}
	return maptable->get(intdual(idx2d.index[0], idx2d.index[1]));
}

