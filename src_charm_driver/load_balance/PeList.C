/*****************************************************************************
 * $Source$
 * $Author$
 * $Date$
 * $Revision$
 *****************************************************************************/

/** \file PeList.C
 *
 */


#include "charm++.h"
#include "PeList.h"
#include <algorithm>
#include "TopoManager.h"
extern TopoManager *topoMgr;
extern Config config;

/**
 * construct the list by iterating through boxes which are sub partitions
 * boxy constructor
 */
PeList::PeList(int boxX, int boxY, int boxZ, int order, int maxX, int maxY, int maxZ, int maxT) 
{
  if(config.torusMap==1)
  {
    CkAssert(topoMgr!=NULL);
    CkAssert(boxX>0);
    CkAssert(boxY>0);
    CkAssert(boxZ>0);
      current=0;
      sorted=false;
      size=config.numPesPerInstance;
      int i=0;
      CkPrintf("Ordering processors along long axis %d in %d X %d X %d\n",order, maxX, maxY, maxZ);
      TheList= new int[size];
      sortIdx= new int[size];
      int numBoxes=0;
      if(order==0)  // long axis along X
	{
	  for(int x=0; x<maxX; x+=boxX) // new box  in X
	    for(int y=0; y<maxY; y+=boxY) // new box in Y 
	      for(int z=0; z<maxZ; z+=boxZ) // new box in Z
		{
		  // fill out this box
		  numBoxes++;

		  for(int bx=0;bx<boxX;bx++)
		    for(int by=0;by<boxY;by++) // make inner planes along X
		      for(int bz=0;bz<boxZ;bz++)
			for(int bt=0;bt<maxT;bt++)
			{
			  sortIdx[i]=i;
			  // CkPrintf("i %d bx %d x %d by %d y bz %d z bt %d\n",i,bx,x,by,y,bz,z,bt);
			  TheList[i++]=topoMgr->coordinatesToRank(bx+x, by+y, bz+z, bt);
			}
		}
	}
      else if(order ==1) // long axis is along Y
	{
	  for(int y=0; y<maxY; y+=boxY) // new box in Y 
	    for(int x=0; x<maxX; x+=boxX) // new box  in X
	      for(int z=0; z<maxZ; z+=boxZ) // new box in Z
		{
		  // fill out this box
		  numBoxes++;


		  for(int by=0;by<boxY;by++)
		    for(int bz=0;bz<boxZ;bz++)
		      for(int bx=0;bx<boxX;bx++)
			for(int bt=0;bt<maxT;bt++)
			{
			  sortIdx[i]=i;
			  TheList[i++]=topoMgr->coordinatesToRank(bx+x, by+y, bz+z, bt);
			}
		}

	}
      else if(order ==2) // long axis is along Z
	{
	  for(int z=0; z<maxZ; z+=boxZ) // new box in Z
	    for(int x=0; x<maxX; x+=boxX) // new box  in X
	      for(int y=0; y<maxY; y+=boxY) // new box in Y 
		{
		  // fill out this box
		  numBoxes++;
		  for(int bz=0;bz<boxZ;bz++)
		    for(int by=0;by<boxY;by++)
		      for(int bx=0;bx<boxX;bx++)
			for(int bt=0;bt<maxT;bt++)
			{
			  sortIdx[i]=i;
			  TheList[i++]=topoMgr->coordinatesToRank(bx+x, by+y, bz+z, bt);
			}
		}

	}
      else
	{
	  CkAbort("unknown order");
	}
      // size is actually boxsize times the number of whole boxes
      // that fit in the mesh
      int end=numBoxes*boxX*boxY*boxZ;
      // fill out remainder 
      if(i<config.numPesPerInstance)
	{
	  PeList remainder(config.numPesPerInstance);
	  for(int i=0; i< size;i++)
	    {
	      int j=0;
	      while(j< size)
		if(remainder.TheList[j]==TheList[i])
		  {
		    CkPrintf("IF %d %d %d %d %d\n", i, j, size, TheList[j], TheList[i]);
		    remainder.remove(j);
		  }
		else
		  {
		    CkPrintf("ELSE %d %d %d %d %d\n", i, j, size, TheList[j], TheList[i]);
		    j++;
		  }
	    }
	  remainder.reindex();
	  // now we just plunk these at the end
	  int i=0;
	  for(; i< remainder.size ; i++)
	    {
	      TheList[i+end]=remainder.TheList[i];
	      sortIdx[i+end]=i+end;
	    }
	}
  }
  else
    CkAbort("You shouldn't be calling this function\n");
}


PeList::PeList() // default constructor
    {
      sorted=true;
      current=0;
      size=config.numPesPerInstance;
      TheList= new int[size+1];
      sortIdx= new int[size+1];
      for(int i=0;i<size;i++)
	{
	  TheList[i]=i;
	  sortIdx[i]=i;
	}
      sortSource(0);
    }

void PeList::rebuild()
{
  sorted=true;
  current=0;
  size=config.numPesPerInstance;
  for(int i=0;i<size;i++)
    {
      TheList[i]=i;
      sortIdx[i]=i;
    }

}


int *PeList::pelower_bound(int pe)
{
  return(std::lower_bound(&TheList[0],&TheList[size],pe));
}
bool PeList::binsearch(int pe)
{
  return(std::binary_search(&TheList[0],&TheList[size],pe));
}

void PeList::sortSource(int srcPe)
{
  if(config.torusMap==1) 
    {
      // sort it using TopoManager 
      //  CkPrintf("PRE: sortIndexByHops\n");
      CkAssert(srcPe>=0);
      CkAssert(srcPe<config.numPes);
      topoMgr->sortRanksByHops(srcPe, TheList, sortIdx, size);
      //  CkPrintf("POST sortIndexByHops\n");
    }
  else 
    {
      // alternate form in non BG/L case (just an ordered list)
      // sort it using CkVec quicksort
      CkVec <int> sortme(size);
      CmiMemcpy(sortme.getVec(), TheList,size*sizeof(int));
      sortme.quickSort();
      CmiMemcpy(TheList, sortme.getVec(), size*sizeof(int));
    }
}

int PeList::minDist(int srcPe)
{
  if(config.torusMap==1)
    return(TheList[topoMgr->pickClosestRank(srcPe, TheList, size)]);
  else
    return(TheList[0]);
  // if not a torus distance is meaningless we just pick the first element
}

int PeList::findNext()        // return next available, increment liststart
{
    int value; 
    if(config.torusMap==1) { 
    if(current>=size)
      {
	//CkPrintf("hey why is current %d >= size %d\n",current, size);
	current=0;
      }
    CkAssert(current<size);
    CkAssert(sortIdx[current]<size);
    value=TheList[sortIdx[current]]; 
    //    TheList.remove(sortIdx[0]);
    //    sortIdx.remove(0);
    }
    else {
    value=TheList[current]; 
    //    TheList.remove(0);
    }
    current++;
    return(value); 
}

