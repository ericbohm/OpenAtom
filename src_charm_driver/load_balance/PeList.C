#include "charm++.h"
#include "PeList.h"
#ifdef CMK_VERSION_BLUEGENE
#include "TopoManager.h"
extern TopoManager *bgltm;
#endif

#ifdef CMK_VERSION_BLUEGENE
//! construct the list by iterating through boxes which are sub
//! partitions 
PeList::PeList(int boxX, int boxY, int boxZ, int order) // boxy constructor
    {

      current=0;

      size=CkNumPes();
      int i=0;
      int maxX=bgltm->getDimX();
      int maxY=bgltm->getDimY();
      int maxZ=bgltm->getDimZ();
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
		    for(int bz=0;bz<boxZ;bz++)
		      for(int by=0;by<boxY;by++) // make inner planes along X
			{
			  sortIdx[i]=i;
			  TheList[i++]=bgltm->coordinatesToRank(bx+x, by+y, bz+z);
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
			{
			  sortIdx[i]=i;
			  TheList[i++]=bgltm->coordinatesToRank(bx+x,by+y, bz+z);
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
			{
			  sortIdx[i]=i;
			  TheList[i++]=bgltm->coordinatesToRank(bx+x,by+y, bz+z);
			}
		}

	}
      else
	{
	  CkAbort("unknown order");
	}
      // size is actually boxsize times the number of whole boxes
      // that fit in the mesh
      int end=numBoxes* boxX*boxY*boxZ;
      // fill out remainder 
      if(i<CkNumPes())
	{
	  PeList remainder(CkNumPes());
	  for(int i=0; i< size;i++)
	    {
	      int j=0;
	      while(j< size)
		if(remainder.TheList[j]==TheList[i])
		  {
		    remainder.remove(j);
		  }
		else
		  {
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

#endif



PeList::PeList() // default constructor
    {
      current=0;
      size=CkNumPes();
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
  current=0;
  size=CkNumPes();
  for(int i=0;i<size;i++)
    {
      TheList[i]=i;
      sortIdx[i]=i;
    }

}


#ifdef CMK_VERSION_BLUEGENE
// BG/L specific PeList implementations

void PeList::sortSource(int srcPe)
{
  // sort it using TopoManager 
  //  CkPrintf("PRE: sortIndexByHops\n");
  CkAssert(srcPe>=0);
  CkAssert(srcPe<CkNumPes());
  bgltm->sortRanksByHops(srcPe, TheList, sortIdx, size);
  //  CkPrintf("POST sortIndexByHops\n");

}
#else
// alternate form in non BG/L case (just an ordered list)
void PeList::sortSource(int srcPe)
{
  // sort it using CkVec quicksort
  CkVec <int> sortme(size);
  CmiMemcpy(sortme.getVec(), TheList,size*sizeof(int));
  sortme.quickSort();
  CmiMemcpy(TheList, sortme.getVec(), size*sizeof(int));
}
#endif
