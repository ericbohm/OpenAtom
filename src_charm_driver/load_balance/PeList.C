#include "charm++.h"
#include "PeList.h"
#ifdef CMK_VERSION_BLUEGENE
#include "bgltorus.h"
extern 	BGLTorusManager *bgltm;
#endif

#ifdef CMK_VERSION_BLUEGENE
//! construct the list by iterating through boxes which are sub
//! partitions 
PeList::PeList(int boxX, int boxY, int boxZ) // boxy constructor
    {

      current=0;
      size=CkNumPes();
      int i=0;
      int maxX=bgltm->getXSize();
      int maxY=bgltm->getYSize();
      int maxZ=bgltm->getZSize();
      TheList= new int[size];
      sortIdx= new int[size];
      for(int x=0; x<maxX; x+=boxX) // new box  in X
	for(int y=0; y<maxY; y+=boxY) // new box in Y 
	  for(int z=0; z<maxZ; z+=boxZ) // new box in Z
	    {
	      // fill out this box
	      for(int bx=0;bx<boxX;bx++)
		for(int by=0;by<boxY;by++)
		  for(int bz=0;bz<boxZ;bz++)
		    {
		      sortIdx[i]=i;
		      TheList[i++]=bgltm->coords2rank(bx+x,by+y, bz+z);
		    }
	    }
    }

#endif

PeList::PeList() // default constructor
    {
      current=0;
      size=CkNumPes();
      TheList= new int[size];
      sortIdx= new int[size];
      for(int i=0;i<size;i++)
	{
	  TheList[i]=i;
	  sortIdx[i]=i;
	}
      sortSource(0);
    }

void PeList::rebuild() // default constructor
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
  // sort it using bgltm
  //  CkPrintf("PRE: sortIndexByHops\n");
  bgltm->sortIndexByHops(srcPe,TheList,sortIdx,size);
  //  CkPrintf("POST sortIndexByHops\n");

}
#else
// alternate form in non BG/L case (just an ordered list)
void PeList::sortSource(int srcPe)
{
  // sort it using CkVec quicksort
  CkVec <int> sortme(size);
  memcpy(sortme.getVec(), TheList,size*sizeof(int));
  sortme.quickSort();
  memcpy(TheList, sortme.getVec(), size*sizeof(int));
}
#endif
