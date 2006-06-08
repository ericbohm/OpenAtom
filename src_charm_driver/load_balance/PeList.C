#include "charm++.h"
#include "PeList.h"
#ifdef CMK_VERSION_BLUEGENE
#include "bgltorus.h"
extern 	BGLTorusManager *bgltm;
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
