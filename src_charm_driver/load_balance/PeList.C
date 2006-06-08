#include "charm++.h"
#include "PeList.h"
#ifdef CMK_VERSION_BLUEGENE
#include "bgltorus.h"
extern 	BGLTorusManager *bgltm;
#endif


PeList::PeList() // default constructor
    {
      current=0;
      TheList=CkVec<int>(CkNumPes());
      for(int i=0;i<CkNumPes();i++)
	TheList[i]=i;
      sortIdx=CkVec<int>(CkNumPes());
      sortSource(TheList[0]);
    }

void PeList::rebuild() // default constructor
    {

      TheList.reserve(CkNumPes());
      TheList.resize(CkNumPes());
      for(int i=0;i<CkNumPes();i++)
	TheList[i]=i;
      sortIdx=CkVec<int>(CkNumPes());
      current=0;
    }

#ifdef CMK_VERSION_BLUEGENE
// BG/L specific PeList implementations

void PeList::sortSource(int srcPe)
{
  // sort it using bgltm
  //  CkPrintf("PRE: sortIndexByHops\n");
  bgltm->sortIndexByHops(srcPe,TheList.getVec(),sortIdx.getVec(),TheList.size());
  //  CkPrintf("POST sortIndexByHops\n");

}
#else
// alternate form in non BG/L case (just an ordered list)
void PeList::sortSource(int srcPe)
{
  // sort it using CkVec quicksort
  TheList.quickSort();
}
#endif
