/** \file PeList.h
 *  Author: Eric J Bohm
 *  Date Created: May 31, 2006
 *
 *  Supports lists of processors which can be sorted by various
 *  criteria, such as distance or latency.  + and - operators are
 *  supported provide union and disjoint operations.  
 * 
 *  Implemented on ckvec.
 *
 *  Used to find the nearest available processor from a given processor.
 *  
 *  Uses bgltorus manager for hopcount sorting if available.
 */

#ifndef _PELIST_H
#define _PELIST_H

#include <math.h>

#ifdef CMK_VERSION_BLUEGENE
#include "bgltorus.h"
#endif
#include "cklists.h"

class PeList 
{
 public:
  CkVec <int> TheList;
  CkVec <int> sortIdx;
  int current;
  PeList();  //default constructor
  PeList(CkVec <int> inlist)
    {
      current=0;
      for(int i=0;i<inlist.size();i++)
	{
	  TheList.push_back(inlist[i]);
	  sortIdx.push_back(i);
	}
      sortSource(TheList[0]);
    }

  PeList(int size, int *a)
    {
      current=0;
      TheList=CkVec<int>(size);
      sortIdx=CkVec<int>(size);
      memcpy(TheList.getVec(),a,size*sizeof(int));
      sortSource(TheList[0]);
    };	 // use an array to construct
  inline bool noPes() 
    {
      return(TheList.size()-current<1);
    }
  
  inline int count() { return(TheList.size()-current);  }
  

  void rebuild(); 

  void reset(){current=0;} 

  void resort(){   sortSource(TheList[0]);   }; 

  inline int findNext()        // return next available, increment liststart
  {
    
#ifdef CMK_VERSION_BLUEGENE
    int value=TheList[sortIdx[current]]; 
    //    TheList.remove(sortIdx[0]);
    //    sortIdx.remove(0);
#else
    int value=TheList[current]; 
    //    TheList.remove(0);
#endif
    current++;
    return(value); 
  };						


  void sortSource(int srcPe); // implementation depends on BG/L or not

  // need to rebuild your sortIdx
  PeList &operator=(PeList &inlist) {TheList=inlist.TheList; sortIdx=inlist.sortIdx;return *this;}
  PeList &operator+(PeList &inlist) {
    for(int i=0; i< inlist.TheList.size();i++)
	TheList.push_back(inlist.TheList[i]);
    return *this;
  }

  // need to rebuild your sortIdx
  PeList &operator-(PeList &inlist) {
    for(int i=0; i< inlist.TheList.size();i++)
      for(int j=0; j< TheList.size();j++)
	if(TheList[j]==inlist.TheList[i])
	  TheList.remove(j);
    return *this;
  }
  
  void dump()
    {
      CkPrintf("PeList\n");
      for(int j=0; j< TheList.size();j++)
	CkPrintf("%d list %d sortidx %d \n",j,TheList[j], sortIdx[j]);
    }
};

#endif
