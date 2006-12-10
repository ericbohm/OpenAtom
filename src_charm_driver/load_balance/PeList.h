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
  int *TheList;
  int *sortIdx;
  int current;
  int size;
  PeList();  //default constructor

#ifdef CMK_VERSION_BLUEGENE
 //! boxy constructor for BG/L
  PeList(int boxX, int boxY, int boxZ, int order);
#endif
  PeList(int _size): size(_size)
    {
      TheList = new int [size];
      sortIdx = new int [size];
      for(int i=0;i<size;i++)
	{
	  TheList[i]=i;
	  sortIdx[i]=i;
	}
      current=0;
    }
  PeList(CkVec <int> inlist)
    {
      current=0;
      size=inlist.size();
      TheList = new int [size];
      sortIdx = new int [size];
      for(int i=0;i<inlist.size();i++)
	{
	  TheList[i]=inlist[i];
	  sortIdx[i]=i;
	}
    }

  PeList(PeList &inlist)
    {
      size=inlist.size;
      TheList = new int [size];
      sortIdx = new int [size];
      for(int i=0;i<inlist.size;i++)
	{
	  TheList[i]=inlist.TheList[i];
	  sortIdx[i]=inlist.sortIdx[i];
	}
      current=inlist.current;
    }


  PeList(int _size, int *a)
    {
      current=0;
      size=_size;
      TheList=new int [size];
      sortIdx=new int [size];
      for(int i=0;i<size;i++)
	{
	  TheList[i]=a[i];
	  sortIdx[i]=i;
	}
    }	 // use an array to construct

  PeList(PeList *a, int start, int _size)
    {
      CkAssert(start<a->size);
      CkAssert(_size<=a->size);
      size=_size;
      current=0;
      TheList=new int [size];
      sortIdx=new int [size];
      for(int i=0;i<size;i++)
	{
	  TheList[i]=a->TheList[a->sortIdx[i+start]];
	}
      reindex();
    };	 // make a copy of a sublist
  
  // given the max grid/torus dimenions and the volume of a desired subpartition

  inline bool noPes() 
    {
      return(size-current<1);
    }
  
  inline int count() { return(size-current);  }
  
  void reindex(){
      for(int i=0;i<size;i++)
	{
	  sortIdx[i]=i;
	}
  } 

  void rebuild(); 

  void reset(){current=0;} 

  void resort(){   sortSource(TheList[0]);   }; 

  void append(PeList &inlist)
  {
    // make array large enough for both, paste together
    int newsize=inlist.size+size;
    int *newlist= new int [newsize];
    int *newIndex= new int [newsize];
    int i=0;
    for(; i< size ; i++)
      {
	newlist[i]=TheList[i];
	newIndex[i]=sortIdx[i];
      }
    for(; i< newsize ; i++)
      {
	newlist[i]=inlist.TheList[i];
	newIndex[i]=i;
      }
    size=newsize;
    delete TheList;
    delete sortIdx;
    TheList=newlist;
    sortIdx=newIndex;
  }

  void mergeOne(int pe)
  {
    // make array large enough for both, paste together
    int i=0;
    bool found=false;
    for(; i< size ; i++)
      {
	if(TheList[i]==pe)
	  found=true;
      }
    if(!found)
      {
	int newsize=size+1;
	int *newlist= new int [newsize];
	int *newIndex= new int [newsize];
	i=0;
	for(; i< size ; i++)
	  {
	    newlist[i]=TheList[i];
	    newIndex[i]=sortIdx[i];
	  }
	for(; i< newsize ; i++)
	  {
	    newlist[i]=pe;
	    newIndex[i]=i;
	  }
	size=newsize;
	delete TheList;
	delete sortIdx;
	TheList=newlist;
	sortIdx=newIndex;
      }
  }


  inline int findNext()        // return next available, increment liststart
    //int findNext()        // return next available, increment liststart 
  {
    
#ifdef CMK_VERSION_BLUEGENE
    if(current>=size)
      {
	//CkPrintf("hey why is current %d >= size %d\n",current, size);
	current=0;
      }
    CkAssert(current<size);
    CkAssert(sortIdx[current]<size);
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


  PeList &operator=(PeList &inlist) {TheList=inlist.TheList; sortIdx=inlist.sortIdx;return *this;}   


  // need to rebuild your sortIdx
  PeList &operator+(PeList &inlist) { 
    // make array large enough for both, paste together
    int newsize=inlist.size+size;
    int *newlist= new int [newsize];
    int *newIndex= new int [newsize];
    int i=0;
    for(; i< size ; i++)
      {
	newlist[i]=TheList[i];
	newIndex[i]=sortIdx[i];
      }
    for(; i< newsize ; i++)
      {
	newlist[i]=inlist.TheList[i];
	newIndex[i]=i;
      }
    size=newsize;
    delete TheList;
    delete sortIdx;
    TheList=newlist;
    sortIdx=newIndex;
    return *this; 
  }
  
  // need to rebuild your sortIdx
  PeList &operator-(PeList &inlist) {
    for(int i=0; i< inlist.size;i++)
      {
	int j=0;
	while(j< size)
	  if(TheList[j]==inlist.TheList[i])
	    {
	      remove(j);
	      //	      CkPrintf("removed %d containing %d leaving %d\n",i,inlist.TheList[i], size);

	    }
	  else
	    {
	      j++;
	    }
      }
    current=0;  // your index is useless now anyway
    return *this;
  }

  void remove(int pos)
    {
      if(pos < size)
	{
	  for (int i=pos; i<size-1; i++)
	    {
	      TheList[i] = TheList[i+1];
	    }
	  size--;
	}
    }

  void trimUsed()
    {
      // build a new array that doesn't include used elements
      // theoretically this could be done in place, but haste makes waste
      if(current>0 && size>0)
	{
	  if(noPes())
	    {
	      size=0;
	      current=0;
	    }
	  else
	    {
	      int newSize= size-current;
	      int *newList = new int[newSize];
	      int *newSortIdx = new int[newSize];
	      for(int i=0; i<newSize;i++)
		{
		  // if we just copy in sorted order we end up with a sorted list
		  newSortIdx[i]=i;
		  newList[i]=TheList[sortIdx[i+current]];
		}
	      delete [] TheList;
	      delete [] sortIdx;
	      TheList=newList;
	      sortIdx=newSortIdx;
	      current=0;
	      size=newSize;
	    }
	}
    }

  void removeIdx(int pos)
    {
      if(pos < size)
	{
	  for (int i=pos; i<size-1; i++)
	    {
	      sortIdx[i] = sortIdx[i+1];
	    }
	  size--;
	}
    }
  void dump()
    {
      CkPrintf("PeList\n");
      for(int j=0; j< size;j++)
	CkPrintf("%d list %d sortidx %d \n",j,TheList[j], sortIdx[j]);
    }
  ~PeList(){
    if(TheList !=NULL)
      delete [] TheList;
    if(sortIdx !=NULL)
      delete [] sortIdx;

  }


};

#endif
