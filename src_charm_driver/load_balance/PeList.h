/*****************************************************************************
 * $Source$
 * $Author$
 * $Date$
 * $Revision$
 *****************************************************************************/

/** \file PeList.h
 *  Author: Eric J Bohm
 *  Date Created: May 31, 2006
 *
 *  Supports lists of processors which can be sorted by various
 *  criteria, such as distance or latency.  + and - operators are
 *  supported provide union and disjoint operations.  
 * 
 *  Implemented on ckvec.
 *  Used to find the nearest available processor from a given processor.
 *  Uses TopoManager for hopcount sorting if available.
 */

/** \ingroup mapping
 *
 */
//@{
#ifndef _PELIST_H
#define _PELIST_H

#include <math.h>
#include "configure.h"
#include "TopoManager.h"
#include "cklists.h"

class PeList 
{
 public:
  int *TheList;
  int *sortIdx;
  int current;
  int size;
  bool sorted;  // if pe list is kept in sorted order we can bin search it
  PeList();  //default constructor

  // boxy constructor for Torus machines
  PeList(int boxX, int boxY, int boxZ, int order, int maxX, int maxY, int maxZ, int maxT);

  PeList(int _size): size(_size)
    {
      sorted=true;
      TheList = new int [size+1];
      sortIdx = new int [size+1];
      for(int i=0;i<size;i++)
	{
	  TheList[i]=i;
	  sortIdx[i]=i;
	}
      current=0;
    }


  PeList(CkVec <int> inlist)
    {
      inlist.quickSort();
      sorted=true;
      current=0;
      size=inlist.size();
      TheList = new int [size+1];
      sortIdx = new int [size+1];
      TheList[0]=inlist[0];
      int count=0;
      for(int i=1;i<inlist.size();i++)
	{ //remove duplicates from sorted list
	  if(inlist[i]!=TheList[count])
	    {
	      count++;
	      TheList[count]=inlist[i];
	      sortIdx[count]=count;
	    }
	  
	}
      size=count;
    }

  PeList(PeList &inlist)
    {
      size=inlist.size;
      sorted=inlist.sorted;
      TheList = new int [size+1];
      sortIdx = new int [size+1];
      for(int i=0;i<inlist.size;i++)
	{
	  TheList[i]=inlist.TheList[i];
	  sortIdx[i]=inlist.sortIdx[i];
	}
      current=inlist.current;
    }


  PeList(int _size, int *a)
    {
      sorted=false;
      current=0;
      size=_size;
      TheList=new int [size+1];
      sortIdx=new int [size+1];
      for(int i=0;i<size;i++)
	{
	  TheList[i]=a[i];
	  sortIdx[i]=i;
	}
    }	 // use an array to construct

  PeList(PeList *a, int start, int _size)
    {
      sorted=false;
      CkAssert(start<a->size);
      CkAssert(_size<=a->size);
      size=_size;
      current=0;
      TheList=new int [size+1];
      sortIdx=new int [size+1];
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
  
  inline void reindex(){
      for(int ri=0;ri<size;ri++)
	sortIdx[ri]=ri;
      current=0;
  } 
  inline int exists(int target)
    {
      for(int i=0;i<size;i++)
	if(TheList[i]==target)
	  return i;
      return -1;
    }
  void rebuild(); 

  void reset(){current=0;} 

  //! return a processor at minumum distance from source
  int  minDist(int source);

  void resort(){   sortSource(TheList[0]);   }; 

  void append(PeList &inlist)
  {
    // make array large enough for both, paste together
    int newsize=inlist.size+size;
    int *newlist= new int [newsize+1];
    int *newIndex= new int [newsize+1];
    int i=0;
    for(; i< size ; i++)
      {
	newlist[i]=TheList[i];
	newIndex[i]=sortIdx[i];
      }
    for(int j=0; (i< newsize && j<inlist.size) ; i++,j++)
      {
	newlist[i]=inlist.TheList[j];
	newIndex[i]=i;
      }
    size=newsize;
    delete [] TheList;
    delete [] sortIdx;
    TheList=newlist;
    sortIdx=newIndex;
  }
  
  bool binsearch(int pe);

  inline bool find(int pe)
    {
      bool found=false;
      if(sorted)
	{ // bin search
	  found=binsearch(pe);
	}
      else // brute
	{
	  int i=0;
	  while(!found && i< size)
	    {
	      if(TheList[i]==pe)
		{
		  found=true;
		}
	      i++;
	    }
	}
      return(found);
    }

  int *pelower_bound(int pe);

  inline void mergeOne(int pe)
  {
    
    // make array large enough for both, paste together
    if(sorted)
      { // maintain sort order by insertion
	if(size==1 && pe>TheList[0])
	  {
	    appendOne(pe);
	  }
	else{
	  int *loc=pelower_bound(pe);
	  int location=loc-TheList;

	  if(location>=size || loc==NULL)
	    { // 
	      appendOne(pe);
	    }
	  else if ((loc>TheList) && (loc[-1] ==pe))
	    {
	      // do nothing
	    }
	  else if(TheList[location] ==pe)
	    {
	      // do nothing
	    }
	  
	  else if((size>1) &&(location<size-1) && (TheList[location+1]==pe))
	    {
	      // do nothing
	    }
	  else
	    { // not already present

	      int newsize=size+1;
	      int *newlist= new int [newsize+1];
	      int *newIndex= new int [newsize+1];
	      if(loc!=TheList) // there are things before location
		CmiMemcpy(newlist,TheList,location*sizeof(int));
	      newlist[location]=pe;
	      CmiMemcpy(&newlist[location+1],loc,(size-(location))*sizeof(int));
	      // your index is shot
	      bzero(newIndex,(newsize+1)*sizeof(int));
	      delete [] TheList;
	      delete [] sortIdx;
	      TheList=newlist;
	      sortIdx=newIndex;
	      size=newsize;
	      current=0;
	    }
	}
      }
    else
      {
	bool found=find(pe);
	if(!found)
	  {
	    appendOne(pe);
	  }
      }
  }

  inline void appendOne(int pe)
    {
	int newsize=size+1;
	int *newlist= new int [newsize+1];
	int *newIndex= new int [newsize+1];
	CmiMemcpy(&newlist[0],&TheList[0],sizeof(int)*size);
	CmiMemcpy(&newIndex[0],&sortIdx[0],sizeof(int)*size);
	newlist[size]=pe;
	newIndex[size]=size;
	size=newsize;
	delete [] TheList;
	delete [] sortIdx;
	TheList=newlist;
	sortIdx=newIndex;
    }

  int findNext();       // return next available, increment liststart

  void sortSource(int srcPe); // implementation depends on BG/L or not

  PeList &operator=(PeList &inlist) {TheList=inlist.TheList; sortIdx=inlist.sortIdx;return *this;}   


  // need to rebuild your sortIdx after set union
  inline PeList &operator+(PeList &inlist) { 
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
    for(int j=0; (i< newsize &j<inlist.size); i++)
      {
	newlist[i]=inlist.TheList[j];
	newIndex[i]=i;
      }
    size=newsize;
    delete [] TheList;
    delete [] sortIdx;
    TheList=newlist;
    sortIdx=newIndex;
    return *this; 
  }
  
  // need to rebuild your sortIdx after unary set difference
  inline PeList &operator-(PeList &inlist) {
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

  inline void remove(int pos)
    {
      if(pos < size)
	{

	  for (int i=pos; i<size-1; i++)
	    {
	      TheList[i] = TheList[i+1];
	    }

	  //	  CmiMemcpy(&TheList[pos],&TheList[pos+1],(size-(pos+1))*sizeof(int));
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
	      int *newList = new int[newSize+1];
	      int *newSortIdx = new int[newSize+1];
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

//@}
#endif
