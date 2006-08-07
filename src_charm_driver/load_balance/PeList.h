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
  PeList::PeList(int boxX, int boxY, int boxZ) // boxy constructor
#endif
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
      sortSource(TheList[0]);
    }

  PeList(int _size, int *a)
    {
      current=0;
      size=_size;
      TheList=new int [size];
      sortIdx=new int [size];
      memcpy(TheList,a,size*sizeof(int));
      sortSource(TheList[0]);
    };	 // use an array to construct
  
  // given the max grid/torus dimenions and the volume of a desired subpartition
  // return the cuboid x,y,z of a subpartition exactly matching that volume
  bool cuboid(int &x, int &y, int &z, int maxX, int maxY, int maxZ, int volume)
    {
      double cubert= cbrt((double) volume);
      int cubetrunc= (int) cubert;
      x=y=z=cubetrunc;
      if(cubetrunc>maxX)
	cubetrunc=maxX;
      if(cubetrunc>maxY)
	cubetrunc=maxY;
      if(cubetrunc>maxZ)
	cubetrunc=maxZ;
      if(volume==x*y*z)
	return true;
      int start=cubetrunc-1;
      for(x=cubetrunc; x<=maxX;x++)
	{
	  for(y=start; y<=maxY;y++)
	    {
	      for(z=start; z<=maxZ;z++)
		{
		  if(volume==x*y*z)
		    return true;
		}
	    }
	}
      return false;
    }


  inline bool noPes() 
    {
      return(size-current<1);
    }
  
  inline int count() { return(size-current);  }
  

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
  //  PeList &operator+(PeList &inlist) { CkAbort("not implemented");}
  /*    for(int i=0; i< inlist.TheList.size();i++)
	TheList.push_back(inlist.TheList[i]);
	return *this; }*/


  // need to rebuild your sortIdx
  PeList &operator-(PeList &inlist) {
    for(int i=0; i< inlist.size;i++)
      for(int j=0; j< size;j++)
	if(TheList[j]==inlist.TheList[i])
	  {
	    remove(j);
	    j--; //j is moving target
	  }
    return *this;
  }

  void remove(int pos)
    {
      if(pos < size)
	{
	  for (int i=pos; i<size-1; i++)
	    TheList[i] = TheList[i+1];
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
    delete [] TheList;
    delete [] sortIdx;

  }
};

#endif
