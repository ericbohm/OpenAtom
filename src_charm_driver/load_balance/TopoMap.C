/** \file TopoMap.C
 *  Author: Abhinav S Bhatele
 *  Date Created: January 04th, 2006 
 *  The functions in this file generate a topology sensitive
 *  mapping scheme for the GSpace, PairCalc and RealSpace objects.
 *  and also for the density objects - RhoR, RhoG, RhoGHart.
 *
 *  Heavily refactored by EJB 2006/5/31 to accelerate map creation
 *
 *  !!!!!!   THIS CODE IS NOW DEFUNCT  !!!!!!!
 *  !!!!!!   THIS CODE IS NOW DEFUNCT  !!!!!!!
 *  !!!!!!   THIS CODE IS NOW DEFUNCT  !!!!!!!

 *  Look in MapTable.[Ch] for actual maps

 */

#include "charm++.h"
#include "main/cpaimd.h"
#include "utility/util.h"
#include "PeList.h"

/** 
 * Function for GSpace objects
 * Does a hybrid chunking of the GSpace grid. Depending on the 
 * no. of processors and the input parameter, Gstates_per_pe (l),
 * the other dimension 'm' is decided. And then this chunk of
 * size l by m is placed on a processor.   
 */

int GSMap::procNum(int handle, const CkArrayIndex &iIndex)
{
#if CMK_TRACE_ENABLED
  double StartTime=CmiWallTimer();
#endif

  int *index=(int *) iIndex.data();
#if 0
  if(maptable==NULL)
  {
    CkPrintf("Warning! GSMap::procnum had to assign maptable on pe %d!\n",CkMyPe());
#ifdef USE_INT_MAP
    maptable= &GSImaptable;
#else
    maptable= &GSmaptable;
#endif
  }
#endif

#ifdef USE_INT_MAP
  int retval=maptable->get(index[0],index[1]);
#else
  int retval=maptable->get(intdual(index[0],index[1]));
#endif
#if CMK_TRACE_ENABLED
  traceUserBracketEvent(10000, StartTime, CmiWallTimer());
#endif

  return retval;
}

/** Function for PairCalc objects
 *
 */

/**
 *
 */
int SCalcMap::procNum(int handle, const CkArrayIndex &iIndex)
{
#if CMK_TRACE_ENABLED
  double StartTime=CmiWallTimer();
#endif


#if 0
  /* this works due to 4D being stored as 2 ints */
  if(maptable==NULL)
  {
    CkPrintf("Warning! SCalc::Procnum had to assign maptable on %d !\n",CkMyPe());
#ifdef USE_INT_MAP
    if(symmetric)
      maptable= &SymScalcImaptable;
    else
      maptable= &AsymScalcImaptable;
#else
    if(symmetric)
      maptable= &SymScalcmaptable;
    else
      maptable= &AsymScalcmaptable;
#endif
  }
#endif
#ifdef USE_INT_MAP
  short *sindex=(short *) iIndex.data();
  int retval=maptable->get(sindex[0], sindex[1], sindex[2], sindex[3]);
#else
  int *index=(int *) iIndex.data();
  int retval=maptable->get(intdual(index[0], index[1]));
#endif
#if CMK_TRACE_ENABLED
  traceUserBracketEvent(30000, StartTime, CmiWallTimer());
#endif
  return retval;
}

/** Function for RealSpace objects
 *
 */

int RSMap::procNum(int handle, const CkArrayIndex &iIndex)
{
#if CMK_TRACE_ENABLED
  double StartTime=CmiWallTimer();
#endif

  int *index=(int *) iIndex.data();
#if 0
  if(maptable==NULL)
  {
    CkPrintf("Warning! RSMap::Procnum had to assign maptable on %d!\n",CkMyPe());
#ifdef USE_INT_MAP
    maptable= &RSImaptable;
#else
    maptable= &RSmaptable;
#endif
  }
#endif

#ifdef USE_INT_MAP
  int retval=maptable->get(index[0], index[1]);
#else
  int retval=maptable->get(intdual(index[0], index[1]));
#endif
#if CMK_TRACE_ENABLED
  traceUserBracketEvent(20000, StartTime, CmiWallTimer());
#endif

  return retval;

}

int RSPMap::procNum(int handle, const CkArrayIndex &index)
{
  CkArrayIndex2D idx2d = *(CkArrayIndex2D *) &index;

#if 0
  if(maptable==NULL)
  {
    CkPrintf("Warning! RSMap::Procnum had to assign maptable on %d!\n",CkMyPe());
#ifdef USE_INT_MAP
    maptable= &RSPImaptable;
#else
    maptable= &RSPmaptable;
#endif
  }
#endif 

#ifdef USE_INT_MAP
  int retval=maptable->get(idx2d.index[0], idx2d.index[1]);
#else
  int retval=maptable->get(intdual(idx2d.index[0], idx2d.index[1]));
#endif
  return retval;

}

/** 
 * New functions beind added for the topology mapping of the
 * density objects - RhoR, RhoG, RhoGHartExt
 */

int RhoRSMap::procNum(int arrayHdl, const CkArrayIndex &iIndex)
{
  int *index=(int *) iIndex.data();

  if(maptable==NULL)
  {
    CkPrintf("Warning! RhoRSMap::Procnum had to assign maptable on %d!\n",CkMyPe() );
#ifdef USE_INT_MAP
    maptable= &RhoRSImaptable;
#else
    maptable= &RhoRSmaptable;
#endif

  }  
#ifdef USE_INT_MAP
  int retval=maptable->get(index[0], 0);
#else
  int retval=maptable->get(intdual(index[0], 0));
#endif
  return retval;

}


int RhoGSMap::procNum(int arrayHdl, const CkArrayIndex &iIndex)
{
  int *index=(int *) iIndex.data();
  if(maptable==NULL)
  {
    CkPrintf("Warning! RhoGSMap::Procnum had to assign maptable on %d!\n",CkMyPe() );
#ifdef USE_INT_MAP
    maptable= &RhoGSImaptable;
#else
    maptable= &RhoGSmaptable;
#endif
  }  
#ifdef USE_INT_MAP
  int retval=maptable->get(index[0], 0);
#else
  int retval=maptable->get(intdual(index[0], 0));
#endif
  CkAssert(retval>=0);
  CkAssert(retval<CkNumPes());
  return retval;

}


int RhoGHartMap::procNum(int arrayHdl, const CkArrayIndex &iIndex)
{
  int *index=(int *) iIndex.data();
  if(maptable==NULL)
  {
    CkPrintf("Warning! RhoGHartMap::Procnum had to assign maptable on %d!\n",CkMyPe() );
#ifdef USE_INT_MAP
    maptable= &RhoGHartImaptable;
#else
    maptable= &RhoGHartmaptable;
#endif

  }  
#ifdef USE_INT_MAP
  int retval=maptable->get(index[0], 0);
#else
  int retval=maptable->get(intdual(index[0], 0));
#endif
  CkAssert(retval>=0);
  CkAssert(retval<CkNumPes());
  return retval;

}    

int RhoRHartMap::procNum(int arrayHdl, const CkArrayIndex &idx)
{
  CkArrayIndex2D idx2d = *(CkArrayIndex2D *) &idx;
  if(maptable==NULL)
  {
    CkPrintf("Warning! RhoRHartMap::Procnum had to assign maptable on %d!\n",CkMyPe() );
#ifdef USE_INT_MAP
    maptable= &RhoRHartImaptable;
#else
    maptable= &RhoRHartmaptable;
#endif

  }  
#ifdef USE_INT_MAP
  int retval=maptable->get(idx2d.index[0], 0);
#else
  int retval=maptable->get(intdual(idx2d.index[0], 0));
#endif
  CkAssert(retval>=0);
  CkAssert(retval<CkNumPes());
  return retval;

}    
