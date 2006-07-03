/** \file TopoMap.C
 *  Author: Abhinav S Bhatele
 *  Date Created: January 04th, 2006 
 *  The functions in this file generate a topology sensitive
 *  mapping scheme for the GSpace, PairCalc and RealSpace objects.
 *  and also for the density objects - RhoR, RhoG, RhoGHart.
 *
 *  Heavily refactored by EJB 2006/5/31 to accelerate map creation
 */

#include "charm++.h"
#include "cpaimd.h"
#include "util.h"
#include "PeList.h"


#ifdef CMK_VERSION_BLUEGENE
#include "bgltorus.h"
#endif

/** 
 * Function for GSpace objects
 * Does a hybrid chunking of the GSpace grid. Depending on the 
 * no. of processors and the input parameter, Gstates_per_pe (l),
 * the other dimension 'm' is decided. And then this chunk of
 * size l by m is placed on a processor.   
 */

int GSMap::procNum(int handle, const CkArrayIndex &iIndex)
{
	int *index=(int *) iIndex.data();
	if(maptable==NULL)
	{
	  CkPrintf("Warning! GSMap::procnum had to assign maptable on pe %d!\n",CkMyPe());
	  maptable= &GSmaptable;
	}
	int retval=maptable->get(intdual(index[0],index[1]));
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
	int *index=(int *) iIndex.data();
	/* this works due to 4D being stored as 2 ints */

	if(maptable==NULL)
	{
	  CkPrintf("Warning! SCalc::Procnum had to assign maptable on %d !\n",CkMyPe());
	  if(symmetric)
	    maptable= &SymScalcmaptable;
	  else
	    maptable= &AsymScalcmaptable;
	}
	int retval=maptable->get(intdual(index[0], index[1]));
	return retval;
}

/** Function for RealSpace objects
 *
 */

int RSMap::procNum(int handle, const CkArrayIndex &iIndex)
{
	int *index=(int *) iIndex.data();
	if(maptable==NULL)
	{
	  CkPrintf("Warning! RSMap::Procnum had to assign maptable on %d!\n",CkMyPe());
	  maptable= &RSmaptable;
	}
	int retval=maptable->get(intdual(index[0], index[1]));
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
        maptable= &RhoRSmaptable;

      }  
	int retval=maptable->get(intdual(index[0], 0));
	CkAssert(retval>=0);
	CkAssert(retval<CkNumPes());
	return retval;

}

    
int RhoGSMap::procNum(int arrayHdl, const CkArrayIndex &iIndex)
{
  int *index=(int *) iIndex.data();
      if(maptable==NULL)
      {
	CkPrintf("Warning! RhoGSMap::Procnum had to assign maptable on %d!\n",CkMyPe() );
        maptable= &RhoGSmaptable;
      }  
      int retval=maptable->get(intdual(index[0], 0));
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
        maptable= &RhoGHartmaptable;

      }  
      int retval=maptable->get(intdual(index[0], 0));
      CkAssert(retval>=0);
      CkAssert(retval<CkNumPes());
      return retval;

}    

