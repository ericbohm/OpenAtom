//**************************************************************************
 /** \file pairCalculator.C
 * This is a matrix multiply library with extra frills to communicate the
 * results back to gspace or the calling ortho char as directed by the
 * callback.
 *
 * The pairCalculator handles initialization and creation of the
 * ckPairCalculator chare arrays, their reduction group, their multicast
 * manager, and their section proxies
 *
 * Expected usage begins with createPairCalculator(*,PairCalcID *,*)
 * the PairCalcID contains meta information about the calculator.
 * In particular, the various section proxies and array ids necessary to
 * handle the expected communication modalities between a parent array and
 * the ckPairCalculator array.
 *
 * After the main::main proc 0 phase the array sections used to populate
 * and return data from the paircalculator are called.
 * The forward path section reduction (to and from ortho) is initialized
 * via initOneRedSect() the backward path is initialized via the
 * appropriate makeOneResultSection_X() call.  In each case the call
 * should be made by each GSP or Ortho object.  That way each one has its
 * own proxy and the section tree will only include relevant processors.
 *
 * Data from GSP travels to the PairCalculators via an InputDataHandler chare 
 * array of the same dimensions as, and bound to, the PairCalculator array. 
 * Appropriate sections/lists of this input handler array for multicasting 
 * the data from GSP to are built in makeLeftTree() and makeRightTree().
 *
 * Each iteration of the GSP-PC-Ortho-PC-GSP loop is started by GSP calling
 * startPairCalcLeft() and startPairCalcRight(). These are simply #defines
 * that turn into the appropriate function: sendLeftData() and sendRightData()
 * or their RDMA equivalents. The input handler chares then collate all the
 * incoming data and then wake their corresponding PC chares once all the data
 * is available. The PCs then do their thing and the result is returned via 
 * the callback set in the create routine (which happens to be Ortho).
 *
 * The backward path is triggered by:
 * finishPairCalcSection(PairCalcID, datasize, data *)
 * Its result is returned via the callback entry point passed in
 * during creation
 *
 * The results of the backward path are returned in a set of section
 * reductions.  The reduction contributes its slice of its matrix of
 * doubles with the offset=thisIndex.z.  The client then returns the
 * sum of each slice to the GSP element that created the section with
 * the offset so it can be copied into the correct place in the points
 * array.
 *
 * The chunk decomposition changes the setup substantially.  In chunk
 * decomposition we send a piece of the nonzero points of gspace to
 * the paircalculators.  It is not a multicast.  Each GSP[P,S] will
 * send its ith chunk to PC[P,0,0,i].  Once nstate chunks arrive at a
 * PC the multiply can proceed.
 *
 * In the hybrid case this becomes a multicast of chunks to the
 * appropriate state decomposition destination as before.
 *
 */
//*************************************************************************
#include "pairCalculator.h"
#include "InputDataHandler.h"
#include "cp_state_ctrl/pcCommManager.h" 
#include <algorithm>

#ifdef USE_COMLIB
extern ComlibInstanceHandle mcastInstanceCP;
extern ComlibInstanceHandle mcastInstanceACP;
extern ComlibInstanceHandle gAsymInstance;
extern ComlibInstanceHandle gSymInstance;
#endif


void createPairCalculator(const pc::pcConfig pcCfg, PairCalcID* pcid, int comlib_flag, CkGroupID *mapid, int priority, CkGroupID mCastGrpId)
{

  traceRegisterUserEvent("calcpairDGEMM", 210);
  traceRegisterUserEvent("calcpairContrib", 220);
  traceRegisterUserEvent("multiplyResultDGEMM1", 230);
  traceRegisterUserEvent("multiplyResultDGEMM2", 240);
  traceRegisterUserEvent("multiplyResultDGEMM1R", 250);

  CkArrayOptions paircalcOpts,handlerOpts;
  CProxy_PairCalculator pairCalculatorProxy;
  CProxy_InputDataHandler<CollatorType,CollatorType> inputHandlerProxy;


  // If a chare mapping is not available, create an empty array
  if(!mapid) 
  {
    pairCalculatorProxy = CProxy_PairCalculator::ckNew();
  }
  // else, create the array with element locations as specified by the map 
  else 
  {
    paircalcOpts.setMap(*mapid);
    pairCalculatorProxy = CProxy_PairCalculator::ckNew(inputHandlerProxy, pcCfg, paircalcOpts);
  }

#ifdef DEBUG_CP_PAIRCALC_CREATION
	CkPrintf("createPairCalculator: Creating empty PairCalculator and InputDataHandler chare arrays for %d loop: asymm(0)/symm(1)\n",pcCfg.isSymmetric);
#endif 
  /// Create an empty input handler chare array that will accept all incoming messages from GSpace
  handlerOpts.bindTo(pairCalculatorProxy);
  inputHandlerProxy = CProxy_InputDataHandler<CollatorType,CollatorType> ::ckNew(pairCalculatorProxy,handlerOpts);

  int proc = 0;
  // Initialize the PairCalcID instance
  pcid->Init(pairCalculatorProxy.ckGetArrayID(), inputHandlerProxy.ckGetArrayID(), pcCfg.grainSize, pcCfg.numChunks, pcCfg.numStates, pcCfg.isSymmetric, comlib_flag, pcCfg.isDoublePackOn, pcCfg.conserveMemory, pcCfg.isLBon, priority, !pcCfg.isInputMulticast);
  pcid->mCastGrpId=mCastGrpId;

#ifdef USE_COMLIB

  // Setup the appropriate multicast strategy
#ifdef CMK_BLUEGENEL
  //  CharmStrategy *multistrat = new RectMulticastStrategy(pairCalculatorProxy.ckGetArrayID());
  Strategy *multistrat = new DirectMulticastStrategy();
#else
  Strategy *multistrat = new DirectMulticastStrategy();
#endif
  
#endif

  // 
  int maxpcstateindex=(pcid->nstates/pcid->GrainSize-1)*pcid->GrainSize;
  // 

#ifdef USE_COMLIB
  if(pcCfg.isSymmetric)
    mcastInstanceCP=ComlibRegister(multistrat);
  else
    mcastInstanceACP=ComlibRegister(multistrat);
#endif

  CkAssert(mapid);
  // If the symmetric loop PC instances are being created
  if(pcCfg.isSymmetric)
	for(int numX = 0; numX < pcCfg.numPlanes; numX ++)
	{
	  for (int s1 = 0; s1 <= maxpcstateindex; s1 += pcCfg.grainSize) 
      {
      	// If phantomSym is turned on
		int s2start=(pcCfg.arePhantomsOn) ? 0 : s1;
		for (int s2 = s2start; s2 <= maxpcstateindex; s2 += pcCfg.grainSize) 
		{
			for (int c = 0; c < pcCfg.numChunks; c++) 
			{
				if(mapid) 
				{
					#ifdef DEBUG_CP_PAIRCALC_CREATION
						CkPrintf("Inserting PC element [%d %d %d %d %d]\n",numX,s1,s2,c,pcCfg.isSymmetric);
					#endif
					pairCalculatorProxy(numX,s1,s2,c).insert(inputHandlerProxy, pcCfg);
				}
				else
				{
					#ifdef DEBUG_CP_PAIRCALC_CREATION
						CkPrintf("Inserting PC element [%d %d %d %d %d] at PE %d\n",numX,s1,s2,c,pcCfg.isSymmetric,proc);
					#endif
					pairCalculatorProxy(numX,s1,s2,c).insert(inputHandlerProxy, pcCfg, proc);
					proc++;
					if (proc >= CkNumPes()) proc = 0;
				}
			}
		}
	  }
	}
  // else, if the asymmetric loop PC instances are being created
  else
  {
	for(int numX = 0; numX < pcCfg.numPlanes; numX ++)
	{
		for (int s1 = 0; s1 <= maxpcstateindex; s1 += pcCfg.grainSize)
		{
			for (int s2 = 0; s2 <= maxpcstateindex; s2 += pcCfg.grainSize)
			{
				for (int c = 0; c < pcCfg.numChunks; c++)
				{
					if(mapid)
					{
						#ifdef DEBUG_CP_PAIRCALC_CREATION
							CkPrintf("Inserting PC element [%d %d %d %d %d]\n",numX,s1,s2,c,pcCfg.isSymmetric);
						#endif
						pairCalculatorProxy(numX,s1,s2,c).insert(inputHandlerProxy, pcCfg);
					}
					else
					{
						#ifdef DEBUG_CP_PAIRCALC_CREATION
							CkPrintf("Inserting PC element [%d %d %d %d %d] on PE %d\n",numX,s1,s2,c,pcCfg.isSymmetric,proc);
						#endif
						pairCalculatorProxy(numX,s1,s2,c).insert(inputHandlerProxy, pcCfg, proc);
						proc++;
						if (proc >= CkNumPes()) proc = 0;
					}
				}
			}
		}
	}
  }
  /// Notify the runtime that we're done inserting all the PC elements
  pairCalculatorProxy.doneInserting();
#ifdef _PAIRCALC_DEBUG_
  CkPrintf("    Finished init {grain=%d, sym=%d, blk=%d, Z=%d, S=%d}\n", pcCfg.grainSize, pcCfg.isSymmetric, pcCfg.numChunks, pcCfg.numPlanes, pcCfg.numStates);
#endif
}




void dumpMatrixDouble(const char *infilename, double *matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric)
{
  char fmt[1000];
  char filename[1000];
  strncpy(fmt,infilename,999);
  strncat(fmt,"_%d_%d_%d_%d_%d.out",999);
  sprintf(filename,fmt, w, x, y, z, symmetric);
  FILE *loutfile = fopen(filename, "w");
  for(int i=0;i<xdim;i++)
    for(int j=0;j<ydim;j++)
      fprintf(loutfile,"%d %d %.12g\n",x+i,y+j,matrix[i*ydim+j]);
  fclose(loutfile);
}



void loadMatrixDouble(const char *infilename, double *matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric)
{
  char fmt[1000];
  char filename[1000];
  strncpy(fmt,infilename,999);
  strncat(fmt,"_%d_%d_%d_%d_%d.out",999);
  sprintf(filename,fmt, w, x, y, z, symmetric);
  FILE *loutfile = fopen(filename, "r");
  if(loutfile!=NULL)
    {
      int junk1,junk2;
      for(int i=0;i<xdim;i++)
	for(int j=0;j<ydim;j++)
	  fscanf(loutfile,"%d %d %lf\n",&junk1,&junk2,&(matrix[i*ydim+j]));
      fclose(loutfile);
    }
  else
    {
      CkAbort(filename);
    }
}

//! NOTE: this uses the evil piny convention
void dumpMatrix2DDouble(const char *infilename, double **matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric)
{
  char fmt[1000];
  char filename[1000];
  strncpy(fmt,infilename,999);
  strncat(fmt,"_%d_%d_%d_%d_%d.out",999);
  sprintf(filename,fmt, w, x, y, z, symmetric);
  FILE *loutfile = fopen(filename, "w");
  for(int i=0;i<xdim;i++)
    for(int j=1;j<=ydim;j++)
      fprintf(loutfile,"%d %d %.12g\n",i,j,matrix[i][j]);
  fclose(loutfile);
}
void loadMatrix2DDouble(const char *infilename, double **matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric)
{
  char fmt[1000];
  char filename[1000];
  strncpy(fmt,infilename,999);
  strncat(fmt,"_%d_%d_%d_%d_%d.out",999);
  sprintf(filename,fmt, w, x, y, z, symmetric);
  FILE *loutfile = fopen(filename, "r");
  if(loutfile!=NULL)
    {
      int junk1,junk2;
      for(int i=0;i<xdim;i++)
	for(int j=1;j<=ydim;j++)
	  fscanf(loutfile,"%d %d %lf\n",&junk1,&junk2,&(matrix[i][j]));
      fclose(loutfile);
    }
  else
    {
      CkAbort(filename);
    }
}


//! NOTE: this uses the evil piny convention
void dumpMatrix2DInt(const char *infilename, int **matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric)
{
  char fmt[1000];
  char filename[1000];
  strncpy(fmt,infilename,999);
  strncat(fmt,"_%d_%d_%d_%d_%d.out",999);
  sprintf(filename,fmt, w, x, y, z, symmetric);
  FILE *loutfile = fopen(filename, "w");
  for(int i=0;i<xdim;i++)
    for(int j=1;j<=ydim;j++)
      fprintf(loutfile,"%d %d %d\n",i,j,matrix[i][j]);
  fclose(loutfile);
}
void loadMatrix2DInt(const char *infilename, int **matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric)
{
  char fmt[1000];
  char filename[1000];
  strncpy(fmt,infilename,999);
  strncat(fmt,"_%d_%d_%d_%d_%d.out",999);
  sprintf(filename,fmt, w, x, y, z, symmetric);
  FILE *loutfile = fopen(filename, "r");
  if(loutfile!=NULL)
    {
      int junk1,junk2;
      for(int i=0;i<xdim;i++)
	for(int j=1;j<=ydim;j++)
	  fscanf(loutfile,"%d %d %d\n",&junk1,&junk2,&(matrix[i][j]));
      fclose(loutfile);
    }
  else
    {
      CkAbort(filename);
    }
}


/**
 * synchronize for migration
 */
void isAtSyncPairCalc(PairCalcID* pcid){
#ifdef _PAIRCALC_DEBUG_
  CkPrintf("     lbsync symm=%d\n", pcid->Symmetric);
#endif
  CkArrayID pairCalculatorID = (CkArrayID)pcid->Aid;
  CProxy_PairCalculator pairCalculatorProxy(pairCalculatorID);
  pairCalculatorProxy.lbsync();
}




#ifdef ROTATE_LIST
bool reorder_elem_list(CkArrayIndexMax *elems, int numelems, int newstart)
{
  if(newstart>numelems)
    return(false);
  CkArrayIndexMax  swap;
  int where=newstart;
  for(int i=0;i<newstart;i++,where++)
    {
      if(where>=numelems)
	where=0;
      swap=elems[i];
      elems[i]=elems[where];
      elems[where]=swap;
    }
  return(true);
}

#else

bool reorder_elem_list(CkArrayIndexMax *elems, int numelems, int newstart)
{
  //  CkPrintf("reordering list of %d elems\n", numelems);
  std::random_shuffle(elems,elems+numelems);
  return(true);
}

bool reorder_elem_list_4D(CkArrayIndex4D *elems, int numelems, int newstart)
{
  //  CkPrintf("reordering list of %d elems\n", numelems);
  std::random_shuffle(elems,elems+numelems);
  return(true);
}

bool reorder_elem_list_max(CkArrayIndexMax *elems, int numelems, int newstart)
{
  //  CkPrintf("reordering list of %d elems\n", numelems);
  std::random_shuffle(elems,elems+numelems);
  return(true);
}
#endif

#include "RDMAMessages.def.h"

