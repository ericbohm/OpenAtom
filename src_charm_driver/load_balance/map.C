
#include <math.h>
#include "charm++.h"
#include "ckarray.h"
#include "util.h"
#include "groups.h"
#include "cpaimd.h"
#include "sim_subroutines.h"
#include "CP_State_Plane.h"

extern int sizeX;
extern Config config;
extern CProxy_CPcharmParaInfoGrp scProxy;
void GSpacePlaneLoad(int , double *, double *);
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/*
 * The maps are intended to spread load around evenly and keep plane
 * and state objects which communicate across phases close together.

 * It basically uses this key idea pe = cumload /  average_load so
 * cumload depends on the object list traversal.  For  plane wise we
 * traverse it plane wise and for state wise do it states in the
 * first loop and planes in the second loop.  For hybrid map, it
 * partitions states and processors in to K partitions  
 * k = numpes/pesperstate so apply the cumulative load logic in each 
 * partition, and add the start pe of the partition to it

 * cumload is the cumulative load till that object average load
 * is the total_load/numpes

 * Assuming object s1,p1 is being inserted s1 is a state and p1 is a
 * plane the map will use a certain  traversal through object space
 * for state map it will traverse all  objects state by state hence
 * cumload = load_sum { s1 * numplanesperstate + p1} for gspace it
 * will be a loop as each plane * has a different load
*/

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================

void GSMap::GSpacePlaneLoad(int idx, double *line_load, double *pt_load){

//============================================================================


  if( (idx < 0) || (idx >= nplane_x) ){
   CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
   CkPrintf("GSpace plane index %d out of range: 0 < idx > %d\n",idx,nplane_x);
   CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
   CkExit();
  }//endif

  line_load[0] = lines_per_plane[idx];
  pt_load[0]   =   pts_per_plane[idx];

//============================================================================
  }//end routine


//============================================================================
void SCalcMap::GSpacePlaneLoad(int idx, double *line_load, double *pt_load){

//============================================================================


  if( (idx < 0) || (idx >= nplane_x) ){
   CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
   CkPrintf("GSpace plane index %d out of range: 0 < idx > %d\n",idx,nplane_x);
   CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
   CkExit();
  }//endif

  line_load[0] = lines_per_plane[idx];
  pt_load[0]   =   pts_per_plane[idx];

//============================================================================
  }//end routine
//============================================================================


//Void allocate each state to more than one processor. Hybrid between 
//state allocation and plane allocation
int basicMap(CkArrayIndex2D &idx2d, int numPlanes) {
  int pes_per_state = config.RpesPerState;
  int np = CkNumPes()/pes_per_state;
 
  if(np < 1)
    np = 1;

  int partition_nstates = config.nstates / np;
  if(config.nstates % np != 0)
    partition_nstates ++;
  
  double total_planes = config.nstates * numPlanes;
  double numPerProc = total_planes / CkNumPes();
    
  //My state in the current partition
  int s = idx2d.index[0]%partition_nstates; 

  int pe = (int) floor((idx2d.index[1]*partition_nstates + s)/numPerProc); 
    
  pe = pe % pes_per_state;
  pe +=  pes_per_state *  (idx2d.index[0]/partition_nstates);

  if(pe >= CkNumPes())
    pe = pe % CkNumPes();

  return pe; 
}

/*
 * map definitions
 * The next two functions are needed to implement a mapping strategy different
 * from the default strategy.
 */
void GSMap::makemap()
{
    int numPlanes=0;
    if(config.doublePack) 
	numPlanes = nplane_x;
    else
	CkAbort("not doublepack broken!");

    for(int state = 0; state < config.nstates; state++){
      for (int plane = 0; plane < numPlanes; plane++) {
	    CkArrayIndex2D idx2d(state,plane);
//	    maptable->put(intdual(state,plane))=slowprocNum(0,idx2d);
//	    CkPrintf("mapping [%d %d] as to pe %d based on slowprocnums %d\n",state,plane,maptable->get(intdual(state,plane)),slowprocNum(0,idx2d));
      }
    }
}

// this one uses a lookup table built by calling the slow version
/*
int GSMap::procNum(int hdl, const CkArrayIndex &idx)
{
    CkArrayIndex2D &idx2d = *(CkArrayIndex2D *)&idx;
    if(maptable==NULL)
    {
    int numPlanes=0;
	if(config.doublePack) 
	numPlanes = nplane_x;
	else
	CkAbort("not doublepack broken!");
	maptable= new CkHashtableT<intdual,int> (numPlanes*config.nstates); 
	makemap();
    }
    return maptable->get(intdual(idx2d.index[0], idx2d.index[1]));
}
*/

//FOOBAR : Breaking planes per chare. They never worked anyway and a
//cleanup was really necessary (-Sameer 04/05)

//int GSMap::slowprocNum(int arrayHdl, const CkArrayIndex2D &idx2d)
int GSMap::procNum(int arrayHdl, const CkArrayIndex &idx)
{
  CkArrayIndex2D idx2d = *(CkArrayIndex2D *) &idx;

  int pe = 0;
  int numPlanes = 0;
  
  if(config.doublePack) 
      numPlanes = nplane_x;
  else
      CkAbort("not doublepack broken\n");
  
  //pe = basicMap(idx2d, numPlanes);
  //return pe;

  if(state_load <= 0.0 ) {
    for(int x = 0; x  < numPlanes; x++) {
      double curload = 0.0;
      double sload = 0.0;
      
      GSpacePlaneLoad(x, &curload, &sload);
      
      state_load += curload;
    }
  }


  int pes_per_state = config.GpesPerState;
  int np = CkNumPes()/pes_per_state;
  
  if(np < 1)
    np = 1;
  
  int partition_nstates = config.nstates / np;
  if(config.nstates % np != 0)
    partition_nstates ++;
  
  int start_pe = (idx2d.index[0]/partition_nstates) * pes_per_state;
  int start_state = (idx2d.index[0]/partition_nstates) * partition_nstates;
  
  double cum_load = 0.0;
  double average_load = state_load * config.nstates / CkNumPes();

  for(int x = 0; x < numPlanes; x++) 
    for(int s = 0; s < partition_nstates; s++) {
      double curload = 0.0;
      double sload = 0.0;
      
      GSpacePlaneLoad(x, &curload, &sload);

      cum_load += curload;

      if((idx2d.index[0] == s + start_state) && idx2d.index[1] == x) {
      
	//if(CkMyPe() == 0)
	//  CkPrintf("Load[%d] = %g\n", pe, load[pe]);
	
	double dpe = 0.0;
	dpe = cum_load / average_load;

	pe = (int)dpe;
	pe = pe % pes_per_state;
	pe += start_pe;

	return pe % CkNumPes();
      }
    }
  
  //  CkPrintf("Warning pe not found for index [%d, %d]\n",idx2d.index[0],idx2d.index[1]);
  return (idx2d.index[0]*1037+idx2d.index[1])%CkNumPes();
}

int RSMap::procNum(int arrayHdl, const CkArrayIndex &idx)
{
    CkArrayIndex2D idx2d = *(CkArrayIndex2D *) &idx;

    int numPlanes = 0;
    numPlanes = sizeX;
    
    int pe  = 0;
    
    double total_planes = config.nstates * numPlanes;
    double planes_per_proc = total_planes / CkNumPes();
    
    pe = basicMap(idx2d, numPlanes);
    return pe;

    //pe = (int)((idx2d.index[1] * config.nstates +
    //idx2d.index[0])/planes_per_proc); 
    //return pe%CkNumPes();
}

inline double scalc_load(int x, int numplanes, double *load) {
  *load = numplanes * numplanes - x*x;
  return *load;
} 


// this one uses a lookup table built by calling the slow version
int SCalcMap::procNum(int hdl, const CkArrayIndex &idx)
{
    CkArrayIndex4D &idx4d = *(CkArrayIndex4D *)&idx;
    int intidx[2];
    CmiMemcpy(intidx,idx4d.index,2*sizeof(int));  // our 4 shorts are now 2 ints
    if(maptable==NULL)
    {
	int calcchunks=max_states/gs;
	maptable= new CkHashtableT<intdual,int> (nplane_x*calcchunks*calcchunks); // times blkSize, but its always 1
	/* that blkSize comment is there just in case someone does something
	 * with blksize ever they can find this reference in a search */
	makemap();
    }
    return maptable->get(intdual(intidx[0], intidx[1]));
}

void SCalcMap::makemap()
{

  if(symmetric)
    for(int numX = 0; numX < nplane_x; numX++){
      for (int s1 = 0; s1 < max_states; s1 += gs) {
	for (int s2 = s1; s2 < max_states; s2 += gs) {
	    CkArrayIndex4D idx4d(numX,s1,s2,0);
	    int intidx[2];
	    CmiMemcpy(intidx,idx4d.index,2*sizeof(int));  // our 4 shorts are now 2 ints
	    maptable->put(intdual(intidx[0],intidx[1]))=slowprocNum(0,idx4d);
//	    CkPrintf("mapping [%d %d %d %d %d] as [%d %d] to pe %d based on slowprocnums %d\n",numX,s1,s2,0,symmetric,intidx[0],intidx[1],maptable->get(intdual(intidx[0],intidx[1])),slowprocNum(0,idx4d));
	}
      }
    }
  else
      for(int numX = 0; numX < nplane_x; numX++){
	  for (int s1 = 0; s1 < max_states; s1 += gs) {
	      for (int s2 = 0; s2 < max_states; s2 += gs) {
		  CkArrayIndex4D idx4d(numX,s1,s2,0);
		  int intidx[2];
		  CmiMemcpy(intidx,idx4d.index,2*sizeof(int));  // our 4 shorts are now 2 ints
		  maptable->put(intdual(intidx[0],intidx[1]))=slowprocNum(0,idx4d);
//		  CkPrintf("mapping [%d %d %d %d %d] as [%d %d] to pe %d based on slowprocnums %d\n",numX,s1,s2,0,symmetric,intidx[0],intidx[1],maptable->get(intdual(intidx[0],intidx[1])),slowprocNum(0,idx4d));
	      }
	  }
      }
}


int SCalcMap::slowprocNum(int hdl, const CkArrayIndex4D &idx4d)
{
//    CkArrayIndex4D &idx4d = *(CkArrayIndex4D *)&idx;
//    CkPrintf("scalc map call for [%d %d %d %d] on pe %d\n", idx4d.index[0],idx4d.index[1],idx4d.index[2],idx4d.index[3],CkMyPe());
#ifndef CMK_OPTIMIZE
      double StartTime=CmiWallTimer();
#endif
  
    //Here maxY is the max number of planes;
    int planeid = idx4d.index[0];
    int numPlanes = 0;

    if(config.doublePack) 
        numPlanes = max_planes/4;
    else
        numPlanes = max_planes/2;

    //for asymetric
    double *load = new double[CkNumPes()];
    memset(load, 0, CkNumPes() * sizeof(double));

    int w=0, x=0, y = 0; 

    if(totalload <= 0.0) { 
      for(w = 0; w < numPlanes; w ++) 
	for(x = 0; x < max_states; x += gs) {
	  if (symmetric)
	    y = x;
	  else
	    y = 0;
	  
	  for(; y < max_states; y += gs) {
	    
	    double curload = 0.0;
	    double gload = 0.0;
	    
	    //scalc_load(w, numPlanes, &curload);
	    GSpacePlaneLoad(w, &gload, &curload);
	    
	    totalload += curload;
	  }
	}
    }
    
    int pe = 0;
    
    for(w = 0; w < numPlanes; w ++)  
      for(x = 0; x < max_states; x += gs) {
	if (symmetric)
	  y = x;
	else
	  y = 0;
	
	for(; y < max_states; y += gs) {
	  double curload = 0.0;
	  double gload = 0.0;
	  
	  //scalc_load(w, numPlanes, &curload);
	  GSpacePlaneLoad(w, &gload, &curload);
	  
	  curload /=  totalload;
          
          if(load[pe] + curload > 0.3/CkNumPes()) 
              pe ++;
	  
	  if(pe >= CkNumPes())
              pe = 0;
          
	  load[pe] += curload;
	  
	  if((w == idx4d.index[0]) && (x == idx4d.index[1]) &&
	     (y == idx4d.index[2])) {

              //if(CkMyPe() == 0)
              //  CkPrintf ("scalc %d %d %d %d assigned to pe %d and curload = %f, load = %f\n", w, x ,y, symmetric, pe, curload, load[pe]);
#ifndef CMK_OPTIMIZE
      traceUserBracketEvent(Scalcmap_, StartTime, CmiWallTimer());    
#endif
              delete [] load;
	    return pe;

	  }
        }
      }
    
    delete [] load;
#ifndef CMK_OPTIMIZE
      traceUserBracketEvent(Scalcmap_, StartTime, CmiWallTimer());    
#endif
    return (idx4d.index[0]*197+idx4d.index[1]*23+idx4d.index[2]*7+idx4d.index[3])%CkNumPes();    

}

/*
int SCalcMap::procNum(int hdl, const CkArrayIndex &idx)
{
    CkArrayIndex4D &idx4d = *(CkArrayIndex4D *)&idx;
  
    //Here maxY is the max number of planes;
    int planeid = idx4d.index[0];
    int numPlanes = 0;

    if(config.doublePack) 
        numPlanes = max_planes/4;
    else
        numPlanes = max_planes/2;

    int w=0, x=0, y = 0; 

    if(totalload <= 0.0) { 

      for(w = 0; w < numPlanes; w ++) 
	for(x = 0; x < max_states; x += gs) {
	  if (symmetric)
	    y = x;
	  else
	    y = 0;
	  
	  for(; y < max_states; y += gs) {
	    
	    double curload = 0.0;
	    double gload = 0.0;
	    
	    //scalc_load(w, numPlanes, &curload);
	    GSpacePlaneLoad(w, &gload, &curload);
	    
	    totalload += curload;
	  }
	}
    }
    
    int pe = 0;
    double cum_load = 0.0;
    double average_load = totalload / CkNumPes();

    for(w = 0; w < numPlanes; w ++)  
      for(x = 0; x < max_states; x += gs) {
	if (symmetric)
	  y = x;
	else
	  y = 0;
	
	for(; y < max_states; y += gs) {
	  double curload = 0.0;
	  double gload = 0.0;
	  
	  //scalc_load(w, numPlanes, &curload);
	  GSpacePlaneLoad(w, &gload, &curload);

	  cum_load += curload;
	  
	  if((w == idx4d.index[0]) && (x == idx4d.index[1]) &&
	     (y == idx4d.index[2])) {
	    double dest_pe = cum_load/average_load;
	    
	    pe = (int) dest_pe;
	    return pe % CkNumPes();
	  }
        }
      }
    
    return (idx4d.index[0]*197+idx4d.index[1]*23+idx4d.index[2]*7+idx4d.index[3])%CkNumPes();    
}


*/
