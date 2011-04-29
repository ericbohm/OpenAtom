//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \file map.C
 *
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
 * will be a loop as each plane has a different load.
 */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================

#include <cmath>
#include "charm++.h"
#include "ckarray.h"
#include "util.h"
#include "main/cpaimd.h"
#include "main/groups.h"
#include "fft_slab_ctrl/fftCacheSlab.h"
#include "cp_state_plane/CP_State_Plane.h"

//============================================================================

extern int sizeX;
extern Config config;
extern CProxy_CPcharmParaInfoGrp scProxy;

//void hackGSpacePlaneLoad(CPcharmParaInfo *sim,int idx, double *line_load, 
//                         double *pt_load);
void GSpacePlaneLoad(int , double *, double *);

//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void GSMap::GSpacePlaneLoad(int idx, double *line_load, double *pt_load){
//============================================================================

  if( (idx < 0) || (idx >= nchareG) ){
   CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
   CkPrintf("GSpace plane index %d out of range: 0 < idx > %d\n",idx,nchareG);
   CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
   CkExit();
  }//endif

  line_load[0] = lines_per_chareG[idx];
  pt_load[0]   =   pts_per_chareG[idx];

//============================================================================
  }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void SCalcMap::GSpacePlaneLoad(int idx, double *line_load, double *pt_load){
//============================================================================


  if( (idx < 0) || (idx >= nchareG) ){
   CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
   CkPrintf("GSpace plane index %d out of range: 0 < idx > %d\n",idx,nchareG);
   CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
   CkExit();
  }//endif

  line_load[0] = lines_per_chareG[idx];
  pt_load[0]   =   pts_per_chareG[idx];

//============================================================================
  }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** 
 * Void allocate each state to more than one processor. Hybrid between 
 * state allocation and plane allocation
 */
//============================================================================
int basicMap(CkArrayIndex2D &idx2d, int numPlanes) {
//============================================================================
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
//============================================================================
   }// end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
 * The next two functions are needed to implement a mapping strategy different
 * from the default strategy.
 */
//============================================================================
void GSMap::makemap(){
//============================================================================
    int numChareG=0;
    if(config.doublePack) 
	numChareG = nchareG;
    else
	CkAbort("not doublepack broken!");

    for(int state = 0; state < config.nstates; state++){
      for (int plane = 0; plane < numChareG; plane++) {
	    CkArrayIndex2D idx2d(state,plane);
      }
    }
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//int GSMap::slowprocNum(int arrayHdl, const CkArrayIndex2D &idx2d)
//============================================================================
int GSMap::procNum(int arrayHdl, const CkArrayIndex &idx){
//============================================================================

  CkArrayIndex2D idx2d = *(CkArrayIndex2D *) &idx;

  int pe        = 0;
  int numChareG = 0;
  
  if(config.doublePack) {
     numChareG = nchareG;
  }else{
     CkAbort("not doublepack broken\n");
  }//endif
  
  if(state_load <= 0.0 ) {
    for(int x = 0; x  < numChareG; x++) {
      double curload = 0.0;
      double sload = 0.0;
      
      GSpacePlaneLoad(x, &curload, &sload);
      
      state_load += curload;
    }//endfor
  }//endif

  int pes_per_state = config.GpesPerState;
  int np = CkNumPes()/pes_per_state;
  
  if(np < 1){np = 1;}
  
  int partition_nstates = config.nstates / np;
  if(config.nstates % np != 0){partition_nstates ++;}
  
  int start_pe = (idx2d.index[0]/partition_nstates) * pes_per_state;
  int start_state = (idx2d.index[0]/partition_nstates) * partition_nstates;
  
  double cum_load = 0.0;
  double average_load = state_load * config.nstates / CkNumPes();

  for(int x = 0; x < numChareG; x++) {
    for(int s = 0; s < partition_nstates; s++) {
      double curload = 0.0;
      double sload = 0.0;
      
      GSpacePlaneLoad(x, &curload, &sload);

      cum_load += curload;

      if((idx2d.index[0] == s + start_state) && idx2d.index[1] == x) {
	double dpe = 0.0;
	dpe = cum_load / average_load;

	pe = (int)dpe;
	pe = pe % pes_per_state;
	pe += start_pe;

	return (pe % CkNumPes());
      }//endif
    }//endfor
  }//endfor
  return (idx2d.index[0]*1037+idx2d.index[1])%CkNumPes();

//============================================================================
   }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
int RSMap::procNum(int arrayHdl, const CkArrayIndex &idx){
    CkArrayIndex2D idx2d = *(CkArrayIndex2D *) &idx;

    int numPlanes = 0;
    numPlanes = sizeX;
    
    int pe  = 0;
    
    double total_planes = config.nstates * numPlanes;
    double planes_per_proc = total_planes / CkNumPes();
    
    pe = basicMap(idx2d, numPlanes);
    return pe;
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
inline double scalc_load(int x, int numplanes, double *load) {
  *load = numplanes * numplanes - x*x;
  return *load;
} 
//============================================================================


//============================================================================
/**
 * this one uses a lookup table built by calling the slow version
 */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
int SCalcMap::procNum(int hdl, const CkArrayIndex &idx){
    CkArrayIndex4D &idx4d = *(CkArrayIndex4D *)&idx;
    int intidx[2];
    CmiMemcpy(intidx,idx4d.index,2*sizeof(int));  // our 4 shorts are now 2 ints
    if(maptable==NULL){
	int calcchunks=max_states/gs;
	maptable= new CkHashtableT<intdual,int> (nchareG*calcchunks*calcchunks); 
                  // times blkSize, but its always 1
	makemap();
    }//endif
    return maptable->get(intdual(intidx[0], intidx[1]));
}
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void SCalcMap::makemap(){

  if(symmetric){

    for(int numX = 0; numX < nchareG; numX++){
      for (int s1 = 0; s1 < max_states; s1 += gs) {
	for (int s2 = s1; s2 < max_states; s2 += gs) {
	    CkArrayIndex4D idx4d(numX,s1,s2,0);
	    int intidx[2];
	    CmiMemcpy(intidx,idx4d.index,2*sizeof(int));  // our 4 shorts are now 2 ints
	    maptable->put(intdual(intidx[0],intidx[1]))=slowprocNum(0,idx4d);
	}//endfor
      }//endfor
    }//endfor

  }else{

      for(int numX = 0; numX < nchareG; numX++){
	  for (int s1 = 0; s1 < max_states; s1 += gs) {
	      for (int s2 = 0; s2 < max_states; s2 += gs) {
		  CkArrayIndex4D idx4d(numX,s1,s2,0);
		  int intidx[2];
		  CmiMemcpy(intidx,idx4d.index,2*sizeof(int));  // our 4 shorts now 2 ints
		  maptable->put(intdual(intidx[0],intidx[1]))=slowprocNum(0,idx4d);
	      }//endfor
	  }//endfor
      }//endfor

  }//endif

//============================================================================
   }//end routine
//============================================================================

int SCalcMap::slowprocNum2(int hdl, const CkArrayIndex4D &idx4d){
  short *idx = reinterpret_cast<short*> ( idx4d.data() );
  // Just use gspace as our guide for (w,x,y,z) use gsp(y+x/grainsize,w);
  return cheesyhackgsprocNum(scProxy.ckLocalBranch()->cpcharmParaInfo, idx[2] + idx[3]/gs, idx[0]);
}

/*
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
int SCalcMap::slowprocNum(int hdl, const CkArrayIndex4D &idx4d){
//============================================================================
//    CkArrayIndex4D &idx4d = *(CkArrayIndex4D *)&idx;
//    CkPrintf("scalc map call for [%d %d %d %d] on pe %d\n", 
//       idx4d.index[0],idx4d.index[1],idx4d.index[2],idx4d.index[3],CkMyPe());
//============================================================================
#ifndef CMK_OPTIMIZE
      double StartTime=CmiWallTimer();
#endif
    //Here maxY is the max number of planes;
    int planeid = idx4d.index[0];
    int numChareG = 0;

    if(config.doublePack){
        numChareG = config.nchareG;
    }else{
        numChareG = max_planes/2;
    }//endif

    //for asymetric
    double *load = new double[CkNumPes()];
    memset(load, 0, CkNumPes() * sizeof(double));

    int w=0, x=0, y = 0; 

    if(totalload <= 0.0) { 
      for(w = 0; w < numChareG; w ++) 
	for(x = 0; x < max_states; x += gs) {
	  if (symmetric){
	    y = x;
	  }else{
	    y = 0;
	  }//endif
	  
	  for(; y < max_states; y += gs) {
	    
	    double curload = 0.0;
	    double gload = 0.0;
	    
	    //scalc_load(w, numChareG, &curload);
	    GSpacePlaneLoad(w, &gload, &curload);
	    
	    totalload += curload;
	  }//endfor
	}//endfor
    }//endif
    
    int pe = 0;
    
    for(w = 0; w < numChareG; w ++) {
      for(x = 0; x < max_states; x += gs){
	if (symmetric){
	  y = x;
	}else{
	  y = 0;
	}//endif

	for(; y < max_states; y += gs) {
	  double curload = 0.0;
	  double gload = 0.0;
	  
	  //scalc_load(w, numChareG, &curload);
	  GSpacePlaneLoad(w, &gload, &curload);
	  
	  curload /=  totalload;
          
          if(load[pe] + curload > 0.3/CkNumPes()){pe ++;}
	  
	  if(pe >= CkNumPes()){pe = 0;}
          
	  load[pe] += curload;
	  
	  if((w == idx4d.index[0]) && (x == idx4d.index[1]) &&
	     (y == idx4d.index[2])) {
#ifndef CMK_OPTIMIZE
             traceUserBracketEvent(Scalcmap_, StartTime, CmiWallTimer());    
#endif
             delete [] load;
             return pe;
	  }//endif
        }//endfor
      }//endfor
    }//endfor : w

    delete [] load;
#ifndef CMK_OPTIMIZE
    traceUserBracketEvent(Scalcmap_, StartTime, CmiWallTimer());    
#endif

//============================================================================
    return (idx4d.index[0]*197+idx4d.index[1]*23+idx4d.index[2]*7
            +idx4d.index[3])%CkNumPes();    
//============================================================================
   }//end routine
//============================================================================
*/


int SCalcMap::slowprocNum(int hdl, const CkArrayIndex4D &idx4d)
{
    short *idx = reinterpret_cast<short*> ( idx4d.data() );

    //Here maxY is the max number of planes;
    int planeid = idx[0];
    int numChareG = 0;


    if(config.doublePack){
        numChareG = config.nchareG;
    }else{
        numChareG = max_planes/2;
    }//endif

    //for asymetric
    double *load = new double[CkNumPes()];
    memset(load, 0, CkNumPes() * sizeof(double));

    int w=0, x=0, y = 0; 

    if(totalload <= 0.0) { 
      for(w = 0; w < numChareG; w ++) 
	for(x = 0; x < max_states; x += gs) {
	  if (symmetric)
	    y = x;
	  else
	    y = 0;
	  
	  for(; y < max_states; y += gs) {
	    
	    double curload = 0.0;
	    double gload = 0.0;
	    
	    //scalc_load(w, numChareG, &curload);
	    GSpacePlaneLoad(w, &gload, &curload);
	    
	    totalload += curload;
	  }
	}
    }
    
    int pe = 0;
    
    for(w = 0; w < numChareG; w ++)  
      for(x = 0; x < max_states; x += gs) {
	if (symmetric)
	  y = x;
	else
	  y = 0;
	
	for(; y < max_states; y += gs) {
	  double curload = 0.0;
	  double gload = 0.0;
	  
	  //scalc_load(w, numChareG, &curload);
	  GSpacePlaneLoad(w, &gload, &curload);
	  
	  curload /=  totalload;
          
          if(load[pe] + curload > 0.3/CkNumPes()) 
              pe ++;
	  
	  if(pe >= CkNumPes())
              pe = 0;
          
	  load[pe] += curload;
	  
	  if((w == idx[0]) && (x == idx[1]) && (y == idx[2])) {

              //if(CkMyPe() == 0)
              //  CkPrintf ("scalc %d %d %d %d assigned to pe %d and curload = %f, load = %f\n", w, x ,y, symmetric, pe, curload, load[pe]);
              
              delete [] load;
	    return pe;
	  }
        }
      }
    
    delete [] load;
    return (idx[0]*197+idx[1]*23+idx[2]*7+idx[3])%CkNumPes();    

}


