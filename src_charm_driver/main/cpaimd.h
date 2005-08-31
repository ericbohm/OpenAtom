//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/*
 * Some basic data structures and the array map classes are defined
 * here. These maps are used to map the array elements to the correct
 * processors.  Please read ../doc/documentation.ps for a detailed
 * documentation on the structure of the program.
 */
//============================================================================

#ifndef _CPAIMD_H
#define _CPAIMD_H

#include "EachToManyMulticastStrategy.h"
#include "RingMulticastStrategy.h"
#include "StreamingStrategy.h"
#include "pairCalculator.h"
#include "fftlib.h"
#include "ckhashtable.h"

#define LOAD_BALANCE_STEP 10

#define PRE_BALANCE_STEP 2

#define FIRST_BALANCE_STEP 10

#define TRACE_ON_STEP 4
#define TRACE_OFF_STEP 7

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class main : public Chare {
 public:
    main(CkMigrateMessage *m) {}
    main(CkArgMsg *);

};
//============================================================================


//helper class for map hashtables copied from femrefine.C,  thanks Sayantan
class intdual{
 private:
    int x,y;
 public:
    intdual(int _x,int _y){
	if(_x <= _y){
	    x = _x; y=_y;
	}else{
	    x = _y; y= _x;
	}
    }
    inline int getx(){return x;};
    inline int gety(){return y;};
    inline CkHashCode hash() const {
	return (CkHashCode)(x+y);
    }
    static CkHashCode staticHash(const void *k,size_t){
	return ((intdual *)k)->hash();
    }
    inline int compare(intdual &t) const{
	return (t.getx() == x && t.gety() == y);
    }
    static int staticCompare(const void *a,const void *b,size_t){
	return ((intdual *)a)->compare((*(intdual *)b));
    }
};



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class GSMap: public CkArrayMap {
  

  int nplane_x;
  double *lines_per_plane;
  double *pts_per_plane;
  double state_load;
  CkHashtableT<intdual, int> *maptable;
 public:
  GSMap() { state_load = 0.0; }
  GSMap(int _nplane_x, double *_lines_per_plane, double *_pts_per_plane): 
        nplane_x(_nplane_x)   
      { 
	  state_load = 0.0; 
	  lines_per_plane= new double[nplane_x];
	  pts_per_plane= new double[nplane_x];
	  memcpy(lines_per_plane,_lines_per_plane,nplane_x*sizeof(double));
	  memcpy(pts_per_plane,_pts_per_plane,nplane_x*sizeof(double));
	  maptable=NULL;
      }
  int procNum(int, const CkArrayIndex &);
//  int slowprocNum(int, const CkArrayIndex2D &);
  void makemap();
  void GSpacePlaneLoad(int idx, double *line_load, double *pt_load);
    void pup(PUP::er &p)
	{
	    CkArrayMap::pup(p);
	    p|nplane_x;
	    p|state_load;
	    maptable=NULL;
	    if (p.isUnpacking()) {
		lines_per_plane= new double[nplane_x];
		pts_per_plane= new double[nplane_x];
	    }	    
	    p(lines_per_plane,nplane_x);
	    p(pts_per_plane,nplane_x);
	}

  ~GSMap(){
    delete [] lines_per_plane;
    delete [] pts_per_plane;
  }
};
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================

class RSMap: public CkArrayMap {

 public:

  RSMap() {}
    int procNum(int, const CkArrayIndex &);
};
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class SCalcMap : public CkArrayMap {
  int nplane_x;
  double *lines_per_plane;
  double *pts_per_plane;
  CkHashtableT<intdual, int> *maptable;
  int max_states, max_planes, gs;
  CmiBool symmetric;
  double totalload;
    
 public:
    
    SCalcMap(int nstates, int nplanes,  int gs, CmiBool flag, int _nplane_x, 
             double *_lines_per_plane, double *_pts_per_plane) { 
        this->gs = gs;
        max_states = nstates;
        max_planes = nplanes;
        symmetric = flag;
	totalload = 0.0;
	nplane_x=_nplane_x;
	lines_per_plane= new double[nplane_x];
	pts_per_plane= new double[nplane_x];
	memcpy(lines_per_plane,_lines_per_plane,nplane_x*sizeof(double));
	memcpy(pts_per_plane,_pts_per_plane,nplane_x*sizeof(double));
	maptable=NULL;
    }
    void GSpacePlaneLoad(int idx, double *line_load, double *pt_load);
    ~SCalcMap(){
      delete [] lines_per_plane;
      delete [] pts_per_plane;
    }
    void pup(PUP::er &p)
	{
	    CkArrayMap::pup(p);
	    p|nplane_x;
	    p|max_states;
	    p|max_planes;
	    p|gs;
	    p|symmetric;
	    p|totalload;
	    maptable=NULL;
	    if (p.isUnpacking()) {
		lines_per_plane= new double[nplane_x];
		pts_per_plane= new double[nplane_x];
	    }	    
	    p(lines_per_plane,nplane_x);
	    p(pts_per_plane,nplane_x);
	}
    void makemap();
    int procNum(int, const CkArrayIndex &);
//    int slowprocNum(int, const CkArrayIndex &);
    int slowprocNum(int, const CkArrayIndex4D &);
};
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
#include "../../include/CPcharmParaInfo.h"
class CPcharmParaInfoGrp: public Group {
 public:
    CPcharmParaInfoGrp(CkMigrateMessage *m) {}
    CPcharmParaInfoGrp(CPcharmParaInfo &s);
    ~CPcharmParaInfoGrp();
    CPcharmParaInfo *cpcharmParaInfo;
};
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
// Using non-symm to do the symm case if defined : Should get rid of 
//#define _CP_TESTSYMM_

// Not doing full CPAIMD computation if undefined : Should get rid of
//#define _FULL_CPAIMD_ 

// Not supporting the full interface reading : Should get rid of
//#define _FIX_FOR_LOCAL_ATOMS_ONLY_ 

/******** General VERBOSE flags ********/
//#define _CP_GS_VERBOSE_
//#define _CP_RS_VERBOSE_
//#define _CP_UTIL_VERBOSE_

/******** Debugging flags ********/
//#define _CP_DEBUG_OLDFORCE_
//#define _CP_DEBUG_NEWFORCE_
//#define _CP_DEBUG_COEFFFINAL_
//#define _CP_DEBUG_LMAT_
//#define _CP_DEBUG_SMAT_
//#define _CP_DEBUG_TMAT_
//#define  _CP_DEBUG_NLMAT_

/******** Optimization flags ********/
//#define _CP_USE_BLAS_GATHER_
//#define _CP_GROUP_REDUCTION_

/******** User Trace Events *********/
#define DoFFTContribute_  100
#define doRealFwFFT_  110
#define doRealBwFFT_  120
#define GspaceFwFFT_  130
#define GspaceBwFFT_  140
#define RhoRtoGxzFFT_  150
#define RhoRtoGyFFT_  151
#define RhoDivRhoXFFT_  152
#define RhoDivRhoYFFT_  153
#define RhoDivRhoZFFT_  154
#define VksofGFFT_  160
#define VksofRFFT_  161
#define AcceptStructFact_  162
//200-300 reserved for paircalculator
#define IntegrateModForces_  1000
#define Scalcmap_            2000
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/*
 * Initialization routines to initialize the GSpace, Real Space,
 * Particle space planes, CP_Rho_GSpacePlane & CP_Rho_RealSpacePlane arrays 
   and the S_Calculators
 */
//============================================================================

void lst_sort_clean(int , int *, int *);
void init_planes(size2d,int, int, int,int,int,int,CPcharmParaInfo *);
void init_rho(size2d , int, int,int);
void control_physics_to_driver();
void get_grp_params(int natm_nl, int numSfGrps, int indexSfGrp, int planeIndex,
		    int *n_ret, int *istrt_ret, int *iend_ret);
int atmGrpMap(int istart, int nsend, int listsize, int *listpe, int AtmGrp, 
              int dup, int planeIndex);
int cheesyhackgsprocNum(CPcharmParaInfo *sim,int state, int plane);
void hackGSpacePlaneLoad(CPcharmParaInfo *sim,int , double *, double *);

//============================================================================

#include "cpaimd.decl.h"

#endif
