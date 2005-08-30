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


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class GSMap: public CkArrayMap {
  

  int nplane_x;
  double *lines_per_plane;
  double *pts_per_plane;
  double state_load;
 public:
  GSMap() { state_load = 0.0; }
  GSMap(int _nplane_x, double *_lines_per_plane, double *_pts_per_plane): 
        nplane_x(_nplane_x)   
    { state_load = 0.0; 
    lines_per_plane= new double[nplane_x];
    pts_per_plane= new double[nplane_x];
    memcpy(lines_per_plane,_lines_per_plane,nplane_x*sizeof(double));
    memcpy(pts_per_plane,_pts_per_plane,nplane_x*sizeof(double));
    }
  int procNum(int, const CkArrayIndex &);
  void GSpacePlaneLoad(int idx, double *line_load, double *pt_load);
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
    int max_states, max_planes, gs;
    CmiBool symmetric;
    double totalload;
    
 public:
    
    SCalcMap(int nstates, int nplanes, int gs, CmiBool flag, int _nplane_x, 
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
    }
    void GSpacePlaneLoad(int idx, double *line_load, double *pt_load);
    ~SCalcMap(){
      delete [] lines_per_plane;
      delete [] pts_per_plane;
    }
    int procNum(int, const CkArrayIndex &);
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
