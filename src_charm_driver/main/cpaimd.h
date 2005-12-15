//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \file cpaimd.h 
 * Some basic data structures and the array map classes are defined
 * here. These maps are used to map the array elements to the correct
 * processors.  Please read ../doc/README for a detailed
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

#ifndef CMK_OPTIMIZE
#define TRACE_ON_STEP 4
#define TRACE_OFF_STEP 7
#endif

#ifndef CmiMemcpy
#define CmiMemcpy(dest, src, size) memcpy((dest), (src), (size))
#endif

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \brief The class which creates the main chare. 
 * 
 *
 *
 */
//============================================================================
 
class main : public Chare {
 public:
    main(CkMigrateMessage *m) { }
    main(CkArgMsg *);
    void doneInit(CkReductionMsg *msg);

};
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \brief Helper class for map hashtables copied from femrefine.C.
 *
 */
//============================================================================
 
class intdual {
 private:
    int x, y;
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


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \brief Class used for instantiation of G-space group objects.
 *
 *
 */
//============================================================================
 
class GSMap: public CkArrayMap {
  
  int nchareG;
  double *lines_per_chareG;
  double *pts_per_chareG;
  double state_load;
  CkHashtableT<intdual, int> *maptable;
 public:
  GSMap() { state_load = 0.0; }
  GSMap(int _nchareG,double *_lines_per_chareG, double *_pts_per_chareG): 
        nchareG(_nchareG)   
      { 
	  state_load = 0.0; 
	  lines_per_chareG= new double[nchareG];
	  pts_per_chareG= new double[nchareG];
	  CmiMemcpy(lines_per_chareG,_lines_per_chareG,nchareG*sizeof(double));
	  CmiMemcpy(pts_per_chareG,_pts_per_chareG,nchareG*sizeof(double));
	  maptable=NULL;
      }
  int procNum(int, const CkArrayIndex &);
// int slowprocNum(int, const CkArrayIndex2D &);
  void makemap();
  void GSpacePlaneLoad(int idx, double *line_load, double *pt_load);
    void pup(PUP::er &p)
	{
	    CkArrayMap::pup(p);
	    p|nchareG;
	    p|state_load;
	    maptable=NULL;
	    if (p.isUnpacking()) {
		lines_per_chareG= new double[nchareG];
		pts_per_chareG= new double[nchareG];
	    }	    
	    p(lines_per_chareG,nchareG);
	    p(pts_per_chareG,nchareG);
	}

  ~GSMap(){
    delete [] lines_per_chareG;
    delete [] pts_per_chareG;
  }
};
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \brief Class used for instantiation of real-space group objects.
 *
 *
 */
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
//============================================================================
/** \brief Class used for instantiation of pair-calculator group objects.
 *
 *
 */
//============================================================================

class SCalcMap : public CkArrayMap {
  int nchareG;
  double *lines_per_chareG;
  double *pts_per_chareG;
  CkHashtableT<intdual, int> *maptable;
  int max_states, max_planes, gs;
  CmiBool symmetric;
  double totalload;
    
 public:
    

    SCalcMap(int _nstates, int _nchareG,  int gs, CmiBool _flag, int _nplanes, 
             double *_lines_per_chareG, double *_pts_per_chareG) { 
        this->gs   = gs;
        nchareG    = _nchareG;
        max_states = _nstates;
        max_planes = _nplanes;
        symmetric  = _flag;
	totalload  = 0.0;
	lines_per_chareG = new double[nchareG];
	pts_per_chareG   = new double[nchareG];
	CmiMemcpy(lines_per_chareG,_lines_per_chareG,nchareG*sizeof(double));
	CmiMemcpy(pts_per_chareG,_pts_per_chareG,nchareG*sizeof(double));
	maptable=NULL;
    }
    void GSpacePlaneLoad(int idx, double *line_load, double *pt_load);
    ~SCalcMap(){
      delete [] lines_per_chareG;
      delete [] pts_per_chareG;
    }
    void pup(PUP::er &p)
	{
	    CkArrayMap::pup(p);
	    p|nchareG;
	    p|max_states;
	    p|max_planes;
	    p|gs;
	    p|symmetric;
	    p|totalload;
	    maptable=NULL;
	    if (p.isUnpacking()) {
		lines_per_chareG= new double[nchareG];
		pts_per_chareG= new double[nchareG];
	    }	    
	    p(lines_per_chareG,nchareG);
	    p(pts_per_chareG,nchareG);
	}
    void makemap();
    int procNum(int, const CkArrayIndex &);
//    int slowprocNum(int, const CkArrayIndex &);
    int slowprocNum(int, const CkArrayIndex4D &);
    int slowprocNum2(int, const CkArrayIndex4D &);
};
//============================================================================

/**
 * provide procnum mapping for RhoR
 */
class RhoRSMap : public CkArrayMap {
  public:
    RhoRSMap(int NN,int ioff):N(NN), off(ioff) {}
    int procNum(int arrayHdl, const CkArrayIndex &idx){
      CkArrayIndex2D idx2d = *(CkArrayIndex2D *) &idx;
      return (((N * idx2d.index[0]) + idx2d.index[1] + off) % CkNumPes());
    }
    void pup(PUP::er &p)
      {
	CkArrayMap::pup(p);
	p|N;
	p|off;
      }
    
  private:
    int N;
    int off;
};


/**
 * provide procnum mapping for RhoG
 */
class RhoGSMap : public CkArrayMap {
  public:
    RhoGSMap(int NN, int ioff, int iavoid, int iavoid_off):N(NN), off(ioff), avoid(iavoid), avoid_off(iavoid_off) {}
    int procNum(int arrayHdl, const CkArrayIndex &idx){
      CkArrayIndex2D idx2d = *(CkArrayIndex2D *) &idx;
//      return (((N * idx2d.index[0]) + idx2d.index[1] + off) % CkNumPes());
      int pe=(((N * idx2d.index[0])  + off) % CkNumPes());
      // avoid PEs favored by the array characterized by the avoid and avoid_off parms
      if(avoid>1 && (pe-avoid_off)%avoid==0)
      {
	  pe+=1;
      }
      return pe;
    }
    void pup(PUP::er &p)
      {
	CkArrayMap::pup(p);
	p|N;
	p|avoid;
	p|off;
	p|avoid_off;
      }

  private:
    int N;
    int off;
    int avoid;
    int avoid_off;
};

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//============================================================================
/**
 * \brief class CPcharmParaInfoGrp.
 *
 *
 */
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


/******** Optimization flags ********/
 //#define _CP_USE_BLAS_GATHER_

/******** User Trace Events *********/
#define DoFFTContribute_     100
#define doRealFwFFT_         110
#define doRealBwFFT_         120
#define GspaceFwFFT_         130
#define GspaceBwFFT_         140
#define RhoRtoGFFT_          150

#define PostByrdfwFFTGtoR_   152
#define BwFFTRtoG_           153
#define ByrdanddoFwFFTGtoR_  154
#define HartExcVksG_         160
#define divRhoVksGspace_     161
#define AcceptStructFact_    162
#define fwFFTGtoR0_          170
#define fwFFTGtoRnot0_       171
#define GradCorrGGA_         180
#define WhiteByrdFFTX_       181
#define WhiteByrdFFTY_       182
#define WhiteByrdFFTZ_       183

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
void init_pair_calculators(int nstates, int indexSize, int *indexZ, int doublePack, CPcharmParaInfo *sim);
void init_ortho_chares(int, int, int *);

void init_commlib_strategies(int, int);
void lst_sort_clean(int , int *, int *);
void init_state_chares(size2d,int,int,int,int,CPcharmParaInfo *);
void init_rho_chares(size2d ,CPcharmParaInfo*);
void control_physics_to_driver();
void get_grp_params(int natm_nl, int numSfGrps, int indexSfGrp, int planeIndex,
		    int *n_ret, int *istrt_ret, int *iend_ret);
int atmGrpMap(int istart, int nsend, int listsize, int *listpe, int AtmGrp, 
              int dup, int planeIndex);
int cheesyhackgsprocNum(CPcharmParaInfo *sim,int state, int plane);
void hackGSpacePlaneLoad(CPcharmParaInfo *sim,int , double *, double *);

//============================================================================

// stuff to be include before the decl or else
#include "../../include/Atoms.h"
#include "energy.h"
#include "cpaimd.decl.h"

#endif
