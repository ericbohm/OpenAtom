/*****************************************************************************
 * $Source$
 * $Author$
 * $Date$
 * $Revision$
 *****************************************************************************/

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \file ortho.h
 *
 *  Ortho is currently decomposed by sGrainSize.  Which makes sense for
 *  the old statewise only decomposition of the paircalculator.  For
 *  the new 4d decomposition it may make sense to decompose the
 *  paircalculator more along the numPoints axis than the states axis.
 *  Which may result in comparitively fewer Ortho objects.  Thereby
 *  further motivating the need to chop up Ortho by some unit smaller
 *  than grainsize.  So as to avoid Amdaling on the orthonormalization
 *  process.
 *
 *  If we allow orthograin to be entirely distinct from sGrainSize we
 *  have arbitrary overlap situations between orthograins and
 *  scalcgrains.  Supporting that requires a rather complicated (and
 *  therefore bug prone) contiguousReducer reduction/multicast
 *  stitcher/splitter implementation to reduce data from scalc->ortho
 *  and multicast data from ortho->scalc.  
 * 
 *  We don't want to do that if we don't have to.  If we restrict
 *  orthograin to be a factor of sGrainsize then we have no section
 *  overlap issues.  Thereby leaving us with ortho sections that need a
 *  simple tiling split of the sgrain sections.  Mirrored by a
 *  stitching of the submatrix inputs for the backward path.  
 *
 *  This can be accomplished manually within the current codebase with
 *  some waste in data replication and computation replication to
 *  handle the splitting/stiching operations.  
 *
 *  A more efficient implementation would adopt the multicast manager
 *  group model of building a tree of participants for these
 *  operations.  The reduction side from the PC would be broken up into
 *  multiple reductions, one for each orthograin within the sgrain.
 *  With a separate contribution for each orthograin.  The multicast
 *  requires us to stitch together the input matrices into one per
 *  sgrain section.  This might be accomplished in two stages, one in
 *  which the stitching is done, and a second in which the stitched
 *  sgrainsize matrices are multicast.  The alternative is to just
 *  multicast the orthograin submatrices where needed and have each
 *  scalc do its strided copying stitching.  As stitching is not
 *  computationally intensive, this may be the simplest and fastest
 *  solution.  The second approach allows you to simply use the
 *  reductions and multicasts as mirror uses of the tree.  Where each
 *  little ortho can run once it gets its input, while the scalcs would
 *  have to assemble their inputs from multiple multicasts.  
 *
 *  Implementation details for this require that each ortho object
 *  participate in a section which has a section multicast client
 *  directed to the sGrainSize PC section.  The converse PC sGrainSize
 *  elements will have an array of section cookies, one for each of the
 *  subsections for all orthograin elements within the sGrain.  The
 *  forward path of the PC will contribute its orthograin tile (via a
 *  strided contribute) which will end up at the correct ortho object.
 *
 *  Note: these PC sections must include all 4th dim blocks.  
 *
 *  OrthoHelper can be used to perform the 2nd of the multiplies in the
 *  3 step S->T process in parallel with the 3rd multiply.  If used,
 *  the results of multiply 1 are sent from ortho[x,y] to
 *  orthoHelper[x,y].  The results are then returned to ortho[x,y].
 *  The last of step2 or step3 will then trigger step4.  Due to the
 *  copy and communication overhead this is only worth doing if the
 *  number of processors is greater than 2 * the number of ortho
 *  chares.
 * 
 *  Allowing sgrainsize choices which are nstates % sgrainsize != 0
 *  forces us to handle remainder logic.  To avoid overlap/straddle
 *  issues between ortho and PC, we still enforce sgrainssize %
 *  orthograinsize ==0.  Complexity cost here comes in two forms.
 *  1. Now ortho tiles are not guaranteed to be of uniform size.  The
 *  remainder states which will reside in the last row and column will
 *  result in tiles larger than orthograinsize*orthograinsize.
 *  2. Ortho tiles are not guaranteed to be square.  Ortho tiles for
 *  the last row and column of PC will have M x N size where M != N.
 *
 *  The total multiply itself will still of course be nstates X
 *  nstates.
 ******************************************************************************/

#include "ortho.decl.h"
#include "paircalc/pairCalculator.h" 
#include "uber/Uber.h"
#include "CLA_Matrix.h"
#include "cpaimd.h"

#ifndef _ortho_h_
#define _ortho_h_


#define INVSQR_TOLERANCE	1.0e-15
#define INVSQR_MAX_ITER		10

extern CkVec <MapType2> OrthoImaptable;
extern CkHashtableT <intdual, int> Orthomaptable;
extern bool fakeTorus;
extern int numPes;
class initCookieMsg : public CkMcastBaseMsg, public CMessage_initCookieMsg {
};

class orthoMtrigger : public CkMcastBaseMsg, public CMessage_initCookieMsg {
};

class Ortho : public CBase_Ortho{
 public:
  Ortho() {}
  Ortho(CkMigrateMessage *m){}
  ~Ortho(){
    delete [] A;
    delete [] B;
    delete [] C;
    delete [] tmp_arr;
    if(thisIndex.x==0 && thisIndex.y==0)
      {
	delete [] ortho;
	delete [] orthoT;
	delete [] wallTimeArr;
      }
  }
  // catch lambda for later non_minimization use
  void acceptAllLambda(CkReductionMsg *msg); 

  void step_2_send(void);

  void recvStep2(double *step2result, int size);

  void sendOrthoTtoAsymm();

  // get our copy of the pcproxy
  void setPCproxy(CProxySection_PairCalculator inproxy);

  // catch lambda for later non_minimization use
  void acceptSectionLambda(CkReductionMsg *msg); 

  void maxCheck(CkReductionMsg *msg);

  void resumeV(CkReductionMsg *msg);

  void resume();

  void multiplyForGamma(double *orthoT, double *lambda, double *gamma,int n);

  /* start step 1 on proxy 1, the callback will be for step 2
   * S2 = 3 * I - T * S1
   * currently A has T, B has S1, need to construct 3*I in C
   */
  void do_iteration(void);

  void do_iteration(orthoMtrigger *m){
  //============================================================================
  // Do not delete msg. Its a nokeep.
  //============================================================================
 
    do_iteration();
  }

  void makeSections(int indexSize, int *indexZ);

  void print_results(void){
    char outname[80];
    snprintf(outname,80,"tmatrix_t:%d_%d_%d.out",numGlobalIter,thisIndex.x,thisIndex.y);
    FILE *outfile = fopen(outname, "w");
	for(int i=0; i<m; i++){
	  for(int j=0; j<n; j++){
           fprintf(outfile, "%d %d %.12g \n",i+thisIndex.x*n+1,j+thisIndex.y*n+1, A[i*n+j]);
	  }
	}
	fclose(outfile);
  }

  void collect_results(orthoMtrigger *m){
  //============================================================================
  // Do not delete msg. Its a nokeep.
  //============================================================================
 
    collect_results();
  }


  void collect_results(void);

  virtual void pup(PUP::er &p);

  void collect_error(CkReductionMsg *msg);

  void orthoCookieinit(initCookieMsg *msg){CkGetSectionInfo(orthoCookie,msg); delete msg;}

  void start_calc(CkReductionMsg *msg);

  /* called from each CLA_Matrix array (3 per multiplication, 3 mults) */
  void all_ready(){
    num_ready++;
    if(num_ready == 9)
      thisProxy.ready();
  }

  void ready(){ // startup initialization synchronization
      // got_start comes from upstream PC reduction the last of got_start and
      // CLA_Matrix readies triggers computation for first iteration
      num_ready = 10;
      if(got_start)
      {
	  do_iteration();
      }
  }

/**
 * OrthoT tolerance check util return max value
 */
  inline double array_diag_max(int sizem, int sizen, double *array){
      double absval, max_ret;
      if(thisIndex.x!=thisIndex.y){ //not diagonal
          max_ret=fabs(array[0]);          
	  for(int i=0;i<sizem;i++){
          for(int j=0;j<sizen;j++){
 	    absval=fabs(array[i*sizen+j]);
	    max_ret = (max_ret>absval) ? max_ret : absval;
	  }}//endfor
      }else{ //on diagonal 
          absval=fabs(fabs(array[0]-2.0));
          max_ret = absval;
	  for(int i=0;i<sizem;i++){
	  for(int j=0;j<sizen;j++){
	    absval=fabs(array[i*sizen+j]);
	    if(i!=j){
 	       max_ret = (max_ret>absval) ? max_ret : absval;
	     }else{
	       absval=fabs(absval-2.0);
	       max_ret = (max_ret>absval) ? max_ret : absval;
	     }//endif
	  }}//endfor
      }//endif
      return max_ret;
  }//end routine

  Ortho(int m, int n, CLA_Matrix_interface matA1,
   CLA_Matrix_interface matB1, CLA_Matrix_interface matC1,
   CLA_Matrix_interface matA2, CLA_Matrix_interface matB2,
   CLA_Matrix_interface matC2, CLA_Matrix_interface matA3,
   CLA_Matrix_interface matB3, CLA_Matrix_interface matC3, 
	int timeKeep, UberCollection );

  void tolerance_check(void);
  void step_2();
  void step_3();
  void gamma_done();

  static inline void step_2_cb(void *obj){
    ((Ortho*) obj)->step_2();
  }

  static inline void step_3_cb(void *obj){
    ((Ortho*) obj)->step_3();
  }

  static inline void tol_cb(void *obj)
    {
      ((Ortho*) obj)->step3done=true;
      if(((Ortho*) obj)->parallelStep2)
	{ 
	  if(((Ortho*) obj)->step2done)
	    {	//if step2 is done do this now, otherwise step2 will trigger
	      ((Ortho*) obj)->tolerance_check();
	    }
	}
      else
	{
	  ((Ortho*) obj)->tolerance_check();
	}
    }

  static inline void gamma_done_cb(void *obj){
    ((Ortho*) obj)->gamma_done();
  }
  bool parallelStep2;
  bool step2done;
  bool step3done;
 private:
  const UberCollection thisInstance;
  int timeKeep;
  double *orthoT; // only used on [0,0]
  double *ortho; //only used on [0,0]
  double *wallTimeArr;//only used on [0,0]
  int numGlobalIter; // global leanCP iterations

  // used in each element
  int iterations; //local inv_sq iterations
  CProxySection_Ortho multiproxy;
  CkSectionInfo orthoCookie;
  int num_ready;
  bool got_start;
  int lbcaught;
  PairCalcID oPairCalcID1;
  PairCalcID oPairCalcID2;
  bool toleranceCheckOrthoT; //trigger tolerance failure PsiV conditions
  double *A, *B, *C, *tmp_arr;
  int step;
  int m, n;
  double invsqr_tolerance;
  int invsqr_max_iter;
  CLA_Matrix_interface matA1, matB1, matC1, matA2, matB2, matC2, matA3,
   matB3, matC3;
#ifdef _CP_ORTHO_DEBUG_COMPARE_TMAT_
  double *savedtmat;
#endif
#ifdef _CP_ORTHO_DEBUG_COMPARE_LMAT_
  double *savedlmat;
#endif

#ifdef _CP_ORTHO_DEBUG_COMPARE_SMAT_
  double *savedsmat;
#endif
#ifdef _CP_ORTHO_DEBUG_COMPARE_GMAT_
  double *savedgmat;
#endif

};

/**
 * provide procnum mapping for Ortho
 * this old class is now defunct: Abhinav
 *
 * class OrthoMap : public CkArrayMap {
 *  public:
 *   OrthoMap(int NN,int _nOrtho, int _stride):N(NN), nOrtho(_nOrtho), stride(_stride)
    {
      offset=0;
      if(nOrtho<CkNumPes())
        offset=1;  //skip proc 0
    }
#ifndef TOPO_ORTHO
    virtual int procNum(int arrayHdl, const CkArrayIndex &iIndex)
    {
      int *index=(int *) iIndex.data();
      return (stride*(N * index[0] + index[1]) + offset) % CkNumPes();
    }
#else
    virtual int procNum(int arrayHdl, const CkArrayIndex &iIndex)
    {
      int *index=(int *) iIndex.data();
      return (stride*(N * index[0] + index[1]) + offset) % CkNumPes();
    }
#endif
  private:
 *   int N;
 *   int nOrtho;
 *   int offset;
 *   int stride;
 * };
 **/

/*
 * new centroid based ortho map
 * actual map creation in MapTable.C
 */
class OrthoMap : public CkArrayMapTable2 {
  public:
    OrthoMap(UberCollection _instance)
    {
      thisInstance=_instance;
#ifdef USE_INT_MAP
      maptable= &OrthoImaptable[thisInstance.getPO()];
#else
      maptable= &Orthomaptable;
#endif
    }

    ~OrthoMap() { }
    
    void pup(PUP::er &p)
    {
      CkArrayMapTable2::pup(p);
#ifdef USE_INT_MAP
      maptable= &OrthoImaptable[thisInstance.getPO()];
#else
      maptable= &Orthomaptable;
#endif
    }
    
    inline int procNum(int, const CkArrayIndex &iIndex)
    {
      int *index=(int *) iIndex.data();
      int proc;
#ifdef USE_INT_MAP
      proc=maptable->get(index[0],index[1]);
#else
      proc=maptable->get(intdual(index[0],index[1]));
#endif
      CkAssert(proc>=0);
      if(numPes != CkNumPes())
	return(proc%CkNumPes());
      else
	return(proc);
	

    }
};

#endif // #ifndef _ortho_h_
