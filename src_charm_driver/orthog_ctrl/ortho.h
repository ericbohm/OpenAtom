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
 *  Ortho is decomposed by orthoGrainSize.  
 * 
 *  We restrict orthograin to be a factor of sGrainsize then we have
 *  no section overlap issues.  Thereby leaving us with ortho sections
 *  that need a simple tiling split of the sgrain sections.  Mirrored
 *  by a stitching of the submatrix inputs for the backward path.
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

#include "debug_flags.h"
#include "orthoConfig.h"
#include "ortho.decl.h"
#include "pcSectionManager.h"
#include "CLA_Matrix.h"
#include "mpi-interoperate.h"
using namespace cp::ortho; ///< @todo: Temporary, till Ortho classes live within namespace ortho

#ifndef _ortho_h_
#define _ortho_h_


#define INVSQR_TOLERANCE	1.0e-15
#define INVSQR_MAX_ITER		10

#ifdef CP_PAIRCALC_USES_COMPLEX_MATH
#define myabs abs
#else
#define myabs std::abs
#endif

extern bool fakeTorus;
extern int numPes;
class initCookieMsg : public CkMcastBaseMsg, public CMessage_initCookieMsg {
};

class orthoMtrigger : public CkMcastBaseMsg, public CMessage_initCookieMsg {
};

class ExtendedOrtho : public CBase_ExtendedOrtho {
  public:
    ExtendedOrtho();
    ExtendedOrtho(CkMigrateMessage *m) {}
    void lambdaSentToDiagonalizer();
};

/** @addtogroup Ortho
  @{
 */
class Ortho : public CBase_Ortho
{
  public:
    /// Default empty constructor. For?
    Ortho() {}
    Ortho(CkMigrateMessage *m){}
    ~Ortho();
    Ortho(int m, int n, 
        CLA_Matrix_interface matA1, CLA_Matrix_interface matB1, CLA_Matrix_interface matC1,
        CLA_Matrix_interface matA2, CLA_Matrix_interface matB2, CLA_Matrix_interface matC2, 
        CLA_Matrix_interface matA3, CLA_Matrix_interface matB3, CLA_Matrix_interface matC3,
        orthoConfig &_cfg,
        CkArrayID step2Helper,
	  int timeKeep, CkGroupID _oMCastGID, CkGroupID _oRedGID, bool lsda);

    /// Trigger the creation of appropriate sections of paircalcs to talk to. Also setup internal comm sections
    void makeSections(const pc::pcConfig &cfgSymmPC, const pc::pcConfig &cfgAsymmPC, CkArrayID symAID, CkArrayID asymAID);

    /// Symmetric PCs contribute data that is summed via this reduction to deposit a portion of the S matrix with this ortho, triggering S->T
    void start_calc(CkReductionMsg *msg);
    /// Triggers the matrix multiplies in step 1 of an ortho iteration.
    void do_iteration(void);
    /// Used when array broadcasts in ortho are delegated to comlib/CkMulticast so as to not involve all PEs in bcast
    void do_iteration(orthoMtrigger *m) { do_iteration(); /* do not delete nokeep msg */ }
    /// Triggers step 2, and optionally step 3 (if ortho helpers are being used)
    void step_2();
    /// Receives the results of the call to OrthoHelper from step_2 
    void recvStep2(CkDataMsg *msg); //double *step2result, int size);
    /// Triggers step 3 in the S->T process
    void step_3();
    /// Computes square of the residuals and contributes to a reduction rooted at Ortho(0,0)::collect_error()
    void tolerance_check();
    /// Computes the RMS error and either launches the next ortho iteration (if needed) or calls collect_results
    void collect_error(CkReductionMsg *msg);
    /// Used when ortho redn/bcasts are delegated to comlib/CkMulticast because charm array broadcasts involve all PEs
    void collect_results(orthoMtrigger *m) { collect_results(); /* do not delete nokeep msg */ }
    /// Computes walltimes, prints simulation status msgs and calls resume() if more openatom iterations are needed
    void collect_results(void);
    /// Sends results (T matrix) to the symm PCs (and also the asymms if this is dynamics)
    void resume();
    /// Used in dynamics, to send the results of S->T to the asymm paircalc instance which will use them
    void sendOrthoTtoAsymm();

    // Accepts lamda reduced from the asymm PC instance. In min, acts as via point and mcasts lambda back to the asymm PCs. In dynamics, triggers computation of gamma = lambda x T
    void acceptSectionLambda(CkReductionMsg *msg);
    void acceptDiagonalizedLambda(int nn, internalType* rmat);
    void transferControlToMPI();
    void waitForQuiescence(int n);
    /// Used in dynamics, to accept computed gamma and send it to the asymm PC instance. Also sends T if it hasnt yet been sent
    void gamma_done();

    /// S should be equal to 2I. This returns max value of deviation from that in this ortho's portion of the S matrix. 
    inline double array_diag_max(int sizem, int sizen, internalType *array);
    /// Called on ortho(0,0). Checks if PsiV update is needed based on the max deviation in S received via a redn across all orthos. Notifies GSpaceDriver if so. Called periodically in start_calc only for dynamics 
    void maxCheck(CkReductionMsg *msg);
    /// Once all GSpaceDriver chares are notified, they resume Ortho execution via a redn broadcast to this method
    void resumeV(CkReductionMsg *msg);
    //output eigen vectors when diagonalizer in use
    void output_eigen(int xoff, int yoff, int size, double *data);
    /// Dumps the T matrix to an appropriately named file
    void print_results(void);
    /// pack/unpack method
    virtual void pup(PUP::er &p);
    void orthoCookieinit(initCookieMsg *msg) { CkGetSectionInfo(orthoCookie,msg); delete msg; }
    /// called from each CLA_Matrix array (3 per multiplication, 3 mults)
    void all_ready() { if(++num_ready == 9) thisProxy.ready(); }


    /// Startup/Init synchronization. When all elements (PC, CLA_Matrix etc) are ready, first iteration is triggered

    void ready()
    { 
      // got_start comes from upstream PC reduction the last of got_start and
      // when all the CLA_Matrix arrays are ready, computation for first iteration is triggered
      num_ready = 10;
      if(got_start)
        do_iteration();
    }

    /// Static methods used as callbacks. Could be replaced by CkCallbacks
    static inline void step_2_cb(void *obj) { ((Ortho*) obj)->step_2(); }
    static inline void step_3_cb(void *obj) { ((Ortho*) obj)->step_3(); }
    static inline void gamma_done_cb(void *obj) { ((Ortho*) obj)->gamma_done(); }
    static inline void tol_cb(void *obj) 
    {
      ((Ortho*) obj)->step3done=true;
      if(((Ortho*) obj)->parallelStep2)
      { 
        //if step2 is done do this now, otherwise step2 will trigger
        if(((Ortho*) obj)->step2done)
          ((Ortho*) obj)->tolerance_check();
      }
      else
        ((Ortho*) obj)->tolerance_check();
    }

    bool parallelStep2;
    bool step2done;
    bool step3done;

  private:
    bool lsda;
    orthoConfig cfg;
    int timeKeep;
    internalType *orthoT; // only used on [0,0]
    internalType *ortho; //only used on [0,0]
    internalType *templambda;
    int numGlobalIter; // global leanCP iterations
    // used in each element
    int iterations; //local inv_sq iterations
    CProxySection_Ortho multiproxy;
    CkSectionInfo orthoCookie;
    int num_ready;
    bool got_start;
    int lbcaught;
    /// Section of symmetric PC chare array used by an Ortho chare
    PCSectionManager symmSectionMgr;
    /// Section of asymmetric PC chare array used by an Ortho chare
    PCSectionManager asymmSectionMgr;
    /// Group IDs for the multicast manager groups
    CkGroupID oMCastGID, oRedGID;
    /// The proxy of the step 2 helper chare array
    CProxy_OrthoHelper step2Helper;
    bool toleranceCheckOrthoT; //trigger tolerance failure PsiV conditions
    internalType *A, *B, *C, *tmp_arr;
    int step;
    int m, n;
    double invsqr_tolerance;
    int invsqr_max_iter;
    CLA_Matrix_interface matA1, matB1, matC1, matA2, matB2, matC2, matA3, matB3, matC3;
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

    int eigenRecv;
    double **eigenData;
};




#include "ckcallback-ccs.h" ///< For definition of CkDataMsg

/** 
 * Result of 0.5 * S3 * S2 arrives from helper
 */
inline void Ortho::recvStep2(CkDataMsg *msg)//double *step2result, int size)
{
  // copy our data into the tmp_arr  
  CmiMemcpy(tmp_arr, msg->getData(), m * n * sizeof(internalType));
  step2done=true;
  delete msg;
  //End of iteration check
  if(step3done)
    tolerance_check();
}

/**
 * OrthoT tolerance check util return max value
 */
inline double Ortho::array_diag_max(int sizem, int sizen, internalType *array)
{
  double absval, max_ret;
  if(thisIndex.x!=thisIndex.y)
  { //not diagonal
    max_ret = myabs(array[0]);
    for(int i=0;i<sizem;i++)
    {
      for(int j=0;j<sizen;j++)
      {
        absval = myabs(array[i*sizen+j]);
        max_ret = (max_ret>absval) ? max_ret : absval;
      }
    }//endfor
  }
  else
  { //on diagonal 
    double diag= (lsda) ? 1.0 : 2.0;
    max_ret = myabs(array[0]-diag);
    for(int i=0;i<sizem;i++)
    {
      for(int j=0;j<sizen;j++)
      {
        absval = myabs(array[i*sizen+j]);
        if(i == j)
          absval = myabs(absval - diag);
        max_ret = (max_ret>absval) ? max_ret : absval;
      }
    }//endfor
  }//endif
  return max_ret;
}//end routine

#endif // #ifndef _ortho_h_

