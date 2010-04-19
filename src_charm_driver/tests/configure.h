/*****************************************************************************
 * $Source$
 * $Author$
 * $Date$
 * $Revision$
 *****************************************************************************/

/** \file configure.h
 *
 */

#ifndef _Configure_
#define _Configure_

#include "dictionary.h"
#include "pup.h"

//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
class Config {
//===================================================================================
 public:

  //==================================
  // input or derived parameters
  //----------------------------------

    int numPes;           // num of procs             (driver)
    int numPesPerInstance;// num of procs for mapping
	int test;		  // of one instance          (driver)
    int nstates;          // num of states            (piny)
    int maxIter;          // num of iterations        (piny)
//    int fftopt;           // fft type[fftw/essl]      (piny)
//    int natm_nl;          // num of nonlocal atoms    (piny)
//    int natm_typ;         // number of atm types      (piny)
    int low_x_size;       // g-grid points along x    (piny)
    int high_x_size;      // g-grid points along x    (piny)
//    int numFFTPoints;     // size number of r-grid    (piny)
    int numData;          // size number of g-grid    (piny)
//    int ees_eext_opt;     // ees eext option on/off   (piny)
    int gen_wave;         // generate initial states  (piny)
    int nchareG;          // num state g-space chares (driver)
//    int nchareRhoG;       // num rho g-space chares   (driver)
    int nchareVdW;        // num VanderWalls chares   (driver)
    int scalc_per_plane;  // num of scalcs g-plane    (driver)

    int UberImax;
    int UberJmax;
    int UberKmax;
    int numInstances;
  //==================================

  //==================================
  // density control flags and values
  //----------------------------------
    int nchareHartAtmT;
    int rhoLineOrder;
    int rhorHartpriority;
    int rhogHartpriority;
    int useGHartInsRhoRP;
    int useGHartInsRHart;
    int useRHartInsGHart;
    int rhorpriority;
    int rhogpriority;
    double gExpandFactRho;
    int lbdensity;
    int rhoGHelpers;
    int rhoRsubplanes;
    int rhoSubPlaneBalance;
    int rhoGToRhoRMsgComb;
    int prioEextFFTMsg;
    //==================================
    // density commlib flags
    //----------------------------------
    int useGIns0RhoRP;
    int useGIns1RhoRP;
    int useGIns2RhoRP;
    int useGIns3RhoRP;
    int useGByrdInsRhoRBP;
    int useRInsRhoGP;
    int useRInsIGXRhoGP;
    int useRInsIGYRhoGP;
    int useRInsIGZRhoGP;
  //==================================

  //==================================
  // state control flags and values
  //----------------------------------
   char dataPath[1024];
   int gBucketSize;
   int rBucketSize;
   double gStreamPeriod;
   double rStreamPeriod;
   double gExpandFact;
   int stateOutput;
   int psipriority;
   int prioFFTMsg;
   int rsfftpriority;
   int gsfftpriority;
   int rsifftpriority;
   int gsifftpriority;
   int conserveMemory;
   int lbgspace;
   int doublePack;
   int useGssInsRealP;
   int useMssInsGP;
   char dataPathOut[1024];
  //==================================

  //==================================
  // PC control flags and values
  //----------------------------------
   int usePairDirectSend;
   int PCCollectTiles;
   int PCdelayBWSend;
   int PCstreamBWout;
   int PCstreamFWblock;
   double invsqr_tolerance;
   int invsqr_max_iter;
   int useOrthoDirect;
   int useOrthoHelpers;
   int useOrthoSection;
   int useOrthoSectionRed;
   int lambdaGrainSize;
   int PCSpanFactor;
   int OrthoRedSpanFactor;
   int OrthoMcastSpanFactor;
   int sGrainSize;
   int gemmSplitFWk;
   int gemmSplitFWm;
   int gemmSplitBW;
   int gemmSplitOrtho;
   int orthoGrainSize;
   int useBWBarrier;
   int phantomSym;
   int lbpaircalc;
   int lambdapriority;
   int toleranceInterval;
   int gSpaceSum;
   int numChunks;
   int numChunksSym;
   int numChunksAsym;
   int prioBW;
   int usePairEtoM;
  //==================================

  //==================================
  // NL control flags and values
  //----------------------------------
   int sfpriority;
   int prioNLFFTMsg;
   int rsNLfftpriority;
   int gsNLfftpriority;
   int numSfGrps;
   int numSfDups;
   int launchNLeesFromRho;
   int useGssInsRealPP;
   int useMssInsGPP;
  //==================================

  //==================================
  // General control flags and values
  //----------------------------------
   int useTimeKeeper;
   int fftprogresssplit;
   int fftprogresssplitReal;
   int useCommlib;
   int useGMulticast;
   int numMulticastMsgs;
   int useCommlibMulticast;
   int atmOutput;
  //==================================

  //==================================
  // Mapping control flags and values
  //----------------------------------
   int torusMap;
   int forceMappingAxis;
   int useCuboidMap;
   int useCuboidMapRS;
   int useStrictCuboid;
   int useCentroidMap;
   int useCentroidMapRho;
   int Gstates_per_pe;
   int Rstates_per_pe;
   int loadMapFiles;
   int dumpMapFiles;
   int dumpMapCoordFiles;
   int useRhoExclusionMap;
   int excludePE0;
   int useReductionExclusionMap;
   int fakeTorus;
   int torusDimNX;	// use these to do torus logic testing
   int torusDimNY;
   int torusDimNZ;
   int torusDimNT;

  //==================================

  //==================================
  // Network progress checking frequencies (useful only on BGL where you have to periodically let call CmiNetworkProgress)
  //----------------------------------
  //  nfreq_classname_methodname (Class and method names may have been shortened a bit)
  int nfreq_cpintegrate;        ///< CPINTEGRATE::CP_integrate_min_STD, CPINTEGRATE::CP_integrate_min_CG
  int nfreq_cplocal_hartext;    ///< CPLOCAL::CP_hart_eext_calc
  int nfreq_cplocal_eeshart;    ///< CPLOCAL::eesHartEextGchare
  int nfreq_cplocal_eesewald;   ///< CPLOCAL::eesEwaldGchare
  int nfreq_cpnonlocal_eke;     ///< CPNONLOCAL::CP_eke_calc
  int nfreq_cpnonlocal_eesfwd;  ///< CPNONLOCAL::eesProjGchare, CPNONLOCAL::eesYlmOnD 
  int nfreq_cpnonlocal_eesbk;   ///< CPNONLOCAL::eesPsiForcGspace
  int nfreq_xcfnctl;            ///< CPXCFNCTS::CP_exc_calc, CPXCFNCTS::CP_getGGAFunctional
  //==================================

  //==================================
  // Class Functions
  //----------------------------------
  Config(){};
  ~Config(){};
   void readConfig(char* input_name,int nstates_in, int nplanes_in, int maxIter_in,int numPes_in);
   void readStateInfo(int &, int &, int &, int &, int &, int &,
                      const char *, int);
   void simpleRangeCheck();
   void rangeExit(int, const char *, int);
   void Finale(int, int);

   void set_config_dict_fun    (int *, DICT_WORD **);
   void set_config_dict_gen    (int *, DICT_WORD **);
   void set_config_dict_rho    (int *, DICT_WORD **);
   void set_config_dict_state  (int *, DICT_WORD **);
   void set_config_dict_pc     (int *, DICT_WORD **);
   void set_config_dict_nl     (int *, DICT_WORD **);
   void set_config_dict_map    (int *, DICT_WORD **);
   void set_config_dict_nfreq  (int *, DICT_WORD **);

   void set_config_params_gen  (DICT_WORD *, char *, char *);
   void set_config_params_rho  (DICT_WORD *, char *, char *, int);
   void set_config_params_state(DICT_WORD *, char *, char *, int);
   void set_config_params_pc   (DICT_WORD *, char *, char *);
   void set_config_params_nl   (DICT_WORD *, char *, char *, int);
   void set_config_params_map  (DICT_WORD *, char *, char *);
   void set_config_params_nfreq(DICT_WORD *, char *, char *);
   void guesstimateParmsConfig (DICT_WORD *, DICT_WORD *, DICT_WORD *,
				     DICT_WORD *, DICT_WORD *,DICT_WORD *);
   int approxFactor(int nstates,int &sGrainSize, int &oGrainSize,int numPes);
   void write_cpaimd_config    (FILE *, DICT_WORD *, int, char *);
   void load_cpaimd_config (DICT_WORD *, int, PINY_NAME *, PINY_NAME *, int, int *);
   bool isPow2(int input){
     int y=0;
     for(int x=0;x<32;x++){
       y = 1<<x;
       if(y==input){return true;}
     }//endfor
     return false;
   }
//============================================================================
  //==================================

//-----------------------------------------------------------------------------------
   }; // end class
//===================================================================================

#ifndef _COOL_CONVERSION_ON_
PUPbytes(Config)
#endif

//===================================================================================

#endif
