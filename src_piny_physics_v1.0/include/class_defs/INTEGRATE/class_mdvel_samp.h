//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                                                                          
//                      class_mdvel_samp.h                                  
//                                                                          
//    Class definition for atom and thermostat velocity resampling          
//                                                                          
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

#ifndef _MDVEL_SAMP_
#define _MDVEL_SAMP_

class MDVEL_SAMP {
 public:
  int ivel_smpl_on;            /* Opt: Periodic atm vel resampl opt   */
  int ivel_scale_on;           /* Opt: Periodic atm vel resampl opt   */
  int nvx_smpl;                /* Num: Freq of atm vel resampling     */
  int nvx_scale;               /* Num: Freq of atm vel rescaling      */
  int nvnhc_smpl;              /* Num: Freq of atm NHC vel resamp     */
  long iseed,iseed2;           /* Num: Random seeds                   */
  double qseed;                /* Num: Real seed for essl ran()       */

//=============================================================================
// Default constructor/destructor

   MDVEL_SAMP(){}
  ~MDVEL_SAMP(){}

//============================================================================
// Pack/Unpack

#ifdef PUP_ON
  void pup(PUP::er &p){

    // PUP ints

    p | ivel_smpl_on;
    p | ivel_scale_on;
    p | nvx_smpl;     
    p | nvx_scale;    
    p | nvnhc_smpl;   
    p | iseed;
    p | iseed2;       

    // PUP doubles
    
    p | qseed;
#ifdef _PARALLEL_DEBUG_        
    if (p.isUnpacking())
     state_class_out ();
#endif  
  } // end pack/unpack
#endif

//==============================================================================
// Print out state of the class

  void state_class_out(){
    
     char fileName [255];
     sprintf (fileName, "%d_mdvel_samp.state", CkMyPe());
     FILE *fp;  fp = fopen(fileName,"w");

     fprintf(fp,"mdvel_samp:  ivel_samp_on %d\n",ivel_smpl_on);
     fprintf(fp,"mdvel_samp:  ivel_scale_on %d\n",ivel_scale_on);
     fprintf(fp,"mdvel_samp:  nvx_samp %d\n",nvx_smpl);
     fprintf(fp,"mdvel_samp:  nvx_scale %d\n",nvx_scale);
     fprintf(fp,"mdvel_samp:  nvnhc_smpl %d\n",nvnhc_smpl);
     fprintf(fp,"mdvel_samp:  iseed1 %d\n",iseed);
     fprintf(fp,"mdvel_samp:  iseed2 %d\n",iseed2);

     fprintf(fp,"mdvel_samp:  qseed %.12g\n",qseed);

     fclose(fp);
  }// end member function
  
}; // end class definition
//==============================================================================

#ifdef PUP_ON
PUPmarshall(MDVEL_SAMP);
#endif

#endif
//==============================================================================

