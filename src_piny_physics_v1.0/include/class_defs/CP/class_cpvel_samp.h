//==========================================================================
//                 Velocity Resampling info                                 
//             {Variables needed for mem allocation:                        
//                                                  }                       
//                                                                          

#ifndef _CPVEL_SAMP_
#define _CPVEL_SAMP_

class CPVEL_SAMP{

  //----------------
  public:
    int ivelc_smpl_on;          // Opt: Periodic PW coef vel resampl    
    int ivelc_scal_on;          // Opt: Periodic PW coef vel rescale    
    int nvc_smpl;               // Num: Freq of PW coef vel resampl    
    int nvcnhc_smpl;            // Num: Freq of PW coef NHC vel resamp 
    int nvc_scal;               // Num: Freq of PW coef vel recale     
    int iauto_vc_scal_opt;      // Opt: auto rescale options on/off    
    long iseed,iseed2;          // Num: Random seeds                   
    double qseed;               // Num: Real seed for essl ran()       
    double vc_scal_tol;         // Num: tol of auto rscale of coef vel 
    double div_scal;

    //----------------
    //con-destruct:
    CPVEL_SAMP(){
      ivelc_smpl_on = 0;
      ivelc_scal_on = 0;
      nvc_smpl      = 0;
      nvcnhc_smpl   = 0;
      nvc_scal      = 0;
      iseed         = 0;
      iseed2        = 0;
      iauto_vc_scal_opt = 0;
      qseed         = 0;
      vc_scal_tol   = 0;
      div_scal      = 0;
    };
    ~CPVEL_SAMP(){};

#ifdef PUP_ON
    //----------------
    //pupping
    void pup(PUP::er &p){
      //pupping ints
      p | ivelc_smpl_on;
      p | ivelc_scal_on;
      p | nvc_smpl;
      p | nvcnhc_smpl;
      p | nvc_scal;
      p | iseed;
      p | iseed2;
      p | iauto_vc_scal_opt;
      //pupping dbles
      p | qseed;
      p | vc_scal_tol;
      p | div_scal;
#ifdef _PARALLEL_DEBUG_        
      if (p.isUnpacking())
        state_class_out ();
#endif           
    } // end pup
#endif

    void state_class_out(){
      char fileName [255];
      sprintf (fileName, "%d_cpvel_samp.state", CkMyPe());
      FILE *fp;  fp = fopen(fileName,"w");

      fprintf(fp,"ivelc_smpl_on %d\n",ivelc_smpl_on);
      fprintf(fp,"ivelc_scal_on %d\n",ivelc_scal_on);
      fprintf(fp,"nvc_smpl %d\n",nvc_smpl);
      fprintf(fp,"nvcnhc_smpl %d\n",nvcnhc_smpl);
      fprintf(fp,"nvc_scal %d\n",nvc_scal);
      fprintf(fp,"iseed %ld\n",iseed);
      fprintf(fp,"iseed2 %ld\n",iseed2);
      fprintf(fp,"iauto_vc_scal_opt %d\n",iauto_vc_scal_opt);
      fprintf(fp,"qseed %g\n",qseed);
      fprintf(fp,"vc_scal_tol %g\n",vc_scal_tol);
      fprintf(fp,"div_scal %g \n",div_scal);
      fclose(fp);
    } // end routine 


}; //CPVEL_SAMP

#ifdef PUP_ON
PUPmarshall(CPVEL_SAMP);
#endif

#endif

//==========================================================================
