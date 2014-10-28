//==========================================================================
//                  Minimization options                                    
//             {Variables needed for mem allocation:}                       

#ifndef _GENMINOPTS_
#define _GENMINOPTS_

class GENMINOPTS {

  //----------------
  public:
    int min_std;                // Opt: Steepest descent min           
    int min_cg;                 // Opt: Conjugate gradiant min         
    int min_diis;               // Opt: DIIS min                       
    int diis_hist_len;          // Num: Length of DIIS history         
    int cp_min_std;             // Opt: Steepest descent min           
    int cp_min_cg;              // Opt: Conjugate gradiant min         
    int cp_min_diis;            // Opt: DIIS min                       
    int cp_diis_hist_len;       // Num: Length of DIIS history         
    int cp_cg_line_min_len;     // Num: COMMENT ME 
    int min_atm_com_fix_opt;    // Opt: Keep the com fixed             

    double tol_coef;            // Num: Tol on PW(plane wave) coef forces   
    double tol_atom;            // Num: Tol on atm forces              

    //----------------
    //con-destruct:
    GENMINOPTS(){
      min_std       = 0;               
      min_cg        = 0;                
      min_diis      = 0;              
      diis_hist_len = 0;         
      cp_min_std    = 0;            
      cp_min_cg     = 0;             
      cp_min_diis   = 0;           
      cp_diis_hist_len    = 0;      
      cp_cg_line_min_len  = 0;    
      min_atm_com_fix_opt = 0;   
    };
    ~GENMINOPTS(){};

#ifdef PUP_ON
    //----------------
    //pupping
    void pup(PUP::er &p){
      //pupping ints
      p | min_std;
      p | min_cg;
      p | min_diis;
      p | diis_hist_len;
      p | cp_min_std;
      p | cp_min_cg;
      p | cp_min_diis;
      p | cp_diis_hist_len;
      p | cp_cg_line_min_len;
      p | min_atm_com_fix_opt;
      //pupping dbles
      p | tol_coef;
      p | tol_atom;
#ifdef _PARALLEL_DEBUG_        
      if (p.isUnpacking())
        state_class_out ();
#endif        
    } // end pup
#endif

    void state_class_out(){

      char fileName [255];
      sprintf (fileName, "%d_genminopts.state", CkMyPe());
      FILE *fp;  fp = fopen(fileName,"w");

      //int
      fprintf(fp,"min_std %d\n",min_std);
      fprintf(fp,"min_cg %d\n",min_cg);
      fprintf(fp,"min_diis %d\n",min_diis);
      fprintf(fp,"diis_hist_len %d\n",diis_hist_len);
      fprintf(fp,"cp_min_std %d\n",cp_min_std);
      fprintf(fp,"cp_min_cg %d\n",cp_min_cg);
      fprintf(fp,"cp_min_diis %d\n",cp_min_diis);
      fprintf(fp,"cp_diis_hist_len %d\n",cp_diis_hist_len);
      fprintf(fp,"cp_cg_line_min_len %d\n",cp_cg_line_min_len);
      fprintf(fp,"min_atm_com_fix_opt %d\n",min_atm_com_fix_opt);
      //dbles
      fprintf(fp,"tol_coef %g\n",tol_coef);
      fprintf(fp,"tol_atom %g\n",tol_atom);
      fclose(fp);
    } // end routine


}; // GENMINOPTS;

#ifdef PUP_ON
PUPmarshall(GENMINOPTS);
#endif

#endif

//==========================================================================
