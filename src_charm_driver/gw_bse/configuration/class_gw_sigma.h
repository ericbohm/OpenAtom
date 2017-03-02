//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//
//                  GWBSE simulation options                                   
//
//                class definition for GW_SIGMA
//                                                                          
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
#ifndef _GW_SIGMA_
#define _GW_SIGMA_


class GW_SIGMA{

  public:
    int PP_nmode;                // Num: number of eigenmodes included in plasmon-pole model

    char nnpFileName[200];       // No need to pup this, contents stored in next three variables
    int num_sig_matels;  // number of <n|Sigma|n'> elements to calculate in GW
    int *n_list_sig_matels, *np_list_sig_matels;  // n an n' lists

    double screened_coulomb_cutoff;     // Num: Cutoff for screened coulomb term (must be less than Ecuteps (eV)
    double bare_coulomb_cutoff;         // Num: cutoff for bare coulomb (eV)

    //----------------
    //con-destruct:
    GW_SIGMA(){
      PP_nmode        = 0;
      screened_coulomb_cutoff = 0;
      bare_coulomb_cutoff     = 0;
      num_sig_matels = 0;
      n_list_sig_matels = 0;
      np_list_sig_matels = 0;
    };
    ~GW_SIGMA(){};


#ifdef PUP_ON
    //-------------
    //pupping
    void pup(PUP::er &p){
      //pupping ints
      p | PP_nmode;
      p | num_sig_matels;
      //puppting dbles
      p | screened_coulomb_cutoff;
      p | bare_coulomb_cutoff;      
      if(p.isUnpacking()) {
	        n_list_sig_matels = new int [num_sig_matels];
	        np_list_sig_matels = new int [num_sig_matels];
      }
      // pup 1d integer arrays
      PUParray(p,n_list_sig_matels,num_sig_matels);
      PUParray(p,np_list_sig_matels,num_sig_matels);
#ifdef _PARALLEL_DEBUG_        
      if (p.isUnpacking())
        state_class_out ();
#endif
    } // end pup
#endif

    void state_class_out(){
      char fileName [255];
      sprintf ( fileName, "%d_gw_sigma.out", CKMYPE());
      FILE *fp; fp = fopen(fileName,"w");
      //ints
      fprintf(fp,"PP_nmode %d\n",PP_nmode);
      fprintf(fp,"n and n' list for sigma matrix elements is %d long\n",num_sig_matels);
      for (int i=0; i < num_sig_matels; i++) {
        	fprintf(fp,"%d %d\n",n_list_sig_matels[i],np_list_sig_matels[i]);
      }
      //dbles
      fprintf(fp,"screened_coulomb_cutoff %lg\n",screened_coulomb_cutoff);
      fprintf(fp,"bare_coulomb_cutoff %lg\n",bare_coulomb_cutoff);
      fclose(fp);
    }// end routine
    
}; // GW_SIGMA

#ifdef PUP_ON
PUPmarshall(GW_SIGMA);
#endif

#endif
//==========================================================================
