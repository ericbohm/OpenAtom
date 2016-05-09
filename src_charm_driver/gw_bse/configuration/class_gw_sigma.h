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
    int band_index_min;          // Band index for sigma
    int band_index_max;          // Band index for sigma

    double screened_coulomb_cutoff;     // Num: Cutoff for screened coulomb term (must be less than Ecuteps (eV)
    double bare_coulomb_cutoff;         // Num: cutoff for bare coulomb (eV)

    //----------------
    //con-destruct:
    GW_SIGMA(){
      PP_nmode        = 0;
      band_index_min  = 0;
      band_index_max  = 0;
      screened_coulomb_cutoff = 0;
      bare_coulomb_cutoff     = 0;
    };
    ~GW_SIGMA(){};


#ifdef PUP_ON
    //-------------
    //pupping
    void pup(PUP::er &p){
      //pupping ints
      p | PP_nmode;
      p | band_index_min;
      p | band_index_max;
      //puppting dbles
      p | screened_coulomb_cutoff;
      p | bare_coulomb_cutoff;      
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
      fprintf(fp,"band_index_min %d\n",band_index_min);
      fprintf(fp,"band_index_max %d\n",band_index_max);
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
