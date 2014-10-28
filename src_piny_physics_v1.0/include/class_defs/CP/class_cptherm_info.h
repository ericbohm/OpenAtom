//==========================================================================
//                 CP coeff Nose'-Hoover chains                             
//             {Variables needed for mem allocation:                        
//                   num_c_nhc,len_c_nhc }                                  
//                                                                          

#ifndef _CPTHERM_INFO_
#define _CPTHERM_INFO_

class CPTHERM_INFO{

  //----------------
  public:
    int cp_any_on;              // Opt: Is cp of any type ``on''
    int nstate;
    int nck_c_nhc;              // Num: # of PW coeff NHC's per plane
    int num_c_nhc;              // Num: # of PW coeff NHC's            
    int num_c_nhc_iso;          // Num: # of PW coeff NHC's            
    int len_c_nhc;              // Num: Lnth of PW coeff NHC's         
    int nres_c_nhc;             // Num: # of PW coeff RESPA NHC steps   
    int nyosh_c_nhc;            // Num: # of PW coeff Yosh NHC's steps  

    //----------------------------------------------------------------------------
    //con-destruct:
    CPTHERM_INFO(){
      cp_any_on      = 0;
      nstate         = 0;        
      nck_c_nhc      = 0;     
      num_c_nhc      = 0;     
      num_c_nhc_iso  = 0;     
      len_c_nhc      = 0;     
      nres_c_nhc     = 0;    
      nyosh_c_nhc    = 0;   
    }
    ~CPTHERM_INFO(){};

    //----------------------------------------------------------------------------
#ifdef PUP_ON
    //pupping
    void pup(PUP::er &p){
      //pupping ints
      p | cp_any_on;      
      p | nstate;      
      p | nck_c_nhc;
      p | num_c_nhc;
      p | num_c_nhc_iso;
      p | len_c_nhc;
      p | nres_c_nhc;
      p | nyosh_c_nhc;
#ifdef _PARALLEL_DEBUG_        
      if (p.isUnpacking())
        state_class_out ();
#endif           
    } // end pup
#endif

    //----------------------------------------------------------------------------
    void state_class_out(){

      int i,j;
      char fileName [255];
      sprintf (fileName, "%d_cptherm_info.state", CkMyPe());
      FILE *fp;  fp = fopen(fileName,"w");

      // ints
      fprintf(fp,"nck_c_nhc %d\n",nck_c_nhc);
      fprintf(fp,"num_c_nhc %d\n",num_c_nhc);
      fprintf(fp,"num_c_nhc_iso %d\n",num_c_nhc_iso);
      fprintf(fp,"len_c_nhc %d\n",len_c_nhc);
      fprintf(fp,"nres_c_nhc %d\n",nres_c_nhc);
      fprintf(fp,"nyosh_c_nhc %d\n",nyosh_c_nhc);
      fclose(fp);
    } // end routine 

    //----------------------------------------------------------------------------
}; //CPTHERM_INFO
//==========================================================================

#ifdef PUP_ON
PUPmarshall(CPTHERM_INFO);
#endif


#endif
//==========================================================================
