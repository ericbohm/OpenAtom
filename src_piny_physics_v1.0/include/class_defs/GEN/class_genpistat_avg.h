//=======================================================================
//
//                      PIMD Static averages                          
//             {Variables needed for mem allocation:}
//                                                                       
//=======================================================================


#ifndef _GENPISTAT_AVG_
#define _GENPISTAT_AVG_

class GENPISTAT_AVG {

 //---------------------------------------------------------------------
 public:

  double kinet_nhc_bead;             //Num:Inst and avg bead NHC
  double aikinet_nhc_bead;
  double akinet_nhc_bead;            
  double pi_ke_prim,pi_ke_vir;       // Num: Quantum KE estimators     
  double api_ke_prim,api_ke_vir;     // Num: Average quantum KE est.   
  double aipi_ke_prim,aipi_ke_vir;   // Num: Inst. Avg. quantum KE est.
  double kin_harm,akin_harm,aikin_harm;// Num: Harmonic KE             

  GENPISTAT_AVG(){
    kinet_nhc_bead   = 0.0;
    aikinet_nhc_bead = 0.0;
    akinet_nhc_bead  = 0.0;
    pi_ke_prim       = 0.0;
    pi_ke_vir        = 0.0;
    api_ke_prim      = 0.0;
    api_ke_vir       = 0.0;
    aipi_ke_prim     = 0.0;
    aipi_ke_vir      = 0.0;
    kin_harm         = 0.0;
    akin_harm        = 0.0;
    aikin_harm       = 0.0;
  }
  ~GENPISTAT_AVG(){}


 //---------------------------------------------------------------------
#ifdef PUP_ON
 //pupping
  void pup(PUP::er &p){

        p | pi_ke_prim; p | api_ke_prim;  p | aipi_ke_prim;
        p | pi_ke_vir;  p | api_ke_vir;   p | aipi_ke_vir;
        p | kin_harm;   p | akin_harm;    p | aikin_harm;
        p | kinet_nhc_bead; p | aikinet_nhc_bead;  p | akinet_nhc_bead;
#ifdef _PARALLEL_DEBUG_        
    if (p.isUnpacking())
     state_class_out ();
#endif
  }// End of PUP
#endif

//---------------------------------------------------------------------
// Output state of the class

  void state_class_out(){


     char fileName [255];
     sprintf (fileName, "%d_genpisstat_avg.state", CkMyPe());
     FILE *fp;  fp = fopen(fileName,"w");


     fprintf(fp,"pi_ke_prim %g\n",pi_ke_prim);
     fprintf(fp,"api_ke_prim %g\n",api_ke_prim);
     fprintf(fp,"aipi_ke_prim %g\n",aipi_ke_prim);
     fprintf(fp,"pi_ke_vir %g\n",pi_ke_vir);
     fprintf(fp,"api_ke_vir %g\n",api_ke_vir);
     fprintf(fp,"aipi_ke_vir %g\n",aipi_ke_vir);
     fprintf(fp,"kin_harm %g\n",kin_harm);
     fprintf(fp,"akin_harm %g\n",akin_harm);
     fprintf(fp,"aikin_harm %g\n",aikin_harm);
     fprintf(fp,"kinet_nhc_bead %g\n",kinet_nhc_bead);
     fprintf(fp,"aikinet_nhc_bead %g\n",aikinet_nhc_bead);
     fprintf(fp,"akinet_nhc_bead %g\n",akinet_nhc_bead);
     fclose(fp);
  }// end output state of class

}; // GENPISTAT_AVG;
//---------------------------------------------------------------------

#ifdef PUP_ON
PUPmarshall(GENPISTAT_AVG);
#endif

#endif

//==========================================================================
