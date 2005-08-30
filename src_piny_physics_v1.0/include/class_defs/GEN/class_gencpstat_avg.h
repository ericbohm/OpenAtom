//=======================================================================
//
//                      CP Static averages                          
//             {Variables needed for mem allocation:}
//                                                                       
//=======================================================================


#ifndef _GENCPSTAT_AVG_
#define _GENCPSTAT_AVG_

class GENCPSTAT_AVG {

 //---------------------------------------------------------------------
 public:
  int iter_shake_cp,iter_ratl_cp;     // Num: CP Shake/Rattle iterations

  double cp_kconv0,cp_kconv;          // Num: CP kinetic energy conv
  double kinet_cp,aikinet_cp,akinet_cp;// Num:Inst and avg CP-Class KE 
  double kinet_nhc_cp,aikinet_nhc_cp,akinet_nhc_cp;// Num: CP-NHC KE   
  double vpotnhc_cp;                  // Num:Inst CP-NHC PE             
  double cp_ehart,aicp_ehart,acp_ehart;// Num:Inst, CP Hartree E      
  double cp_eext,aicp_eext,acp_eext;  // Num:Inst  and avg CP Vext      
  double cp_exc,aicp_exc,acp_exc;     // Num:Inst and avg CP exc  E     
  double cp_muxc;                     // Num:Integral of xc pot times rho 
  double cp_eke,aicp_eke,acp_eke;     // Num:Inst and avg CP KE         
  double cp_enl,aicp_enl,acp_enl;     // Num:Inst and avg CP V_nonlocal 
  double aiter_shake_cp,aiter_ratl_cp;// Num: Inst and avg # Shake/Ratl 
                                      //         iterations for CP         
  double maxfc,maxf;                  // Num: MAX component of the force
                                      //      maxfc=coefs,maxf=atms      
  double fc_mag_up,fc_mag_dn;         // Num: avg mag of coef force     
  double fc_max_up,fc_max_dn;         // Num: Max mag of coef force     
  double fatm_max,fatm_mag;           // Num: Max mag of atm  force     
  double count_diag_srot;             // Num: Number of rotations to 
                                      //      diagonal ovlap basis
  double max_off_diag;                // Num: maximum off-diagonal overlap
                                      //         matrix element  
  double max_diag;                    // Num: maximum diagonal overlap
                                      //      matrix element                 
  double econv0,econv;               // Num: Energy conservation       


  GENCPSTAT_AVG(){
   iter_shake_cp = 0;
   iter_ratl_cp  = 0;

   cp_kconv0  = 0;
   cp_kconv   = 0;
   kinet_cp   = 0;
   aikinet_cp = 0;
   akinet_cp  = 0;
   kinet_nhc_cp   = 0;
   aikinet_nhc_cp = 0;
   akinet_nhc_cp  = 0;
   vpotnhc_cp = 0;
   cp_ehart   = 0;
   aicp_ehart = 0;
   acp_ehart  = 0;
   cp_eext    = 0;
   aicp_eext  = 0;
   acp_eext   = 0; 
   cp_exc     = 0;
   aicp_exc   = 0;
   acp_exc    = 0;
   cp_muxc    = 0;
   cp_eke     = 0;
   aicp_eke   = 0;
   acp_eke    = 0;
   cp_enl     = 0;
   aicp_enl   = 0;
   acp_enl    = 0;
   aiter_shake_cp = 0;
   aiter_ratl_cp = 0;
   maxfc     = 0;
   maxf      = 0;
   fc_mag_up = 0;
   fc_mag_dn = 0;
   fc_max_up = 0;
   fc_max_dn = 0;
   fatm_max  = 0;
   fatm_mag  = 0;
   count_diag_srot = 0;
   max_off_diag = 0;
   max_diag = 0;
   econv0   = 0;
   econv    = 0;
  }
  ~GENCPSTAT_AVG(){}

 //---------------------------------------------------------------------
#ifdef PUP_ON
 //pupping
  void pup(PUP::er &p){

    //pupping ints and doubles

        p | econv;      p | econv0;
        p | cp_kconv0;  p | cp_kconv;
        p | kinet_cp;   p | aikinet_cp;   p | akinet_cp;
        p | kinet_nhc_cp; p | aikinet_nhc_cp; p | akinet_nhc_cp;
        p | vpotnhc_cp;
        p | cp_ehart;   p | aicp_ehart;   p | acp_ehart;
        p | cp_eext;    p | aicp_eext;    p | acp_eext;
        p | cp_exc;     p | aicp_exc;     p | acp_exc;
        p | cp_muxc;
        p | cp_eke;     p | aicp_eke;     p | acp_eke;
        p | cp_enl;     p | aicp_enl;     p | acp_enl;
        p | aiter_shake_cp;p | aiter_ratl_cp;
        p | iter_shake_cp;
        p | iter_ratl_cp;
        p | maxfc;      p | maxf;
        p | fc_mag_up;  p | fc_mag_dn;
        p | fc_max_up;  p | fc_max_dn;
        p | fatm_max;   p | fatm_mag;
        p | count_diag_srot;
        p | max_diag;   p | max_off_diag; 
#ifdef _PARALLEL_DEBUG_        
    if (p.isUnpacking())
     state_class_out ();
#endif
  }// This is the tail of the PUP
#endif

//---------------------------------------------------------------------
// Output state of the class

  void state_class_out(){

     char fileName [255];
     sprintf (fileName, "%d_gencpstat_avg.state", CkMyPe());
     FILE *fp;  fp = fopen(fileName,"w");


     fprintf(fp,"cp_kconv0 %g\n",cp_kconv0);
     fprintf(fp,"cp_kconv %g\n",cp_kconv);
     fprintf(fp,"kinet_cp %g\n",kinet_cp);
     fprintf(fp,"aikinet_cp %g\n",aikinet_cp);
     fprintf(fp,"akinet_cp %g\n",akinet_cp);
     fprintf(fp,"kinet_nhc_cp %g\n",kinet_nhc_cp);
     fprintf(fp,"aikinet_nhc_cp %g\n",aikinet_nhc_cp); 
     fprintf(fp,"akinet_nhc_cp %g\n",akinet_nhc_cp);
     fprintf(fp,"vpotnhc_cp %g\n",vpotnhc_cp);
     fprintf(fp,"cp_ehart %g\n",cp_ehart);
     fprintf(fp,"aicp_ehart %g\n",aicp_ehart);
     fprintf(fp,"acp_ehart %g\n",acp_ehart);
     fprintf(fp,"cp_eext %g\n",cp_eext);
     fprintf(fp,"aicp_eext %g\n",aicp_eext);
     fprintf(fp,"acp_eext %g\n",acp_eext);
     fprintf(fp,"cp_exc %g\n",cp_exc);
     fprintf(fp,"aicp_exc %g\n",aicp_exc);
     fprintf(fp,"acp_exc %g\n",acp_exc);
     fprintf(fp,"cp_muxc %g\n",cp_muxc);
     fprintf(fp,"cp_eke %g\n",cp_eke);
     fprintf(fp,"aicp_eke %g\n",aicp_eke);
     fprintf(fp,"acp_eke %g\n",acp_eke);
     fprintf(fp,"cp_enl %g\n",cp_enl);
     fprintf(fp,"aicp_enl %g\n",aicp_enl);
     fprintf(fp,"acp_enl %g\n",acp_enl);
     fprintf(fp,"aiter_shake_cp %g\n",aiter_shake_cp);
     fprintf(fp,"aiter_ratl_cp %g\n",aiter_ratl_cp);
     fprintf(fp,"maxfc %g\n",maxfc);
     fprintf(fp,"maxf %g\n",maxf);
     fprintf(fp,"fc_mag_up %g\n",fc_mag_up);
     fprintf(fp,"fc_mag_dn %g\n",fc_mag_dn);
     fprintf(fp,"fc_max_up %g\n",fc_max_up);
     fprintf(fp,"fc_max_dn %g\n",fc_max_dn);
     fprintf(fp,"fatm_max %g\n",fatm_max);
     fprintf(fp,"fatm_mag %g\n",fatm_mag);
     fprintf(fp,"count_diag_srot %g\n",count_diag_srot);
     fprintf(fp,"max_diag %g\n",max_diag);
     fprintf(fp,"max_off_diag %g\n",max_off_diag);

  }// end output state of class

//---------------------------------------------------------------------
   }; // GENCPSTAT_AVG;
//=======================================================================

#ifdef PUP_ON
PUPmarshall(GENCPSTAT_AVG);
#endif

#endif
//=======================================================================
