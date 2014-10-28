//=======================================================================
//               Static averages                                         
//             {Variables needed for mem allocation:                     
//                                                                       

#ifndef _GENSTAT_AVG_
#define _GENSTAT_AVG_

class GENSTAT_AVG {

  //----------------
  public:
    int nexclt;                        // Num: Tot # of exclusions       
    int iter_shake, iter_ratl;         // Num: Shake/Rattle iterations
    int iter_shake_cp,iter_ratl_cp;    // Num: CP Shake/Rattle iterations
    int iswit_vdw;                     // Opt: Switch inter shift on/off 
    int itime_update,itime_update_w;   // Num: TIme of updates           
    int write_cp_atm_flag;             // Opt : write cp atm flag

    double vintert,aivintert,avintert; // Num:Inst, and avgs inter PE    
    double vintrat,aivintrat,avintrat; // Num:Inst, and avgs intra PE    
    double vsurft;                     // Num: surface energy pot        
    double vbondt,vbendt,vtorst,vonfot;// Num: Bond,bend,tors, and onfo  
    double vbend_bndt;                 // NUM: Uri Bradleys              
    double vbend_bnd_bond,vbend_bnd_bend;// NUM: Uri Bradleys decomp   
    double vbondt_watts,vbendt_watts,vtot_watts;// NUM: Watts dcomp  
    double vreal,vrecip;               // Num: Inter mol PE              
    double vvdw,vcoul;                 // Num: Van der Waals and Coulomb 
    double vlong;                      // Num: Long range correction to LJ
    double vbond_free,vbend_free,vtors_free;// Num: Free intra pots     
    double vbar_free;                  // Num: more free energy pots     
    double kinet,aikinet,akinet;       // Num:Inst, and avg atm KE       
    double kinet_v,aikinet_v,akinet_v; // Num:Inst, and avg volume KE    
    double vol,aivol,avol;             // Num:Inst, and avg volume       
    double kinet_nhc,aikinet_nhc,akinet_nhc;//Num:Inst and avg Atm NHC KE
    double kinet_nhc_bead;             //Num:Inst and avg bead NHC
    double aikinet_nhc_bead;
    double akinet_nhc_bead;            
    double vpot_v;                     // Num: Volume PE                 
    double vpotnhc;                    // Num: NHC PE                    
    double aiter_shake,aiter_ratl;     // Num:Inst and avg # Shake/Ratl   
    double aiter_23,aiter_33, aiter_46;// Num:Inst and avg # grp Shake/Ratl   
    double iter_23, iter_33, iter_46, iter_43, iter_21; 
    double aiter_21,aiter_43;          // Num:Inst and avg # grp Shake/Ratl   
    double iter_23r, iter_33r, iter_46r, iter_43r, iter_21r; 
    double aiter_23r,aiter_33r, aiter_46r;// Num:Inst and avg # grp Shake/Ratl
    double aiter_21r,aiter_43r;        // Num:Inst and avg # grp Shake/Ratl
    double acella,acellb,acellc;
    double aicella,aicellb,aicellc;    // Num:Inst and avg cell lngth    
    double acellab,acellbc,acellac;
    double aicellab,aicellbc,aicellac; // Num:Inst and avg cell angles   
    double apress,aipress;             // Num: Avg, inst avg, pressure   
    double press_inter,press_intra;    // Num: Inter, intra pressure     
    double apress_inter,aipress_inter; // Num: Avg, inst avg, inter pressure   
    double apress_intra,aipress_intra; // Num: Avg, inst avg, intra pressure   
    double press_kin;
    double apress_kin,aipress_kin;     // Num: Avg, inst avg, intra pressure   
    double econv0,econv;               // Num: Energy conservation       
    double cp_kconv0,cp_kconv;         // Num: CP kinetic energy conv
    double cpu1,cpu2,acpu,cpu_now;     // Num: Cpu time                  
    double updates,updates_w;          // Num: Ver-list,NL-list updates  
    double kinet_cp,aikinet_cp,akinet_cp;// Num:Inst and avg CP-Class KE 
    double kinet_nhc_cp,aikinet_nhc_cp,akinet_nhc_cp;// Num: CP-NHC KE   
    double vpotnhc_cp;                 // Num:Inst CP-NHC PE             
    double cp_ehart,aicp_ehart,acp_ehart;// Num:Inst, CP Hartree E      
    double cp_eext,aicp_eext,acp_eext; // Num:Inst  and avg CP Vext      
    double cp_exc,aicp_exc,acp_exc;    // Num:Inst and avg CP exc  E     
    double cp_muxc;                    // Num:Integral of xc pot times rho 
    double cp_eke,aicp_eke,acp_eke;    // Num:Inst and avg CP KE         
    double cp_enl,aicp_enl,acp_enl;    // Num:Inst and avg CP V_nonlocal 
    double aiter_shake_cp,aiter_ratl_cp;// Num: Inst and avg # Shake/Ratl 
    //         iterations for CP         
    double maxfc,maxf;                 // Num: MAX component of the force
    //      maxfc=coefs,maxf=atms      
    double pi_ke_prim,pi_ke_vir;       // Num: Quantum KE estimators     
    double api_ke_prim,api_ke_vir;     // Num: Average quantum KE est.   
    double aipi_ke_prim,aipi_ke_vir;   // Num: Inst. Avg. quantum KE est.
    double kin_harm,akin_harm,aikin_harm;// Num: Harmonic KE             

    double fc_mag_up,fc_mag_dn;        // Num: avg mag of coef force     
    double fc_max_up,fc_max_dn;        // Num: Max mag of coef force     
    double fatm_max,fatm_mag;          // Num: Max mag of atm  force     
    double count_diag_srot;            // Num: Number of rotations to 
    //      diagonal ovlap basis      
    double econv_now;                  // Num: The conserved quantity    

    double *apten, *aipten, *apten_out;// Lst: Avg pressure tensor; Lth:9

    //----------------
    //con-destruct:
    GENSTAT_AVG(){
      nexclt            = 0;
      iter_shake        = 0;
      iter_ratl         = 0;        
      iter_shake_cp     = 0;
      iter_ratl_cp      = 0;   
      iswit_vdw         = 0;                    
      itime_update      = 0;
      itime_update_w    = 0;  
      write_cp_atm_flag = 0;            
      apten      = (double *)cmalloc(9*sizeof(double),"class_genpstat_cons")-1;
      aipten     = (double *)cmalloc(9*sizeof(double),"class_genptens_cons")-1;
      apten_out  = (double *)cmalloc(9*sizeof(double),"class_genptens_cons")-1;
    };
    ~GENSTAT_AVG(){};

#ifdef PUP_ON
    //----------------
    //pupping
    void pup(PUP::er &p){
      //pupping ints
      p | nexclt;
      p | iter_shake;
      p | iter_ratl; 
      p | iswit_vdw;
      p | itime_update;
      p | itime_update_w;
      p | iter_shake_cp;
      p | iter_ratl_cp;
      p | write_cp_atm_flag;
      //pupping dbles
      p | vintert;    p | aivintert;    p | avintert;
      p | vintrat;    p | aivintrat;    p | avintrat;
      p | vsurft;
      p | vbondt;
      p | vbendt;
      p | vtorst;
      p | vonfot;
      p | vbend_bndt;
      p | vbend_bnd_bond;
      p | vbend_bnd_bend;
      p | vbondt_watts;
      p | vbendt_watts;
      p | vtot_watts;
      p | vreal;
      p | vrecip;
      p | vvdw;
      p | vcoul;
      p | vlong;
      p | vbond_free;
      p | vbend_free;
      p | vtors_free;
      p | vbar_free;
      p | kinet;      p | aikinet;      p | akinet;
      p | kinet_v;    p | aikinet_v;    p | akinet_v;
      p | vol;        p | aivol;        p | avol;
      p | kinet_nhc;  p | aikinet_nhc;  p | akinet_nhc;
      p | kinet_nhc_bead; p | aikinet_nhc_bead;  p | akinet_nhc_bead;
      p | vpot_v;
      p | vpotnhc;
      p | aiter_shake;p | aiter_ratl;
      p | aiter_23;   p | aiter_33;     p | aiter_46;   p | aiter_43;
      p | iter_23;    p | iter_33;      p | iter_46;
      p | iter_43;    p | iter_21;      p | aiter_21;
      p | iter_23r;   p | iter_33r;     p | iter_46r;
      p | iter_43r;   p | iter_21r;     p | aiter_23r;
      p | aiter_33r;  p | aiter_46r;    p | aiter_21r;   p | aiter_43r;
      p | acella;     p | acellb;       p | acellc;
      p | aicella;    p | aicellb;      p | aicellc;
      p | acellab;    p | acellbc;      p | acellac;
      p | aicellab;   p | aicellbc;     p | aicellac;
      p | apress;     p | aipress;
      p | press_inter;p | apress_inter; p | aipress_inter;
      p | press_intra;p | apress_intra; p | aipress_intra;
      p | press_kin;  p | apress_kin;   p | aipress_kin;
      p | econv0;     p | econv;
      p | cp_kconv0;  p | cp_kconv;
      p | cpu1;       p | cpu2;         p | acpu;
      p | cpu_now;
      p | updates;    p | updates_w;
      p | kinet_cp;   p | aikinet_cp;   p | akinet_cp;
      p | kinet_nhc_cp;p | aikinet_nhc_cp; p | akinet_nhc_cp;
      p | vpotnhc_cp;
      p | cp_ehart;   p | aicp_ehart;   p | acp_ehart;
      p | cp_eext;    p | aicp_eext;    p | acp_eext;
      p | cp_exc;     p | aicp_exc;     p | acp_exc;
      p | cp_muxc;
      p | cp_eke;     p | aicp_eke;     p | acp_eke;
      p | cp_enl;     p | aicp_enl;     p | acp_enl;
      p | aiter_shake_cp;p | aiter_ratl_cp;
      p | maxfc;      p | maxf;
      p | pi_ke_prim; p | api_ke_prim;  p | aipi_ke_prim;
      p | pi_ke_vir;  p | api_ke_vir;   p | aipi_ke_vir;
      p | kin_harm;   p | akin_harm;    p | aikin_harm;
      p | fc_mag_up;  p | fc_mag_dn;
      p | fc_max_up;  p | fc_max_dn;
      p | fatm_max;   p | fatm_mag;
      p | count_diag_srot;
      //pupping dbl  arrays
      pup1d_dbl(p,&apten,9);
      pup1d_dbl(p,&aipten,9);
      pup1d_dbl(p,&apten_out,9);
#ifdef _PARALLEL_DEBUG_        
      if (p.isUnpacking())
        state_class_out ();
#endif        
    } // end pup
#endif

    void state_class_out(){

      char fileName [255];
      sprintf (fileName, "%d_genstat_avg.state", CkMyPe());
      FILE *fp;  fp = fopen(fileName,"w");

      fprintf(fp,"nexclt %d\n",nexclt);
      fprintf(fp,"iter_shake %d\n",iter_shake);
      fprintf(fp,"iter_ratl %d\n",iter_ratl); 
      fprintf(fp,"iswit_vdw %d\n",iswit_vdw);
      fprintf(fp,"itime_update %d\n",itime_update);
      fprintf(fp,"itime_update_w %d\n",itime_update_w);
      fprintf(fp,"iter_shake_cp %d\n",iter_shake_cp);
      fprintf(fp,"iter_ratl_cp %d\n",iter_ratl_cp);
      fprintf(fp,"write_cp_atm_flag %d\n",write_cp_atm_flag);
      // dbles
      fprintf(fp,"vintert %g\n",vintert);
      fprintf(fp,"aivintert %g\n",aivintert);
      fprintf(fp,"avintert %g\n",avintert);
      fprintf(fp,"vintrat %g\n",vintrat);
      fprintf(fp,"aivintrat %g\n",aivintrat);
      fprintf(fp,"avintrat %g\n",avintrat);
      fprintf(fp,"vsurft %g\n",vsurft);
      fprintf(fp,"vbondt %g\n",vbondt);
      fprintf(fp,"vbendt %g\n",vbendt);
      fprintf(fp,"vtorst %g\n",vtorst);
      fprintf(fp,"vonfot %g\n",vonfot);
      fprintf(fp,"vbend_bndt %g\n",vbend_bndt);
      fprintf(fp,"vbend_bnd_bond %g\n",vbend_bnd_bond);
      fprintf(fp,"vbend_bnd_bend %g\n",vbend_bnd_bend);
      fprintf(fp,"vbondt_watts %g\n",vbondt_watts);
      fprintf(fp,"vbendt_watts %g\n",vbendt_watts);
      fprintf(fp,"vtot_watts %g\n",vtot_watts);
      fprintf(fp,"vreal %g\n",vreal);
      fprintf(fp,"vrecip %g\n",vrecip);
      fprintf(fp,"vvdw %g\n",vvdw);
      fprintf(fp,"vcoul %g\n",vcoul);
      fprintf(fp,"vlong %g\n",vlong);
      fprintf(fp,"vbond_free %g\n",vbond_free);
      fprintf(fp,"vbend_free %g\n",vbend_free);
      fprintf(fp,"vtors_free %g\n",vtors_free);
      fprintf(fp,"vbar_free %g\n",vbar_free);
      fprintf(fp,"kinet %g\n",kinet);
      fprintf(fp,"aikinet %g\n",aikinet);
      fprintf(fp,"akinet %g\n",akinet);
      fprintf(fp,"kinet_v %g\n",kinet_v);
      fprintf(fp,"aikinet_v %g\n",aikinet_v);
      fprintf(fp,"akinet_v %g\n",akinet_v);
      fprintf(fp,"vol %g\n",vol);
      fprintf(fp,"aivol %g\n",aivol);
      fprintf(fp,"avol %g\n",avol);
      fprintf(fp,"kinet_nhc %g\n",kinet_nhc);
      fprintf(fp,"aikinet_nhc %g\n",aikinet_nhc);
      fprintf(fp,"akinet_nhc %g\n",akinet_nhc);
      fprintf(fp,"kinet_nhc_bead %g\n",kinet_nhc_bead);
      fprintf(fp,"aikinet_nhc_bead %g\n",aikinet_nhc_bead);
      fprintf(fp,"akinet_nhc_bead %g\n",akinet_nhc_bead);
      fprintf(fp,"vpot_v %g\n",vpot_v);
      fprintf(fp,"vpotnhc %g\n",vpotnhc);
      fprintf(fp,"aiter_shake %g\n",aiter_shake);
      fprintf(fp,"aiter_ratl %g\n",aiter_ratl);
      fprintf(fp,"aiter_23 %g\n",aiter_23);
      fprintf(fp,"aiter_33 %g\n",aiter_33);
      fprintf(fp,"aiter_46 %g\n",aiter_46);
      fprintf(fp,"aiter_43 %g\n",aiter_43);
      fprintf(fp,"iter_23 %g\n",iter_23);
      fprintf(fp,"iter_33 %g\n",iter_33);
      fprintf(fp,"iter_46 %g\n",iter_46);
      fprintf(fp,"iter_43 %g\n",iter_43);
      fprintf(fp,"iter_21 %g\n",iter_21);
      fprintf(fp,"aiter_21 %g\n",aiter_21);
      fprintf(fp,"iter_23r %g\n",iter_23r);
      fprintf(fp,"iter_33r %g\n",iter_33r);
      fprintf(fp,"iter_46r %g\n",iter_46r);
      fprintf(fp,"iter_43r %g\n",iter_43r);
      fprintf(fp,"iter_21r %g\n",iter_21r);
      fprintf(fp,"aiter_23r %g\n",aiter_23r);
      fprintf(fp,"aiter_33r %g\n",aiter_33r);
      fprintf(fp,"aiter_46r %g\n",aiter_46r);
      fprintf(fp,"aiter_21r %g\n",aiter_21r);
      fprintf(fp,"aiter_43r %g\n",aiter_43r);
      fprintf(fp,"acella %g\n",acella);
      fprintf(fp,"acellb %g\n",acellb);
      fprintf(fp,"acellc %g\n",acellc);
      fprintf(fp,"aicella %g\n",aicella);
      fprintf(fp,"aicellb %g\n",aicellb);
      fprintf(fp,"aicellc %g\n",aicellc);
      fprintf(fp,"acellab %g\n",acellab);
      fprintf(fp,"acellbc %g\n",acellbc);
      fprintf(fp,"acellac %g\n",acellac);
      fprintf(fp,"aicellab %g\n",aicellab);
      fprintf(fp,"aicellbc %g\n",aicellbc);
      fprintf(fp,"aicellac %g\n",aicellac);
      fprintf(fp,"apress %g\n",apress);
      fprintf(fp,"aipress %g\n",aipress);
      fprintf(fp,"press_inter %g\n",press_inter);
      fprintf(fp,"apress_inter %g\n",apress_inter);
      fprintf(fp,"aipress_inter %g\n",aipress_inter);
      fprintf(fp,"press_intra %g\n",press_intra);
      fprintf(fp,"apress_intra %g\n",apress_intra);
      fprintf(fp,"aipress_intra %g\n",aipress_intra);
      fprintf(fp,"press_kin %g\n",press_kin);
      fprintf(fp,"apress_kin %g\n",apress_kin);
      fprintf(fp,"aipress_kin %g\n",aipress_kin);
      fprintf(fp,"econv0 %g\n",econv0);
      fprintf(fp,"econv %g\n",econv);
      fprintf(fp,"cp_kconv0 %g\n",cp_kconv0);
      fprintf(fp,"cp_kconv %g\n",cp_kconv);
      fprintf(fp,"cpu1 %g\n",cpu1);
      fprintf(fp,"cpu2 %g\n",cpu2);
      fprintf(fp,"acpu %g\n",acpu);
      fprintf(fp,"cpu_now %g\n",cpu_now);
      fprintf(fp,"updates %g\n",updates);
      fprintf(fp,"updates_w %g\n",updates_w);
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
      fprintf(fp,"pi_ke_prim %g\n",pi_ke_prim);
      fprintf(fp,"api_ke_prim %g\n",api_ke_prim);
      fprintf(fp,"aipi_ke_prim %g\n",aipi_ke_prim);
      fprintf(fp,"pi_ke_vir %g\n",pi_ke_vir);
      fprintf(fp,"api_ke_vir %g\n",api_ke_vir);
      fprintf(fp,"aipi_ke_vir %g\n",aipi_ke_vir);
      fprintf(fp,"kin_harm %g\n",kin_harm);
      fprintf(fp,"akin_harm %g\n",akin_harm);
      fprintf(fp,"aikin_harm %g\n",aikin_harm);
      fprintf(fp,"fc_mag_up %g\n",fc_mag_up);
      fprintf(fp,"fc_mag_dn %g\n",fc_mag_dn);
      fprintf(fp,"fc_max_up %g\n",fc_max_up);
      fprintf(fp,"fc_max_dn %g\n",fc_max_dn);
      fprintf(fp,"fatm_max %g\n",fatm_max);
      fprintf(fp,"fatm_mag %g\n",fatm_mag);
      fprintf(fp,"count_diag_srot %g\n",count_diag_srot);
      fclose(fp);
    }// end routine


}; // GENSTAT_AVG;

#ifdef PUP_ON
PUPmarshall(GENSTAT_AVG);
#endif

#endif

//==========================================================================

