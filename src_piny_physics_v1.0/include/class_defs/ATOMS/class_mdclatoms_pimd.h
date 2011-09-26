//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                                                                          
//                      class_mdclatoms_pimd.h                              
//                                                                          
//         Class definition for classical atom 
//         currently used for path integrals                                
//                                                                          
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

#ifndef _MDCLATOMS_PIMD_
#define _MDCLATOMS_PIMD_

class MDCLATOMS_PIMD{
 public:
  int pimd_on;                  // Opt: Path integral are on
  int pi_beads;                 // Num: Number of path integral beads
  int natm_tot;                 // Num: Total number of atoms     
  int nfree_pimd;               // Num: # degrees of freedom     
  int pimd_freez_typ;           // Opt: type of freezing 
  int pi_beads_full_ter;        // Num: full beads per break-up          
  int pi_beads_res_ter;         // Num: res_ter beads per break-up       
  int pi_beads_full_tra;        // Num: full_tra beads per res_ter break-up
  int pi_beads_res_tra;         // Num: res_tra beads per full_tra break-up

  double gamma_adb;             // Num: Adiabaticity parameter  
  double pi_temperature;        // Num: Path integral temperature 
                                //       (used in harmonic bead interaction)
  double wght_pimd;             // Num: bead harmonic RESPA wght 
  double pi_beads_full_ter_wght;// Num: full beads wght          
  double pi_beads_res_ter_wght; // Num: res_ter beads wght       
  double pi_beads_full_tra_wght;// Num: full_tra beads wght      
  double pi_beads_res_tra_wght; // Num: res_tra beads wght       
  double rcut_spread;

  int *ip_lab;                  // Lst: bead types          Lth: pi_beads 
  double *rat1_stag,*rat2_stag; // Lst: bead convert factor Lth: pi_beads
  double *path_eig;             // Lst: List of bead eigenvalues Lth: pi_beads
  double *prekf;                // Lst : natm_tot 

//---------------------------------------------------------------------------
// Default constructor/destructor

   MDCLATOMS_PIMD(){
     pimd_on = 0;
     pi_beads = 0;
     natm_tot = 0;
     nfree_pimd = 0;
     pimd_freez_typ = 0;
     pi_beads_full_ter = 0;
     pi_beads_res_ter = 0;
     pi_beads_full_tra = 0;
     pi_beads_res_tra = 0;
   }
  ~MDCLATOMS_PIMD(){}

//---------------------------------------------------------------------------
// Pack/Unpack

#ifdef PUP_ON
  void pup(PUP::er &p){
    // PUP the ints
    p | pimd_on;
    p | pi_beads;
    p | natm_tot;
    p | nfree_pimd;
    p | pimd_freez_typ;
    p | pi_beads_full_ter;
    p | pi_beads_res_ter;
    p | pi_beads_full_tra;
    p | pi_beads_res_tra;

    // PUP the doubles
    p | gamma_adb;
    p | pi_temperature;
    p | wght_pimd;
    p | pi_beads_full_ter_wght;
    p | pi_beads_res_ter_wght;
    p | pi_beads_full_tra_wght;
    p | pi_beads_res_tra_wght;
    p | rcut_spread;

    // PUP the arrays (For now, but later needs the path integral on option)
#ifdef _PARALLEL_DEBUG_        
    if (p.isUnpacking())
     state_class_out ();
#endif         
  }// end pack/unpack 
#endif

//---------------------------------------------------------------------------
// Print out the state of the class
  void state_class_out(){
     char fileName [255];
     sprintf (fileName, "%d_mdclatom_pimd.state", CkMyPe());
     FILE *fp;  fp = fopen(fileName,"w");

     fprintf(fp,"mdclatom_pimd:  pimd_on %d\n",pimd_on);
     fprintf(fp,"mdclatom_pimd:  pi_beads %d\n",pi_beads);
     fprintf(fp,"mdclatom_pimd:  natm_tot %d\n",natm_tot);
     fprintf(fp,"mdclatom_pimd:  nfree_pimd %d \n",nfree_pimd);
     fprintf(fp,"mdclatom_pimd:  pimd_freez_typ %d \n",pimd_freez_typ);
     fprintf(fp,"mdclatom_pimd:  pi_beads_full_ter %d \n",pi_beads_full_ter);
     fprintf(fp,"mdclatom_pimd:  pi_beads_res_ter %d \n",pi_beads_res_ter);
     fprintf(fp,"mdclatom_pimd:  pi_beads_full_tra %d \n",pi_beads_full_tra);
     fprintf(fp,"mdclatom_pimd:  pi_beads_res_tra %d \n",pi_beads_res_tra);

     fprintf(fp,"mdclatom_pimd: gamma_adb %lg \n",gamma_adb);
     fprintf(fp,"mdclatom_pimd: pi_temperature %lg \n",pi_temperature);
     fprintf(fp,"mdclatom_pimd: wght_pimd %lg \n",wght_pimd);
     fprintf(fp,"mdclatom_pimd: pi_beads_full_ter_wght %lg \n",
                      pi_beads_full_ter_wght);
     fprintf(fp,"mdclatom_pimd: pi_beads_res_ter_wght %lg \n",
                      pi_beads_res_ter_wght);
     fprintf(fp,"mdclatom_pimd: pi_beads_full_tra_wght %lg \n",
                      pi_beads_full_tra_wght);
     fprintf(fp,"mdclatom_pimd: pi_beads_res_tra_wght %lg \n",
                      pi_beads_res_tra_wght);
     fprintf(fp,"mdclatom_pimd: rcut_spread %lg \n",rcut_spread);

     fclose(fp);

  }// end member function 
//---------------------------------------------------------------------------

}; // end class definition 

#ifdef PUP_ON
PUPmarshall(MDCLATOMS_PIMD);
#endif

#endif
//==========================================================================
