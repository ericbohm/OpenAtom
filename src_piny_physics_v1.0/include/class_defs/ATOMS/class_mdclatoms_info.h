//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                                                                          
//                      class_mdclatoms_info.h                              
//                                                                          
//                Class definition for classical atom information           
//                                                                          
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

#ifndef _MDCLATOMS_INFO_
#define _MDCLATOMS_INFO_

class MDCLATOMS_INFO{
  public:
    int natm_tot;                /* Num: # of atm                          */
    int nmol_typ;                /* Num: Number of molecule types          */
    int pimd_on;                 /* Opt : path integrals on                */
    int pi_beads;                /* Num: # of path integral beads          */
    int nfree;                   /* Num: # degrees of freedom              */
    int nchrg;                   /* Num: # chrged atms       Lth:<natm_tot */
    int cp_grimme;               /* Opt: Grimme Vdw correction */

    double mass_sc_fact;         /* Num: classical mass scaling factor     */
    double rheal;                /* num: cutoff healing length */
    double s6Grimme;             /* Num: Grimme scaling factor*/
    double dGrimme;              /* Num: Grimme Fermi functoin parameter   */

    int *ichrg;                  /* lst: indx of chrged atms lth:<natm_tot */

    double *mass;                /* Lst: atm masses          Lth: natm_tot */
    double *q;                   /* Lst: atm chrges          Lth: natm_tot */
    double *alp_pol;             /* Lst: atm polarizability  Lth: natm_tot */
    double *b_neut;              /* Lst: atm scattering fact Lth: natm_tot */ 
    double *text_atm;            /* Lst: Atomic temperaures  Lth: natm_tot */
    double *c6Grimme;            /* Lst: Grimme C6           Lth: natm_tot */
    double *r0Grimme;            /* Lst: Grimme R0           Lth: natm_tot */
    double *text_mol;            /* Lst: mol  temperaures    Lth: nmol_typ */

    //-----------------------------------------------------------------------------
    // Default constructor/destructor

    MDCLATOMS_INFO(){
      natm_tot   = 0;    
      nmol_typ   = 0;    
      pimd_on    = 0;    
      pi_beads   = 0;    
      nfree      = 0;       
      nchrg      = 0;               
      cp_grimme  = 0;
      s6Grimme   = 0.0; 
      dGrimme    = 20.0;  
      rheal      = 0.5/0.529177; 
    }
    ~MDCLATOMS_INFO(){}

    //-----------------------------------------------------------------------------
    // Pack/Unpack

#ifdef PUP_ON
    void pup(PUP::er &p){
      //PUP the ints
      p | natm_tot;            
      p | nmol_typ;
      p | pimd_on;            
      p | pi_beads;            
      p | nfree;               
      p | nchrg;               
      p | cp_grimme;            

      // PUP the doubles
      p | mass_sc_fact;          
      p | s6Grimme;
      p | dGrimme;
      p | rheal;
      // PUP the arrays

      pup1d_int(p,&ichrg,nchrg);    

      pup1d_dbl(p,&mass,natm_tot);
      pup1d_dbl(p,&q,natm_tot);
      pup1d_dbl(p,&alp_pol,natm_tot);
      pup1d_dbl(p,&b_neut,natm_tot);
      pup1d_dbl(p,&text_atm,natm_tot);
      if(cp_grimme==1){
        pup1d_dbl(p,&c6Grimme,natm_tot);
        pup1d_dbl(p,&r0Grimme,natm_tot);
      }//endif
      pup1d_dbl(p,&text_mol,nmol_typ);
#ifdef _PARALLEL_DEBUG_        
      if (p.isUnpacking())
        state_class_out ();
#endif     
    }// end pack/unpack
#endif

    //-------------------------------------------------------------------------------
    // Write out the state of the class

    void state_class_out(){

      int i;
      char fileName [255];
      sprintf (fileName, "%d_mdclatom_info.state", CkMyPe());
      FILE *fp;  fp = fopen(fileName,"w");

      //Print the ints
      fprintf(fp,"mdclatoms_info:  natm_tot %d\n",natm_tot);
      fprintf(fp,"mdclatoms_info:  nmol_typ %d\n",nmol_typ);
      fprintf(fp,"mdclatoms_info:  pimd_on %d\n",pimd_on);
      fprintf(fp,"mdclatoms_info:  pi_beads %d\n",pi_beads);
      fprintf(fp,"mdclatoms_info:  nfree %d\n",nfree);
      fprintf(fp,"mdclatoms_info:  nchrg %d\n",nchrg);

      // Print the doubles
      fprintf(fp,"mdclatoms_info:  mass_sc_fact %.12g\n",mass_sc_fact);

      // Arrays
      for(i=1;i<=nchrg;i++){fprintf(fp,"mdclatoms_info:  ichrg[%d] %d\n",
          i,ichrg[i]);}

      for(i=1;i<=natm_tot;i++){fprintf(fp,"mdclatoms_info:  mass[%d] %.12g\n",
          i,mass[i]);}
      for(i=1;i<=natm_tot;i++){fprintf(fp,"mdclatoms_info:  q[%d] %.12g\n",
          i,q[i]);}
      for(i=1;i<=natm_tot;i++){fprintf(fp,"mdclatoms_info:  alp_pol[%d] %.12g\n",
          i,alp_pol[i]);}
      for(i=1;i<=natm_tot;i++){fprintf(fp,"mdclatoms_info:  b_neut[%d] %.12g\n",
          i,b_neut[i]);}

      for(i=1;i<=natm_tot;i++){fprintf(fp,"mdclatoms_info:  text_atm[%d] %.12g\n",
          i,text_atm[i]);}
      for(i=1;i<=nmol_typ;i++){fprintf(fp,"mdclatoms_info:  text_mol[%d] %.12g\n",
          i,text_mol[i]);}

      fclose(fp);
    } // end member function 
    //-------------------------------------------------------------------------------

}; // end class definition

#ifdef PUP_ON
PUPmarshall(MDCLATOMS_INFO);
#endif

#endif

//==========================================================================
