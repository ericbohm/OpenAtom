//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                                                                          
//                      class_mdclatoms_tran.h                              
//                                                                          
//         Class definition for classical atom variable transformations     
//         currently used for path integrals                                
//                                                                          
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

#ifndef _MDCLATOMS_TRAN_
#define _MDCLATOMS_TRAN_

class MDCLATOMS_TRAN{
 public:
  int pimd_on;                 // Opt: Path integral are on
  int pi_beads;                // Num: Number of path integral beads
  int natm_tot;                // Num: Total number of atoms     
  double *rat1_stag,*rat2_stag;// Lst: bead convert factor Lth: pi_beads
  double *path_eig;            // Lst: List of bead eigenvalues Lth: pi_beads
  double *x_trans,*y_trans,*z_trans; //  Lth:  pi_beads 
  double *prekf;               // Lst : natm_tot 

//---------------------------------------------------------------------------
// Default constructor/destructor

   MDCLATOMS_TRAN(){}
  ~MDCLATOMS_TRAN(){}

//---------------------------------------------------------------------------
// Pack/Unpack

#ifdef PUP_ON
  void pup(PUP::er &p){
    // PUP the ints
    p | pimd_on;
    p | pi_beads;
    p | natm_tot;
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
     sprintf (fileName, "%d_mdclatom_tran.state", CkMyPe());
     FILE *fp;  fp = fopen(fileName,"w");
       
     fprintf(fp,"mdclatom_tran:  pimd_on %d\n",pimd_on);
     fprintf(fp,"mdclatom_tran:  pi_beads %d\n",pi_beads);
     fprintf(fp,"mdclatom_tran:  natm_tot %d\n",natm_tot);

  }// end member function 
//---------------------------------------------------------------------------

}; // end class definition 

#ifdef PUP_ON
PUPmarshall(MDCLATOMS_TRAN);
#endif

#endif
//==========================================================================
