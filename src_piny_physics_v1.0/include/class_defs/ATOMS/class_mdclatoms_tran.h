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
    if(pimd_on == 1){
      pup1d_dbl(p,&rat1_stag,pi_beads);
      pup1d_dbl(p,&rat2_stag,pi_beads); 
      pup1d_dbl(p,&path_eig,pi_beads);
      pup1d_dbl(p,&x_trans,pi_beads);
      pup1d_dbl(p,&y_trans,pi_beads);
      pup1d_dbl(p,&z_trans,pi_beads);
      pup1d_dbl(p,&prekf,natm_tot);
    }// endif 
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

     if(pimd_on == 1){
       for(int i=1;i<=pi_beads; i++){
           fprintf(fp,"mdclatom_tran:  rat1_stag[%d] %.12g\n",i,rat1_stag[i]);}
       for(int i=1;i<=pi_beads; i++){
           fprintf(fp,"mdclatom_tran:  rat2_stag[%d] %.12g\n",i,rat2_stag[i]);}
       for(int i=1;i<=pi_beads; i++){
           fprintf(fp,"mdclatom_tran:   path_eig[%d] %.12g\n",i,path_eig[i]);}
       for(int i=1;i<=pi_beads; i++){
           fprintf(fp,"mdclatom_tran:    x_trans[%d] %.12g\n",i,x_trans[i]);}
       for(int i=1;i<=pi_beads; i++){
           fprintf(fp,"mdclatom_tran:    y_trans[%d] %.12g\n",i,y_trans[i]);}
       for(int i=1;i<=pi_beads; i++){
           fprintf(fp,"mdclatom_tran:    z_trans[%d] %.12g\n",i,z_trans[i]);}
       for(int i=1;i<=natm_tot; i++){
           fprintf(fp,"mdclatom_tran:      prekf[%d] %.12g\n",i,prekf[i]);}
     } // endif
     fclose(fp);

  }// end member function 
//---------------------------------------------------------------------------

}; // end class definition 

#ifdef PUP_ON
PUPmarshall(MDCLATOMS_TRAN);
#endif

#endif
//==========================================================================
