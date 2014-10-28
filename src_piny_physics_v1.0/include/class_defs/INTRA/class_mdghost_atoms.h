//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                                                                          
//                      class_mdghost_atoms.h                               
//                                                                          
//         Class definition for ghost atom information                      
//                                                                          
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

#ifndef _MDGHOST_ATOMS_
#define _MDGHOST_ATOMS_

class MDGHOST_ATOMS {
  public:
    int natm_tot;
    int nghost_tot;              /* Num: # of ghost atoms                  */
    int natm_comp_max;           /* Num: Max # of atoms comprising a ghost */

    int *ighost_flag;            /* Lst: 1/0 if atoms is/is not a ghost
Lth: natm_tot                     */
    int *ighost_map;             /* Map: Lst of indicies of ghost atoms
                                    in atoms list. Ex: The 3rd ghost
                                    is atom number 27. Lth:nghost_tot */
    int *natm_comp;              /* Lst: # of atoms comprising each ghost  */

    int **iatm_comp;             /* Lst: Atoms comprising each ghost       */

    double **coef;               /* Lst: Coefficients of atoms comprising
                                    the ghost                         */

    //==========================================================================
    // Default constructor/destructor

    MDGHOST_ATOMS(){
      natm_tot      = 0;           
      nghost_tot    = 0;           
      natm_comp_max = 0;        
    }
    ~MDGHOST_ATOMS(){}

    //==========================================================================
#ifdef PUP_ON
    void pup(PUP::er &p){

      // PUP ints

      p | natm_tot;        
      p | nghost_tot;        
      p | natm_comp_max;    

      // PUP arrays

      if(nghost_tot>0){
        pup1d_int(p,&ighost_map,nghost_tot);
        pup1d_int(p,&natm_comp,nghost_tot);
        pup2d_int(p,&iatm_comp,natm_comp_max,nghost_tot);
        pup2d_dbl(p,&coef,natm_comp_max,nghost_tot,"mdghost_atoms");
      }//endif

      if(natm_tot>0){
        pup1d_int(p,&ighost_flag,natm_tot);
      }//endif
#ifdef _PARALLEL_DEBUG_        
      if (p.isUnpacking())
        state_class_out ();
#endif 
    }/* end pack/unpack */
#endif

    //==========================================================================
    void state_class_out(){
      int i;
      char fileName [255];
      sprintf (fileName, "%d_mdghost_atoms.state", CkMyPe());
      FILE *fp;  fp = fopen(fileName,"w");

      fprintf(fp,"nghost_tot %d\n",nghost_tot);        
      fprintf(fp,"natm_comp_max %d\n",natm_comp_max);    
      //  arrays
      for(i=1;i<=nghost_tot;i++){fprintf(fp,"ighost_map[%d] %d ",i,ighost_map[i]);}
      for(i=1;i<=nghost_tot;i++){fprintf(fp,"natm_comp[%d] %d ",i,natm_comp[i]);}
      fclose(fp);
    }// end routine

}; // end class definition
//==========================================================================

#ifdef PUP_ON
PUPmarshall(MDGHOST_ATOMS);
#endif

#endif
//==========================================================================
