//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                                                                          
//                      class_mdbrnch_root.                                 
//                                                                          
//    Class definition for branch root list scheme                          
//                                                                          
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

#ifndef _MDBRNCH_ROOT_
#define _MDBRNCH_ROOT_

class MDBRNCH_ROOT {
 public:
 int brnch_root_list_opt;      /* Opt: Option to shave brnch interactions */
 int natm_tot;                 /* Num: Total number of atoms */
 int nbrnch_tot;               /* Num: Number of brnchs*/
 int nroot_tot;                /* Num: Number of roots*/
 int natm_add_tot;             /* Num: Number of additions*/
 int nbrnch_root_max;          /* Num: Max Number of brnch's off a root*/

 double brnch_root_dist_max;   /* Num: Max brnch root distance*/
 double cutdif_min_brnch_root; /* Num: MIN(root_root_cut - root_brnch_cut)*/
 double cutdif_min_brnch_brnch;/* Num: MIN(root_root_cut - brnch_brnch_cut)*/
 double cutdif_min_brnch_root_res;/* Num: MIN(root_root_cut - root_brnch_cut)*/
 double cutdif_min_brnch_brnch_res;/* Num: MIN(root_root_cut-brnch_brnch_cut)*/
 double r2_nocheck_dist;       /* Num: pare down no check distance*/
 double r2_nocheck_dist_res;   /* Num: pare down no check distance*/
 double cutskin_bb_min;        /* Num: minimum branch-branch cutskin*/
 double cutskin_bb_min_res;    /* Num: minimum branch-branch cutskin*/

 int *brnch_atm_list; /* Lst: List of brnch atoms
                         Exp: brnch_atm_list[3] = 25; 
                              The third brnch is atm 25 
                         Lth: nbrnch_tot */
 int *brnch_atm_root; /* Lst: List of brnch atom root
                         Exp: brnch_atm_roots[3] = 25; 
                              The third brnch's root is atm 25 
                         Lth: nbrnch_tot */
 int *root_atm_list;  /* Lst: List of root atoms
                         Exp: root_atm_list[3] = 25; 
                              The third root is atm 25 
                          Lth: nroot_tot */
 int *root_atm_map;   /* Map: Map of atoms into the root list 
                          Exp: root_atm_map[25] = 3; 
                              The 25th atom is the third root
                          Lth: natm_tot */
 int *brnch_atm_map;  /* Map: Map of atoms of big list into the branch list 
                          Exp: branch_atm_map[25] = 3; 
                              The 25th atom is the third branch
                          Lth: natm_tot */
 int *nbrnch_of_root_big;/* Lst: The number of brnchs of each atom in big
                                 list including self. Brnch atoms have 0. 
                          Exp: nbrnch_of_root[3] = 2; 
                              The 3rd atom has two brnchs.
                          Lth: nroot_tot */
 int *nbrnch_of_root;/* Lst: The number of brnchs of each root
                             doesn't not included itself.
                          Exp: nbrnch_of_root[3] = 2; 
                              The 3rd root has two brnchs
                          Lth: nroot_tot */
 int **ibrnch_of_root;/* Lst: The list of brnchs of each root
                          Exp: ibrnch_of_root[3][1] = 30; 
                              The first brnch of the 3rd root is atm 30.
                          Lth: nroot_tot x nbrnch_max */
 int *natm_add;       /* Lst: The number of atm's to be added to each atm's
                              neighbor list to correct for the root-brnch 
                              method
                          Exp: natm_add[3] = 2; 
                              The 3rd atm needs 2 atm's added to its list
                          Lth: natm_tot */
 int *iatm_add;       /* Lst: The list of atm's to be added to each atm's 
                              list
                          Exp: iatm_add[2+iatm_add_off[3]] = 27; 
                              The 2nd atm added to the list of the third 
                              atom is atm number 27
                          Lth: natm_add_tot */
 int *iatm_add_off;     /* Lst: The offset list
                          Exp: iatm_add[2+iatm_add_off[3]] = 27; 
                              The 2nd atm added to the list of the third 
                              atom is atm number 27
                          Lth: natm_add_tot */


//----------------------------------------------------------------------------
//  Default constructor/destructor

   MDBRNCH_ROOT(){
     brnch_root_list_opt = 0;
     natm_tot   = 0;              
     nbrnch_tot = 0;            
     nroot_tot  = 0;             
     natm_add_tot    = 0;          
     nbrnch_root_max = 0;       
   }
  ~MDBRNCH_ROOT(){}

//----------------------------------------------------------------------------
// Pack/Unpack

#ifdef PUP_ON
  void pup(PUP::er &p){
    // PUP ints
    p | brnch_root_list_opt;
    p | natm_tot;
    p | nbrnch_tot;    
    p | nroot_tot;     
    p | natm_add_tot;  
    p | nbrnch_root_max;

    // PUP doubles 
    p | brnch_root_dist_max;   
    p | cutdif_min_brnch_root; 
    p | cutdif_min_brnch_brnch;
    p | cutdif_min_brnch_root_res;
    p | cutdif_min_brnch_brnch_res;
    p | r2_nocheck_dist;       
    p | r2_nocheck_dist_res;   
    p | cutskin_bb_min;        
    p | cutskin_bb_min_res;    

    // PUP arrays
#ifdef ANYLIST_IMPLEMENTED
    if(brnch_root_list_opt > 0){
      pup1d_int(p,&brnch_atm_list,nbrnch_tot);
      pup1d_int(p,&brnch_atm_root,nbrnch_tot);
      pup1d_int(p,&root_atm_list,nroot_tot);
      pup1d_int(p,&root_atm_map,natm_tot);
      pup1d_int(p,&brnch_atm_map,natm_tot);
      pup1d_int(p,&nbrnch_of_root_big,nroot_tot);
      pup1d_int(p,&nbrnch_of_root,nroot_tot);
      pup2d_int(p,&ibrnch_of_root,nroot_tot,nbrnch_root_max);
      pup1d_int(p,&natm_add,natm_tot);
      pup1d_int(p,&iatm_add,natm_add_tot);
      pup1d_int(p,&iatm_add_off,natm_add_tot);
    }// endif
#endif
#ifdef _PARALLEL_DEBUG_        
    if (p.isUnpacking())
     state_class_out ();
#endif
  } // end pack unpack 
#endif

//----------------------------------------------------------------------------
// Print out state of the class
  void state_class_out(){
     char fileName [255];
     sprintf (fileName, "%d_mdbrnch_root.state", CkMyPe());
     FILE *fp;  fp = fopen(fileName,"w");

    // Print ints

    fprintf(fp,"mdbrnch_root: brnch_root_list_opt %d\n",brnch_root_list_opt);
    fprintf(fp,"mdbrnch_root: natm_tot %d\n",natm_tot);
    fprintf(fp,"mdbrnch_root: nbrnch_tot %d\n",nbrnch_tot);    
    fprintf(fp,"mdbrnch_root: nroot_tot %d\n",nroot_tot);     
    fprintf(fp,"mdbrnch_root: natm_add_tot %d\n",natm_add_tot);  
    fprintf(fp,"mdbrnch_root: nbrnch_root_max %d\n",nbrnch_root_max);

    // Print doubles 
    fprintf(fp,"mdbrnch_root: brnch_root_dist_max %.12g\n",brnch_root_dist_max);
    fprintf(fp,"mdbrnch_root: cutdif_min_brnch_root %.12g\n",
                              cutdif_min_brnch_root); 
    fprintf(fp,"mdbrnch_root: cutdif_min_brnch_brnch %.12g\n",
                              cutdif_min_brnch_brnch);
    fprintf(fp,"mdbrnch_root: cutdif_min_brnch_root_res %.12g\n",
                              cutdif_min_brnch_root_res);
    fprintf(fp,"mdbrnch_root: cutdif_min_brnch_brnch_res %.12g\n",
                              cutdif_min_brnch_brnch_res);
    fprintf(fp,"mdbrnch_root: r2_nocheck_dist %.12g\n",
                              r2_nocheck_dist);       
    fprintf(fp,"mdbrnch_root: r2_nocheck_dist_res %.12g\n",
                              r2_nocheck_dist_res);   
    fprintf(fp,"mdbrnch_root: cutskin_bb_min %.12g\n",cutskin_bb_min);        
    fprintf(fp,"mdbrnch_root: cutskin_bb_min_res %.12g\n",cutskin_bb_min_res);

    // Print arrays
#ifdef ANYLIST_IMPLEMENTED
   int i,j;
    if(brnch_root_list_opt > 0){
      for(i=1;i<=nbrnch_tot;i++){
       fprintf(fp,"mdbrnch_root: brnch_atm_list[%d] %d\n",i,brnch_atm_list[i]);}
      for(i=1;i<=nbrnch_tot;i++){
       fprintf(fp,"mdbrnch_root: brnch_atm_root[%d] %d\n",i,brnch_atm_root[i]);}
      for(i=1;i<=nroot_tot;i++){
       fprintf(fp,"mdbrnch_root: root_atm_list[%d] %d\n",i,root_atm_list[i]);}
      for(i=1;i<=natm_tot;i++){
       fprintf(fp,"mdbrnch_root: root_atm_map[%d] %d\n",i,root_atm_map[i]);}
      for(i=1;i<=natm_tot;i++){
       fprintf(fp,"mdbrnch_root: brnch_atm_map[%d] %d\n",i,brnch_atm_map[i]);}
      for(i=1;i<=nroot_tot;i++){  
       fprintf(fp,"mdbrnch_root: nbrnch_of_root_big[%d] %d\n",
                              i,nbrnch_of_root_big[i]);}
      for(i=1;i<=nroot_tot;i++){
       fprintf(fp,"mdbrnch_root: nbrnch_of_root[%d] %d\n",i,nbrnch_of_root[i]);}
      for(i=1;i<=nroot_tot;i++){
       for(j=1;j<=nbrnch_root_max;j++){
        fprintf(fp,"mdbrnch_root:  ibrnch_of_root[%d][%d] %d\n",
                           i,j,ibrnch_of_root[i][j]);
       }// endfor j
      }// endfor i 
      for(i=1;i<=natm_tot;i++){
       fprintf(fp,"mdbrnch_root: natm_add[%d] %d\n",i,natm_add[i]);}
      for(i=1;i<=natm_add_tot;i++){
       fprintf(fp,"mdbrnch_root: iatm_add[%d] %d\n",i,iatm_add[i]);}
      for(i=1;i<=natm_add_tot;i++){
       fprintf(fp,"mdbrnch_root: iatm_add_off[%d] %d\n",i,iatm_add_off[i]);}
    }//endif : opt on
#endif
    fclose(fp);
  }/* end member function */
//----------------------------------------------------------------------------

}; /* end class definition */

#ifdef PUP_ON
PUPmarshall(MDBRNCH_ROOT);
#endif

#endif
