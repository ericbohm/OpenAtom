//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                                                                          
//                          class_mdexcl.h                                  
//                                                                          
//            Class definition for exclusion lists                          
//                                                                          
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

#ifndef _MDEXCL_
#define _MDEXCL_

class MDEXCL {
  public:
    int natm_tot;                /* Num: Total number of atoms */
    int nlst;                    /* Num: Total # atm pairs in excl list */
    int nlst_root;               /* Num: Total # atm pairs in excl list */
    int num_cp;
    int brnch_root_list_opt;     /* Num: Option to shave brnch interactions */
    int nroot_tot;

    int *num,*num_root;          /* Lst: # of excls of each atm: 
Exp: Atm 27 has 3 exclusions all of
which have  indices < 27;
Lth: natm_tot                        
Can be eliminated using joff        */
    int *j,*j_root;              /* Lst: Indices of excluded atms.   
Lth: num                            */
    int *j_off,*j_off_root;      /* Lst: Starting pts in *j of excls of 
                                    atms 2-natm_tot;
Exp: Excls of atm  27 start at index 
102 in *j
Lth: natm_tot                       */
    int *j1_cp,*j2_cp;  /*Lst: Indices of atms needed for ewald summation */
    /*  between ab initio and classical atoms         */

    //-----------------------------------------------------------------------------
    // Default constructor/destructor

    MDEXCL(){
      natm_tot           = 0;
      nlst               = 0;                
      nlst_root          = 0;           
      num_cp             = 0;
      brnch_root_list_opt=0;
      nroot_tot          = 0;
    }

    ~MDEXCL(){}

    //-----------------------------------------------------------------------------
    // Pack/Unpack

#ifdef PUP_ON
    void pup(PUP::er &p){

      // PUP integers
      p | natm_tot;
      p | nlst;
      p | nlst_root;
      p | num_cp;
      p | brnch_root_list_opt;
      p | nroot_tot;

      // PUP arrays
      pup1d_int(p,&num,natm_tot);
      pup1d_int(p,&j_off,natm_tot);

      if(nroot_tot>0 && brnch_root_list_opt>0){
        pup1d_int(p,&num_root,nroot_tot);
        pup1d_int(p,&j_off_root,nroot_tot);
      }// endif

      if(nlst>0){pup1d_int(p,&j,nlst);}

      if(nlst_root>0 && brnch_root_list_opt>0){pup1d_int(p,&j_root,nlst_root);}

      if(num_cp>0){pup1d_int(p,&j1_cp,num_cp);}

      if(num_cp>0){pup1d_int(p,&j2_cp,num_cp);}

#ifdef _PARALLEL_DEBUG_        
      if (p.isUnpacking())
        state_class_out ();
#endif     
    } // end pack/unpack
#endif

    //-----------------------------------------------------------------------------
    // Print out state of function 
    void state_class_out(){
      int i;
      char fileName [255];
      sprintf (fileName, "%d_mdexcl.state", CkMyPe());
      FILE *fp;  fp = fopen(fileName,"w");


      // Print integers
      fprintf(fp,"mdexcl:  natm_tot %d\n",natm_tot);
      fprintf(fp,"mdexcl:  nlst %d\n",nlst);  
      fprintf(fp,"mdexcl:  nlst_root %d\n",nlst_root);    
      fprintf(fp,"mdexcl:  brnch_root_list_opt %d\n",brnch_root_list_opt);
      fprintf(fp,"mdexcl:  nroot_tot %d\n",nroot_tot);
      fprintf(fp,"mdexcl:  num_cp %d\n",num_cp);
      // Print Integer Arrays 
      for(i=1;i<=nlst;i++){fprintf(fp,"mdexcl: j[%d] %d\n",i,j[i]);}
      for(i=1;i<=natm_tot;i++){fprintf(fp,"mdexcl: j_off[%d] %d\n",i,j_off[i]);}
      for(i=1;i<=natm_tot;i++){fprintf(fp,"mdexcl: num[%d] %d\n",i,num[i]);}
      if(brnch_root_list_opt>0){
        for(i=1;i<=nlst_root;i++){fprintf(fp,"mdexcl: j_root[%d] %d\n",
            i,j_root[i]);}
        for(i=1;i<=nroot_tot;i++){fprintf(fp,"mdexcl: j_off_root[%d] %d\n",
            i,j_off_root[i]);}
        for(i=1;i<=nroot_tot;i++){fprintf(fp,"mdexcl: num_root[%d] %d\n",
            i,num_root[i]);}
      }// endif
      for(i=1;i<=num_cp;i++){fprintf(fp,"mdexcl: j1_cp[%d] %d\n",i,j1_cp[i]);}
      for(i=1;i<=num_cp;i++){fprintf(fp,"mdexcl: j2_cp[%d] %d\n",i,j2_cp[i]);}

      fclose(fp);

    }// end member function 
    //-----------------------------------------------------------------------------

}; // end class definition 

#ifdef PUP_ON
PUPmarshall(MDEXCL);
#endif

#endif
//==========================================================================
