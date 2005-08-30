//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                                                                          
//                      class_mdlnk_list.h                                   
//                                                                          
//    Class definition for intermolecular link list                         
//                                                                          
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

#ifndef _MDLNKLIST_
#define _MDLNKLIST_

class MDLNKLIST {
 public:
    int ilnk;                    // Opt: Link list option 
    int lnk_ver_update;          // Opt: Link list update Verlist option 
    int ilnk_init;               // Opt: 1st lnklst fill                
    int lnk_for_odd;             // Opt: Force the # of divisions along 
                                 //      each cell edge to be odd       
    int ncell_div_avg;           // Num: Input about how many divisions
                                 //     of each cell edge to use        
    int lnk_excl_safe;           // Num: 
    int lnk_vol_safe;            // Num: 
    int ncell_a,ncell_b,ncell_c; // Num: Cell divisions along a,b,c     
    int natm_cell_max;           // Num: Max # of atms in any cell      
    int lnk_list_dim;            // Num: Size of lnkcell                
    int nshft_lnk;               // Num: Total # of shifts needed to
                                 //       find the inter-mol interacts   
    int nshft_lnk_res;           // Num:  RESPA nshft_lnk               


    double  vol_lnk_map;         // Lst: Lnk volume                     
    double rcut_max,rcut_max_res,rexl_max; //Num: Lnk cell cutoffs      

    int *lnk_list;                 // Lst: List of atms in each cell
                                   // Lth: lnk_list_dim                   
    int *ishft_a,*ishft_b,*ishft_c;// Lst: Vector of each shift (a,b,c)
                                   // Lth: nshft_lnk                    
    int *iexl_chk;                 // Lst: Indicator of if a shift may 
                                   //     contains exlcs;
                                   // Lth: nshft_lnk                      
    int *ishft_a_res,*ishft_b_res,*ishft_c_res; 
                                   // Lst: RESPA shifts
                                   // Num: nshft_lnk_res                  
    int *iexl_chk_res;             // Lst: Excl indicator   
                                   // Lth: nshft_lnk_res                  
    int *natm_cell;                // Lst: # atms in each cell            

    double *shft_wght;             // Lst: Weight of each shift;    
                                   // Lth: nshft_lnk_res                  
    double *shft_wght_res;         // Lst: Weight of each RESPA shift;
                                   // Lth: nshft_lnk_res                  
    double *hmat_lnk,*hmati_lnk;   // Lst: Lnk cell box shape             

//===============================================================================
// Default constructor/destructor

    MDLNKLIST(){
      ilnk          = 0;
      lnk_ver_update =0;
      ilnk_init     = 0;              
      lnk_for_odd   = 0;           
      ncell_div_avg = 0;         
      lnk_excl_safe = 0;         
      lnk_vol_safe  = 0;          
      ncell_a = 0;
      ncell_b = 0;
      ncell_c = 0; 
      natm_cell_max = 0;           
      lnk_list_dim  = 0;            
      nshft_lnk     = 0;               
      nshft_lnk_res = 0;           

      hmat_lnk =(double *) cmalloc(9*sizeof(double),"constructor mdlnklist ")-1;
      hmati_lnk=(double *) cmalloc(9*sizeof(double),"constructor mdlnklist ")-1;
    }
    ~MDLNKLIST(){}

//===============================================================================
// Pack/Unpack

#ifdef PUP_ON
    void pup(PUP::er &p){

      // PUP ints

      p | ilnk;
      p | lnk_ver_update;
      p | ilnk_init;         
      p | lnk_for_odd;        
                              
      p | ncell_div_avg;      
                              
      p | lnk_excl_safe;      
      p | lnk_vol_safe;       

      p | ncell_a;
      p | ncell_b;
      p | ncell_c; 
      p | natm_cell_max;   
      p | lnk_list_dim;    
      p | nshft_lnk;       
                           
      p | nshft_lnk_res;   

      // PUP doubles 

      p | vol_lnk_map;       
      p | rcut_max;
      p | rcut_max_res;
      p | rexl_max;

      // PUP arrays

#ifdef ANYLIST_IMPLEMENTED
      if(ilnk==1 || ilnk_ver==1){
       pup1d_int(p,&lnk_list,lnk_list_dim);
       pup1d_int(p,&ishft_a,nshft_lnk);
       pup1d_int(p,&ishft_b,nshft_lnk);
       pup1d_int(p,&ishft_c,nshft_lnk);
       pup1d_int(p,&iexl_chk,nshft_lnk);
       pup1d_int(p,&ishft_a_res,nshft_lnk_res);
       pup1d_int(p,&ishft_b_res,nshft_lnk_res);
       pup1d_int(p,&ishft_c_res,nshft_lnk_res); 
       pup1d_int(p,&iexl_chk_res,nshft_lnk_res);
       pup1d_int(p,&natm_cell,natm_cell_max);  

       pup1d_dbl(p,&shft_wght,nshft_lnk_res);
       pup1d_dbl(p,&shft_wght_res,nshft_lnk_res); 
       pup1d_dbl(p,&hmat_lnk,9);
       pup1d_dbl(p,&hmati_lnk,9);
      }//endif
#endif
#ifdef _PARALLEL_DEBUG_        
    if (p.isUnpacking())
     state_class_out ();
#endif
  } // end pack unpack 
#endif
//===============================================================================
// Print out state of class

  void state_class_out(){
    

     char fileName [255];
     sprintf (fileName, "%d_mdlnklist.state", CkMyPe());
     FILE *fp;  fp = fopen(fileName,"w");

     //Print Ints

     fprintf(fp,"mdlnklist: ilnk_init %d\n",ilnk_init);         
     fprintf(fp,"mdlnklist: lnk_for_odd %d\n",lnk_for_odd);        
     fprintf(fp,"mdlnklist: ncell_div_avg %d\n",ncell_div_avg);      
     fprintf(fp,"mdlnklist: lnk_excl_safe %d\n",lnk_excl_safe);      
     fprintf(fp,"mdlnklist: lnk_vol_safe %d\n",lnk_vol_safe);       

     fprintf(fp,"mdlnklist: ncell_a %d\n",ncell_a);
     fprintf(fp,"mdlnklist: ncell_b %d\n",ncell_b);
     fprintf(fp,"mdlnklist: ncell_c %d\n",ncell_c); 
     fprintf(fp,"mdlnklist: natm_cell_max %d\n",natm_cell_max);   
     fprintf(fp,"mdlnklist: lnk_list_dim %d\n",lnk_list_dim);    
     fprintf(fp,"mdlnklist: nshft_lnk %d\n",nshft_lnk);       
                          
     fprintf(fp,"mdlnklist: nshft_lnk_res %d\n",nshft_lnk_res); 

      // Print doubles

     fprintf(fp,"mdlnklist: vol_lnk_map %.12g\n",vol_lnk_map);       
     fprintf(fp,"mdlnklist: rcut_max %.12g\n",rcut_max);
     fprintf(fp,"mdlnklist: rcut_max_res %.12g\n",rcut_max_res);
     fprintf(fp,"mdlnklist: rexl_max %.12g\n",rexl_max);

#ifdef ANYLIST_IMPLEMENTED
     if(ilnk==1 || ilnk_ver==1){
     int i;
      // Print arrays
      for(=1;i<=lnk_list_dim;i++){fprintf(fp,"mdlnklist:lnk_list[%d] %d\n"
                                 ,i,lnk_list[i]);}
      for(=1;i<=nshft_lnk;i++){fprintf(fp,"mdlnklist:ishft_a[%d] %d\n"
                                 ,i,ishft_a[i]);}
      for(=1;i<=nshft_lnk;i++){fprintf(fp,"mdlnklist:ishft_b[%d] %d\n"
                                 ,i,ishft_b[i]);}
      for(=1;i<=nshft_lnk;i++){fprintf(fp,"mdlnklist:ishft_c[%d] %d\n"
                                 ,i,ishft_c[i]);}
      for(=1;i<=nshft_lnk;i++){fprintf(fp,"mdlnklist:iexl_chk[%d] %d\n"
                                 ,i,iexl_chk[i]);}
      for(=1;i<=nshft_lnk_res;i++){fprintf(fp,"mdlnklist:ishft_a_res[%d] %d\n"
                                 ,i,ishft_a_res[i]);}
      for(=1;i<=nshft_lnk_res;i++){fprintf(fp,"mdlnklist:ishft_b_res[%d] %d\n"
                                 ,i,ishft_b_res[i]);}
      for(=1;i<=nshft_lnk_res;i++){fprintf(fp,"mdlnklist:ishft_c_res[%d] %d\n"
                                 ,i,ishft_c_res[i]);}
      for(=1;i<=nshft_lnk_res;i++){fprintf(fp,"mdlnklist:iexl_chk_res[%d] %d\n"
                                 ,i,iexl_chk_res[i]);}
      for(=1;i<=natm_cell_max;i++){fprintf(fp,"mdlnklist:natm_cell[%d] %d\n"
                                 ,i,natm_cell[i]);}
      for(=1;i<=nshft_lnk_res;i++){fprintf(fp,"mdlnklist:shft_wght[%d] %.12g\n"
                                 ,i,shft_wght[i]);}
      for(=1;i<=nshft_lnk_res;i++){fprintf(fp,"mdlnklist:shft_wght_res[%d] %g\n"
                                 ,i,shft_wght_res[i]);}
      for(=1;i<=9;i++){fprintf(fp,"mdlnklist:hmat_lnk[%d] %.9g\n"
                                 ,i,hmat_lnk[i]);}
      for(=1;i<=9;i++){fprintf(fp,"mdlnklist:hmati_lnk[%d] %.9g\n"
                                 ,i,hmati_lnk[i]);}
     }//endif
#endif
     fclose(fp);
  }// end member function 


}; // end class definition 

#ifdef PUP_ON
PUPmarshall(MDLNKLIST);
#endif

#endif

