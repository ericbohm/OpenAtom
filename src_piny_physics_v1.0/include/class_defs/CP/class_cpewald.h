//==========================================================================
//  More CP/Ewald stuff 


#ifndef _CPEWALD_
#define _CPEWALD_

class CPEWALD {

 //----------------
 public:
  int cp_any_on;                   // Opt: Is cp of any type ``on''
  int nktot_sm;                    // Num: # of PW coeffon small sphere
                                   //  cutoff g-space grid=ncoef  dens_cp_box
  int nktot_lg;                    // Num: # of PW coeffon large sphere
                                   //  cutoff g-space grid=ncoef  dens_cp_box
  int nktot_dens_cp_box;           // When have 2 boxes  dens_cp_box
                                   // Num: # of PW coeffon large sphere
                                   //      cutoff g-space cp_grid=ncoef_cp_l  
  int box_rat;                     // Box ratio for dual gridding 
                                    


  double dbox_rat;                 // Double value for box ratio 
  double gw_gmin,gw_gmax;
  double gw_gmin_dens_cp_box,gw_gmax_dens_cp_box;


  int *kmax_cp;                    // Lst: Int cutoff in a,b,c directions 
  int *kmax_cp_dens_cp_box;        // Lst: cutoff for the small box       
                                   //  dual=0 kmax_cp_dens_cp_box=kmax_cp 
                                   //  dual=1 kmax_cp_dens_cp_box=kmax_cp 
                                   //         but kmax_cp requires box_rat
                                   //  dual=2 kmax_cp has true values     
  int *nfft;
  int *nfft_dens;
 //----------------
 //con-destruct:
   CPEWALD(){
    nktot_sm          = 0;            
    nktot_lg          = 0;            
    nktot_dens_cp_box = 0;   
    box_rat           = 0;   
    kmax_cp = (int *) cmalloc(3*sizeof(int),"CPEWALD constructor")-1;
    kmax_cp_dens_cp_box = (int *) cmalloc(3*sizeof(int),"CPEWALD constructor")-1;
    nfft = (int *) cmalloc(3*sizeof(int),"CPEWALD constructor")-1;
    nfft_dens = (int *) cmalloc(3*sizeof(int),"CPEWALD constructor")-1;
   };
  ~CPEWALD(){};

#ifdef PUP_ON
 //----------------
 //pupping
  void pup(PUP::er &p){
    //pupping ints
      p | cp_any_on;
      p | nktot_sm;
      p | nktot_lg;
      p | nktot_dens_cp_box;
      p | box_rat;          
    //pupping dbles
      p | dbox_rat;
      p | gw_gmin;
      p | gw_gmax;
      p | gw_gmin_dens_cp_box;
      p | gw_gmax_dens_cp_box;
    //pupping int arrays
    if(cp_any_on==1){
       pup1d_int(p,&kmax_cp,3);
       pup1d_int(p,&kmax_cp_dens_cp_box,3);
       pup1d_int(p,&nfft,3);
       pup1d_int(p,&nfft_dens,3);
    }//endif
#ifdef _PARALLEL_DEBUG_        
    if (p.isUnpacking())
     state_class_out ();
#endif         
  } // end pup
#endif

  void state_class_out(){
     int i;
     char fileName [255];
     sprintf (fileName, "%d_cpewald.state", CkMyPe());
     FILE *fp;  fp = fopen(fileName,"w");
    //ints
      fprintf(fp,"cp_any_on %d\n",cp_any_on);
      fprintf(fp,"nktot_sm %d\n",nktot_sm);
      fprintf(fp,"nktot_lg %d\n",nktot_lg);
      fprintf(fp,"nktot_dens_cp_box %d\n",nktot_dens_cp_box);
      fprintf(fp,"box_rat %d\n",box_rat);          
    //dbles
      fprintf(fp,"dbox_rat %g\n",dbox_rat);
      fprintf(fp,"gw_gmin %g\n",gw_gmin);
      fprintf(fp,"gw_gmax %g\n",gw_gmax);
      fprintf(fp,"gw_gmin_dens_cp_box %g\n",gw_gmin_dens_cp_box);
      fprintf(fp,"gw_gmax_dens_cp_box %g\n",gw_gmax_dens_cp_box);
    //int arrays
      if(cp_any_on==1){
       for(i=1;i<=3;i++){fprintf(fp,"kmax_cp[%d] %d\n",i,kmax_cp[i]);}
       for(i=1;i<=3;i++){fprintf(fp,"kmax_cp_dens_cp_box[%d] %d\n",
                                     i,kmax_cp_dens_cp_box[i]);}
      }//endif
   fclose(fp);
  }// end routine

}; //CPEWALD

#ifdef PUP_ON
PUPmarshall(CPEWALD);
#endif


#endif

//==========================================================================
