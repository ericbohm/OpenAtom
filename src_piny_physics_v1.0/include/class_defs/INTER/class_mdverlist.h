/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                      class_mdverlist.h                                   */
/*                                                                          */
/*    Class definition for intermolecular Verlet list                       */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#ifndef _MDVERLIST_
#define _MDVERLIST_

class MDVERLIST {
  public:
    int iver;                    /* Opt: Verlet on option */
    int natm_tot;                /* Num: number of atoms */
    int iver_init;               /* Opt: 1st verlist fill               */
    int iver_fill,iver_count;    /* Opt: Fill and count options         */
    int nolst_ver_update,lnk_ver_update; 
    /* Opt: Update using nolst or lnk list */
    int jver_pad;                /* Num: Padding used in lnk_ver_update */
    int nmem_min_lst;            /* Num: */
    int nver_lst;                /* Num: # of atms in Verlet list       */
    int nver_lst_res;            /* Num: # of atms in RESPA Verlet list */ 
    int nver_lst_now;         /* Num: actual # of atms in Verlet list       */
    int nver_lst_now_res;     /* Num: actual # of atms in RESPA Verlet list */ 
    int *nter;                   /* Lst: # of neighbors of each atm: 
Lth: natm_tot                         
Can be eliminated using joff        */
    int *jter;           /* Lst: Indices of neighboring atms.  
Lth: nver_lst                       */
    int *jver_off;               /* Lst: Starting pts of nghbors in jter 
Lth: natm_tot                       */
    int *nter_res;               /* Lst: RESPA nter                        
Lth: natm_tot                          
Can be eliminated using joff_res    */
    int *jter_res;       /* Lst: RESPA jter
Lth: nver_lst_res                   */
    int *jver_off_res;           /* Lst: RESPA jver      
Lth: natm_tot                       */
    double mem_safe;                /* Num: */


    //=====================================================================================
    // Default constructor/destructor

    MDVERLIST(){
      iver       = 0;
      natm_tot   = 0;              
      iver_init  = 0;             
      iver_fill  = 0;
      iver_count = 0;  
      nolst_ver_update = 0;
      lnk_ver_update   = 0; 
      jver_pad     = 0;              
      nmem_min_lst = 0;          
      nver_lst     = 0;              
      nver_lst_res = 0;          
      nver_lst_now = 0;        
      nver_lst_now_res = 0;    
    }
    ~MDVERLIST(){}

    //=====================================================================================
    // Pack/Unpack

#ifdef PUP_ON
    void pup(PUP::er &p){

      // PUP ints

      p | iver;
      p | natm_tot;
      p | iver_init;    
      p | iver_fill;
      p | iver_count;  
      p | nolst_ver_update;
      p | lnk_ver_update; 

      p | jver_pad;              
      p | nmem_min_lst;          
      p | nver_lst;              
      p | nver_lst_res;          
      p | nver_lst_now;         
      p | nver_lst_now_res;     

      // PUP doubles 

      p | mem_safe;

      // PUP arrays

#ifdef ANYLIST_IMPLEMENTED
      if(iver == 1){
        pup1d_int(p,&nter,natm_tot);
        pup1d_int(p,&jter,nver_lst);  
        pup1d_int(p,&jver_off,natm_tot);  
        pup1d_int(p,&nter_res,natm_tot);
        pup1d_int(p,&jter_res,nver_lst_res);  
        pup1d_int(p,&jver_off_res,natm_tot);
      }/* endif */
#endif
#ifdef _PARALLEL_DEBUG_        
      if (p.isUnpacking())
        state_class_out ();
#endif  
    } /* end pack unpack */
#endif

    //=====================================================================================
    // Print out state of class

    void state_class_out(){

      char fileName [255];
      sprintf (fileName, "%d_mdverlist.state", CkMyPe());
      FILE *fp;  fp = fopen(fileName,"w");

      // Print ints

      fprintf(fp,"mdverlist:   natm_tot %d\n",natm_tot);
      fprintf(fp,"mdverlist:   iver_init %d\n",iver_init);    
      fprintf(fp,"mdverlist:   iver_fill %d\n",iver_fill);
      fprintf(fp,"mdverlist:   iver_count %d\n",iver_count);  
      fprintf(fp,"mdverlist:   nolst_ver_update %d\n",nolst_ver_update);
      fprintf(fp,"mdverlist:   lnk_ver_update %d\n",lnk_ver_update); 

      fprintf(fp,"mdverlist:   jver_pad %d\n",jver_pad);              
      fprintf(fp,"mdverlist:   nmem_min_lst %d\n",nmem_min_lst);          
      fprintf(fp,"mdverlist:   nver_lst %d\n",nver_lst);              
      fprintf(fp,"mdverlist:   nver_lst_res %d\n",nver_lst_res);          
      fprintf(fp,"mdverlist:   nver_lst_now %d\n",nver_lst_now);         
      fprintf(fp,"mdverlist:   nver_lst_now_res %d\n",nver_lst_now_res);     

      // Print doubles 

      fprintf(fp,"mdverlist:   mem_safe %.12g\n",mem_safe);

      fclose(fp);
    }/* end member function */



};  /* end class definition */

#ifdef PUP_ON
PUPmarshall(MDVERLIST);
#endif


#endif




