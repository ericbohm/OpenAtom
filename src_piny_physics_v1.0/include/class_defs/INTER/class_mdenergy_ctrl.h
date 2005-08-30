/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                      class_mdenergy_ctrl.h                               */
/*                                                                          */
/*    Class definition for intermolecular force calculation options         */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#ifndef _MDENERGY_CTRL_
#define _MDENERGY_CTRL_

class MDENERGY_CTRL {
 public:
   int pme_on;                        /* Opt: Particle mesh ewald on    */
   int lj_coul_explicit_opt;          /* Opt: grind out those flops baby */
   int int_res_ter,int_res_tra;       /* Opt: Inter/Intra res on/off    */
   int isep_vvdw;                     /* Opt: Separate out Vander Waals */
                                      /*      and Coulomb energies      */
   int iswit_vdw;                     /* Opt: VdW  smoothly switched    */
   int iget_pe_real_inter_freq;       /* Opt: Freq calc of inter PE     */
   int iget_pe_real_inter,iget_pv_real_inter; 
                                      /* Opt: Calculate inter PE and PV */

   int itime;                         /* Num: Present time step;        */
   int iget_full_inter,iget_res_inter;/* Opt: Full or respa inter       */
   int iget_full_intra,iget_res_intra;/* Opt: Full or respa intra       */
   int iget_full_pimd,iget_res_pimd;  /* Opt: Full or respa intra       */


//============================================================================
// Default constructor/destructor

   MDENERGY_CTRL(){
    pme_on      = 0;                        
    lj_coul_explicit_opt = 0;
    int_res_ter = 0;
    int_res_tra = 0;       
    isep_vvdw   = 0;                     
    iswit_vdw   = 0;                     
    iget_pe_real_inter_freq = 0;       
    iget_pe_real_inter = 0;
    iget_pv_real_inter = 0; 
    itime = 0;                         
    iget_full_inter = 0;
    iget_res_inter  = 0;
    iget_full_intra = 0;
    iget_res_intra  = 0;
    iget_full_pimd  = 0;
    iget_res_pimd   = 0;  
   }
   ~MDENERGY_CTRL(){}

//=============================================================================
// Pack/Unpack

#ifdef PUP_ON
   void pup(PUP::er &p){
     // PUP ints
     p | pme_on;  
     p | lj_coul_explicit_opt;
     p | int_res_ter;
     p | int_res_tra; 

     p | isep_vvdw;               
     p | iswit_vdw;               
     p | iget_pe_real_inter_freq; 
     p | iget_pe_real_inter;

     p | iget_pv_real_inter; 
     p | itime;                   
     p | iget_full_inter;
     p | iget_res_inter;

     p | iget_full_intra; 
     p | iget_res_intra;
     p | iget_full_pimd;
     p | iget_res_pimd;  
#ifdef _PARALLEL_DEBUG_        
    if (p.isUnpacking())
     state_class_out ();
#endif  
   } /* end pack unpack */
#endif

//===========================================================================
// Print out state of the class

  void state_class_out(){
    
     char fileName [255];
     sprintf (fileName, "%d_mdenergy_ctrl.state", CkMyPe());
     FILE *fp;  fp = fopen(fileName,"w");
     fprintf(fp,"class_mdenergy_ctrl:  pme_on %d\n",pme_on);  
     fprintf(fp,"class_mdenergy_ctrl:  lj_coul_explicit_opt %d\n",lj_coul_explicit_opt);  
     fprintf(fp,"class_mdenergy_ctrl:  int_res_ter %d\n",int_res_ter);
     fprintf(fp,"class_mdenergy_ctrl:  int_res_tra %d\n",int_res_tra); 

     fprintf(fp,"class_mdenergy_ctrl:  isep_vvdw %d\n",isep_vvdw);               
     fprintf(fp,"class_mdenergy_ctrl:  iswit_vdw %d\n",iswit_vdw);               
     fprintf(fp,"class_mdenergy_ctrl:  iget_pe_real_inter_freq %d\n",iget_pe_real_inter_freq); 
     fprintf(fp,"class_mdenergy_ctrl:  iget_pe_real_inter %d\n",iget_pe_real_inter);

     fprintf(fp,"class_mdenergy_ctrl:  iget_pv_real_inter %d\n",iget_pv_real_inter);
     fprintf(fp,"class_mdenergy_ctrl:  itime %d\n",itime);                   
     fprintf(fp,"class_mdenergy_ctrl:  iget_full_inter %d\n",iget_full_inter);
     fprintf(fp,"class_mdenergy_ctrl:  iget_res_inter %d\n",iget_res_inter);

     fprintf(fp,"class_mdenergy_ctrl:  iget_full_intra %d\n",iget_full_intra); 
     fprintf(fp,"class_mdenergy_ctrl:  iget_res_intra %d\n",iget_res_intra);
     fprintf(fp,"class_mdenergy_ctrl:  iget_full_pimd %d\n",iget_full_pimd);
     fprintf(fp,"class_mdenergy_ctrl:  iget_res_pimd %d\n",iget_res_pimd);  

     fclose(fp);
  }/* end member function */


 }; /* end class definition */
//---------------------------------------------------------------------------

#ifdef PUP_ON
PUPmarshall(MDENERGY_CTRL);
#endif

#endif
//============================================================================
