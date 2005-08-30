//==========================================================================
//                 CP coeff Nose'-Hoover chains                             
//             {Variables needed for mem allocation:                        
//                   num_c_nhc,len_c_nhc }                                  
//                                                                          

#ifndef _CPTHERM_INFO_
#define _CPTHERM_INFO_

class CPTHERM_INFO{

 //----------------
 public:
  int cp_any_on;              // Opt: Is cp of any type ``on''
  int nstate;                 // Num: # of states
  int num_c_nhc;              // Num: # of PW coeff NHC's            
  int len_c_nhc;              // Num: Lnth of PW coeff NHC's         
  int nres_c_nhc;             // Num: # of PW coeff RESPA NHC steps   
  int nyosh_c_nhc;            // Num: # of PW coeff Yosh NHC's steps  
  int massiv_flag;            // opt: Massive thermo flag         
  int istate_nhc_opt;         // opt: NHC option                     


  double cp_therm_heat_fact;  // opt: Sample or Rescale Hot NHCs     
  double dt_nhc,dti_nhc,wght; // Num: NHC respa time steps           
  double c_gkt_massiv;        // Num: gkt variable under massive  
  double cmass_nhc_massiv;    // Num: nhc mass under massive      


  int *icmapup_nhc;           // Map: Coeff index -> coeff NHC index;
                              // Lth: nstate_up               
  int *icmapdn_nhc;           // Map: Coeff index -> coeff NHC index; 
                              // Lth: nstate_dn                        

  double *wdti,*wdti2,*wdti4,*wdti8,*wdti16;// Lst: Yosh steps;  Lth:25

  double **c_gkt;             // Lst: PW coeff NHC deg free;
                              // Lth: num_c_nhc x len_c_nhc          
  double **cmass_nhc;         // Lst: PW coeff NHC mass; 
                              // Lth: num_c_nhc x len_c_nhc            

//----------------------------------------------------------------------------
//con-destruct:
   CPTHERM_INFO(){
    cp_any_on      = 0;
    nstate         = 0;        
    num_c_nhc      = 0;     
    len_c_nhc      = 0;     
    nres_c_nhc     = 0;    
    nyosh_c_nhc    = 0;   
    massiv_flag    = 0;   
    istate_nhc_opt = 0;
    wdti   = (double *)cmalloc(25*sizeof(double),"CPTHERM_INFO constructor")-1;
    wdti2  = (double *)cmalloc(25*sizeof(double),"CPTHERM_INFO constructor")-1;
    wdti4  = (double *)cmalloc(25*sizeof(double),"CPTHERM_INFO constructor")-1;
    wdti8  = (double *)cmalloc(25*sizeof(double),"CPTHERM_INFO constructor")-1;
    wdti16 = (double *)cmalloc(25*sizeof(double),"CPTHERM_INFO constructor")-1;
   }
  ~CPTHERM_INFO(){};

//----------------------------------------------------------------------------
#ifdef PUP_ON
 //pupping
  void pup(PUP::er &p){
    //pupping ints
      p | cp_any_on;      
      p | nstate;      
      p | num_c_nhc;
      p | len_c_nhc;
      p | nres_c_nhc;
      p | nyosh_c_nhc;
      p | massiv_flag;
      p | istate_nhc_opt;
    //pupping dbles
      p | cp_therm_heat_fact;
      p | dt_nhc;
      p | dti_nhc;
      p | wght;
      p | c_gkt_massiv;
      p | cmass_nhc_massiv;
      if(cp_any_on==1 && num_c_nhc> 0 && len_c_nhc > 0){
        pup1d_dbl(p,&wdti,25);
        pup1d_dbl(p,&wdti2,25);
        pup1d_dbl(p,&wdti4,25);
        pup1d_dbl(p,&wdti8,25);
        pup1d_dbl(p,&wdti16,25);
        if(massiv_flag==0){
          pup1d_int(p,&icmapup_nhc,nstate);
          pup1d_int(p,&icmapdn_nhc,nstate);
          pup2d_dbl(p,&c_gkt,num_c_nhc,len_c_nhc);
          pup2d_dbl(p,&cmass_nhc,num_c_nhc,len_c_nhc);
	}//endif
      }//endif
#ifdef _PARALLEL_DEBUG_        
    if (p.isUnpacking())
     state_class_out ();
#endif           
  } // end pup
#endif

//----------------------------------------------------------------------------
  void state_class_out(){

      int i,j;
      char fileName [255];
      sprintf (fileName, "%d_cptherm_info.state", CkMyPe());
      FILE *fp;  fp = fopen(fileName,"w");

    // ints
      fprintf(fp,"num_c_nhc %d\n",num_c_nhc);
      fprintf(fp,"len_c_nhc %d\n",len_c_nhc);
      fprintf(fp,"nres_c_nhc %d\n",nres_c_nhc);
      fprintf(fp,"nyosh_c_nhc %d\n",nyosh_c_nhc);
      fprintf(fp,"massiv_flag %d\n",massiv_flag);
      fprintf(fp,"istate_nhc_opt %d\n",istate_nhc_opt);
    // dbles
      fprintf(fp,"cp_therm_heat_fact %g\n",cp_therm_heat_fact);
      fprintf(fp,"dt_nhc %g\n",dt_nhc);
      fprintf(fp,"dti_nhc %g\n",dti_nhc);
      fprintf(fp,"wght %g\n",wght);
      fprintf(fp,"c_gkt_massiv %g\n",c_gkt_massiv);
      fprintf(fp,"cmass_nhc_massiv %g\n",cmass_nhc_massiv);
    //dble arrays
      if(cp_any_on==1){
       for(i=1;i<=25;i++){fprintf(fp,"wdti[%d] %g\n",i,wdti[i]);}
       for(i=1;i<=25;i++){fprintf(fp,"wdti2[%d] %g\n",i,wdti2[i]);}
       for(i=1;i<=25;i++){fprintf(fp,"wdti4[%d] %g\n",i,wdti4[i]);}
       for(i=1;i<=25;i++){fprintf(fp,"wdti8[%d] %g\n",i,wdti8[i]);}
       for(i=1;i<=25;i++){fprintf(fp,"wdti16[%d] %g\n",i,wdti16[i]);}
       for(i=1;i<=num_c_nhc;i++){
       for(j=1;j<=len_c_nhc;j++){
         fprintf(fp,"c_gkt[%d][%d] %g\n",i,j,c_gkt[i][j]);
         fprintf(fp,"cmass_nhc[%d][%d] %g\n",i,j,cmass_nhc[i][j]);
       }}// end :i,j
      }//endif
    fclose(fp);
  } // end routine 

//----------------------------------------------------------------------------
  }; //CPTHERM_INFO
//==========================================================================

#ifdef PUP_ON
PUPmarshall(CPTHERM_INFO);
#endif


#endif
//==========================================================================
