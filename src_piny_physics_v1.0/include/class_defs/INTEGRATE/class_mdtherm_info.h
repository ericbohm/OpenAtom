//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                                                                          
//                      class_mdtherm_info.h                                
//                                                                          
//         Class definition for atom thermostat information                 
//                                                                          
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

#ifndef _MDTHERM_INFO_
#define _MDTHERM_INFO_

class MDTHERM_INFO{
 public:
  int my_type;                 /* Opt: 0-- classical 1--- bead */
  int natm_tot;                /* Num: Total number of atoms          */
  int iextended_on;            /* Opt: extended system on             */
  int num_nhc;                 /* Num: # of NHC's                     */
  int len_nhc;                 /* Num: Length of NHC's                */
  int nres_nhc;                /* Num: # of RESPA NHC steps           */  
  int nyosh_nhc;               /* Num: # of Yosh NHC steps            */  
  int therm_typ;               /* Num: 1=NHC, 2=GGMT                  */
  int isokin_opt;
  int num_nhc_iso;
  double dt_nhc,dti_nhc,wght;  /* Num: NHC respa time steps           */

  int *inhc_x,*inhc_y,*inhc_z; /* Map: Atm index -> to atm NHC;
                                  Lth: natm_tot                       */
  double *text_nhc;            /* Lst: T_ext of NHC                   */
                               /*   Lth: num_nhc                      */ 
  double *wdti,*wdti2,*wdti4,*wdti8,*wdti16;/* Lst: Yosh steps;  Lth:9*/
  double **mass_nhc,**gkt;     /* Lst: Mass,degs free of NHC's;       */
                               /* Lth: num_nhc x len_nhc              */

//--------------------------------------------------------------------------
// Default constructor/destructor

   MDTHERM_INFO(){
    my_type      = 0;
    natm_tot     = 0;  
    iextended_on = 0;
    num_nhc      = 0;   
    len_nhc      = 0;   
    nres_nhc     = 0;  
    nyosh_nhc    = 0; 
    therm_typ    = 0; 
    isokin_opt   = 0;
    num_nhc_iso  = 0;

    wdti   = ((double *)cmalloc(25*sizeof(double),"MDTHERM_INFO constr"))-1;
    wdti2  = ((double *)cmalloc(25*sizeof(double),"MDTHERM_INFO constr"))-1;
    wdti4  = ((double *)cmalloc(25*sizeof(double),"MDTHERM_INFO constr"))-1;
    wdti8  = ((double *)cmalloc(25*sizeof(double),"MDTHERM_INFO constr"))-1;
    wdti16 = ((double *)cmalloc(25*sizeof(double),"MDTHERM_INFO constrr"))-1;
   }
  ~MDTHERM_INFO(){
    cfree(&wdti[1],"MDTHERM_INFO destruct");
    cfree(&wdti2[1],"MDTHERM_INFO destruct");
    cfree(&wdti4[1],"MDTHERM_INFO destruct");
    cfree(&wdti8[1],"MDTHERM_INFO destruct");
    cfree(&wdti16[1],"MDTHERM_INFO destruct");
   }

//--------------------------------------------------------------------------
// Pack/Unpack

#ifdef PUP_ON
  void pup(PUP::er &p){
    // PUP ints
    p | iextended_on;
    p | my_type;
    p | natm_tot;
    p | num_nhc;   
    p | len_nhc;   
    p | nres_nhc;  
    p | nyosh_nhc; 
    p | therm_typ;   
    p | isokin_opt;
    p | num_nhc_iso;

    // PUP doubles 
    p | dt_nhc;
    p | dti_nhc;
    p | wght;
    // PUP arrays
    if(iextended_on > 0) {
      if(natm_tot>0){
        pup1d_int(p,&inhc_x,natm_tot);
        pup1d_int(p,&inhc_y,natm_tot);
        pup1d_int(p,&inhc_z,natm_tot);
      }//endif
      pup1d_dbl(p,&wdti,25);
      pup1d_dbl(p,&wdti2,25);
      pup1d_dbl(p,&wdti4,25);
      pup1d_dbl(p,&wdti8,25);
      pup1d_dbl(p,&wdti16,25);
      if(num_nhc>0 && len_nhc>0){
        pup1d_dbl(p,&text_nhc,num_nhc);
        pup2d_dbl(p,&mass_nhc,len_nhc,num_nhc,"mdtherm_info");
        pup2d_dbl(p,&gkt,len_nhc,num_nhc,"mdtherm_info");
      }//endif
    }// endif
#ifdef _PARALLEL_DEBUG_        
    if (p.isUnpacking())
     state_class_out ();
#endif
  } // pack/unpack
#endif

//--------------------------------------------------------------------------
// Print out the state of the class

  void state_class_out(){  
     int i,j;

     char fileName [255];
     sprintf (fileName, "%d_mdtherm_info.state", CkMyPe());
     FILE *fp;  fp = fopen(fileName,"w"); 


     fprintf(fp,"mdtherm_info: num_nhc %d\n",num_nhc);   
     fprintf(fp,"mdtherm_info: len_nhc %d\n",len_nhc);   
     fprintf(fp,"mdtherm_info: nres_nhc %d\n",nres_nhc);   
     fprintf(fp,"mdtherm_info: nyosh_nhc %d\n",nyosh_nhc);   
     fprintf(fp,"mdtherm_info: therm_typ %d\n",therm_typ);   
     fprintf(fp,"mdtherm_info: isokin_opt %d\n",isokin_opt);
     fprintf(fp,"mdtherm_info: num_nhc_iso %d\n",num_nhc_iso);

     fprintf(fp,"mdtherm_info: dt_nhc %.12g\n",dt_nhc);   
     fprintf(fp,"mdtherm_info: dti_nhc %.12g\n",dti_nhc);   
     fprintf(fp,"mdtherm_info: wght %.12g\n",wght);   

     if(iextended_on > 0){
      for(i=1;i<=natm_tot;i++){fprintf(fp,"mdtherm_info: inhc_x[%d] %d\n",
                               i,inhc_x[i]);}
      for(i=1;i<=natm_tot;i++){fprintf(fp,"mdtherm_info: inhc_y[%d] %d\n",
                               i,inhc_y[i]);}
      for(i=1;i<=natm_tot;i++){fprintf(fp,"mdtherm_info: inhc_z[%d] %d\n",
                               i,inhc_z[i]);}
      for(i=1;i<=num_nhc;i++){fprintf(fp,"mdtherm_info: text_nhc[%d] %g\n",
                               i,text_nhc[i]);}

      for(i=1;i<=9;i++){fprintf(fp,"mdtherm_info:   wdti[%d] %g\n",i,wdti[i]);}
      for(i=1;i<=9;i++){fprintf(fp,"mdtherm_info:  wdti2[%d] %g\n",i,wdti2[i]);}
      for(i=1;i<=9;i++){fprintf(fp,"mdtherm_info:  wdti4[%d] %g\n",i,wdti4[i]);}
      for(i=1;i<=9;i++){fprintf(fp,"mdtherm_info:  wdti8[%d] %g\n",i,wdti8[i]);}
      for(i=1;i<=9;i++){fprintf(fp,"mdtherm_info: wdti16[%d] %g\n",i,wdti16[i]);}

      for(i=1;i<=len_nhc;i++){
       for(j=1;j<=len_nhc;j++){
        fprintf(fp,"mdtherm_info: mass_nhc[%d][%d] %g\n",i,j,mass_nhc[i][j]);
       }//endfor : j
      }// endfor : i
      for(i=1;i<=len_nhc;i++){
       for(j=1; j<=len_nhc;j++){
         fprintf(fp,"mdtherm_info:  gkt[%d][%d] %g\n",i,j,gkt[i][j]);
       }//endfor : j
      }// endfor : i

     }//endif : iextended_on
     fclose(fp);
  }// end member function

//--------------------------------------------------------------------------
   }; // end class definition
//==========================================================================

#ifdef PUP_ON
PUPmarshall(MDTHERM_INFO);
#endif

#endif
//==========================================================================

