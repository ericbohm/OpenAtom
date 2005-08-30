//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                                                                          
//                            class_mdecor.h                                  
//                                                                          
//            Class definition for classical Ewald corrections              
//                                                                          
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

#ifndef _MDECOR_
#define _MDECOR_

class MDECOR{
 public:
  int num;                     /* Num: # ewald corrections            */
  int nsplin,nsplin_m2;        /* Num: spline coeffs                  */
  int nktot_res;

  double alp_ewd;              /* Num: Ewald convergence parameter    */
  double ecut,ecut_res;
  double rmin_spl;
  double rmax_spl;
  double dr_spl,dri_spl;

  int *j1,*j2;                 /* Lst: Indices of atms in ecors;
                                  Lth: num                            */
  double *cv0,*cdv0,*cdv0_res;           /* Lst: spline coefs */

//-----------------------------------------------------------------------------
// Default Constructor/destructor

   MDECOR(){
    num       = 0;       
    nsplin    = 0;
    nsplin_m2 = 0; 
    nktot_res = 0;
   }
  ~MDECOR(){}

//-----------------------------------------------------------------------------
// Pack/Unpack

#ifdef PUP_ON
  void pup(PUP::er &p){

  // PUP integers
   p | num;    
   p | nsplin;
   p | nsplin_m2;  
   p | nktot_res;
  // PUP doubles
   p | alp_ewd;  
   p | ecut;
   p | ecut_res;
   p | rmin_spl;
   p | rmax_spl;
   p | dr_spl;
   p | dri_spl;
  // PUP Arrays
   if(num>0){
    pup1d_int(p,&j1,num);
    pup1d_int(p,&j2,num);
    pup1d_dbl(p,&cv0,nsplin);
    pup1d_dbl(p,&cdv0,nsplin);
    if(nktot_res>0){pup1d_dbl(p,&cdv0_res,nsplin);}
   }/*endif*/
#ifdef _PARALLEL_DEBUG_        
    if (p.isUnpacking())
     state_class_out ();
#endif 
  }/* end pack/unpack */
#endif

//-----------------------------------------------------------------------------
// Print out state of class

  void state_class_out(){
     int i;
     char fileName [255];
     sprintf (fileName, "%d_mdecor.state", CkMyPe());
     FILE *fp;  fp = fopen(fileName,"w");

  // Print integers
    fprintf(fp,"mdecor:  nsplin %d\n",nsplin);
    fprintf(fp,"mdecor:  nsplin_m2 %d\n",nsplin_m2);  
    fprintf(fp,"mdecor:  num %d\n",num);    
    fprintf(fp,"mdecor:  nktot_res %d\n",nktot_res);
  // Print doubles
    fprintf(fp,"mdecor:  alp_ewd %.12g\n",alp_ewd);  
    fprintf(fp,"mdecor:  ecut %.12g\n",ecut);
    fprintf(fp,"mdecor:  ecut_res %.12g\n",ecut_res);
    fprintf(fp,"mdecor:  rmin_spl %.12g\n",rmin_spl);
    fprintf(fp,"mdecor:  rmax_spl %.12g\n",rmax_spl);
    fprintf(fp,"mdecor:  dr_spl %.12g\n",dr_spl);
    fprintf(fp,"mdecor:  dri_spl %.12g\n",dri_spl);
   // Print Integer Arrays 
    for(i=1;i<=num;i++){fprintf(fp,"mdonfo:  j1[%d] %d\n",i,j1[i]);}
    for(i=1;i<=num;i++){fprintf(fp,"mdonfo:  j2[%d] %d\n",i,j2[i]);}

    fclose(fp);
  }/* end member function */
//-----------------------------------------------------------------------------

}; /* End class definition */

#ifdef PUP_ON
PUPmarshall(MDECOR);
#endif

#endif
//==========================================================================
