/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                          class_mdbond.h                                  */
/*                                                                          */
/*                Class definition for classical bonds                      */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#ifndef _MDBOND_
#define _MDBOND_

class MDBOND {
 public:
  int npow;                    /* Num: # pow series bonds             */
  int ntyp_pow;                /* Num: # pow series bond typs         */
  int ncon;                    /* Num: # constrained bonds            */
  int ntyp_con;                /* Num: # constrained bond typs        */

  int ncon_tot;

  int *j1_pow,*j2_pow;         /* Lst: Indices of atms in pow bond;
                                  Lth: npow                           */
  int *jtyp_pow;               /* Map: index of bond -> bond type;    
                                  Lth: npow                           */
  int *j1_con,*j2_con;         /* Lst: Indices of atms in cons bond;   
                                  Lth: ncon                           */
  int *jtyp_con;               /* Map: index of bond -> bond type;
                                  Lth: ncon                           */
  double *eq_pow;              /* Lst: List of eq. bond lgths;        
                                  Lth: ntyp_pow                       */
  double *c_0,*c_1,*c_2,*c_3,*c_4,*c_5,*c_6;            
  double *dc_0,*dc_1,*dc_2,*dc_3,*dc_4,*dc_5,*dc_6;        
                               /* Lst: Bond power series coefficients
                                  Lth: ntyp_pow                       */
  double *al_con;              /* Lst: List of cons bond multipliers;    
                                  Lth: ncon                           */
  double *eq_con;              /* Lst: List of eq. bond lgths;  
                                  Lth: ntyp_con                       */

//============================================================================
// Default Constructor/Destructor

   MDBOND(){
    npow     = 0;            
    ntyp_pow = 0;        
    ncon     = 0;            
    ntyp_con = 0;        
    ncon_tot = 0;
   }
  ~MDBOND(){}

//=============================================================================
// Pack/Unpac

#ifdef PUP_ON
  void pup(PUP::er &p){

   // Pupping ints

   p | npow;          
   p | ntyp_pow;      
   p | ncon;          
   p | ntyp_con;      
   p | ncon_tot;


   // Integer Arrays 

   if(npow > 0){
     pup1d_int(p,&j1_pow,npow);
     pup1d_int(p,&j2_pow,npow);
     pup1d_int(p,&jtyp_pow,npow);
   }/* endif */
                              
   if(ncon > 0){
     pup1d_int(p,&j1_con,ncon);
     pup1d_int(p,&j2_con,ncon);
     pup1d_int(p,&jtyp_con,ncon);
     pup1d_dbl(p,&al_con,ncon);     
   }/* endif */
                              
   // Double Arrays 


   if(ntyp_pow > 0){
     pup1d_dbl(p,&eq_pow,ntyp_pow);    
                                   
     pup1d_dbl(p,&c_0,ntyp_pow);
     pup1d_dbl(p,&c_1,ntyp_pow);
     pup1d_dbl(p,&c_2,ntyp_pow);
     pup1d_dbl(p,&c_3,ntyp_pow);
     pup1d_dbl(p,&c_4,ntyp_pow);
     pup1d_dbl(p,&c_5,ntyp_pow);
     pup1d_dbl(p,&c_6,ntyp_pow);                
     pup1d_dbl(p,&dc_0,ntyp_pow);
     pup1d_dbl(p,&dc_1,ntyp_pow);
     pup1d_dbl(p,&dc_2,ntyp_pow);
     pup1d_dbl(p,&dc_3,ntyp_pow);
     pup1d_dbl(p,&dc_4,ntyp_pow);
     pup1d_dbl(p,&dc_5,ntyp_pow);
     pup1d_dbl(p,&dc_6,ntyp_pow);          
   }/* endif */

   if(ntyp_con > 0){
      pup1d_dbl(p,&eq_con,ntyp_con);
   }/* endif */
#ifdef _PARALLEL_DEBUG_        
    if (p.isUnpacking())
     state_class_out ();
#endif 
  } /* Pack/Unpack */
#endif

//=====================================================================================
// Print out state of class

  void state_class_out(){
    
     char fileName [255];
     sprintf (fileName, "%d_mdbond.state", CkMyPe());
     FILE *fp;  fp = fopen(fileName,"w");

     fprintf(fp,"mdbond:  npow %d\n",npow);
     fprintf(fp,"mdbond:  ntyp_pow %d\n",ntyp_pow);
     fprintf(fp,"mdbond:  ncon %d\n",ncon);
     fprintf(fp,"mdbond:  ntyp_con %d\n",ntyp_con);
     fprintf(fp,"mdbond:  ncon_tot %d\n",ncon_tot);

   // Print Integer Arrays 

     for(int i=1;i<=npow;i++){fprintf(fp,"mdbond:  j1_pow[%d] %d\n",i,j1_pow[i]);}
     for(int i=1;i<=npow;i++){fprintf(fp,"mdbond:  j2_pow[%d] %d\n",i,j2_pow[i]);}

     for(int i=1;i<=npow;i++){fprintf(fp,"mdbond:  jtyp_pow[%d] %d\n",i,jtyp_pow[i]);}
                              
     for(int i=1;i<=ncon;i++){fprintf(fp,"mdbond:  j1_con[%d] %d\n",i,j1_con[i]);}
     for(int i=1;i<=ncon;i++){fprintf(fp,"mdbond:  j2_con[%d] %d\n",i,j2_con[i]);}
                            
     for(int i=1;i<=ncon;i++){fprintf(fp,"mdbond:  jtyp_con[%d] %d\n",i,jtyp_con[i]);}
                              
   // Print Double Arrays 

     for(int i=1;i<=ntyp_pow;i++){fprintf(fp,"mdbond:  eq_pow[%d] %.12g\n",i,eq_pow[i]);}
                                   
     for(int i=1;i<=ntyp_pow;i++){fprintf(fp,"mdbond:  c_0[%d] %.12g\n",i,c_0[i]);}
     for(int i=1;i<=ntyp_pow;i++){fprintf(fp,"mdbond:  c_1[%d] %.12g\n",i,c_1[i]);}
     for(int i=1;i<=ntyp_pow;i++){fprintf(fp,"mdbond:  c_2[%d] %.12g\n",i,c_2[i]);}
     for(int i=1;i<=ntyp_pow;i++){fprintf(fp,"mdbond:  c_3[%d] %.12g\n",i,c_3[i]);}
     for(int i=1;i<=ntyp_pow;i++){fprintf(fp,"mdbond:  c_4[%d] %.12g\n",i,c_4[i]);}
     for(int i=1;i<=ntyp_pow;i++){fprintf(fp,"mdbond:  c_5[%d] %.12g\n",i,c_5[i]);}
     for(int i=1;i<=ntyp_pow;i++){fprintf(fp,"mdbond:  c_6[%d] %.12g\n",i,c_6[i]);}
     for(int i=1;i<=ntyp_pow;i++){fprintf(fp,"mdbond:  dc_0[%d] %.12g\n",i,dc_0[i]);}
     for(int i=1;i<=ntyp_pow;i++){fprintf(fp,"mdbond:  dc_1[%d] %.12g\n",i,dc_1[i]);}
     for(int i=1;i<=ntyp_pow;i++){fprintf(fp,"mdbond:  dc_2[%d] %.12g\n",i,dc_2[i]);}
     for(int i=1;i<=ntyp_pow;i++){fprintf(fp,"mdbond:  dc_3[%d] %.12g\n",i,dc_3[i]);}
     for(int i=1;i<=ntyp_pow;i++){fprintf(fp,"mdbond:  dc_4[%d] %.12g\n",i,dc_4[i]);}
     for(int i=1;i<=ntyp_pow;i++){fprintf(fp,"mdbond:  dc_5[%d] %.12g\n",i,dc_5[i]);}
     for(int i=1;i<=ntyp_pow;i++){fprintf(fp,"mdbond:  dc_6[%d] %.12g\n",i,dc_6[i]);}

     for(int i=1;i<=ncon;i++){fprintf(fp,"mdbond:  al_con[%d] %.12g\n",i,al_con[i]);}
     for(int i=1;i<=ntyp_con;i++){fprintf(fp,"mdbond:  eq_con[%d] %.12g\n",i,eq_con[i]);}

     fclose(fp);
  }/* end member function */

};/* end class def for BOND */

#ifdef PUP_ON
PUPmarshall(MDBOND);
#endif

#endif
