/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                          class_mdgrp_bond_con.h                          */
/*                                                                          */
/*                Class definition for group bond constraints               */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#ifndef _MDGRP_BOND_CON_
#define _MDGRP_BOND_CON_

class MDGRP_BOND_CON {
 public:
  int num_33;                  /* Num: # 3 atom - 3bond units     */
  int num_21;                  /* Num: # 3 atom - 3bond units     */
  int num_43;                  /* Num: # 3 atom - 3bond units     */
  int ntyp_33;                 /* Num: 33 types                   */
  int num_23;                  /* Num: # 3 atom - 2bond units     */
  int ntyp_23,ntyp_21,ntyp_43; /* Num: 33 types                   */
  int num_46;                  /* Lst: 4 atom - 6bond units lists */
  int ntyp_46;                 /* Num: 33 types                   */
  int max_iter;                /* Num: maximum iterations             */

  double tol;                  /* Num: bond_con tol               */

  int *j1_33,*j2_33,*j3_33;    /* Lst: 3 atom - 3bond units lists */
  int *j1_21,*j1_43,*j2_21;    /* Lst: 3 atom - 3bond units lists */
  int *j2_43,*j3_43,*j4_43;    /* Lst: 3 atom - 3bond units lists */
  int *jtyp_33,*jtyp_21,*jtyp_43; /* Map: Unit to unit type          */
  int *j1_23,*j2_23,*j3_23;    /* Lst: 3 atom - 3bond units lists */
  int *jtyp_23;                /* Map: Unit to unit type          */
  int *j1_46,*j2_46,*j3_46,*j4_46;
                               /* Lst: 4 atom - 6bond units lists */
  int *jtyp_46;                /* Map: Unit to unit type          */

  double **al_33,**eq_33;      /* Lst: Multipliers and eqbonds    */
  double **al_23,**eq_23;      /* Lst: Multipliers and eqbonds    */
  double **al_21,**al_43;      /* Lst: Multipliers and eqbonds    */
  double **eq_21,**eq_43;      /* Lst: Multipliers and eqbonds    */
  double **al_46,**eq_46;      /* Lst: Multipliers and eqbonds    */

// Constructor/Destructor

  MDGRP_BOND_CON(){
    num_33 = 0;                 
    num_21 = 0;                 
    num_43 = 0;                 
    ntyp_33 = 0;                
    num_23 = 0;                 
    ntyp_23 = 0;
    ntyp_21 = 0;
    ntyp_43 = 0;
    num_46 = 0;                 
    ntyp_46 = 0;                
    max_iter = 0;               
  }
  ~MDGRP_BOND_CON(){}

#ifdef PUP_ON
  void pup(PUP::er &p){
    // PUP the ints 
    p | num_33;  
    p | num_21;  
    p | num_43;  
    p | ntyp_33; 
    p | num_23;  
    p | ntyp_23;
    p | ntyp_21;
    p | ntyp_43;
    p | num_46;                 
    p | ntyp_46;                
    p | max_iter;               

    // PUP the doubles
    p | tol;

    // PUP the Integer arrays

    if(num_33>0){
     pup1d_int(p,&j1_33,num_33);
     pup1d_int(p,&j2_33,num_33);
     pup1d_int(p,&j3_33,num_33);
     pup1d_int(p,&jtyp_33,num_33);
    }
    if(num_21>0){
     pup1d_int(p,&j1_21,num_21);
     pup1d_int(p,&j2_21,num_21);
     pup1d_int(p,&jtyp_21,num_21);
    }
    if(num_43>0){
     pup1d_int(p,&j1_43,num_43);
     pup1d_int(p,&j2_43,num_43);
     pup1d_int(p,&j3_43,num_43);
     pup1d_int(p,&j4_43,num_43);
     pup1d_int(p,&jtyp_43,num_43);
    }
    if(num_23>0){
     pup1d_int(p,&j1_23,num_23);
     pup1d_int(p,&j2_23,num_23);
     pup1d_int(p,&j3_23,num_23);
     pup1d_int(p,&jtyp_23,num_23);
    }
    if(num_46>0){
     pup1d_int(p,&j1_46,num_46);
     pup1d_int(p,&j2_46,num_46);
     pup1d_int(p,&j3_46,num_46);
     pup1d_int(p,&j4_46,num_46);
     pup1d_int(p,&jtyp_46,num_46);
    }

     //PUP the Double arrays

    if(ntyp_33>0){
     pup2d_dbl(p,&al_33,3,ntyp_33,"mdgrpbondcon");
     pup2d_dbl(p,&eq_33,3,ntyp_33,"mdgrpbondcon");
    }
    if(ntyp_23>0){
     pup2d_dbl(p,&al_23,2,ntyp_23,"mdgrpbondcon");
     pup2d_dbl(p,&eq_23,2,ntyp_23,"mdgrpbondcon");
    }
    if(ntyp_21>0){
     pup2d_dbl(p,&al_21,1,ntyp_21,"mdgrpbondcon");
     pup2d_dbl(p,&eq_21,1,ntyp_21,"mdgrpbondcon");
    }
    if(ntyp_43>0){
     pup2d_dbl(p,&al_43,3,ntyp_43,"mdgrpbondcon");
     pup2d_dbl(p,&eq_43,3,ntyp_43,"mdgrpbondcon");
    }
    if(ntyp_46>0){
     pup2d_dbl(p,&al_46,6,ntyp_46,"mdgrpbondcon");
     pup2d_dbl(p,&eq_46,6,ntyp_46,"mdgrpbondcon");
    }
#ifdef _PARALLEL_DEBUG_        
    if (p.isUnpacking())
     state_class_out ();
#endif 
  }/* Pack/Unpack */
#endif  

  void state_class_out(){
   int i;
   char fileName [255];
   sprintf (fileName, "%d_mdgrp_bond_con.state", CkMyPe());
   FILE *fp;  fp = fopen(fileName,"w");

    fprintf(fp,"num_33 %d\n",num_33);  
    fprintf(fp,"num_21 %d\n",num_21);  
    fprintf(fp,"num_43 %d\n",num_43);  
    fprintf(fp,"ntyp_33 %d\n",ntyp_33); 
    fprintf(fp,"num_23 %d\n",num_23);  
    fprintf(fp,"ntyp_23 %d\n",ntyp_23);
    fprintf(fp,"ntyp_21 %d\n",ntyp_21);
    fprintf(fp,"ntyp_43 %d\n",ntyp_43);
    fprintf(fp,"num_46 %d\n",num_46);                 
    fprintf(fp,"ntyp_46 %d\n",ntyp_46);                
    fprintf(fp,"max_iter %d\n",max_iter);               

    // doubles
    fprintf(fp,"tol %g\n",tol);
    // Integer arrays
    for(i=1;i<=num_33;i++){fprintf(fp,"j1_33[%d] %d\n",i,j1_33[i]);}
    for(i=1;i<=num_33;i++){fprintf(fp,"j2_33[%d] %d\n",i,j2_33[i]);}
    for(i=1;i<=num_33;i++){fprintf(fp,"j3_33[%d] %d\n",i,j3_33[i]);}

    for(i=1;i<=num_21;i++){fprintf(fp,"j1_21[%d] %d\n",i,j1_21[i]);}
    for(i=1;i<=num_21;i++){fprintf(fp,"j2_21[%d] %d\n",i,j2_21[i]);}

    for(i=1;i<=num_43;i++){fprintf(fp,"j1_43[%d] %d\n",i,j1_43[i]);}
    for(i=1;i<=num_43;i++){fprintf(fp,"j2_43[%d] %d\n",i,j2_43[i]);}
    for(i=1;i<=num_43;i++){fprintf(fp,"j3_43[%d] %d\n",i,j3_43[i]);}
    for(i=1;i<=num_43;i++){fprintf(fp,"j4_43[%d] %d\n",i,j4_43[i]);}

    for(i=1;i<=num_23;i++){fprintf(fp,"j1_23[%d] %d\n",i,j1_23[i]);}
    for(i=1;i<=num_23;i++){fprintf(fp,"j2_23[%d] %d\n",i,j2_23[i]);}
    for(i=1;i<=num_23;i++){fprintf(fp,"j3_23[%d] %d\n",i,j3_23[i]);}

    for(i=1;i<=num_46;i++){fprintf(fp,"j1_46[%d] %d\n",i,j1_46[i]);}
    for(i=1;i<=num_46;i++){fprintf(fp,"j2_46[%d] %d\n",i,j2_46[i]);}
    for(i=1;i<=num_46;i++){fprintf(fp,"j3_46[%d] %d\n",i,j3_46[i]);}
    for(i=1;i<=num_46;i++){fprintf(fp,"j4_46[%d] %d\n",i,j4_46[i]);}

    for(i=1;i<=ntyp_33;i++){fprintf(fp,"jtyp_33[%d] %d\n",i,jtyp_33[i]);}
    for(i=1;i<=ntyp_21;i++){fprintf(fp,"jtyp_21[%d] %d\n",i,jtyp_21[i]);}
    for(i=1;i<=ntyp_43;i++){fprintf(fp,"jtyp_43[%d] %d\n",i,jtyp_43[i]);}
    for(i=1;i<=ntyp_23;i++){fprintf(fp,"jtyp_23[%d] %d\n",i,jtyp_23[i]);}
    for(i=1;i<=ntyp_46;i++){fprintf(fp,"jtyp_46[%d] %d\n",i,jtyp_46[i]);}

   fclose(fp);

  }// end routine


}; /* end class definition */

#ifdef PUP_ON
PUPmarshall(MDGRP_BOND_CON);
#endif


#endif

