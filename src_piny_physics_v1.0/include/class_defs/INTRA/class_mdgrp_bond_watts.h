//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                                                                          
//                          class_mdgrp_bond_watts.h                        
//                                                                          
//                Class definition for Watts type group bonds               
//                                                                          
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

#ifndef _MDGRP_BOND_WATTS_
#define _MDGRP_BOND_WATTS_

class MDGRP_BOND_WATTS {
 public:
  int num_33;                  // Num: # 3 atom - 3bond units     
  int num_33_tot;              // Num: # 3 atom - 3bond units     
  int ntyp_33;                 // Num: 33 types                   

  int ngrp_33;
  int ngrp_typ_33;

  int *j1_33,*j2_33,*j3_33;    // Lst: 3 atom - 3bond units lists 
  int *jtyp_33;                 // Map: Unit to unit type          

  double *cos_thet0_2,*sin_thet0_2;
  double **eq_33;      // Lst: Multipliers and eqbonds    
  double **c_0_33,**c_1_33,**c_2_33,**c_3_33;
  double **c_4_33,**c_5_33,**c_6_33;
  double **dc_0_33,**dc_1_33,**dc_2_33,**dc_3_33;
  double **dc_4_33,**dc_5_33,**dc_6_33;

//-------------------------------------------------------------------------
// Constructor/Destructor

   MDGRP_BOND_WATTS(){
    num_33 = 0;                 
    num_33_tot = 0;             
    ntyp_33 = 0;                
    ngrp_33 = 0;
    ngrp_typ_33 = 0;
   }
  ~MDGRP_BOND_WATTS(){}

//-------------------------------------------------------------------------
#ifdef PUP_ON
  void pup(PUP::er &p){
    // PUP the ints 

    p | num_33;          
    p | num_33_tot; 
    p | ntyp_33;    
    p | ngrp_33;
    p | ngrp_typ_33;

    // PUP the Integer Arrays

    if(ngrp_33>0){
     pup1d_int(p,&j1_33,ngrp_33);
     pup1d_int(p,&j2_33,ngrp_33);
     pup1d_int(p,&j3_33,ngrp_33);
     pup1d_int(p,&jtyp_33,ngrp_33);

     pup1d_dbl(p,&cos_thet0_2,ngrp_typ_33);
     pup1d_dbl(p,&sin_thet0_2,ngrp_typ_33);
     pup2d_dbl(p,&eq_33, 3,ngrp_typ_33,"mdgrpbondwatts");
     pup2d_dbl(p,&c_0_33,3,ngrp_typ_33,"mdgrpbondwatts");
     pup2d_dbl(p,&c_1_33,3,ngrp_typ_33,"mdgrpbondwatts");
     pup2d_dbl(p,&c_2_33,3,ngrp_typ_33,"mdgrpbondwatts");
     pup2d_dbl(p,&c_3_33,3,ngrp_typ_33,"mdgrpbondwatts");
     pup2d_dbl(p,&c_4_33,3,ngrp_typ_33,"mdgrpbondwatts");
     pup2d_dbl(p,&c_5_33,3,ngrp_typ_33,"mdgrpbondwatts");
     pup2d_dbl(p,&c_6_33,3,ngrp_typ_33,"mdgrpbondwatts");
     pup2d_dbl(p,&dc_0_33,3,ngrp_typ_33,"mdgrpbondwatts");
     pup2d_dbl(p,&dc_1_33,3,ngrp_typ_33,"mdgrpbondwatts");
     pup2d_dbl(p,&dc_2_33,3,ngrp_typ_33,"mdgrpbondwatts");
     pup2d_dbl(p,&dc_3_33,3,ngrp_typ_33,"mdgrpbondwatts");
     pup2d_dbl(p,&dc_4_33,3,ngrp_typ_33,"mdgrpbondwatts");
     pup2d_dbl(p,&dc_5_33,3,ngrp_typ_33,"mdgrpbondwatts");
     pup2d_dbl(p,&dc_6_33,3,ngrp_typ_33,"mdgrpbondwatts");
    }//endif
#ifdef _PARALLEL_DEBUG_        
    if (p.isUnpacking())
     state_class_out ();
#endif 
  }// end Pack/Unpack 
#endif

//-------------------------------------------------------------------------
  void state_class_out(){
   char fileName [255];
   sprintf (fileName, "%d_mdgrp_bond_watts.state", CkMyPe());
   FILE *fp;  fp = fopen(fileName,"w");

    fprintf(fp,"num_33 %d\n",num_33);  
    fprintf(fp,"ntyp_33 %d\n",ntyp_33); 
    fprintf(fp,"ngrp_33 %d\n",ngrp_33);
    fprintf(fp,"ngrp_typ_33 %d\n",ngrp_typ_33);
    // Integer Arrays
    for(int i=1;i<=ngrp_33;i++){fprintf(fp,"j1_33[%d] %d\n",i,j1_33[i]);}
    for(int i=1;i<=ngrp_33;i++){fprintf(fp,"j2_33[%d] %d\n",i,j2_33[i]);}
    for(int i=1;i<=ngrp_33;i++){fprintf(fp,"j3_33[%d] %d\n",i,j3_33[i]);}
    for(int i=1;i<=ngrp_typ_33;i++){fprintf(fp,"jtyp_33[%d] %d\n",
                                     i,jtyp_33[i]);}
   fclose(fp);
  }// end routine

//-------------------------------------------------------------------------
  }; // end class definition 
//==========================================================================


#ifdef PUP_ON
PUPmarshall(MDGRP_BOND_WATTS);
#endif


#endif
//==========================================================================
