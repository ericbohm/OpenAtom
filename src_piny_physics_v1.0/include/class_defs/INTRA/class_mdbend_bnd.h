/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                          class_mdbend_bnd.h                              */
/*                                                                          */
/*            Class definition for classical bend-bonds (Urey-Bradley)       */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#ifndef _MDBEND_BND_
#define _MDBEND_BND_

class MDBEND_BND {
  public:
    int num;                     /* Num: # pow series bonds             */
    int ntyp;                    /* Num: # pow series bond typs         */
    int nbend_bnd;    
    int nbend_bnd_typ;

    int *j1,*j2,*j3;             /* Lst: Indices of atms in pow bond;
Lth: num                            */
    int *jtyp;                   /* Map: index of bond -> bond type;    
Lth: num                            */
    double *eq_bond;             /* Lst: List of eq. bond lgths;        
Lth: ntyp                           */
    double *eq_bend;             /* Lst: List of eq. bond lgths;        
Lth: ntyp                           */
    double *cbond_0,*cbond_1,*cbond_2,*cbond_3,*cbond_4,*cbond_5,*cbond_6;   
    double *dcbond_0,*dcbond_1,*dcbond_2,*dcbond_3,*dcbond_4,*dcbond_5,*dcbond_6;
    /* Lst: Bond power series coefficients
Lth: ntyp                           */
    double *cbend_0,*cbend_1,*cbend_2,*cbend_3,*cbend_4,*cbend_5,*cbend_6;
    double *sbend_0,*sbend_1,*sbend_2,*sbend_3,*sbend_4,*sbend_5,*sbend_6;
    double *dcbend_0,*dcbend_1,*dcbend_2,*dcbend_3,*dcbend_4,*dcbend_5,*dcbend_6;
    double *dsbend_0,*dsbend_1,*dsbend_2,*dsbend_3,*dsbend_4,*dsbend_5,*dsbend_6;
    /* Lst: Bend power series coefficients
Lth: ntyp                           */

    //============================================================================
    // Default Constructor/Destructor

    MDBEND_BND(){
      num       = 0;         
      ntyp      = 0;        
      nbend_bnd = 0;    
      nbend_bnd_typ = 0;
    }
    ~MDBEND_BND(){}

    //=============================================================================
    // Pack/Unpack

#ifdef PUP_ON
    void pup(PUP::er &p){

      // Pupping ints

      p | num;
      p | ntyp;
      p | nbend_bnd;
      p | nbend_bnd_typ;

      // PUP Arrays

      if(num>0){
        pup1d_int(p,&j1,num);
        pup1d_int(p,&j2,num);
        pup1d_int(p,&j3,num);
        pup1d_int(p,&jtyp,num);
      }/* endif */
      if(ntyp > 0){
        pup1d_dbl(p,&eq_bond,ntyp);
        pup1d_dbl(p,&eq_bend,ntyp);

        pup1d_dbl(p,&cbond_0,ntyp);
        pup1d_dbl(p,&cbond_1,ntyp);
        pup1d_dbl(p,&cbond_2,ntyp);
        pup1d_dbl(p,&cbond_3,ntyp);
        pup1d_dbl(p,&cbond_4,ntyp);
        pup1d_dbl(p,&cbond_5,ntyp);
        pup1d_dbl(p,&cbond_6,ntyp);

        pup1d_dbl(p,&dcbond_0,ntyp);
        pup1d_dbl(p,&dcbond_1,ntyp);
        pup1d_dbl(p,&dcbond_2,ntyp);
        pup1d_dbl(p,&dcbond_3,ntyp);
        pup1d_dbl(p,&dcbond_4,ntyp);
        pup1d_dbl(p,&dcbond_5,ntyp);
        pup1d_dbl(p,&dcbond_6,ntyp);

        pup1d_dbl(p,&cbend_0,ntyp);
        pup1d_dbl(p,&cbend_1,ntyp);
        pup1d_dbl(p,&cbend_2,ntyp);
        pup1d_dbl(p,&cbend_3,ntyp);
        pup1d_dbl(p,&cbend_4,ntyp);
        pup1d_dbl(p,&cbend_5,ntyp);
        pup1d_dbl(p,&cbend_6,ntyp);

        pup1d_dbl(p,&sbend_0,ntyp);
        pup1d_dbl(p,&sbend_1,ntyp);
        pup1d_dbl(p,&sbend_2,ntyp);
        pup1d_dbl(p,&sbend_3,ntyp);
        pup1d_dbl(p,&sbend_4,ntyp);
        pup1d_dbl(p,&sbend_5,ntyp);
        pup1d_dbl(p,&sbend_6,ntyp);

        pup1d_dbl(p,&dcbend_0,ntyp);
        pup1d_dbl(p,&dcbend_1,ntyp);
        pup1d_dbl(p,&dcbend_2,ntyp);
        pup1d_dbl(p,&dcbend_3,ntyp);
        pup1d_dbl(p,&dcbend_4,ntyp);
        pup1d_dbl(p,&dcbend_5,ntyp);
        pup1d_dbl(p,&dcbend_6,ntyp);

        pup1d_dbl(p,&dsbend_0,ntyp);
        pup1d_dbl(p,&dsbend_1,ntyp);
        pup1d_dbl(p,&dsbend_2,ntyp);
        pup1d_dbl(p,&dsbend_3,ntyp);
        pup1d_dbl(p,&dsbend_4,ntyp);
        pup1d_dbl(p,&dsbend_5,ntyp);
        pup1d_dbl(p,&dsbend_6,ntyp);
      }/* endif */
#ifdef _PARALLEL_DEBUG_        
      if (p.isUnpacking())
        state_class_out ();
#endif 
    }/* end Pack/Unpack */
#endif

    //=====================================================================================
    // Print out state of class

    void state_class_out(){

      char fileName [255];
      sprintf (fileName, "%d_mdbend_bnd.state", CkMyPe());
      FILE *fp;  fp = fopen(fileName,"w");

      fprintf(fp,"mdbend_bnd:  num %d\n",num);
      fprintf(fp,"mdbend_bnd:  ntyp %d\n",ntyp);
      fprintf(fp,"mdbend_bnd:  nbend_bnd %d\n",nbend_bnd);
      fprintf(fp,"mdbend_bnd:  nbend_bnd_typ %d\n",nbend_bnd_typ);

      // Print Integer Arrays 
      int i;
      for(i=1;i<=num;i++){fprintf(fp,"mdbend_bnd:  j1[%d] %d\n",i,j1[i]);}
      for(i=1;i<=num;i++){fprintf(fp,"mdbend_bnd:  j2[%d] %d\n",i,j2[i]);}
      for(i=1;i<=num;i++){fprintf(fp,"mdbend_bnd:  j3[%d] %d\n",i,j3[i]);}

      for(i=1;i<=num;i++){fprintf(fp,"mdbend_bnd:  jtyp[%d] %d\n",i,jtyp[i]);}

      // Print Double Arrays 

      for(i=1;i<=ntyp;i++){fprintf(fp,"mdbend_bnd:  eq_bond[%d] %.12g\n",i,eq_bond[i]);}
      for(i=1;i<=ntyp;i++){fprintf(fp,"mdbend_bnd:  eq_bend[%d] %.12g\n",i,eq_bend[i]);}

      for(i=1;i<=ntyp;i++){fprintf(fp,"mdbend_bnd:  cbond_0[%d] %.12g\n",i,cbond_0[i]);}
      for(i=1;i<=ntyp;i++){fprintf(fp,"mdbend_bnd:  cbond_1[%d] %.12g\n",i,cbond_1[i]);}
      for(i=1;i<=ntyp;i++){fprintf(fp,"mdbend_bnd:  cbond_2[%d] %.12g\n",i,cbond_2[i]);}
      for(i=1;i<=ntyp;i++){fprintf(fp,"mdbend_bnd:  cbond_3[%d] %.12g\n",i,cbond_3[i]);}
      for(i=1;i<=ntyp;i++){fprintf(fp,"mdbend_bnd:  cbond_4[%d] %.12g\n",i,cbond_4[i]);}
      for(i=1;i<=ntyp;i++){fprintf(fp,"mdbend_bnd:  cbond_5[%d] %.12g\n",i,cbond_5[i]);}
      for(i=1;i<=ntyp;i++){fprintf(fp,"mdbend_bnd:  cbond_6[%d] %.12g\n",i,cbond_6[i]);}

      for(i=1;i<=ntyp;i++){fprintf(fp,"mdbend_bnd:  cbend_0[%d] %.12g\n",i,cbend_0[i]);}
      for(i=1;i<=ntyp;i++){fprintf(fp,"mdbend_bnd:  cbend_1[%d] %.12g\n",i,cbend_1[i]);}
      for(i=1;i<=ntyp;i++){fprintf(fp,"mdbend_bnd:  cbend_2[%d] %.12g\n",i,cbend_2[i]);}
      for(i=1;i<=ntyp;i++){fprintf(fp,"mdbend_bnd:  cbend_3[%d] %.12g\n",i,cbend_3[i]);}
      for(i=1;i<=ntyp;i++){fprintf(fp,"mdbend_bnd:  cbend_4[%d] %.12g\n",i,cbend_4[i]);}
      for(i=1;i<=ntyp;i++){fprintf(fp,"mdbend_bnd:  cbend_5[%d] %.12g\n",i,cbend_5[i]);}
      for(i=1;i<=ntyp;i++){fprintf(fp,"mdbend_bnd:  cbend_6[%d] %.12g\n",i,cbend_6[i]);}

      for(i=1;i<=ntyp;i++){fprintf(fp,"mdbend_bnd:  sbend_0[%d] %.12g\n",i,sbend_0[i]);}
      for(i=1;i<=ntyp;i++){fprintf(fp,"mdbend_bnd:  sbend_1[%d] %.12g\n",i,sbend_1[i]);}
      for(i=1;i<=ntyp;i++){fprintf(fp,"mdbend_bnd:  sbend_2[%d] %.12g\n",i,sbend_2[i]);}
      for(i=1;i<=ntyp;i++){fprintf(fp,"mdbend_bnd:  sbend_3[%d] %.12g\n",i,sbend_3[i]);}
      for(i=1;i<=ntyp;i++){fprintf(fp,"mdbend_bnd:  sbend_4[%d] %.12g\n",i,sbend_4[i]);}
      for(i=1;i<=ntyp;i++){fprintf(fp,"mdbend_bnd:  sbend_5[%d] %.12g\n",i,sbend_5[i]);}
      for(i=1;i<=ntyp;i++){fprintf(fp,"mdbend_bnd:  sbend_6[%d] %.12g\n",i,sbend_6[i]);}

      for(i=1;i<=ntyp;i++){fprintf(fp,"mdbend_bnd:  dcbond_0[%d] %.12g\n",i,dcbond_0[i]);}
      for(i=1;i<=ntyp;i++){fprintf(fp,"mdbend_bnd:  dcbond_1[%d] %.12g\n",i,dcbond_1[i]);}
      for(i=1;i<=ntyp;i++){fprintf(fp,"mdbend_bnd:  dcbond_2[%d] %.12g\n",i,dcbond_2[i]);}
      for(i=1;i<=ntyp;i++){fprintf(fp,"mdbend_bnd:  dcbond_3[%d] %.12g\n",i,dcbond_3[i]);}
      for(i=1;i<=ntyp;i++){fprintf(fp,"mdbend_bnd:  dcbond_4[%d] %.12g\n",i,dcbond_4[i]);}
      for(i=1;i<=ntyp;i++){fprintf(fp,"mdbend_bnd:  dcbond_5[%d] %.12g\n",i,dcbond_5[i]);}
      for(i=1;i<=ntyp;i++){fprintf(fp,"mdbend_bnd:  dcbond_6[%d] %.12g\n",i,dcbond_6[i]);}

      for(i=1;i<=ntyp;i++){fprintf(fp,"mdbend_bnd:  dcbend_0[%d] %.12g\n",i,dcbend_0[i]);}
      for(i=1;i<=ntyp;i++){fprintf(fp,"mdbend_bnd:  dcbend_1[%d] %.12g\n",i,dcbend_1[i]);}
      for(i=1;i<=ntyp;i++){fprintf(fp,"mdbend_bnd:  dcbend_2[%d] %.12g\n",i,dcbend_2[i]);}
      for(i=1;i<=ntyp;i++){fprintf(fp,"mdbend_bnd:  dcbend_3[%d] %.12g\n",i,dcbend_3[i]);}
      for(i=1;i<=ntyp;i++){fprintf(fp,"mdbend_bnd:  dcbend_4[%d] %.12g\n",i,dcbend_4[i]);}
      for(i=1;i<=ntyp;i++){fprintf(fp,"mdbend_bnd:  dcbend_5[%d] %.12g\n",i,dcbend_5[i]);}
      for(i=1;i<=ntyp;i++){fprintf(fp,"mdbend_bnd:  dcbend_6[%d] %.12g\n",i,dcbend_6[i]);}

      for(i=1;i<=ntyp;i++){fprintf(fp,"mdbend_bnd:  dsbend_0[%d] %.12g\n",i,dsbend_0[i]);}
      for(i=1;i<=ntyp;i++){fprintf(fp,"mdbend_bnd:  dsbend_1[%d] %.12g\n",i,dsbend_1[i]);}
      for(i=1;i<=ntyp;i++){fprintf(fp,"mdbend_bnd:  dsbend_2[%d] %.12g\n",i,dsbend_2[i]);}
      for(i=1;i<=ntyp;i++){fprintf(fp,"mdbend_bnd:  dsbend_3[%d] %.12g\n",i,dsbend_3[i]);}
      for(i=1;i<=ntyp;i++){fprintf(fp,"mdbend_bnd:  dsbend_4[%d] %.12g\n",i,dsbend_4[i]);}
      for(i=1;i<=ntyp;i++){fprintf(fp,"mdbend_bnd:  dsbend_5[%d] %.12g\n",i,dsbend_5[i]);}
      for(i=1;i<=ntyp;i++){fprintf(fp,"mdbend_bnd:  dsbend_6[%d] %.12g\n",i,dsbend_6[i]);}

      fclose(fp);
    }/* end member function */


}; /* end class definition */

#ifdef PUP_ON
PUPmarshall(MDBEND_BND);
#endif

#endif
