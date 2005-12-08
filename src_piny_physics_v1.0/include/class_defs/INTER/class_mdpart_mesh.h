/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                      class_mdpart_mesh.h                                 */
/*                                                                          */
/*    Class definition for particle mesh Ewald                              */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#ifndef _MDPART_MESH_
#define _MDPART_MESH_

class MDPART_MESH {
 public:
  int pme_on;                       /*Opt: PME on                         */
  int kmax_pme;                     /*Num: User PME Mesh estimate         */
  int n_interp;                     /*Num: Order of interpolation         */
  int pme_res_on;
  int kmax_pme_res;                 
  int n_interp_res;
  int nlen_pme;                     /*Num: Scr lngth                      */

  int nktot_pme;                    /*Num: equal to ewald->nktot          */
  int ngrid_a,ngrid_b,ngrid_c;      /*Num: PME mesh                       */
  int nktot_pme_res;                /*Num: Equal to ewald->nktot_res      */
  int ngrid_a_res,ngrid_b_res,ngrid_c_res;      

  int ninterp;    
  int ninterp_tot;   

  int *iatemp,*ibtemp,*ictemp;      /*Lst: Lth: nlen_pme                  */
  int *nc,*ioff_c;                  /*Lst: Lth ngrid_c;                   */

  int **igrid_a,**igrid_b,**igrid_c;/*Lst: Lth: ninterp*nlen_pme          */
  int **igrid_now;                  /*Lst: Lth: ninterp*nlen_pme         */

  double ecut;                      /*Num: Cutoff in Hartree              */
  double ecut_res;                  /*Num: Cutoff in Hartree              */
  double *bweight_tot_res;          
  double *bweight_tot;              /*Lst: Lth: nktot                     */
  double *qgrid;                    /*Lst: Interpolation grid 
                                            Lth:2*ngrid_a*_b*_c           */
  double *qgrid_scr;                /*Lst: Interpolation grid 
                                            Lth:2*ngrid_a*_b*_c           */
  double *qgrid_tmp_real;           /*Lst: Lth: nktot                     */
  double *qgrid_tmp_imag;           /*Lst: Lth: nktot                     */
  double *frac_a,*frac_b,*frac_c;   /*Lst: Lth:nlen_pme                   */
  double *aj,*rn,*rn1;              /*Lst: Lth: ninterp*nlen_pme          */

  double **ua,**ub,**uc;            /*Lst: Lth:ninterp*nlen_pme           */
  double **mn_a,**mn_b,**mn_c;      /*Lst: Lth: ninterp*nlen_pme          */
  double **dmn_a,**dmn_b,**dmn_c;   /*Lst: Lth: ninterp*nlen_pme          */
  double **qgrid_now;               /*Lst: Lth: ninterp*nlen_pme         */

//=====================================================================================
// Default constructor/destructor

   MDPART_MESH(){
    pme_on   = 0;                      
    kmax_pme = 0;                    
    n_interp = 0;                    
    pme_res_on   = 0;
    kmax_pme_res = 0;                 
    n_interp_res = 0;
    nlen_pme     = 0;                    
    nktot_pme    = 0;                   
    ngrid_a = 0;
    ngrid_b = 0;
    ngrid_c = 0;     
    nktot_pme_res = 0;               
    ngrid_a_res   = 0;
    ngrid_b_res   = 0;
    ngrid_c_res   = 0;      
    ninterp = 0;    
    ninterp_tot = 0;   
   }
  ~MDPART_MESH(){}

//=====================================================================================
// Pack/Unpack

#ifdef PUP_ON
  void pup(PUP::er &p){

    // PUP ints

    p | pme_on;            
    p | kmax_pme;          
    p | n_interp;          
    p | pme_res_on;
    p | kmax_pme_res;                 
    p | n_interp_res;
    p | nlen_pme;          

    p | nktot_pme;         
    p | ngrid_a;
    p | ngrid_b;
    p | ngrid_c; 
    p | nktot_pme_res;           
    p | ngrid_a_res;
    p | ngrid_b_res;
    p | ngrid_c_res;      

    p | ninterp;    
    p | ninterp_tot;   

    // PUP doubles 

    p | ecut; 
    p | ecut_res;  

    // PUP arrays
#ifdef PME_IMPLEMENTED
    pup1d_int(p,&iatemp,nlen_pme);
    pup1d_int(p,&ibtemp,nlen_pme);
    pup1d_int(p,&ictemp,nlen_pme); 
    pup1d_int(p,&nc,ngrid_c);
    pup1d_int(p,&ioff_c,ngrid_c); 

    pup2d_int(p,&igrid_a,ninterp,nlen_pme);
    pup2d_int(p,&igrid_b,ninterp,nlen_pme);
    pup2d_int(p,&igrid_c,ninterp,nlen_pme);
    pup2d_int(p,&igrid_now,ninterp,nlen_pme);

    pup1d_dbl(p,&bweight_tot_res,nktot_pme);          
    pup1d_dbl(p,&bweight_tot,nktot_pme);
    pup1d_dbl(p,&qgrid,2*ngrid_a*ngrid_b*ngrid_c); 
    pup1d_dbl(p,&qgrid_scr,2*ngrid_a*ngrid_b*ngrid_c);
    pup1d_dbl(p,&qgrid_tmp_real,nktot_pme); 
    pup1d_dbl(p,&qgrid_tmp_imag,nktot_pme); 
    pup1d_dbl(p,&frac_a,nlen_pme);
    pup1d_dbl(p,&frac_b,nlen_pme);
    pup1d_dbl(p,&frac_c,nlen_pme);
    pup1d_dbl(p,&aj,ninterp*nlen_pme);
    pup1d_dbl(p,&rn,ninterp*nlen_pme);
    pup1d_dbl(p,&rn1,ninterp*nlen_pme);

    pup2d_dbl(p,&ua,ninterp,nlen_pme,"mdpart_mesh");
    pup2d_dbl(p,&ub,ninterp,nlen_pme,"mdpart_mesh");
    pup2d_dbl(p,&uc,ninterp,nlen_pme,"mdpart_mesh");
    pup2d_dbl(p,&mn_a,ninterp,nlen_pme,"mdpart_mesh");
    pup2d_dbl(p,&mn_b,ninterp,nlen_pme,"mdpart_mesh");
    pup2d_dbl(p,&mn_c,ninterp,nlen_pme,"mdpart_mesh");
    pup2d_dbl(p,&dmn_a,ninterp,nlen_pme,"mdpart_mesh");
    pup2d_dbl(p,&dmn_b,ninterp,nlen_pme,"mdpart_mesh");
    pup2d_dbl(p,&dmn_c,ninterp,nlen_pme,"mdpart_mesh");
    pup2d_dbl(p,&qgrid_now,ninterp,nlen_pme,"mdpart_mesh");
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
     sprintf (fileName, "%d_mdpart_mesh.state", CkMyPe());
     FILE *fp;  fp = fopen(fileName,"w");

    // Print ints

    fprintf(fp,"mdpart_mesh: pme_on %d\n",pme_on);            
    fprintf(fp,"mdpart_mesh: kmax_pme %d\n",kmax_pme);          
    fprintf(fp,"mdpart_mesh: n_interp %d\n",n_interp);          
    fprintf(fp,"mdpart_mesh: pme_res_on %d\n",pme_res_on);
    fprintf(fp,"mdpart_mesh: kmax_pme_res %d\n",kmax_pme_res);                 
    fprintf(fp,"mdpart_mesh: n_interp_res %d\n",n_interp_res);
    fprintf(fp,"mdpart_mesh: nlen_pme %d\n",nlen_pme);          

    fprintf(fp,"mdpart_mesh: nktot_pme %d\n",nktot_pme);         
    fprintf(fp,"mdpart_mesh: ngrid_a %d\n",ngrid_a);
    fprintf(fp,"mdpart_mesh: ngrid_b %d\n",ngrid_b);
    fprintf(fp,"mdpart_mesh: ngrid_c %d\n",ngrid_c); 
    fprintf(fp,"mdpart_mesh: nktot_pme_res %d\n",nktot_pme_res);           
    fprintf(fp,"mdpart_mesh: ngrid_a_res %d\n",ngrid_a_res);
    fprintf(fp,"mdpart_mesh: ngrid_b_res %d\n",ngrid_b_res);
    fprintf(fp,"mdpart_mesh: ngrid_c_res %d\n",ngrid_c_res);      

    fprintf(fp,"mdpart_mesh: ninterp %d\n",ninterp);    
    fprintf(fp,"mdpart_mesh: ninterp_tot %d\n",ninterp_tot);   

    // Print doubles 

    fprintf(fp,"mdpart_mesh: ecut %.12g\n",ecut); 
    fprintf(fp,"mdpart_mesh: ecut_res %.12g\n",ecut_res);  

    fclose(fp);

  }/* end member function */

}; /* end class definition */

#ifdef PUP_ON
PUPmarshall(MDPART_MESH);
#endif

#endif



