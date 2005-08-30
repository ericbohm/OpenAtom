//==========================================================================
//                  Nuevo-Parrinello-Rahman info                            
//             {Variables needed for mem allocation:                        
//                                                  }                       
//                                                                          

#ifndef _GENPAR_RAHMAN_
#define _GENPAR_RAHMAN_

class GENPAR_RAHMAN{

 //----------------
 public:
  int npt_f;                  // Opt: npt_f on or off
  double vol;                 // Num: volume
  double mass_hm;             // Num: Mass of hh^-1 matrix           
  double c1_hm;               // Num: Useful constant                
  double area;                // Num: area                           

  double *vgmat,*vgmat_g;     // Lst: Velocity of hh^-1 matrix Lth:9 
  double *vgmat_glob;         // Lst: Velocity of hh^-1 matrix Lth:9 
  double *fgmat_p,*fgmat_v;   // Lst: Force on hh^-1 matrix    Lth:9 
  double *vtemps,*veigv,*veig;// Lst: Integration temp         Lth:9 
  double *vexpdt,*vsindt;     // Lst: Integration temp         Lth:9 
  double *vtempx, *vtempv;    // Lst: Integration temp         Lth:9 
  double *hmat_t, *hmato;     // Lst: cell shape temps         Lth:9 
  double *fv1,*fv2;           // Lst: Rs temps                 Lth:3  
  double *vgmat_g_wght;       // Lst: weighted vgmat_g         Lth:9 

 //----------------
 //con-destruct:
  GENPAR_RAHMAN(){
   vol     = 0;
   mass_hm = 0;
   c1_hm   = 0;
   area    = 0;
   vgmat   = (double *) cmalloc(9*sizeof(double),"Constructor GENPAR_RAHMAN")-1;
   vgmat_g = (double *) cmalloc(9*sizeof(double),"Constructor GENPAR_RAHMAN")-1;
   vgmat_glob = (double *) cmalloc(9*sizeof(double),"Constructor GENPAR_RAHMAN")-1;
   fgmat_p = (double *) cmalloc(9*sizeof(double),"Constructor GENPAR_RAHMAN")-1;
   fgmat_v = (double *) cmalloc(9*sizeof(double),"Constructor GENPAR_RAHMAN")-1;
   vtemps  = (double *) cmalloc(9*sizeof(double),"Constructor GENPAR_RAHMAN")-1;
   veigv   = (double *) cmalloc(9*sizeof(double),"Constructor GENPAR_RAHMAN")-1;
   veig    = (double *) cmalloc(9*sizeof(double),"Constructor GENPAR_RAHMAN")-1;
   vexpdt  = (double *) cmalloc(9*sizeof(double),"Constructor GENPAR_RAHMAN")-1;
   vsindt  = (double *) cmalloc(9*sizeof(double),"Constructor GENPAR_RAHMAN")-1;
   vtempx  = (double *) cmalloc(9*sizeof(double),"Constructor GENPAR_RAHMAN")-1;
   vtempv  = (double *) cmalloc(9*sizeof(double),"Constructor GENPAR_RAHMAN")-1;
   hmat_t  = (double *) cmalloc(9*sizeof(double),"Constructor GENPAR_RAHMAN")-1;
   hmato   = (double *) cmalloc(9*sizeof(double),"Constructor GENPAR_RAHMAN")-1;
   fv1     = (double *) cmalloc(3*sizeof(double),"Constructor GENPAR_RAHMAN")-1;
   fv2     = (double *) cmalloc(3*sizeof(double),"Constructor GENPAR_RAHMAN")-1;
   vgmat_g_wght = (double *) cmalloc(9*sizeof(double),"Constructor GENPAR_RAHMAN")-1;
  }
 ~GENPAR_RAHMAN(){}

#ifdef PUP_ON
 //----------------
 //pupping
  void pup(PUP::er &p){
    //pupping ints
        p | npt_f;
    //pupping dbles
        p | vol;
        p | mass_hm;
        p | c1_hm;
        p | area;
    //pupping dbl  arrays
        if(npt_f==1){
         pup1d_dbl(p,&vgmat,9);
         pup1d_dbl(p,&vgmat_g,9);
         pup1d_dbl(p,&vgmat_glob,9);
         pup1d_dbl(p,&fgmat_p,9);
         pup1d_dbl(p,&fgmat_v,9);
         pup1d_dbl(p,&vtemps,9);
         pup1d_dbl(p,&veigv,9);
         pup1d_dbl(p,&veig,9);
         pup1d_dbl(p,&vexpdt,9);
         pup1d_dbl(p,&vsindt,9);
         pup1d_dbl(p,&vtempx,9);
         pup1d_dbl(p,&vtempv,9);
         pup1d_dbl(p,&hmat_t,9);
         pup1d_dbl(p,&hmato,9);
         pup1d_dbl(p,&fv1,9);
         pup1d_dbl(p,&fv2,9);
         pup1d_dbl(p,&vgmat_g_wght,9);
	}// endif
#ifdef _PARALLEL_DEBUG_        
    if (p.isUnpacking())
     state_class_out ();
#endif        
  } // end pup
#endif

  void state_class_out(){

     char fileName [255];
     sprintf (fileName, "%d_genpar_rahman.state", CkMyPe());
     FILE *fp;  fp = fopen(fileName,"w");

    fprintf(fp,"npt_f %d\n",npt_f);
    fprintf(fp,"vol %g\n",vol);
    fprintf(fp,"mass_hm %g\n",mass_hm);
    fprintf(fp,"c1_hm %g\n",c1_hm);
    fprintf(fp,"area %g\n",area);
   fclose(fp);

 } // end routine

}; // GENPAR_RAHMAN;

#ifdef PUP_ON
PUPmarshall(GENPAR_RAHMAN);
#endif


#endif

//==========================================================================
