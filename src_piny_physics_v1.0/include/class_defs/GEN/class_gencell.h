//==========================================================================
//                  Simulation cell data                                    

#ifndef _GENCELL_
#define _GENCELL_

class GENCELL{

 //----------------
 public:
  int cp_any_on;              // Opt: Any cp option on
  int iperd;                  // Opt: Periodicity of system(0-3)     
  int intra_perds;            // Opt: Intramol periodic flag         
  int cubic_box_flag;         // Opt: Cubic box flag for fast imaging 
  int hmat_int_typ;           // Opt: Flag to distinguish between normal and
                              //      upper triangular integration 
  int hmat_cons_typ;          // Opt: Constraint on the box/h_mat 
                              //       0=none,1=orthorhombic,2=monoclinic
  int imov_cp_box;            // Opt: Flag to move center of cp_box 

  double vol;                 // Num: The volume                     
  double vol_cp;              // Num: The volume                     
  double vol0;                // Num: The initial volume            
  double area;                // Num: The area                      
  double alpha_conv_dual;

  double *hmat,*hmati;            // Lst: Cell shape; Lth:9              
  double *hmat_ewd,*hmat_ewd_cp;  // Lst: Cell shape for g-vectors Lth:9 
  double *hmat_cp,*hmati_cp;      // Lst: Cell shape; Lth:9              
  double *cp_box_center;           
  double *cp_box_center_rel;
  double *cp_vbox_center;         //box center velocities 
  double *cp_fbox_center;         //box center forces     

 //----------------
 //con-destruct:
  GENCELL(){
   cp_any_on       = 0;
   iperd           = 0;
   intra_perds     = 0;
   cubic_box_flag  = 0;
   hmat_int_typ    = 0;
   hmat_cons_typ   = 0;
   imov_cp_box     = 0;
   vol             = 0;
   vol_cp          = 0;
   vol0            = 0;
   area            = 0;
   alpha_conv_dual = 0;
   hmat       = (double *)cmalloc(9*sizeof(double),"class_gencell_constructor")-1;
   hmati      = (double *)cmalloc(9*sizeof(double),"class_gencell_constructor")-1;
   hmat_ewd   = (double *)cmalloc(9*sizeof(double),"class_gencell_constructor")-1;
   hmat_ewd_cp= (double *)cmalloc(9*sizeof(double),"class_gencell_constructor")-1;
   hmat_cp    = (double *)cmalloc(9*sizeof(double),"class_gencell_constructor")-1;
   hmati_cp   = (double *)cmalloc(9*sizeof(double),"class_gencell_constructor")-1;
   cp_box_center = (double *)
                   cmalloc(3*sizeof(double),"class_gencell_constructor")-1;
   cp_box_center_rel = (double *)
                   cmalloc(3*sizeof(double),"class_gencell_constructor")-1;
   cp_vbox_center  = (double *)
                   cmalloc(3*sizeof(double),"class_gencell_constructor")-1;
   cp_fbox_center  = (double *)
                   cmalloc(3*sizeof(double),"class_gencell_constructor")-1;
  };

  ~GENCELL(){
    cfree(&hmat[1],"class_gencell_destructor");
    cfree(&hmati[1],"class_gencell_destructor");
    cfree(&hmat_ewd[1],"class_gencell_destructor");
    cfree(&hmat_ewd_cp[1],"class_gencell_destructor");
    cfree(&hmat_cp[1],"class_gencell_destructor");
    cfree(&hmati_cp[1],"class_gencell_destructor");
    cfree(&cp_box_center[1],"class_gencell_destructor");
    cfree(&cp_box_center_rel[1],"class_gencell_destructor");
    cfree(&cp_vbox_center[1],"class_gencell_destructor");
    cfree(&cp_fbox_center[1],"class_gencell_destructor");
   };

#ifdef PUP_ON
 //----------------
 //pupping
  void pup(PUP::er &p){
    //pupping ints
      p | cp_any_on;
      p | iperd;
      p | intra_perds;
      p | cubic_box_flag;
      p | hmat_int_typ;
      p | hmat_cons_typ;
      p | imov_cp_box;
    //pupping dbles
      p | vol;
      p | vol_cp;
      p | vol0;
      p | area;
      p | alpha_conv_dual;
    //pupping dbl  arrays
      pup1d_dbl(p,&hmat,9);
      pup1d_dbl(p,&hmati,9);
      pup1d_dbl(p,&hmat_ewd,9);
      pup1d_dbl(p,&hmat_ewd_cp,9);
      pup1d_dbl(p,&hmat_cp,9);
      pup1d_dbl(p,&hmati_cp,9);
      if(cp_any_on==1){
       pup1d_dbl(p,&cp_box_center,3);
       pup1d_dbl(p,&cp_box_center_rel,3);
       pup1d_dbl(p,&cp_vbox_center,3);
       pup1d_dbl(p,&cp_fbox_center,3);
      }//endif
#ifdef _PARALLEL_DEBUG_        
    if (p.isUnpacking())
     state_class_out ();
#endif      
  } // end pup
#endif

  void state_class_out(){

     int i;
     char fileName [255];
     sprintf (fileName, "%d_gencell.state", CkMyPe());
     FILE *fp;  fp = fopen(fileName,"w");

    // ints
     fprintf(fp,"iperd %d\n",iperd);
     fprintf(fp,"intra_perds %d\n",intra_perds);
     fprintf(fp,"cubic_box_flag %d\n",cubic_box_flag);
     fprintf(fp,"hmat_int_typ %d\n",hmat_int_typ);
     fprintf(fp,"hmat_cons_typ %d\n",hmat_cons_typ);
     fprintf(fp,"imov_cp_box %d\n",imov_cp_box);
    // dbles
     fprintf(fp,"vol %g\n",vol);
     fprintf(fp,"vol_cp %g\n",vol_cp);
     fprintf(fp,"vol0 %g\n",vol0);
     fprintf(fp,"area %g\n",area);
     fprintf(fp,"alpha_conv_dual %g\n",alpha_conv_dual);
    // dbl  arrays
      for(i=1;i<=9;i++){fprintf(fp,"hmat[%d], %lg\n",i,hmat[i]);}
      for(i=1;i<=9;i++){fprintf(fp,"hmati[%d], %g\n",i,hmati[i]);}
      for(i=1;i<=9;i++){fprintf(fp,"hmat_ewd[%d], %g\n",i,hmat_ewd[i]);}
      for(i=1;i<=9;i++){fprintf(fp,"hmat_ewd_cp[%d], %g\n",i,hmat_ewd_cp[i]);}
      for(i=1;i<=9;i++){fprintf(fp,"hmat_cp[%d], %g\n",i,hmat_cp[i]);}
      for(i=1;i<=9;i++){fprintf(fp,"hmati_cp[%d], %g\n",i,hmati_cp[i]);}
      if(cp_any_on==1){
       for(i=1;i<=3;i++){fprintf(fp,"cp_box_center[%d] %g\n",i,cp_box_center[i]);}
       for(i=1;i<=3;i++){fprintf(fp,"cp_box_center_rel[%d] %g\n",i,cp_box_center_rel[i]);}
      }//endif
    fclose(fp);
  }// end routine


}; // GENCELL;

#ifdef PUP_ON
PUPmarshall(GENCELL);
#endif

#endif

//==========================================================================
