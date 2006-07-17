
//==========================================================================
// Should be invoked by group members with real space planes ONCE.
// The invoker must provide a flag vec, allowed_planes[k] 1..nfftz, 
// which is 1 if the plane is assigned to the group and 0 otherwise.
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void NonLocalGrp::setup_non_local_atm_Bspline(){
     igrid   = cnew_imat3(natm,n_interp2,n_interp);
     mn
     dmn_x
     dmn_y 
     dmn_z
     plane_flag;
}
//==========================================================================



//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
double **cnew_mat2(int n1, int n2){
  double **m;
  if(n1 >= 1 && n2 >= 1){
    m = new double *[n1];
    for(int i1=0;i1<n1;i1++){
      m[i1] = new double[n2];
    }//endfor
  }else{
   m = (double **)NULL;
  }//endif
  return m;
}// end routine
//===========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
double **cnew_mat3(int n1, int n2, int n3){
  double ***m;
  if(n1 >= 1 && n2 >= 1 && n3 >=1){
     m = new double **[n1];
     for(int i1=0;i1<n1;i1++){
       m[i1] = new double *[n2];
       for(int i2=0;i2<n2;i1++){
         m[i1][i2] = new double *[n2];
       }//endfor
     }//endfor
  }else{
     m = (double ***)NULL;
  }//endif
  return m;
}// end routine
//===========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
int **cnew_imat2(int n1, int n2){
  int **m;
  if(n1 >= 1 && n2 >= 1){
    m = new int *[n1];
    for(int i1=0;i1<n1;i1++){
      m[i1] = new int[n2];
    }//endfor
  }else{
    m = (int **)NULL;
  }//endif
  return m;
}// end routine
//===========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
int **cnew_imat3(int n1, int n2, int n3){
  int ***m;
  if(n1 >= 1 && n2 >= 1 && n3 >=1){
     m = new int **[n1];
     for(int i1=0;i1<n1;i1++){
       m[i1] = new int *[n2];
       for(int i2=0;i2<n2;i1++){
         m[i1][i2] = new int *[n2];
       }//endfor
     }//endfor
  }else{
     m = (int ***)NULL;
  }//endif
  return m;
}// end routine
//===========================================================================
