//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                          
//==========================================================================

#ifndef _PhysicsParamTransfer_
#define _PhysicsParamTransfer_

//==========================================================================
class PhysicsParamTransfer{

 //---------------------------------------------------------------------------
 public:

 //---------------------------------------------------------------------------
 //default con-destruct:
   PhysicsParamTransfer(){};
  ~PhysicsParamTransfer(){};

 //---------------------------------------------------------------------------
 // functions

  static void ParaInfoInit(CPcharmParaInfo *);
  static void control_mapping_function(CPcharmParaInfo *,int);
  static void control_new_mapping_function(CPcharmParaInfo *,int);
  static void compute_lines_per_plane_half_plane(int ,int ,int ,int ,
                             double *, double *,int *, int *, int *);
  static void compute_lines_per_plane_half_sphere(int ,int ,int ,int ,
                            double *, double *,int *, int *, int *);
  static void get_Sfgrp_params(int ,int ,int ,int *, int *, int *);
  static void get_Sfgrp_max(int , int, int *);

};
//==========================================================================

#endif
