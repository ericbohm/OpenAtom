//=========================================================================
//=========================================================================
// This routine computes the energy and forces from                        
// the intramolecular bend potential.                                      
//=========================================================================
//=========================================================================

#include "bend.h"

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

void bend (const Vector& pos1,      // 1st atom's position
    const Vector& pos2,      // 2nd atom's position
    const Vector& pos3,      // 3rd atom's position
    Vector&       force1,    // 1st atom's force
    Vector&       force2,    // 2nd atom's force
    Vector&       force3,    // 3rd atom's force
    double*       pvten,     // Pressure tensor 
    double*       pvten_tot, // Pressure tensor 
    int           ityp,      // bend typ between atoms 
    StepOutput*   out)
  //==========================================================================
{//begin routine
  //==========================================================================

  // Define local pointers                


  // reference to read only class that contains data
  const MDBEND & bend     = readonly_mdintra.mdbend;
  const GENCELL & gencell = readonly_general_data.gencell;
  int iget_pv_real_inter  = readonly_mdinter.mdenergy_ctrl.iget_pv_real_inter;

  double bend_c_0        = bend.c_0[ityp];
  double bend_c_1        = bend.c_1[ityp];
  double bend_c_2        = bend.c_2[ityp];
  double bend_c_3        = bend.c_3[ityp];
  double bend_c_4        = bend.c_4[ityp];
  double bend_c_5        = bend.c_5[ityp];
  double bend_c_6        = bend.c_6[ityp];

  double bend_s_0        = bend.s_0[ityp];
  double bend_s_1        = bend.s_1[ityp];
  double bend_s_2        = bend.s_2[ityp];
  double bend_s_3        = bend.s_3[ityp];
  double bend_s_4        = bend.s_4[ityp];
  double bend_s_5        = bend.s_5[ityp];
  double bend_s_6        = bend.s_6[ityp];

  double bend_dc_0       = bend.dc_0[ityp];
  double bend_dc_1       = bend.dc_1[ityp];
  double bend_dc_2       = bend.dc_2[ityp];
  double bend_dc_3       = bend.dc_3[ityp];
  double bend_dc_4       = bend.dc_4[ityp];
  double bend_dc_5       = bend.dc_5[ityp];
  double bend_dc_6       = bend.dc_6[ityp];

  double bend_ds_0       = bend.ds_0[ityp];
  double bend_ds_1       = bend.ds_1[ityp];
  double bend_ds_2       = bend.ds_2[ityp];
  double bend_ds_3       = bend.ds_3[ityp];
  double bend_ds_4       = bend.ds_4[ityp];
  double bend_ds_5       = bend.ds_5[ityp];
  double bend_ds_6       = bend.ds_6[ityp];

  //----------------------------------------------------------------------
  // Local Variables 

  double dx12,dy12,dz12;
  double dx13,dy13,dz13;
  double dx23,dy23,dz23;
  double r122i,r12i;       // (r-r_0)^2 and (r-r_0)  
  double r322i,r32i;
  double cost,sint;
  double vbendc,vbends;    // Bond potential       
  double dvbendc,dvbends;  // Derivative of bend pot 
  double vpot;
  double sisum;            // sine sum         
  double pre;              // Force prefactor  
  double rpmag;            // 1/(r12*r32)      
  double cos122,cos322;    // cos/r122, cos/r322 
  double fx1,fy1,fz1;
  double fx2,fy2,fz2;
  double fx3,fy3,fz3;
  double wfor;
  double pvten_tmp[10];  

  //=======================================================================
  // 0) Lower cutoff for sin(theta)                                        

  const double seps = 1.0e-8;

  wfor  = 1.0;   //(wght_tra_res) not assigned yet? 

  int i;
  for(i=1;i<=9;i++){(pvten_tmp)[i]=0;}

  //=======================================================================

  dx12 = pos1.x - pos2.x;
  dy12 = pos1.y - pos2.y;
  dz12 = pos1.z - pos2.z;

  dx23 = pos3.x - pos2.x;
  dy23 = pos3.y - pos2.y;
  dz23 = pos3.z - pos2.z;

  //-------------------------------------------------------------------

  r122i  = 1.0/(dx12*dx12 + dy12*dy12 + dz12*dz12);
  r12i   = sqrt(r122i);

  r322i  = 1.0/(dx23*dx23 + dy23*dy23 + dz23*dz23);
  r32i   = sqrt(r322i);
  rpmag  = r12i*r32i;

  //-------------------------------------------------------------------
  //  F) Calculate the cosine and sine of the angle between them       
  //     (cos(theta_{123}), sine(theta_{123}))                         

  cost = (dx12*dx23 + dy12*dy23 + dz12*dz23)*rpmag;

  cost   = (cost <  1.0 ? cost: 1.0);
  cost   = (cost > -1.0 ? cost:-1.0);

  sint   = sqrt(1.0 - cost*cost);
  sint   = (sint > seps ? sint:seps);

  cos122 = cost*r122i;
  cos322 = cost*r322i;

  //----------------------------------------------------------------------
  // H) Get the bending potential energy from the cosine and sine power
  //     series using Horner's method                                   

  vbendc = ((((((bend_c_6
                *cost + bend_c_5)
              *cost + bend_c_4)
            *cost + bend_c_3)
          *cost + bend_c_2)
        *cost + bend_c_1)
      *cost + bend_c_0);

  sisum = (((((bend_s_6
              *sint + bend_s_5)
            *sint + bend_s_4)
          *sint + bend_s_3)
        *sint + bend_s_2)
      *sint);

  vbends = sisum*cost + bend_s_1*sint;

  vpot = vbendc + vbends;

  //-------------------------------------------------------------------
  //  I) Get the force on the atoms using the chain rule and Horner's  

  dvbendc =  (((((bend_dc_6
              *cost + bend_dc_5)
            *cost + bend_dc_4)
          *cost + bend_dc_3)
        *cost + bend_dc_2)
      *cost + bend_dc_1);

  dvbends =  ((((bend_ds_6
            *sint + bend_ds_5)
          *sint + bend_ds_4)
        *sint + bend_ds_3)
      *sint + bend_ds_2);

  dvbends = dvbends*cost + bend_ds_1;

  pre     = -(dvbendc + sisum - dvbends*cost/sint);

  fx1 = (dx23*rpmag - dx12*cos122)*pre;
  fy1 = (dy23*rpmag - dy12*cos122)*pre;
  fz1 = (dz23*rpmag - dz12*cos122)*pre;

  fx3 = (dx12*rpmag - dx23*cos322)*pre;
  fy3 = (dy12*rpmag - dy23*cos322)*pre;
  fz3 = (dz12*rpmag - dz23*cos322)*pre;

  fx2 = -(fx1 + fx3); 
  fy2 = -(fy1 + fy3);
  fz2 = -(fz1 + fz3);

  //----------------------------------------------------------------------
  //  J) Get the pressure tensor                                         

  if(gencell.iperd == 2 || gencell.iperd == 3) {
    pvten_tmp[1]  = dx12*fx1 + dx23*fx3;
    pvten_tmp[5]  = dy12*fy1 + dy23*fy3;
    pvten_tmp[9]  = dz12*fz1 + dz23*fz3;
    pvten_tmp[2]  = dx12*fy1 + dx23*fy3;
  }//endif
  if(gencell.iperd == 3) {
    pvten_tmp[3] = dx12*fz1 + dx23*fz3;
    pvten_tmp[6] = dy12*fz1 + dy23*fz3;
  }//endif

  //---------------------------------------------------------------------
  //  K) Sum the potential energy                                        

  //======================================================================
  // Increment the Pressure tensor 

  force1.x += wfor*fx1;
  force1.y += wfor*fy1;
  force1.z += wfor*fz1;

  force2.x += wfor*fx2;
  force2.y += wfor*fy2;
  force2.z += wfor*fz2;

  force3.x += wfor*fx3;
  force3.y += wfor*fy3;
  force3.z += wfor*fz3;

  pvten_tmp[4] = pvten_tmp[2];
  pvten_tmp[7] = pvten_tmp[3];
  pvten_tmp[8] = pvten_tmp[6];

  for(i=1;i<=9;i++){
    pvten[i]  += pvten_tmp[i]*wfor;
  }//endfor

  if(iget_pv_real_inter==1){
    for(i=1;i<=9;i++){
      pvten_tot[i] += pvten_tmp[i];
    }//endfor
  }//endif

  out->pe    += vpot;
  out->vbend += vpot;

  //----------------------------------------------------------------------
}//end routine
//======================================================================




