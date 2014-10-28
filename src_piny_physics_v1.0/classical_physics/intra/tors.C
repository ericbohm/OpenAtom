//==========================================================================
//==========================================================================
//                                                                          
// This routine computes the energy and forces from                         
// the intramolecular tors potential.                                       
//                                                                          
//==========================================================================
//==========================================================================

#include "tors.h"

//==========================================================================
//==========================================================================

void tors (const Vector& pos1,  // 1st atom's position
    const Vector& pos2,  // 2nd atom's position
    const Vector& pos3,  // 3rd atom's position
    const Vector& pos4,  // 4th atom's position
    Vector&       force1,      // 1st atom's force
    Vector&       force2,      // 2nd atom's force
    Vector&       force3,      // 3rd atom's force
    Vector&       force4,      // 4th atom's force
    double*       pvten,          // Pressure tensor 
    double*       pvten_tot,      // Pressure tensor 
    int           ityp, 
    StepOutput*   out)
  //==========================================================================
{//begin routine
  //=======================================================================

  // read-only class with bond constants 
  const MDTORS & tors      = readonly_mdintra.mdtors;
  const GENCELL & gencell  = readonly_general_data.gencell;
  int iget_pv_real_inter   = readonly_mdinter.mdenergy_ctrl.iget_pv_real_inter;
  double r122,r12;
  double r132,r13;
  double r232,r23;
  double r422,r42;
  double r432,r43;
  double temp1p,temp1q,temp2q,temp2p;
  double temp3p,temp3q,temp4q,temp4p;
  // Distances              
  double dxp,dyp,dzp,rp2i,rpi,rpqi;
  double dxq,dyq,dzq,rq2i,rqi;     
  double dp13,dp42,dq13,dq42;
  double cost,sint,seps,vtorsc,vtorss,dvtorsc,dvtorss,pre;          
  double sisum;
  double costrp2,costrq2,dp13crp2,dq42crq2;
  double dp13rpq,dq13rpq,dp42rpq,dq42rpq;
  double d1323i,palp,palp1;
  double d4223i,qalp,qalp1; 
  double fpx1,fpy1,fpz1;
  double fpx2,fpy2,fpz2;  
  double fqx4,fqy4,fqz4;
  double fqx2,fqy2,fqz2;  
  int i;

  double wfor;
  double pvten_tmp[10];

  double dx12,dy12,dz12;
  double dx13,dy13,dz13;
  double dx42,dy42,dz42;
  double dx43,dy43,dz43;
  double dx23,dy23,dz23;

  double fx1,fy1,fz1;
  double fx2,fy2,fz2;
  double fx3,fy3,fz3;
  double fx4,fy4,fz4;

  double  vpot;

  //--------------------------------------
  // Local Pointers

  double tors_eq_pow     = tors.eq_pow[ityp];

  double tors_c_0        = tors.c_0[ityp];
  double tors_c_1        = tors.c_1[ityp];
  double tors_c_2        = tors.c_2[ityp];
  double tors_c_3        = tors.c_3[ityp];
  double tors_c_4        = tors.c_4[ityp];
  double tors_c_5        = tors.c_5[ityp];
  double tors_c_6        = tors.c_6[ityp];

  double tors_s_0        = tors.s_0[ityp];
  double tors_s_1        = tors.s_1[ityp];
  double tors_s_2        = tors.s_2[ityp];
  double tors_s_3        = tors.s_3[ityp];
  double tors_s_4        = tors.s_4[ityp];
  double tors_s_5        = tors.s_5[ityp];
  double tors_s_6        = tors.s_6[ityp];

  double tors_dc_0       = tors.dc_0[ityp];
  double tors_dc_1       = tors.dc_1[ityp];
  double tors_dc_2       = tors.dc_2[ityp];
  double tors_dc_3       = tors.dc_3[ityp];
  double tors_dc_4       = tors.dc_4[ityp];
  double tors_dc_5       = tors.dc_5[ityp];
  double tors_dc_6       = tors.dc_6[ityp];

  double tors_ds_0       = tors.ds_0[ityp];
  double tors_ds_1       = tors.ds_1[ityp];
  double tors_ds_2       = tors.ds_2[ityp];
  double tors_ds_3       = tors.ds_3[ityp];
  double tors_ds_4       = tors.ds_4[ityp];
  double tors_ds_5       = tors.ds_5[ityp];
  double tors_ds_6       = tors.ds_6[ityp];

  //=======================================================================
  // 0) Lower cutoff for sin(theta)                   

  wfor  = 1.0;   //(wght_tra_res) not assigned yet? 
  seps = 1.0e-8;

  for(i=1;i<=9;i++){(pvten_tmp)[i]=0;}

  //=======================================================================
  // I) loop over all the tors in steps of nlen to save memory             

  dx12 = pos1.x - pos2.x;
  dy12 = pos1.y - pos2.y;
  dz12 = pos1.z - pos2.z;

  dx13 = pos1.x - pos3.x;
  dy13 = pos1.y - pos3.y;
  dz13 = pos1.z - pos3.z;

  dx42 = pos4.x - pos2.x;
  dy42 = pos4.y - pos2.y;
  dz42 = pos4.z - pos2.z;

  dx43 = pos4.x - pos3.x;
  dy43 = pos4.y - pos3.y;
  dz43 = pos4.z - pos3.z;

  dx23 = pos2.x - pos3.x;
  dy23 = pos2.y - pos3.y;
  dz23 = pos2.z - pos3.z;

  r122   = (dx12*dx12 + dy12*dy12 + dz12*dz12);
  r12    = sqrt(r122);
  r132   = (dx13*dx13 + dy13*dy13 + dz13*dz13);
  r13    = sqrt(r132);
  r422   = (dx42*dx42 + dy42*dy42 + dz42*dz42);
  r42    = sqrt(r422);
  r432   = (dx43*dx43 + dy43*dy43 + dz43*dz43);
  r43    = sqrt(r432);
  r232   = (dx23*dx23 + dy23*dy23 + dz23*dz23);
  r23    = sqrt(r232);

  //------------------------------------------------------------------
  // F) Construct vector in plane of 123 perp to r23 call it(dxp,dyp,dzp) 

  d1323i  = 1.0/(dx13*dx23 + dy13*dy23 + dz13*dz23);
  palp   =-(dx12*dx23 + dy12*dy23 + dz12*dz23)*d1323i;
  palp1  = 1.0+palp;
  dxp    = palp*dx13+dx12;
  dyp    = palp*dy13+dy12;
  dzp    = palp*dz13+dz12;
  rp2i   = 1.0/(dxp*dxp+dyp*dyp+dzp*dzp);
  rpi    = sqrt(rp2i);

  //-----------------------------------------------------------------
  //  G) Get dot product of this vector (dxp,dyp,dzp) with r13 and r42 

  dp13  = dxp*dx13 + dyp*dy13 + dzp*dz13;
  dp42  = dxp*dx42 + dyp*dy42 + dzp*dz42;

  //-----------------------------------------------------------------
  //  H) Get derivative of parameter palp used to construct (dxp,dyp,dzp)
  //     with respect to atoms (1-4)                                     

  fpx1   = -palp1*d1323i*dx23;
  fpy1   = -palp1*d1323i*dy23;
  fpz1   = -palp1*d1323i*dz23;
  fpx2   = (dx23-dxp)*d1323i;
  fpy2   = (dy23-dyp)*d1323i;
  fpz2   = (dz23-dzp)*d1323i;

  //------------------------------------------------------------------
  //  I) Construct vector in plane of 234 perp to r23 (dxq,dyq,dzq) 

  d4223i = 1.0/(dx42*dx23 + dy42*dy23 + dz42*dz23);
  qalp   =-(dx43*dx23 + dy43*dy23 + dz43*dz23)*d4223i;
  qalp1  = 1.0+qalp;
  dxq    = qalp*dx42+dx43;
  dyq    = qalp*dy42+dy43;
  dzq    = qalp*dz42+dz43;
  rq2i   = 1.0/(dxq*dxq+dyq*dyq+dzq*dzq);
  rqi    = sqrt(rq2i);

  //------------------------------------------------------------------
  //  J) Get dot product of this vector (dxq,dyq,dzq) with r13 and r42

  dq13  = (dxq*dx13 + dyq*dy13 + dzq*dz13);
  dq42  = (dxq*dx42 + dyq*dy42 + dzq*dz42);

  //-------------------------------------------------------------------
  //  K) Get derivative of parameter qalp used to construct (dxq,dyq,dzq)
  //     with respect to atoms (1-4)                                   

  fqx4   = -qalp1*d4223i*dx23;
  fqy4   = -qalp1*d4223i*dy23;
  fqz4   = -qalp1*d4223i*dz23;

  fqx2   = (qalp*dx23 - dxq)*d4223i;
  fqy2   = (qalp*dy23 - dyq)*d4223i;
  fqz2   = (qalp*dz23 - dzq)*d4223i;

  //--------------------------------------------------------------------
  //  L) Get cosine  and sine of the angle                              

  rpqi = rpi*rqi;
  cost   = (dxp*dxq+dyp*dyq+dzp*dzq)*rpqi;
  cost   = (cost < 1.0 ? cost:1.0);
  cost   = (cost > -1.0 ? cost:-1.0);

  sint = sqrt(1.0 - cost*cost);
  sint = (sint > seps ? sint:seps);

  //---------------------------------------------------------------------
  //  M) Calculate the potential energy from cosine and sine power       
  //     series using Horner's method                                    

  vtorsc = ((((((tors_c_6
                *cost + tors_c_5)
              *cost + tors_c_4)
            *cost + tors_c_3)
          *cost + tors_c_2)
        *cost + tors_c_1)
      *cost + tors_c_0);

  sisum = (((((tors_s_6
              *sint + tors_s_5)
            *sint + tors_s_4)
          *sint + tors_s_3)
        *sint + tors_s_2)
      *sint);

  vtorss = sisum*cost + (tors_s_1)*sint;

  vpot = vtorsc + vtorss;

  //-------------------------------------------------------------------
  //  N) Get the force on the atoms using the chain rule 

  dvtorsc = (((((tors_dc_6
              *cost + tors_dc_5)
            *cost + tors_dc_4)
          *cost + tors_dc_3)
        *cost + tors_dc_2)
      *cost + tors_dc_1);

  dvtorss =  ((((tors_ds_6
            *sint + tors_ds_5)
          *sint + tors_ds_4)
        *sint + tors_ds_3)
      *sint + tors_ds_2);

  dvtorss = dvtorss*cost + (tors_ds_1);

  pre     = -(dvtorsc + sisum - dvtorss*cost/sint);

  costrp2  = cost*rp2i;
  costrq2  = cost*rq2i;
  dp13crp2 = dp13*costrp2;
  dq42crq2 = dq42*costrq2;
  dp13rpq  = dp13*rpqi;
  dq13rpq  = dq13*rpqi;
  dp42rpq  = dp42*rpqi;
  dq42rpq  = dq42*rpqi;
  temp1p   = (dq13rpq-dp13crp2)*pre;
  temp1q   = (dp42rpq-dq42crq2)*pre;
  temp2p   = (costrp2*palp1)*pre;
  temp2q   = (costrq2*qalp1)*pre;
  temp3q   = (rpqi*palp1)*pre;  
  temp3p   = (rpqi*qalp1)*pre;  
  temp4p   = (costrp2-qalp*rpqi)*pre;
  temp4q   = (qalp*costrq2-rpqi)*pre;

  fx1 = -dxp*temp2p+dxq*temp3q+temp1p*fpx1;
  fy1 = -dyp*temp2p+dyq*temp3q+temp1p*fpy1;
  fz1 = -dzp*temp2p+dzq*temp3q+temp1p*fpz1;

  fx4 = -dxq*temp2q+dxp*temp3p+temp1q*fqx4;
  fy4 = -dyq*temp2q+dyp*temp3p+temp1q*fqy4;
  fz4 = -dzq*temp2q+dzp*temp3p+temp1q*fqz4;

  fx2 =  dxp*temp4p+dxq*temp4q+temp1p*fpx2+temp1q*fqx2;
  fy2 =  dyp*temp4p+dyq*temp4q+temp1p*fpy2+temp1q*fqy2;
  fz2 =  dzp*temp4p+dzq*temp4q+temp1p*fpz2+temp1q*fqz2;

  fx3 = -(fx1 + fx2 + fx4);
  fy3 = -(fy1 + fy2 + fy4);
  fz3 = -(fz1 + fz2 + fz4);

  //--------------------------------------------------------------------
  //  O) Get the pressure tensor                                        

  pvten_tmp[1] = dx13*fx1 + dx23*fx2 + dx43*fx4; 
  pvten_tmp[5] = dy13*fy1 + dy23*fy2 + dy43*fy4;
  pvten_tmp[9] = dz13*fz1 + dz23*fz2 + dz43*fz4;

  pvten_tmp[2] = dx13*fy1 + dx23*fy2 + dx43*fy4;
  pvten_tmp[4] = dy13*fx1 + dy23*fx2 + dy43*fx4;
  pvten_tmp[3] = dx13*fz1 + dx23*fz2 + dx43*fz4;

  pvten_tmp[7] = dz13*fx1 + dz23*fx2 + dz43*fx4;
  pvten_tmp[6] = dy13*fz1 + dy23*fz2 + dy43*fz4;
  pvten_tmp[8] = dz13*fy1 + dz23*fy2 + dz43*fy4; 

  //---------------------------------------------------------------------
  //  P) Sum the potential energy                                             

  //----------------------------------------------------------------------

  force1.x += wfor*fx1;
  force1.y += wfor*fy1;
  force1.z += wfor*fz1;

  force2.x += wfor*fx2;
  force2.y += wfor*fy2;
  force2.z += wfor*fz2;

  force3.x += wfor*fx3;
  force3.y += wfor*fy3;
  force3.z += wfor*fz3;

  force4.x += wfor*fx4;
  force4.y += wfor*fy4;
  force4.z += wfor*fz4;

  //==========================================================================
  // Increment the Pressure tensor 

  for(i=1;i<=9;i++){
    pvten[i] += pvten_tmp[i]*wfor;
  }//endfor

  if(iget_pv_real_inter==1){    
    for(i=1;i<=9;i++){
      pvten_tot[i] += pvten_tmp[i];
    }//endfor
  }//endif

  out->pe    += vpot;
  out->vtors += vpot;
  //--------------------------------------------------------------------------
}// end routine
//==========================================================================


