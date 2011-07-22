//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================

#include "standard_include.h"
#include "../../../include/debug_flags.h"
#include "../class_defs/allclass_gen.h"
#include "../class_defs/allclass_mdintegrate.h"
#include "../../../include/Atoms.h"
#include "../class_defs/ATOM_OPERATIONS/class_atomintegrate.h"

//============================================================================
//        NVE Integration : 2nd half
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void ATOMINTEGRATE::integrate_nve_2nd_half(int itime,int natm,Atom *atoms,
                                double *eKinetic,int natmNow,int natmStr,int natmEnd)
//============================================================================
    {//begin routine 
//============================================================================
// Local variables

   GENERAL_DATA *general_data = GENERAL_DATA::get();
#include "../class_defs/allclass_strip_gen.h"
   double dt  = gentimeinfo->dt;
   double dt2 = dt*0.5;

//============================================================================
// Evolve the system : if its not the first time step

   if(itime>0){
     for(int i=natmStr;i<natmEnd;i++){
       atoms[i].vx += dt2*atoms[i].fx/atoms[i].m;
       atoms[i].vy += dt2*atoms[i].fy/atoms[i].m;
       atoms[i].vz += dt2*atoms[i].fz/atoms[i].m;
     }//endfor  
   }//endif
 
   computeEkin(natm,atoms,eKinetic,natmNow,natmStr,natmEnd);

//----------------------------------------------------------------------------
  }//end routine
//============================================================================


//============================================================================
//        NVE Integration : 1st half
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void ATOMINTEGRATE::integrate_nve_1st_half(int natm,Atom *atoms,
                                           int natmNow,int natmStr,int natmEnd)
//============================================================================
    {//begin routine 
//============================================================================
// Local variables

   GENERAL_DATA *general_data = GENERAL_DATA::get();
#include "../class_defs/allclass_strip_gen.h"
   double dt  = gentimeinfo->dt;
   double dt2 = dt*0.5;

//============================================================================
// Evolve the system 

   for(int i=natmStr;i<natmEnd;i++){
     atoms[i].vx += dt2*atoms[i].fx/atoms[i].m;
     atoms[i].vy += dt2*atoms[i].fy/atoms[i].m;
     atoms[i].vz += dt2*atoms[i].fz/atoms[i].m;
   }//endfor  

   for(int i=natmStr;i<natmEnd;i++){
     atoms[i].x += dt*atoms[i].vx;
     atoms[i].y += dt*atoms[i].vy;
     atoms[i].z += dt*atoms[i].vz;
   }//endfor  

//----------------------------------------------------------------------------
  }//end routine
//============================================================================


//============================================================================
//        NVT Integration : 2nd half
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void ATOMINTEGRATE::integrate_nvt_2nd_half(int itime,int natm,int len_nhc,
                           Atom *atoms,AtomNHC *atomsNHC,double *eKinetic,
                           double *eKineticNhc,double *potNhc,
                           int natmNow,int natmStr,int natmEnd)
//============================================================================
    {//begin routine 
//============================================================================
// Local variables

   MDINTEGRATE  *mdintegrate  = MDINTEGRATE::get();
   GENERAL_DATA *general_data = GENERAL_DATA::get();
   MDTHERM_INFO *mdtherm_info = &(mdintegrate->mdtherm_info);
#include "../class_defs/allclass_strip_gen.h"

   double dt  = gentimeinfo->dt;
   double dt2 = dt*0.5;
   int nresp  = mdtherm_info->nres_nhc;
   int nyosh  = mdtherm_info->nyosh_nhc;

//============================================================================
// Evolve the system : if its not the first time step

   if(itime>0){
     for(int i=natmStr;i<natmEnd;i++){
       atoms[i].vx += dt2*atoms[i].fx/atoms[i].m;
       atoms[i].vy += dt2*atoms[i].fy/atoms[i].m;
       atoms[i].vz += dt2*atoms[i].fz/atoms[i].m;
     }//endfor  
     applyNHC(natm,len_nhc,atoms,atomsNHC,dt,nyosh,nresp,natmNow,natmStr,natmEnd);
   }//endif

   computeEkin(natm,atoms,eKinetic,natmNow,natmStr,natmEnd);
   double eKineticNhc1=0.0;
   computeENHC(natm,len_nhc,atomsNHC,eKineticNhc,&eKineticNhc1,potNhc,0,
               natmNow,natmStr,natmEnd);

//----------------------------------------------------------------------------
  }//end routine
//============================================================================


//============================================================================
//        NVT Integration : 1st half
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void ATOMINTEGRATE::integrate_nvt_1st_half(int natm,int len_nhc,
                                           Atom *atoms,AtomNHC *atomsNHC,
                                           int natmNow,int natmStr,int natmEnd)
//============================================================================
    {//begin routine 
//============================================================================
// Local variables

   MDINTEGRATE  *mdintegrate  = MDINTEGRATE::get();
   GENERAL_DATA *general_data = GENERAL_DATA::get();
   MDTHERM_INFO *mdtherm_info = &(mdintegrate->mdtherm_info);
#include "../class_defs/allclass_strip_gen.h"

   double dt  = gentimeinfo->dt;
   double dt2 = dt*0.5;
   int nresp  = mdtherm_info->nres_nhc;
   int nyosh  = mdtherm_info->nyosh_nhc;

//============================================================================
// Evolve the system : 

   applyNHC(natm,len_nhc,atoms,atomsNHC,dt,nyosh,nresp,natmNow,natmStr,natmEnd);
   for(int i=natmStr;i<natmEnd;i++){
     atoms[i].vx += dt2*atoms[i].fx/atoms[i].m;
     atoms[i].vy += dt2*atoms[i].fy/atoms[i].m;
     atoms[i].vz += dt2*atoms[i].fz/atoms[i].m;
   }//endfor  

   for(int i=natmStr;i<natmEnd;i++){
     atoms[i].x += dt*atoms[i].vx;
     atoms[i].y += dt*atoms[i].vy;
     atoms[i].z += dt*atoms[i].vz;
   }//endfor  

//----------------------------------------------------------------------------
  }//end routine
//============================================================================

//============================================================================
//        NVT Integration : 2nd half
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void ATOMINTEGRATE::integrate_isonvt_2nd_half(int itime,int natm,int len_nhc,
                           Atom *atoms,AtomNHC *atomsNHC,double *eKinetic,
                           double *eKineticNhc,double *potNhc,
                           int natmNow,int natmStr,int natmEnd)
//============================================================================
    {//begin routine 
//============================================================================
// Local variables

   MDINTEGRATE  *mdintegrate  = MDINTEGRATE::get();
   GENERAL_DATA *general_data = GENERAL_DATA::get();
   MDTHERM_INFO *mdtherm_info = &(mdintegrate->mdtherm_info);
#include "../class_defs/allclass_strip_gen.h"

   double dt    = gentimeinfo->dt;
   double dt2   = dt*0.5;
   int nresp    = mdtherm_info->nres_nhc;
   int nyosh    = mdtherm_info->nyosh_nhc;
   double gamma = 0.5;

//============================================================================
// Evolve the system : if its not the first time step

   if(itime>0){
     applyIsoVel(natm,atoms,atomsNHC,dt,natmNow,natmStr,natmEnd);
     applyIsoNHC(natm,len_nhc,atoms,atomsNHC,dt,nyosh,nresp,natmNow,natmStr,natmEnd);
   }//endif
   scaleIso(natm,atoms,atomsNHC,natmNow,natmStr,natmEnd); // prevent round-off
                                                          // set initial condition

   double eKineticNhc1=0.0;
   computeEkin(natm,atoms,eKinetic,natmNow,natmStr,natmEnd);
   computeENHC(natm,len_nhc,atomsNHC,eKineticNhc,&eKineticNhc1,potNhc,1,
               natmNow,natmStr,natmEnd);
   (*eKinetic) += (gamma*eKineticNhc1);

//----------------------------------------------------------------------------
  }//end routine
//============================================================================


//============================================================================
//        NVT Integration : 1st half
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void ATOMINTEGRATE::integrate_isonvt_1st_half(int natm,int len_nhc,
                                           Atom *atoms,AtomNHC *atomsNHC,
                                           int natmNow,int natmStr,int natmEnd)
//============================================================================
    {//begin routine 
//============================================================================
// Local variables

   MDINTEGRATE  *mdintegrate  = MDINTEGRATE::get();
   GENERAL_DATA *general_data = GENERAL_DATA::get();
   MDTHERM_INFO *mdtherm_info = &(mdintegrate->mdtherm_info);
#include "../class_defs/allclass_strip_gen.h"

   double dt  = gentimeinfo->dt;
   double dt2 = dt*0.5;
   int nresp  = mdtherm_info->nres_nhc;
   int nyosh  = mdtherm_info->nyosh_nhc;

//============================================================================
// Evolve the system : 

   applyIsoNHC(natm,len_nhc,atoms,atomsNHC,dt,nyosh,nresp,natmNow,natmStr,natmEnd);
   applyIsoVel(natm,atoms,atomsNHC,dt,natmNow,natmStr,natmEnd);

   for(int i=natmStr;i<natmEnd;i++){
     atoms[i].x += dt*atoms[i].vx;
     atoms[i].y += dt*atoms[i].vy;
     atoms[i].z += dt*atoms[i].vz;
   }//endfor  

//----------------------------------------------------------------------------
  }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void ATOMINTEGRATE::computeEkin(int natm, Atom *atoms, double *eKinetic_ret,
                                int natmNow,int natmStr,int natmEnd){

   double eKinetic = 0.0;
   for(int i=natmStr;i<natmEnd;i++){
          eKinetic += atoms[i].m*(atoms[i].vx*atoms[i].vx+
                                  atoms[i].vy*atoms[i].vy+
                                  atoms[i].vz*atoms[i].vz);
   }//endfor
   eKinetic *= 0.5;

   (*eKinetic_ret) = eKinetic;

//----------------------------------------------------------------------------
  }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void ATOMINTEGRATE::computeENHC(int natm, int len_nhc,AtomNHC *atomsNHC, 
                                double *eKinNHC_ret,double *eKinNHC1_ret,
                                double *potNHC_ret,int iopt,
                                int natmNow,int natmStr,int natmEnd)
//============================================================================
  {//begin routine
//============================================================================

   double potNHC = 0.0;
   for(int i=natmStr;i<natmEnd;i++){potNHC += atomsNHC[i].posKT;}
 
   double eKinNHC1 = 0.0;
   double eKinNHC  = 0.0;
   for(int i=natmStr;i<natmEnd;i++){
     int j=0;
     eKinNHC1 += atomsNHC[i].m[j]*(atomsNHC[i].vx[j]*atomsNHC[i].vx[j]
                                  +atomsNHC[i].vy[j]*atomsNHC[i].vy[j]+
                                  +atomsNHC[i].vz[j]*atomsNHC[i].vz[j]);
     for(int j=1;j<len_nhc;j++){
        eKinNHC += atomsNHC[i].m[j]*(atomsNHC[i].vx[j]*atomsNHC[i].vx[j]
                                    +atomsNHC[i].vy[j]*atomsNHC[i].vy[j]+
                                    +atomsNHC[i].vz[j]*atomsNHC[i].vz[j]);
     }//endfor
   }//endfor
   eKinNHC  *= 0.5;
   eKinNHC1 *= 0.5;

   if(iopt==0){
     eKinNHC+=eKinNHC1;
     eKinNHC1=0.0;
   }//endif
   
   eKinNHC_ret[0]  = eKinNHC;
   eKinNHC1_ret[0] = eKinNHC1;
   potNHC_ret[0]   = potNHC;

//----------------------------------------------------------------------------
  }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void ATOMINTEGRATE::applyNHC(int natm,int len_nhc,Atom *atoms,AtomNHC *atomsNHC,
                             double dt,int nyosh,int nresp,
                             int natmNow,int natmStr,int natmEnd)
//============================================================================
  {//begin routine
//============================================================================
// Get the yoshida stuff out

  double wdti[20],wdti2[20],wdti4[20],wdti8[20];
  double dti = dt/((double)nresp);
  set_yosh(nyosh,dti,wdti,wdti2,wdti4,wdti8);

//============================================================================
// (I)Compute the forces : Precompute  m*v^2 of atoms

   for(int i=natmStr;i<natmEnd;i++){
     atoms[i].mvx2 = atoms[i].m*atoms[i].vx*atoms[i].vx;
     atoms[i].mvy2 = atoms[i].m*atoms[i].vy*atoms[i].vy;
     atoms[i].mvz2 = atoms[i].m*atoms[i].vz*atoms[i].vz;
   }//endfor
   get_forc_NHC0(natm,atoms,atomsNHC,natmNow,natmStr,natmEnd);
   for(int ic=1;ic<len_nhc;ic++){get_forc_NHC(natm,ic,atomsNHC,natmNow,natmStr,natmEnd);}

//============================================================================
// (II) Yoshida-Suzuki step yourself to heaven

  for(int iresn=1;iresn<=nresp;iresn++){
  for(int iyosh=1;iyosh<=nyosh;iyosh++){
//--------------------------------------------------------------------------
//  1) Evolve the last therm velocity in each chain                         
    evolve_vNHCM(natm,len_nhc,atomsNHC,wdti4[iyosh],natmNow,natmStr,natmEnd);
//--------------------------------------------------------------------------
//  2) Evolve the last-1 to the first therm velocity in each chain        
    for(int ic=len_nhc-2;ic>=0;ic--){
      evolve_vNHC(natm,ic,atomsNHC,wdti4[iyosh],wdti8[iyosh],natmNow,natmStr,natmEnd);
    }//endfor
//--------------------------------------------------------------------------
//  3) Evolve the particle velocities 
    evolve_vAtmNHC(natm,atoms,atomsNHC,wdti2[iyosh],natmNow,natmStr,natmEnd);
//--------------------------------------------------------------------------
//  4) Evolve the therm positions                                           
    evolve_pNHC(natm,len_nhc,atomsNHC,wdti2[iyosh],natmNow,natmStr,natmEnd);
//--------------------------------------------------------------------------
//  5) Evolve the 1 to last-1 therm velocity in each chain : get forces
    get_forc_NHC0(natm,atoms,atomsNHC,natmNow,natmStr,natmEnd);
    for(int ic=0,icp=1;ic<(len_nhc-1);ic++,icp++){
      evolve_vNHC(natm,ic,atomsNHC,wdti4[iyosh],wdti8[iyosh],natmNow,natmStr,natmEnd);
      get_forc_NHC(natm,icp,atomsNHC,natmNow,natmStr,natmEnd);
    }//endfor
//--------------------------------------------------------------------------
//  6) Evolve the last therm velocotiy in each chain                        
    evolve_vNHCM(natm,len_nhc,atomsNHC,wdti4[iyosh],natmNow,natmStr,natmEnd);
//--------------------------------------------------------------------------
  }}//endfor: iyosh,iresn

//==========================================================================
  }//end routine
//==========================================================================

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void ATOMINTEGRATE::get_forc_NHC0(int natm,Atom *atoms,AtomNHC *atomsNHC,
                                  int natmNow,int natmStr,int natmEnd){
  for(int i=natmStr;i<natmEnd;i++){
    atomsNHC[i].fx[0] = (atoms[i].mvx2-atomsNHC[i].kT)/atomsNHC[i].m[0];
    atomsNHC[i].fy[0] = (atoms[i].mvy2-atomsNHC[i].kT)/atomsNHC[i].m[0];
    atomsNHC[i].fz[0] = (atoms[i].mvz2-atomsNHC[i].kT)/atomsNHC[i].m[0];
  }//endfor
}//end routine
//==========================================================================

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void ATOMINTEGRATE::get_forc_NHC(int natm,int ic,AtomNHC *atomsNHC,
                                 int natmNow,int natmStr,int natmEnd){
  int icm = ic-1;
  for(int i=natmStr;i<natmEnd;i++){
    atomsNHC[i].fx[ic]=(atomsNHC[i].m[icm]*atomsNHC[i].vx[icm]*atomsNHC[i].vx[icm]
                       -atomsNHC[i].kT)/atomsNHC[i].m[ic];
    atomsNHC[i].fy[ic]=(atomsNHC[i].m[icm]*atomsNHC[i].vy[icm]*atomsNHC[i].vy[icm]
                       -atomsNHC[i].kT)/atomsNHC[i].m[ic];
    atomsNHC[i].fz[ic]=(atomsNHC[i].m[icm]*atomsNHC[i].vz[icm]*atomsNHC[i].vz[icm]
                       -atomsNHC[i].kT)/atomsNHC[i].m[ic];
  }//endfor
}//end routine
//==========================================================================

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void ATOMINTEGRATE::evolve_vNHCM(int natm,int len_nhc,AtomNHC *atomsNHC,
                                 double wdti4,int natmNow,int natmStr,int natmEnd)
//==========================================================================
{//begin routine
  int ic = len_nhc-1;
  for(int i=natmStr;i<natmEnd;i++){
    atomsNHC[i].vx[ic] += atomsNHC[i].fx[ic]*wdti4;
    atomsNHC[i].vy[ic] += atomsNHC[i].fy[ic]*wdti4;
    atomsNHC[i].vz[ic] += atomsNHC[i].fz[ic]*wdti4;
  }//endfor
}//end routine
//==========================================================================

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void ATOMINTEGRATE::evolve_vNHC(int natm,int ic,AtomNHC *atomsNHC,
                                double wdti4,double wdti8,
                                int natmNow,int natmStr,int natmEnd)
//==========================================================================
{//begin routine
  double argx,argy,argz;
  double aax,aay,aaz;
  int icp = ic+1;
  for(int i=natmStr;i<natmEnd;i++){
    argx = -wdti8*atomsNHC[i].vx[icp]; aax = exp(argx);  
    argy = -wdti8*atomsNHC[i].vy[icp]; aay = exp(argy);
    argz = -wdti8*atomsNHC[i].vz[icp]; aaz = exp(argz);
    atomsNHC[i].vx[ic] = atomsNHC[i].vx[ic]*aax*aax + wdti4*atomsNHC[i].fx[ic]*aax;
    atomsNHC[i].vy[ic] = atomsNHC[i].vy[ic]*aay*aay + wdti4*atomsNHC[i].fy[ic]*aay;
    atomsNHC[i].vz[ic] = atomsNHC[i].vz[ic]*aaz*aaz + wdti4*atomsNHC[i].fz[ic]*aaz;
  }//endfor
}//end routine
//==========================================================================

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void ATOMINTEGRATE::evolve_vAtmNHC(int natm,Atom *atoms,AtomNHC *atomsNHC,
                                   double wdti2,int natmNow,int natmStr,int natmEnd)
//==========================================================================
{//begin routine
  double argx,argy,argz;
  double aax,aay,aaz;
  for(int i=natmStr;i<natmEnd;i++){
    argx = -wdti2*atomsNHC[i].vx[0]; aax = exp(argx); 
    argy = -wdti2*atomsNHC[i].vy[0]; aay = exp(argy); 
    argz = -wdti2*atomsNHC[i].vz[0]; aaz = exp(argz);
    atoms[i].vx *= aax;              atoms[i].mvx2 *= (aax*aax); 
    atoms[i].vy *= aay;              atoms[i].mvy2 *= (aay*aay); 
    atoms[i].vz *= aaz;              atoms[i].mvz2 *= (aaz*aaz); 
  }//endfor
}//end routine
//==========================================================================

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void ATOMINTEGRATE::evolve_pNHC(int natm,int len_nhc, AtomNHC *atomsNHC,
                                double wdti2,int natmNow,int natmStr,int natmEnd)
//==========================================================================
{//begin routine
  for(int i=natmStr;i<natmEnd;i++){
  double pre = (atomsNHC[i].kT*wdti2);
  for(int ic=0;ic<len_nhc;ic++){
    atomsNHC[i].posKT += (pre*(atomsNHC[i].vx[ic]+atomsNHC[i].vy[ic]
                              +atomsNHC[i].vz[ic]));
             
  }}//endfor
}//end routine
//==========================================================================

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void ATOMINTEGRATE::set_yosh(int nyosh,double dt,double *wdt,double *wdt2,
                             double *wdt4,double *wdt8)
//========================================================================
  {//begin routine
//========================================================================
//             Local variable declarations                                

   double temp,p2,onethird;

//========================================================================
// Find the yoshida steps you want

  switch(nyosh){
    case 1: 
          wdt[1] = 1.0;
         break;
    case 3:
          onethird = 1.0/3.0;
          temp    = pow(2.0,onethird);
          wdt[1] =  1.0/(2.0-temp);
          wdt[2] = -temp/(2.0-temp);
          wdt[3] =  1.0/(2.0-temp);
         break;
    case 5:
          onethird = 1.0/3.0;
          temp = pow(4.0,onethird);
          p2    = 1.0/(4.0-temp);
          wdt[1]  = p2;
          wdt[2]  = p2;
          wdt[3]  = 1.0-4.0*p2;
          wdt[4]  = p2;
          wdt[5]  = p2;
         break;
    case 7:
          wdt[1] =  0.784513610477560;
          wdt[2] =  0.235573213359357;
          wdt[3] = -1.17767998417887;
          wdt[4] =  1.0 - 2.0*(wdt[1]+wdt[2]+wdt[3]);
          wdt[5] = -1.17767998417887;
          wdt[6] =  0.235573213359357;
          wdt[7] =  0.784513610477560;
         break;
    case 9:
          wdt[1] =  0.192;
          wdt[2] =  0.554910818409783619692725006662999;
          wdt[3] =  0.124659619941888644216504240951585;
          wdt[4] = -0.843182063596933505315033808282941;
          wdt[5] =  1.0 - 2.0*(wdt[1]+wdt[2]+wdt[3]+wdt[4]);
          wdt[6] = -0.843182063596933505315033808282941;
          wdt[7] =  0.124659619941888644216504240951585;
          wdt[8] =  0.554910818409783619692725006662999;
          wdt[9] =  0.192;
         break;
  }//switch

//========================================================================
// Scale the yosida steps by the time step

  for(int i=1;i<=nyosh;i++){
    wdt[i]  = wdt[i]*dt;
    wdt2[i] = wdt[i]*0.5;
    wdt4[i] = wdt[i]*0.25;
    wdt8[i] = wdt[i]*0.125;
  }//endfor

//-------------------------------------------------------------------------
   }// end routine
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void ATOMINTEGRATE::scaleIso(int natm,Atom *atoms,AtomNHC *atomsNHC,
                             int natmNow,int natmStr,int natmEnd){
//==========================================================================

   double gamma = 0.5;
   for(int i=natmStr;i<natmEnd;i++){
      double kinAtm,kinNhc,scale;
      double kT = atomsNHC[i].kT;

      kinAtm = atoms[i].m*atoms[i].vx*atoms[i].vx;
      kinNhc = atomsNHC[i].m[0]*atomsNHC[i].vx[0]*atomsNHC[i].vx[0];
      scale  = sqrt(kT/(kinAtm+gamma*kinNhc));
      atoms[i].vx       *= scale;
      atomsNHC[i].vx[0] *= scale;

      kinAtm = atoms[i].m*atoms[i].vy*atoms[i].vy;
      kinNhc = atomsNHC[i].m[0]*atomsNHC[i].vy[0]*atomsNHC[i].vy[0];
      scale  = sqrt(kT/(kinAtm+gamma*kinNhc));
      atoms[i].vy       *= scale;
      atomsNHC[i].vy[0] *= scale;

      kinAtm = atoms[i].m*atoms[i].vz*atoms[i].vz;
      kinNhc = atomsNHC[i].m[0]*atomsNHC[i].vz[0]*atomsNHC[i].vz[0];
      scale  = sqrt(kT/(kinAtm+gamma*kinNhc));
      atoms[i].vz       *= scale;
      atomsNHC[i].vz[0] *= scale;
   }//endfor

//-------------------------------------------------------------------------
   }// end routine
//==========================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void ATOMINTEGRATE::applyIsoVel(int natm,Atom *atoms,AtomNHC *atomsNHC,
                                double dt,int natmNow,int natmStr,int natmEnd)
//============================================================================
  {//begin routine
//============================================================================
// Scale initial conditions for safety 

   double ffdot,pfdot,arg,s,sdot,ffdot_r;
   double gkt = atomsNHC[0].kT;
   double dt2 = dt*0.5;

//============================================================================
// x evolution

    for(int i=natmStr;i<natmEnd;i++){
      ffdot = (atoms[i].fx*atoms[i].fx)/(atoms[i].m*gkt);
      pfdot = (atoms[i].fx*atoms[i].vx)/gkt;
      ffdot_r = sqrt(ffdot);
      arg     = ffdot_r*dt2;
      if(arg>1.0e-5){
        s     = (pfdot/ffdot)*(cosh(arg)-1.0) + (1.0/ffdot_r)*sinh(arg);
        sdot  = (pfdot/ffdot_r)*sinh(arg) + cosh(arg);
      }else{
        s     = (((ffdot*pfdot/24.0*dt2+ffdot/6.0)*dt2+0.50*pfdot)*dt2+1.0)*dt2;
        sdot  = ((ffdot*pfdot/6.0*dt2+ffdot/2.0)*dt2+pfdot)*dt2+1.0;
      }//endif
      atoms[i].vx       = (atoms[i].vx + atoms[i].fx*s/atoms[i].m)/sdot;
      atomsNHC[i].vx[0] = atomsNHC[i].vx[0]/sdot;
    }//endfor

//============================================================================
// y evolution

    for(int i=natmStr;i<natmEnd;i++){
      ffdot   = (atoms[i].fy*atoms[i].fy)/(atoms[i].m*gkt);
      pfdot   = (atoms[i].fy*atoms[i].vy)/gkt;
      ffdot_r = sqrt(ffdot);
      arg     = ffdot_r*dt2;
      if(arg>1.0e-5){
        s     = (pfdot/ffdot)*(cosh(arg)-1.0) + (1.0/ffdot_r)*sinh(arg);
        sdot  = (pfdot/ffdot_r)*sinh(arg) + cosh(arg);
      }else{
        s     = (((ffdot*pfdot/24.0*dt2+ffdot/6.0)*dt2+0.50*pfdot)*dt2+1.0)*dt2;
        sdot  = ((ffdot*pfdot/6.0*dt2+ffdot/2.0)*dt2+pfdot)*dt2+1.0;
      }//endif
      atoms[i].vy       = (atoms[i].vy + atoms[i].fy*s/atoms[i].m)/sdot;
      atomsNHC[i].vy[0] = atomsNHC[i].vy[0]/sdot;
    }//endfor

//============================================================================
// z evolution

    for(int i=natmStr;i<natmEnd;i++){
      ffdot   = (atoms[i].fz*atoms[i].fz)/(atoms[i].m*gkt);
      pfdot   = (atoms[i].fz*atoms[i].vz)/gkt;
      ffdot_r = sqrt(ffdot);
      arg     = ffdot_r*dt2;
      if(arg>1.0e-5){
        s     = (pfdot/ffdot)*(cosh(arg)-1.0) + (1.0/ffdot_r)*sinh(arg);
        sdot  = (pfdot/ffdot_r)*sinh(arg) + cosh(arg);
      }else{
        s     = (((ffdot*pfdot/24.0*dt2+ffdot/6.0)*dt2+0.50*pfdot)*dt2+1.0)*dt2;
        sdot  = ((ffdot*pfdot/6.0*dt2+ffdot/2.0)*dt2+pfdot)*dt2+1.0;
      }//endif
      atoms[i].vz       = (atoms[i].vz + atoms[i].fz*s/atoms[i].m)/sdot;
      atomsNHC[i].vz[0] = atomsNHC[i].vz[0]/sdot;
    }//endfor

//============================================================================
  }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void ATOMINTEGRATE::applyIsoNHC(int natm,int len_nhc,Atom *atoms,AtomNHC *atomsNHC,
                                double dt,int nyosh,int nresp,
                                int natmNow,int natmStr,int natmEnd)
//============================================================================
  {//begin routine
//============================================================================
// Get the yoshida stuff out

  double wdti[20],wdti2[20],wdti4[20],wdti8[20];
  double dti = dt/((double)nresp);
  set_yosh(nyosh,dti,wdti,wdti2,wdti4,wdti8);

//============================================================================
// (I)Compute the forces : Precompute  m*v^2 of atoms

   for(int i=natmStr;i<natmEnd;i++){
     atomsNHC[i].fx[0]=0.0;
     atomsNHC[i].fy[0]=0.0;
     atomsNHC[i].fz[0]=0.0;
   }//endfor
   for(int ic=1;ic<len_nhc;ic++){get_forc_NHC(natm,ic,atomsNHC,natmNow,natmStr,natmEnd);}

//============================================================================
// (II) Yoshida-Suzuki step yourself to heaven

  for(int iresn=1;iresn<=nresp;iresn++){
  for(int iyosh=1;iyosh<=nyosh;iyosh++){
//--------------------------------------------------------------------------
//  1) Evolve the last therm velocity in each chain                         
    evolve_vNHCM(natm,len_nhc,atomsNHC,wdti4[iyosh],natmNow,natmStr,natmEnd);
//--------------------------------------------------------------------------
//  2) Evolve the last-1 to the first therm velocity in each chain        
    for(int ic=len_nhc-2;ic>=1;ic--){
      evolve_vNHC(natm,ic,atomsNHC,wdti4[iyosh],wdti8[iyosh],natmNow,natmStr,natmEnd);
    }//endfor
//--------------------------------------------------------------------------
//  3) Evolve the particle velocities 
    evolve_vAtmIsoNHC(natm,atoms,atomsNHC,wdti2[iyosh],natmNow,natmStr,natmEnd);
//--------------------------------------------------------------------------
//  4) Evolve the therm positions                                           
    evolve_pIsoNHC(natm,len_nhc,atomsNHC,wdti2[iyosh],natmNow,natmStr,natmEnd);
//--------------------------------------------------------------------------
//  5) Evolve the particle velocities 
    evolve_vAtmIsoNHC(natm,atoms,atomsNHC,wdti2[iyosh],natmNow,natmStr,natmEnd);
//--------------------------------------------------------------------------
//  6) Evolve the 1 to last-1 therm velocity in each chain : get forces
    get_forc_NHC(natm,1,atomsNHC,natmNow,natmStr,natmEnd);
    for(int ic=1,icp=2;ic<(len_nhc-1);ic++,icp++){
      evolve_vNHC(natm,ic,atomsNHC,wdti4[iyosh],wdti8[iyosh],natmNow,natmStr,natmEnd);
      get_forc_NHC(natm,icp,atomsNHC,natmNow,natmStr,natmEnd);
    }//endfor
//--------------------------------------------------------------------------
//  7) Evolve the last therm velocotiy in each chain                        
    evolve_vNHCM(natm,len_nhc,atomsNHC,wdti4[iyosh],natmNow,natmStr,natmEnd);
//--------------------------------------------------------------------------
  }}//endfor: iyosh,iresn

//==========================================================================
  }//end routine
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void ATOMINTEGRATE::evolve_pIsoNHC(int natm,int len_nhc,AtomNHC *atomsNHC,double wdti2,
                                   int natmNow,int natmStr,int natmEnd){
//==========================================================================

  for(int i=natmStr;i<natmEnd;i++){
    atomsNHC[i].posKT -= (wdti2*atomsNHC[i].m[0]*(
                          atomsNHC[i].vx[0]*atomsNHC[i].vx[0]*atomsNHC[i].vx[1]
                         +atomsNHC[i].vy[0]*atomsNHC[i].vy[0]*atomsNHC[i].vy[1]
                         +atomsNHC[i].vz[0]*atomsNHC[i].vz[0]*atomsNHC[i].vz[1]));
  }//endfor

  double pre = (atomsNHC[0].kT*wdti2);
  for(int i=natmStr;i<natmEnd;i++){
  for(int ic=1;ic<len_nhc;ic++){
    atomsNHC[i].posKT += pre*(atomsNHC[i].vx[ic]
                             +atomsNHC[i].vy[ic]
                             +atomsNHC[i].vz[ic]);
  }}//endfor

//==========================================================================
  }//end routine
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void ATOMINTEGRATE::evolve_vAtmIsoNHC(int natm,Atom *atoms,AtomNHC *atomsNHC,double wdti2,
                                      int natmNow,int natmStr,int natmEnd){
//==========================================================================

  double wdti4 = wdti2*0.5;
  double gamma = 0.5;
  double gkt   = atomsNHC[0].kT;

//==========================================================================

  for(int i=natmStr;i<natmEnd;i++){
    double aa   = exp(-atomsNHC[i].vx[1]*wdti4);
    double temp = atoms[i].m*atoms[i].vx*atoms[i].vx
                 +atomsNHC[i].m[0]*atomsNHC[i].vx[0]*atomsNHC[i].vx[0]*gamma*aa*aa;
    double s    = sqrt(gkt/temp);
    atomsNHC[i].vx[0] *= (s*aa);
    atoms[i].vx       *= s;
  }//endfor

//==========================================================================

  for(int i=natmStr;i<natmEnd;i++){
    double aa   = exp(-atomsNHC[i].vy[1]*wdti4);
    double temp = atoms[i].m*atoms[i].vy*atoms[i].vy
                 +atomsNHC[i].m[0]*atomsNHC[i].vy[0]*atomsNHC[i].vy[0]*gamma*aa*aa;
    double s    = sqrt(gkt/temp);
    atomsNHC[i].vy[0] *= (s*aa);
    atoms[i].vy       *= s;
  }//endfor

//==========================================================================

  for(int i=natmStr;i<natmEnd;i++){
    double aa   = exp(-atomsNHC[i].vz[1]*wdti4);
    double temp = atoms[i].m*atoms[i].vz*atoms[i].vz
                +atomsNHC[i].m[0]*atomsNHC[i].vz[0]*atomsNHC[i].vz[0]*gamma*aa*aa;
    double s    = sqrt(gkt/temp);
    atomsNHC[i].vz[0] *= (s*aa);
    atoms[i].vz       *= s;
  }//endfor

//==========================================================================
  }//end routine
//==========================================================================
