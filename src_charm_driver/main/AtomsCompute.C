//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
/** \file AtomsCompute.C
 *          Processor group class Functions : Atoms
 */
//==============================================================================
#include "AtomsCache.h"
#include "AtomsCompute.h"
#include "energyGroup.h"
#include "eesCache.h"
#include "cp_state_ctrl/CP_State_GSpacePlane.h"
#include "fft_slab_ctrl/fftCacheSlab.h"
#include "utility/util.h"
#include "CPcharmParaInfoGrp.h"
#include "load_balance/IntMap.h"
#include "charm++.h"
#include "src_piny_physics_v1.0/include/class_defs/Interface_ctrl.h"
#include "PhysScratchCache.h"

#include <cmath>

//#include "CP_State_Plane.h"

//----------------------------------------------------------------------------
#define CHARM_ON
#include "src_piny_physics_v1.0/include/class_defs/piny_constants.h"
#include "src_piny_physics_v1.0/include/class_defs/ATOM_OPERATIONS/class_atomintegrate.h"
#include "src_piny_physics_v1.0/include/class_defs/ATOM_OPERATIONS/class_atomoutput.h"
#include "src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cprspaceion.h"

//----------------------------------------------------------------------------

extern CkVec <CProxy_PIBeadAtoms>       UPIBeadAtomsProxy;
extern CkVec <IntMap2on2> GSImaptable;
extern CkVec <CProxy_EnergyGroup>          UegroupProxy;
extern Config                      config;
extern CkVec <CProxy_CP_State_GSpacePlane> UgSpacePlaneProxy;
extern CkVec <CProxy_GSpaceDriver>         UgSpaceDriverProxy;
extern CkVec <CProxy_AtomsCache>             UatomsCacheProxy;
extern CkVec <CProxy_AtomsCompute>             UatomsComputeProxy;
extern CkVec <CProxy_EnergyGroup>          UegroupProxy;
extern CkVec <CProxy_StructFactCache>      UsfCacheProxy;
extern CkVec <CProxy_eesCache>             UeesCacheProxy;
extern CProxy_TemperController temperControllerProxy;
extern CProxy_InstanceController instControllerProxy;
extern CProxy_PhysScratchCache pScratchProxy;

//----------------------------------------------------------------------------

//#define _CP_DEBUG_PSI_OFF_
//#define _CP_ENERGY_GRP_VERBOSE_
//#define _CP_DEBUG_ATMS_
//#define _CP_DEBUG_ATMS_EXIT_

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
/**
 *
 * @addtogroup Atoms
 * @{ 
 *
 */
//==============================================================================
AtomsCompute::AtomsCompute(int n, int n_nl, int len_nhc_, int iextended_on_,int cp_min_opt_,
			   int cp_wave_opt_, int cp_bomd_opt_, int isokin_opt_,int cp_grimme_,
                           double kT_, Atom* a, AtomNHC *aNHC, int nChareAtoms_,
			   UberCollection _thisInstance) : natm(n), natm_nl(n_nl), len_nhc(len_nhc_), iextended_on(iextended_on_), cp_min_opt(cp_min_opt_), cp_wave_opt(cp_wave_opt_), cp_bomd_opt(cp_bomd_opt_), isokin_opt(isokin_opt_),cp_grimme(cp_grimme_), kT(kT_), nChareAtoms(nChareAtoms_), thisInstance(_thisInstance) 
                                    //==============================================================================
{// begin routine
  //==============================================================================
  // parameters, options and energies
  handleForcesCount=0;
  ktemps          = 0;
  pot_ewd_rs      = 0.0;
  potGrimmeVdw   = 0.0;
  eKinetic        = 0.0;
  eKineticNhc     = 0.0;
  potNhc          = 0.0;    
  potPIMDChain    = 0.0;
  countAtm        = 0;
  acceptCountfu   = 0;
  acceptCountX    = 0;
  acceptCountu    = 0;
  atomsCMrecv=atomsPIMDXrecv=false;
  temperScreenFile = (UatomsCacheProxy[thisInstance.proxyOffset].ckLocalBranch())->temperScreenFile;

  //==============================================================================
  // Initial positions, forces, velocities 
  numPIMDBeads    = config.UberImax;
  PIBeadIndex     = thisInstance.idxU.x;
  TemperIndex     = thisInstance.idxU.z;
  atoms           = new Atom[natm];
  atomsNHC        = new AtomNHC[natm];
  AtomsCache *ag         = UatomsCacheProxy[thisInstance.proxyOffset].ckLocalBranch();
  // we keep our own copy of this annoying thing to avoid data races
  // with the force accumulation in the AtomsCache.
  fastAtoms.q    = new double[natm];
  fastAtoms.x    = new double[natm];
  fastAtoms.y    = new double[natm];
  fastAtoms.z    = new double[natm];
  fastAtoms.fx   = new double[natm];
  fastAtoms.fy   = new double[natm];
  fastAtoms.fz   = new double[natm];
  fastAtoms.fxu   = new double[natm];
  fastAtoms.fyu   = new double[natm];
  fastAtoms.fzu   = new double[natm];
  for(int i=0;i<natm;i++){
    fastAtoms.x[i]  = a[i].x;
    fastAtoms.y[i]  = a[i].y;
    fastAtoms.z[i]  = a[i].z;
    fastAtoms.fx[i] = a[i].fx;
    fastAtoms.fy[i] = a[i].fy;
    fastAtoms.fz[i] = a[i].fz;
  }//endfor
  iteration = &(ag->iteration);

  CmiMemcpy(atoms, a, natm * sizeof(Atom));           // atoms has no vectors
  for(int i=0;i<natm;i++){atomsNHC[i].Init(&aNHC[i]);}
  ftot = (double *)fftw_malloc((3*natm+3)*sizeof(double));

  zeroforces();
  if(iextended_on==1 && (cp_min_opt==0 || cp_bomd_opt==1)){
    zeronhc();
  }//endif

  //==============================================================================
  // A copy of the atoms for fast routines : only need to copy charge once

  copyFastToSlow();
  double *qq = fastAtoms.q;
  for(int i=0;i<natm;i++){qq[i]=atoms[i].q;}

  //==============================================================================
  // Number of messages to be received when atoms are moved : simple parallel

  int nproc = nChareAtoms;
  int div   = (natm / nproc);
  int rem   = (natm % nproc);

  nAtmMsgRecv = 0;
  for(int myid=0;myid<nproc;myid++){
    int natmNow = div;
    if(myid< rem){natmNow++;}
    if(natmNow>0){nAtmMsgRecv++;}
  }//endfor

  //==============================================================================
  // PIMD set up : Even if classical, this can run. It is harmless and
  //               prevents careless bugs

  massPIMDScal = (double *)fftw_malloc(numPIMDBeads*sizeof(double));
  initPIMD();

  //-----------------------------------------------------------------------------
}//end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
/** Destructor
 *
 *
 **/
//==============================================================================
AtomsCompute::~AtomsCompute(){
  delete [] atoms;
  delete [] atomsNHC;
  delete [] fastAtoms.q;
  delete [] fastAtoms.x;
  delete [] fastAtoms.y;
  delete [] fastAtoms.z;
  delete [] fastAtoms.fx;
  delete [] fastAtoms.fy;
  delete [] fastAtoms.fz;
  delete [] fastAtoms.fxu;
  delete [] fastAtoms.fyu;
  delete [] fastAtoms.fzu;
  fftw_free(PIMD_CM_Atoms.x);
  fftw_free(PIMD_CM_Atoms.y);
  fftw_free(PIMD_CM_Atoms.z);
  fftw_free(ftot);
}
//==============================================================================


//==========================================================================
// Bead root is the 0th element.  The AtomsCompute map will lookup the
// GSP map to ensure that the 0th AtomCompute is co-mapped with GSP(0,0) in
// each bead instance.  

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

void AtomsCompute::init()
{
  // Bead root should be the 0th element.  The AtomsCompute map will use
  // the GSP map to ensure that the 0th AtomCompute is co-mapped
  // with GSP(0,0) in each bead instance.
  const int offset=thisInstance.getPO();
  amBeadRoot      = (CkMyPe()==GSImaptable[offset].get(0,0));
  amZerothBead    = (thisInstance.idxU.x==0);

  if(numPIMDBeads>1)
  {
    // make lists of Array IDs and their Commander PE
    CkArrayID *atomsArrayids= new CkArrayID[numPIMDBeads];
    CkArrayIndex **elems  = new CkArrayIndex*[numPIMDBeads];
    CkArrayIndex **elemsAll  = new CkArrayIndex*[numPIMDBeads];
    int *naelems = new int[numPIMDBeads];
    int *naelemsAll = new int[numPIMDBeads];
    for(int i=0;i<numPIMDBeads;i++){
      elems[i]= new CkArrayIndex1D[1];
      elemsAll[i]= new CkArrayIndex1D[nChareAtoms];
      UberCollection instance=thisInstance;
      naelems[i]=1;
      naelemsAll[i]=nChareAtoms;
      instance.idxU.x=i;
      int offset=instance.calcPO();
      elems[i][0]=CkArrayIndex1D(0);
      for(int j=0;j<nChareAtoms;j++)
        elemsAll[i][j]=CkArrayIndex1D(j);
      atomsArrayids[i]=UatomsComputeProxy[offset].ckGetArrayID();
      //CkPrintf("{%d}[%d] AtomsCompute::init elems[%d][0]=%d, atomsgrpids[%d]=%d amBeadRoot=%d\n",thisInstance.proxyOffset, CkMyPe(), i, elems[i][0], i, atomsgrpids[i], amBeadRoot);     
    }//endfor
    //      proxyHeadBeads=CProxySection_AtomsCompute(numPIMDBeads, atomsArrayids, elems, naelems);
    proxyAllBeads=CProxySection_AtomsCompute(numPIMDBeads, atomsArrayids, elemsAll, naelemsAll);
    delete [] naelems;
    delete [] naelemsAll;
    for(int i=0;i<numPIMDBeads;i++){delete [] elems[i]; delete [] elemsAll[i];}
    delete [] elemsAll;
    delete [] atomsArrayids;

  }//endif

  //==============================================================================
}//end routine
//==============================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
/** recvContribute
 * Every Array member is sent all the forces at present.
 */
//==========================================================================
void AtomsCompute::recvContribute(CkReductionMsg *msg) {
  //============================================================
#ifdef _CP_DEBUG_ATMS_
  CkPrintf("{%d}[%d] AtomsCompute::recvContribute count %d\n ", thisInstance.proxyOffset, thisIndex, handleForcesCount);     
#endif
  //==========================================================================
  // Local pointers

  int i,j;
  int myid        = thisIndex;
  contribMsg[handleForcesCount++]=msg;
  EnergyGroup *eg = UegroupProxy[thisInstance.proxyOffset].ckLocalBranch();

  //============================================================
  // Copy out the reduction of energy and forces and nuke msg
  double *ftot    = (double *) msg->getData();

  double pot_ewd_rs_loc   = ftot[3*natm];
  double potPerdCorrLoc   = ftot[3*natm+1];
  double potGrimmeVdwLoc  = ftot[3*natm+2];

  //============================================================
  // Tuck things that can be tucked.

  eg->estruct.eewald_real = pot_ewd_rs_loc;  
  eg->estruct.grimmeVdw   = potGrimmeVdwLoc;  


#ifdef _CP_DEBUG_ATMS_EXIT_
  if(myid==0){CkExit();}
#endif
  if(handleForcesCount==2)
    handleForces();

  //==========================================================================
}//end routine
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
/** recvContributeForces
 * Every Array member has all the forces at present.
 */
//==========================================================================
void AtomsCompute::recvContributeForces(CkReductionMsg *msg) {
  //============================================================
#ifdef _CP_DEBUG_ATMS_
  CkPrintf("{%d}[%d] AtomsCompute::recvContributeForces count %d\n ", thisInstance.proxyOffset, thisIndex, handleForcesCount);     
#endif
  //==========================================================================
  // Local pointers
  int myid        = thisIndex;
  contribMsg[handleForcesCount++]=msg;

#ifdef _CP_DEBUG_ATMS_EXIT_
  if(myid==0){CkExit();}
#endif

  if(handleForcesCount==2)
    handleForces();
  //=================================================================
}//end routine
//=================================================================

//=================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=================================================================
void AtomsCompute::handleForces()  {
  //=================================================================
#ifdef _CP_DEBUG_ATMS_
  CkPrintf("{%d}[%d] AtomsCompute::handleForces \n ", thisInstance.proxyOffset, thisIndex);     
#endif
  //=================================================================
  AtomsCache *ag         = UatomsCacheProxy[thisInstance.proxyOffset].ckLocalBranch();
  int iteration = ag->iteration;
  //=================================================================
  copyFastToSlow();
  zeroforces(); // we're now getting everyone's forces
  double fmag=0.0;
  handleForcesCount=0;
  EnergyGroup *eg = UegroupProxy[thisInstance.proxyOffset].ckLocalBranch();
  // Model forces  : This is fine for path integral checking too

  // sum the contributions
  double *ftot0    = (double *) contribMsg[0]->getData();
  double *ftot1    = (double *) contribMsg[1]->getData();
  // no loop carried dependencies here, so a compiler worth its salt
  // should easily vectorize these ops.
  int i, j;
  for(i=0,j=0;i<natm;i++,j+=3){
    atoms[i].fx = ftot0[j] + ftot1[j]; 
    atoms[i].fy = ftot0[j+1] + ftot1[j+1];
    atoms[i].fz = ftot0[j+2] + ftot1[j+2];
#ifdef _CP_DEBUG_ATMS_
    if(CkMyPe()==0)
    {
      CkPrintf("AtomCompute handleForces forces %d %.5g,%.5g,%.5g\n",i, atoms[i].fx, atoms[i].fy, atoms[i].fz);
    }
#endif
    fmag += atoms[i].fx * atoms[i].fx + atoms[i].fy * atoms[i].fy + atoms[i].fz * atoms[i].fz;
  }// endfor

  fmag /= (double)(3*natm);
  fmag  = sqrt(fmag);
  eg->estruct.fmag_atm    = fmag;
  delete contribMsg[0];
  delete contribMsg[1];

#ifdef  _CP_DEBUG_PSI_OFF_
  double omega    = (0.0241888/15.0); // 15 fs^{-1}
double omega2   = omega*omega;
int npts        = 200;
if(iteration==0){
  px = new double *[natm];
  for(i =0;i<natm;i++){
    px[i] = new double[npts];
    for(j=0;j<npts;j++){px[i][j]=0.0;}
  }//endfor
  for(i=0;i<natm;i++){
    atoms[i].x = 0.0;
    atoms[i].y = 0.0;
    atoms[i].z = 0.0;
  }//endfor
}//endif

for(i=0;i<natm;i++){
  atoms[i].fx = -omega2*atoms[i].x*atoms[i].m;
  atoms[i].fy = -omega2*atoms[i].y*atoms[i].m;;
  atoms[i].fz = -omega2*atoms[i].z*atoms[i].m;;
}//endfor
#endif

//=================================================================
// Check the forces

#ifdef _GLENN_CHECK_FORCES
FILE *fp = fopen("forces.save","w");
for(i=0;i<natm;i++){
  fprintf(fp,"%.12g %.12g %.12g\n",atoms[i].fx,atoms[i].fy,atoms[i].fz);
}
fclose(fp);
CkExit();
#endif


//==========================================================================
// if classical go on to integration, otherwise Fx -> Fu

if(numPIMDBeads>1){
    eg->estruct.totalpotPIMDChain=0;  
  if(iteration==0){
    send_PIMD_Fx_and_x(); // atom integration must wait on both transforms
  }else{
    send_PIMD_Fx();   // atom integration must wait on Fx->Fu transformation
  }//endif
}else{
  integrateAtoms(); // Fx is all you need so go forth and integrate
}//endif

//==========================================================================
}//end routine
//==========================================================================

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
/**
 * Integrate atoms.  This is parallelized so that a subset of the atoms are 
 * computed on each processor and their results sent to AtomCompute->acceptAtoms(). 
 */
//==========================================================================
void AtomsCompute::integrateAtoms(){
  //============================================================
#ifdef _CP_DEBUG_ATMS_
  CkPrintf("GJM_DBG: Before atom integrate %d : %d\n",thisIndex,natm);
#endif
  //============================================================
  // Local pointers

  int  mybead = thisInstance.idxU.x+1;
  int myid    = thisIndex;

  EnergyGroup *eg     = UegroupProxy[thisInstance.proxyOffset].ckLocalBranch();

  //=================================================================
  // Atom parallelization: This decompostion is done dynamically because
  // its a pretty trivial choice and it allows us to support atom
  // migration into and out of the QM region as long that atom migration
  // mechanism updates natm and the atom array accordingly.

  int div     = (natm / nChareAtoms);
  int rem     = (natm % nChareAtoms);
  int natmStr = div*myid;
  int natmNow = div;
  if(myid>=rem){natmStr += rem;}
  if(myid< rem){natmStr += myid;}
  if(myid< rem){natmNow++;}
  int natmEnd = natmNow+natmStr;

  //============================================================
  // Zero some local copies energies and flags

  double eKinetic_loc     = 0.0;
  double eKineticNhc_loc  = 0.0;
  double potNhc_loc       = 0.0;
  double potPIMDChain_loc = 0.0;
  int iwrite_atm          = 0;
  int myoutput_on         = 0;

  //============================================================
  // DEBUGGING : Compute the distribution function for model

#ifdef  _CP_DEBUG_PSI_OFF_
  if(nChareAtoms>1 && numPIMDBeads>1){
    CkPrintf("Harmonic oscillator debug test currently broken for PIMD\n");
    CkPrintf("Need a reduction over beads of the chain energy\n");
    CkExit();
  }//endif
  double omega    = (0.0241888/15.0); // 15 fs^{-1} Ok for PIMD
  double omega2   = omega*omega;
  int npts        = 200;
  double sigma;
  if(numPIMDBeads==1){
    sigma    = 1.0/sqrt(kT*atoms[0].m*omega2);
  }else{
    sigma    = 1.0/sqrt(2.0*atoms[0].m*omega); // Assuming lots of beads
  }//endif

  double dx       = 6.0*sigma/(double)npts;
  for(int i=natmStr;i<natmEnd;i++){
    int i1 = (atoms[i].x+3.0*sigma)/dx;
    int i2 = (atoms[i].y+3.0*sigma)/dx;
    int i3 = (atoms[i].z+3.0*sigma)/dx;
    i1     = (i1 >= 0    ? i1 : 0);
    i2     = (i2 >= 0    ? i2 : 0);
    i3     = (i3 >= 0    ? i3 : 0);
    i1     = (i1 < npts ? i1 : npts-1);
    i2     = (i2 < npts ? i2 : npts-1);
    i3     = (i3 < npts ? i3 : npts-1);
    px[i][i1] += 1.0;
    px[i][i2] += 1.0;
    px[i][i3] += 1.0;
  }//endif

  if(((*iteration+1) % 1000)==0){
    for(int i=natmStr;i<natmEnd;i++){
      char fname[100];
      sprintf(fname,"atmtest.dat.%d",i);
      double anorm = sqrt(2.0/(sigma*sigma*acos(-1.0)));
      FILE *fp = fopen(fname,"w");
      for(int j =1;j<npts-1;j++){
        double x   = ((double)(2*j+1))*0.5*dx - 3.0*sigma;
        double ans = anorm*exp(-0.5*x*x/(sigma*sigma));
        double cnt = 3.0*((double)(*iteration+1));
        fprintf(fp,"%g %g %g\n",x,(px[i][j]*anorm/cnt),ans);
      }//endfor       
      fclose(fp);
    }//endfor
  }//endif

  double pot_harm = 0.0;
  for(int i=natmStr;i<natmEnd;i++){
    pot_harm   += (atoms[i].m*omega2*(atoms[i].x*atoms[i].x+
          atoms[i].y*atoms[i].y+
          atoms[i].z*atoms[i].z));
  }//endfor
  pot_harm *= 0.5;
#endif

  //============================================================
  // Move atoms
  int move_atoms = 0;
  if(cp_min_opt==0 && cp_wave_opt==0) { move_atoms = 1; }
  if(cp_bomd_opt==1 && tol_reached==1) { move_atoms = 1; }

  //============================================================
  // Path integral :  Overwrite variables to keep integrator clean:
  //                  Add the chain force to transformed forces.
  //                  Scale the masses to fict bead masses
  //                  This is PIMD Respa implementation friendly 

  if(numPIMDBeads>1 && move_atoms && natmNow>0){
    switchPIMDBeadForceMass(mybead,natmStr,natmEnd,&potPIMDChain_loc);
  }//endif

  //============================================================
  // Just a little debug, early in the computation, beats a cup of coffee
  /*#define _DEBUG_PIMD_TRANSFORM_*/
#ifdef _DEBUG_PIMD_TRANSFORM_
   if(numPIMDBeads>1 && move_atoms && natmNow>0){
      AtomsCache *ag = UatomsCacheProxy[thisInstance.proxyOffset].ckLocalBranch();
      int iter_now = ag->iteration;
      char fname[100];
      sprintf(fname,"BeadTransformTest_Bead.%d_Atm.%d_%d_Iter.%d",mybead,natmStr,natmEnd,iter_now);
      FILE *fp = fopen(fname,"w");
      for(int ii=natmStr; ii<natmEnd;ii++){
        fprintf(fp,"%g %g %g %g %g %g %g %g %g %g %g %g\n",
                    atoms[ii].x  ,atoms[ii].y  ,atoms[ii].z,
                    atoms[ii].xu ,atoms[ii].yu ,atoms[ii].zu,
                    atoms[ii].fx ,atoms[ii].fy ,atoms[ii].fz,
   		    atoms[ii].fxu,atoms[ii].fyu,atoms[ii].fzu);
      }//endfor
      fclose(fp);
   }//endif
#endif
  //============================================================
  // Path integral :  Overwrite variables to keep integrator clean:
  //                  Add the chain force to transformed forces.
  //                  Scale the masses to fict bead masses
  //                  This is PIMD Respa implementation friendly 

  if(numPIMDBeads>1 && cp_min_opt==0 && cp_wave_opt==0 && natmNow>0){
    switchPIMDBeadForceMass(mybead,natmStr,natmEnd,&potPIMDChain_loc);
  }//endif
  //============================================================
  // Integrate the atoms : Path Integral Ready

#ifndef  _CP_DEBUG_SCALC_ONLY_ 
  ATOMINTEGRATE::ctrl_atom_integrate(*iteration,natm,len_nhc,cp_min_opt,
      cp_wave_opt,cp_bomd_opt,tol_reached,iextended_on,atoms,atomsNHC,myid,
      &eKinetic_loc,&eKineticNhc_loc,&potNhc_loc,&iwrite_atm,
      myoutput_on,natmNow,natmStr,natmEnd,mybead);
#endif

  //============================================================
  // Path integral :  Rescale to the physical masses

  if(numPIMDBeads>1 && move_atoms==1 && natmNow>0){
    unswitchPIMDMass(mybead,natmStr,natmEnd);
  }//endif

  //============================================================
  // Debug output :  Not correct for PIMD which needs chain PE

#ifdef  _CP_DEBUG_PSI_OFF_
  double etot_atm;
  if(isokin_opt==0){
    etot_atm = eKinetic_loc+eKineticNhc_loc+potNhc_loc+pot_harm+potPIMDChain_loc;
    CkPrintf("iteration %d : tot class energy %.12g on %d\n",*iteration,etot_atm,myid);
  }else{
    etot_atm = eKineticNhc_loc+potNhc_loc;
    CkPrintf("iteration %d : tot class energy %.12g %.12g on %d\n",*iteration,etot_atm,
        (eKinetic_loc+pot_harm),myid);
  }//endif
#endif

  //============================================================
  // Get ready for the next iteration : 
  //      zero forces, outputAtmEnergy, atomsDone

  zeroforces();

  if(move_atoms==1){
    if(natmNow>0){
      sendAtoms(eKinetic_loc,eKineticNhc_loc,potNhc_loc,potPIMDChain_loc,natmNow,natmStr,natmEnd);
    }//endif
  }else{
    // in the minimization case the atoms don't move so we don't
    // update the coordinates
    eKinetic                    = 0.0;
    eKineticNhc                 = 0.0;
    potNhc                      = 0.0;
    potPIMDChain                = 0.0;
    eg->estruct.eKinetic_atm    = 0.0;
    eg->estruct.eKineticNhc_atm = 0.0;
    eg->estruct.potNhc_atm      = 0.0;
    eg->estruct.potPIMDChain    = 0.0;

    copySlowToFast();
    outputAtmEnergy();

    int i=0;
    CkCallback cb(CkIndex_AtomsCompute::atomsDone(NULL), CkArrayIndex1D(0), UatomsComputeProxy[thisInstance.proxyOffset]);
    contribute(sizeof(int),&i,CkReduction::sum_int,cb);
  }//endif

  //-------------------------------------------------------------------------
}//end routine
//==========================================================================

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
/** trigger force computation
 * based on real forces available in each processor's chare then contribute to
 * global group reduction of all atom forces (within this bead) -> recvContribute
 **/
//==============================================================================
void AtomsCompute::startRealSpaceForces(int t_reached){
  //==========================================================================
#ifdef _CP_DEBUG_ATMS_
  CkPrintf("{%d}[%d] AtomsCompute::startRealSpaceForces\n ", thisInstance.proxyOffset, thisIndex);     
#endif
  tol_reached = t_reached;
  //==========================================================================
  // Atom parallelization : same as for integrate

  int myid   = thisIndex;
  int nproc  = nChareAtoms;

  //==========================================================================
  // Get the real space atom forces plus these 4 quantities

  pot_ewd_rs   = 0.0;
  potGrimmeVdw = 0.0;
  vself        = 0.0;
  vbgr         = 0.0;
  potPerdCorr  = 0.0;

#ifndef _CP_DEBUG_PSI_OFF_
#ifndef _CP_DEBUG_SCALC_ONLY_ 
  CPRSPACEION::CP_getionforce(natm,&fastAtoms,myid,nproc,&pot_ewd_rs,&vself,&vbgr,&potPerdCorr, pScratchProxy.ckLocalBranch()->psscratch,&potGrimmeVdw);
#endif
#endif

  //==========================================================================
  // Everybody contributes to the reduction (see contribute forces comment)

#ifdef _CP_DEBUG_ATMS_
  CkPrintf("GJM_DBG: calling contribute atm forces %d\n",myid);
#endif
  // get the atomCache working to collect the other force contribution for us
  if(thisIndex==0)
    UatomsCacheProxy[thisInstance.proxyOffset].contributeforces();

  double *ftot           = new double[(3*natm+3)];
  for(int i=0,j=0; i<natm; i++,j+=3){
    ftot[j]   = fastAtoms.fx[i];
    ftot[j+1] = fastAtoms.fy[i];
    ftot[j+2] = fastAtoms.fz[i];
  }//endfor
  ftot[3*natm]  =pot_ewd_rs;
  ftot[3*natm+1]=potPerdCorr;
  ftot[3*natm+2]=potGrimmeVdw;
  CkPrintf("Grimme %d %g\n",myid,potGrimmeVdw);
  CkCallback cb(CkIndex_AtomsCompute::recvContribute(NULL), UatomsComputeProxy[thisInstance.proxyOffset]);
  contribute((3*natm+3)*sizeof(double),ftot,CkReduction::sum_double,cb);
  delete [] ftot;
  //==========================================================================
}//end routine
//==========================================================================


//==========================================================================
// Atom energy output
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void AtomsCompute::outputAtmEnergy() {
  //==========================================================================

  EnergyGroup *eg         = UegroupProxy[thisInstance.proxyOffset].ckLocalBranch();
  double eKinetic         = eg->estruct.eKinetic_atm;
  double eKineticNhc      = eg->estruct.eKineticNhc_atm;
  double fmag             = eg->estruct.fmag_atm;
  double pot_ewd_rs_now   = eg->estruct.eewald_real;
  double potGrimmeVdw_now = eg->estruct.grimmeVdw;
  double potNhc           = eg->estruct.potNhc_atm;
  double free_atm         = 3*((double)natm);

  CPcharmParaInfo *sim  = CPcharmParaInfo::get(); 
  int iperd             = sim->iperd;

  // Atoms cache increments it's iteration counter at the end of it's work
  // where as GSpace increments it at the beginning of an iteration. So we
  // add 1 to keep output consistent with the output from GSpace.
  AtomsCache *ag        = UatomsCacheProxy[thisInstance.proxyOffset].ckLocalBranch();
  int iteration         = ag->iteration + 1;

  int move_atoms=0;
  int do_output = 0;
  if (cp_min_opt==0 && cp_wave_opt==0) { move_atoms = 1; }
  if (cp_bomd_opt==1 && tol_reached==1) { move_atoms = 1; }
  if (cp_bomd_opt==0 || tol_reached==1) { do_output = 1; }

  int myid = CkMyPe();  
  if(myid==0 && do_output){
    fprintf(temperScreenFile,"AtomsCompute printing energies computed in iteration %d\n",iteration);
    if(iperd!=0){
      fprintf(temperScreenFile,"Iter [%d] EWALD_REAL  = %5.8lf\n",iteration, pot_ewd_rs_now);
      fprintf(temperScreenFile,"Iter [%d] EWALD_SELF  = %5.8lf\n",iteration, vself);
      fprintf(temperScreenFile,"Iter [%d] EWALD_BGR   = %5.8lf\n",iteration, vbgr);
      if(iperd!=3){
        fprintf(temperScreenFile,"Iter [%d] EWALD_Perd  = %5.8lf\n",iteration, potPerdCorr);
      }//endif
    }else{
      fprintf(temperScreenFile,"Iter [%d] ATM_COUL    = %5.8lf\n",iteration, pot_ewd_rs_now);
    }//endif
    if(cp_grimme==1){
      fprintf(temperScreenFile,"Iter [%d] Grimme VdW  = %5.8lf\n",iteration,potGrimmeVdw_now);
    }//endif
    if(move_atoms==1){
      fprintf(temperScreenFile,"Iter [%d] atm eKin    = %5.8lf\n",iteration, eKinetic);
      fprintf(temperScreenFile,"Iter [%d] atm Temp    = %5.8lf\n",iteration, (2.0*eKinetic*BOLTZ/free_atm));
      fprintf(temperScreenFile,"Iter [%d] atm fmag    = %5.8lf\n",iteration, fmag);
      if(numPIMDBeads>1){
        fprintf(temperScreenFile,"Iter [%d] Bead Chain  = %5.8lf\n",iteration, potPIMDChain);
      }//endif
      if(iextended_on==1){
        double free_Nhc;
        if(isokin_opt==0){
          free_Nhc    = free_atm*((double)len_nhc);
        }else{
          free_Nhc    = free_atm*((double)(len_nhc-1));
        }//endif
        fprintf(temperScreenFile,"Iter [%d] atm eKinNhc = %5.8lf\n",iteration, eKineticNhc);
        fprintf(temperScreenFile,"Iter [%d] atm TempNHC = %5.8lf\n",iteration, (2.0*eKineticNhc*BOLTZ/free_Nhc));
        fprintf(temperScreenFile,"Iter [%d] atm potNHC  = %5.8lf\n",iteration, potNhc);
      }//endif
    }else{
      fprintf(temperScreenFile,"Iter [%d] atm fmag    = %5.8lf\n",iteration, fmag);
    }//endif
  }//endif

  //-------------------------------------------------------------------------
}//end routine
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void AtomsCompute::send_PIMD_Fx(){ 
  //==========================================================================
  //  CkPrintf("{%d}[%d] AtomsCompute::send_PIMD_fx iteration %d\n ", thisInstance.proxyOffset, CkMyPe(), *iteration);     
  // every BOC has all Fx might as well just bcast from beadroot

  if(amBeadRoot){
    AtomXYZMsg *msg= new (natm, natm, natm, 8*sizeof(int)) AtomXYZMsg;
    for(int atomI=0;atomI<natm;atomI++){
      msg->x[atomI]=atoms[atomI].fx;
      msg->y[atomI]=atoms[atomI].fy;
      msg->z[atomI]=atoms[atomI].fz;
    }//endfor
    msg->index=PIBeadIndex;
    UPIBeadAtomsProxy[thisInstance.proxyOffset].accept_PIMD_Fx(msg);
  }else{ 
    // everyone else should chill out for Fu
  }//endif
  //-------------------------------------------------------------------------
}//end routine
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void AtomsCompute::send_PIMD_Fx_and_x(){ 
  //==========================================================================
  //  CkPrintf("{%d}[%d] AtomsCompute::send_PIMD_fx iteration %d\n ", thisInstance.proxyOffset, CkMyPe(), *iteration);     
  // every BOC has all Fx might as well just bcast from beadroot

  if(amBeadRoot){
    AtomXYZMsg *msg= new (2*natm, 2*natm, 2*natm, 8*sizeof(int)) AtomXYZMsg;
    for(int atomI=0;atomI<natm;atomI++){
      msg->x[atomI]=atoms[atomI].fx;
      msg->y[atomI]=atoms[atomI].fy;
      msg->z[atomI]=atoms[atomI].fz;
    }//endfor
    int ioff = natm;
    for(int atomI=0;atomI<natm;atomI++){
      msg->x[(atomI+ioff)]=atoms[atomI].x;
      msg->y[(atomI+ioff)]=atoms[atomI].y;
      msg->z[(atomI+ioff)]=atoms[atomI].z;
    }//endfor
    msg->index=PIBeadIndex;
    UPIBeadAtomsProxy[thisInstance.proxyOffset].accept_PIMD_Fx_and_x(msg);
  }else{ 
    // everyone else should chill out for Fu
  }//endif
  //-------------------------------------------------------------------------
}//end routine
//==========================================================================



//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
/** 
 * This thing is sending its updated atoms to all elements.
 * Effectively an all-to-all implementation of what is basically an
 * all-gather.  Cost of this isn't too bad for the array version if we
 * use a section proxy for the multicast it constrains the number of
 * destinations fairly managably. Formerly, this degraded to nAtoms
 * global broadcasts in the old group implementation.  If you've ever
 * wondered why network fabrics quake in fear at the approach of an
 * openatom run, this sort of thing should give you an inkling of the
 * abuse we dish out.  The only saving grace is that the number of
 * atoms is relatively small.
 */
void AtomsCompute::sendAtoms(double eKinetic_loc,double eKineticNhc_loc,double potNhc_loc,
    double potPIMDChain_loc,int natmNow,int natmStr,int natmEnd){
  //==========================================================================
  // Malloc the message
  //  CkPrintf("{%d}[%d] AtomsCompute::sendAtoms.\n ", thisInstance.proxyOffset, CkMyPe());     

  int nsize    = 9*natmNow+4;
  AtomMsg *msg = new (nsize,8*sizeof(int)) AtomMsg;
  CkSetQueueing(msg, CK_QUEUEING_IFIFO);
  *(int*)CkPriorityPtr(msg) = config.sfpriority-10;

  double *atmData = msg->data;
  msg->nsize      = nsize;
  msg->natmStr    = natmStr;
  msg->natmEnd    = natmEnd;

  //==========================================================================
  // pack atom positions : new for use : old for output

  for(int i=natmStr,j=0;i<natmEnd;i++,j+=9){
    atmData[(j)  ]=atoms[i].x;  // for PIMD these are the XU
    atmData[(j+1)]=atoms[i].y;
    atmData[(j+2)]=atoms[i].z;
    atmData[(j+3)]=atoms[i].xold; // for PIMD these are the X
    atmData[(j+4)]=atoms[i].yold;
    atmData[(j+5)]=atoms[i].zold;
    atmData[(j+6)]=atoms[i].vxold; // for PIMD these are the Vxu
    atmData[(j+7)]=atoms[i].vyold;
    atmData[(j+8)]=atoms[i].vzold;
  }//endfor

  //==========================================================================
  // pack the 3 energies

  atmData[(nsize-4)] = potPIMDChain_loc;
  atmData[(nsize-3)] = eKinetic_loc;
  atmData[(nsize-2)] = eKineticNhc_loc;
  atmData[(nsize-1)] = potNhc_loc;

  //==========================================================================
  // Send the message to everyone in the group

  // FIX THIS: should use delegated section to minimize stupid anytime
  // migration overhead.  Also, this is really an all-gather.

  UatomsComputeProxy[thisInstance.proxyOffset].acceptAtoms(msg);

  //-------------------------------------------------------------------------
}//end routine
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
/** 
 * For dynamics we have to update all the caches with the new coordinates
 */
void AtomsCompute::bcastAtomsToAtomCache()
{
  //==========================================================================
  // Malloc the message
  //  CkPrintf("{%d}[%d] AtomsCompute::bcastAtomsToAtomCache.\n ", thisInstance.proxyOffset, thisIndex);     

  int nsize    = 9*natm+4;
  AtomMsg *msg = new (nsize,8*sizeof(int)) AtomMsg;
  CkSetQueueing(msg, CK_QUEUEING_IFIFO);
  *(int*)CkPriorityPtr(msg) = config.sfpriority-10;

  double *atmData = msg->data;
  msg->nsize      = nsize;
  msg->natmStr    = 0;
  msg->natmEnd    = natm;

  //==========================================================================
  // pack atom positions : new for use : old for output

  for(int i=0,j=0;i<natm;i++,j+=9){
    atmData[(j)  ]=atoms[i].x;  // for PIMD these are the XU
    atmData[(j+1)]=atoms[i].y;
    atmData[(j+2)]=atoms[i].z;
    atmData[(j+3)]=atoms[i].xold; // for PIMD these are the X
    atmData[(j+4)]=atoms[i].yold;
    atmData[(j+5)]=atoms[i].zold;
    atmData[(j+6)]=atoms[i].vxold; // for PIMD these are the Vxu
    atmData[(j+7)]=atoms[i].vyold;
    atmData[(j+8)]=atoms[i].vzold;
  }//endfor
  EnergyGroup *eg             = UegroupProxy[thisInstance.proxyOffset].ckLocalBranch();
  atmData[nsize-4]=eg->estruct.eKinetic_atm;
  atmData[nsize-3]=eg->estruct.eKineticNhc_atm;
  atmData[nsize-2]=eg->estruct.potNhc_atm;
  atmData[nsize-1]=eg->estruct.potPIMDChain;

  //==========================================================================
  // Send the message to everyone in the group

  // FIX THIS: should use delegated section to minimize stupid anytime
  // migration overhead.  

  UatomsCacheProxy[thisInstance.proxyOffset].acceptAtoms(msg);

  if(amBeadRoot && !amZerothBead){
      UberCollection instance=thisInstance;
      instance.idxU.x=0;
      int offset=instance.calcPO();
      UatomsCacheProxy[offset][0].acceptChainContribution(potPIMDChain);
    }
  //-------------------------------------------------------------------------
}//end routine
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
/**
 * Take packed message of a chunk of the atoms with updated positions. Update
 * local copy of atoms.  Update local energyGroup members.  Print atom energies
 * when we have all of them.  Do file output of atoms if desired.
 */
//==========================================================================
void AtomsCompute::acceptAtoms(AtomMsg *msg) {
  //==========================================================================
  //    CkPrintf("{%d}[%d] AtomsCompute::acceptAtoms.\n ", thisInstance.proxyOffset, CkMyPe());     
  double *atmData = msg->data;
  int    nsize    = msg->nsize;
  int    natmStr  = msg->natmStr;
  int    natmEnd  = msg->natmEnd;

  //==========================================================================
  // unpack atom position and velocity : For PIMD x is reconstructed from xu

  if(numPIMDBeads>1){
    for(int i=natmStr,j=0;i<natmEnd;i++,j+=9){
      atoms[i].xu     = atmData[(j)  ];
      atoms[i].yu     = atmData[(j+1)];
      atoms[i].zu     = atmData[(j+2)];
      atoms[i].xold  = atmData[(j+3)];
      atoms[i].yold  = atmData[(j+4)];
      atoms[i].zold  = atmData[(j+5)];
      atoms[i].vxold = atmData[(j+6)];
      atoms[i].vyold = atmData[(j+7)];
      atoms[i].vzold = atmData[(j+8)];
    }//endfor
  }else{
    for(int i=natmStr,j=0;i<natmEnd;i++,j+=9){
      atoms[i].x     = atmData[(j)  ];
      atoms[i].y     = atmData[(j+1)];
      atoms[i].z     = atmData[(j+2)];
      atoms[i].xold  = atmData[(j+3)];
      atoms[i].yold  = atmData[(j+4)];
      atoms[i].zold  = atmData[(j+5)];
      atoms[i].vxold = atmData[(j+6)];
      atoms[i].vyold = atmData[(j+7)];
      atoms[i].vzold = atmData[(j+8)];
    }//endfor
  }//endif
#ifdef _CP_DEBUG_ATMS_  
  if(CkMyPe()==0)
  {
    for(int i=natmStr,j=0;i<natmEnd;i++,j+=9){

      CkPrintf("AtomCompute acceptAtoms[%d] updated to %.5g,%.5g,%.5g from %.5g,%.5g,%.5g\n",i, atoms[i].x, atoms[i].y, atoms[i].z, atoms[i].xold, atoms[i].yold, atoms[i].zold);
    }
  }
#endif
  //==========================================================================
  // unpack energy

  countAtm++;
  if(countAtm==1){
    eKinetic    = 0;
    eKineticNhc = 0;
    potNhc      = 0;
    potPIMDChain = 0;
  }//endif

  potPIMDChain+= atmData[(nsize-4)];
  eKinetic    += atmData[(nsize-3)];
  eKineticNhc += atmData[(nsize-2)];  
  potNhc      += atmData[(nsize-1)];

  //==========================================================================
  // Delete the message

  delete msg;

  //==========================================================================
  // Copy to the fast vectors and phone home : If PIMD Transform U to X

  if(countAtm==nAtmMsgRecv){
    countAtm = 0;

    //---------------------------------------------------------------------------------
    // Do these quantities need to be updated on energy group chares that are not co-resident with an AtomsCompute?

    EnergyGroup *eg             = UegroupProxy[thisInstance.proxyOffset].ckLocalBranch();
    eg->estruct.eKinetic_atm    = eKinetic;
    eg->estruct.eKineticNhc_atm = eKineticNhc;
    eg->estruct.potNhc_atm      = potNhc;
    eg->estruct.potPIMDChain    = potPIMDChain;

    //---------------------------------------------------------------------------------
    copySlowToFast();  // not complete for PIMD as you have xu not x but that's OK
    outputAtmEnergy();

    //---------------------------------------------------------------------------------
    // Output : Iteration is time of atoms[i].xold
    //          maxIter is 1 more than you need, slightly annoying but livable.

    int output_on = config.atmOutput;
    if(output_on==1 && *iteration<=config.maxIter-1){ 
      int pi_beads   = 1;
      int iwrite_atm = 0;
      int myid       = thisIndex;
      ATOMOUTPUT::ctrl_piny_output(*iteration,natm,len_nhc,pi_beads,myid,atoms,atomsNHC,
          &iwrite_atm,output_on,TemperIndex,PIBeadIndex);
      if(myid==0 && iwrite_atm>0){
        fprintf(temperScreenFile,"-----------------------------------\n");
        fprintf(temperScreenFile, "Writing atoms to disk at time %d\n",*iteration);
        fprintf(temperScreenFile,"-----------------------------------\n");
      }//endif
    }//endif

    //---------------------------------------------------------------------------------
    // Transform PIMD U to X

    if(numPIMDBeads>1){
      // CkPrintf("{%d}[%d] AtomComputep::acceptAtoms. numPIMDBeads >1 transform PIMD U to X iteration %d\n ", 
      //           thisInstance.proxyOffset, CkMyPe(), *iteration);
      if(amBeadRoot){send_PIMD_u();}
      if(amBeadRoot && amZerothBead){//For staging, this is the 1st bead not the CM but that's fine
        AtomXYZMsg *msg = new (natm,natm,natm) AtomXYZMsg;
        for(int atomI=0;atomI<natm;atomI++){
          msg->x[atomI]=atoms[atomI].xu;
          msg->y[atomI]=atoms[atomI].yu;
          msg->z[atomI]=atoms[atomI].zu;
        }//endfor
        //	  proxyHeadBeads.accept_PIMD_CM(msg);
        proxyAllBeads.accept_PIMD_CM(msg);
      }//endif : I am King of the Beads
    }else{
      //non PIMD case we're done
      int i=0;
      CkCallback cb(CkIndex_AtomsCompute::atomsDone(NULL), CkArrayIndex1D(0), UatomsComputeProxy[thisInstance.proxyOffset]);
      contribute(sizeof(int),&i,CkReduction::sum_int,cb);
    }//endif : 

  }//endif : I have received all my messages 

  //-------------------------------------------------------------------------
}//end routine
//==========================================================================

void AtomsCompute::atomsDone(CkReductionMsg *msg)
{
  int move_atoms = 0;
  if(cp_min_opt==0 && cp_wave_opt==0) { move_atoms = 1; }
  if(cp_bomd_opt==1 && tol_reached==1) { move_atoms = 1; }
  // This means that all the atomCompute chares are done integrating
  // in minimization.  Tell the atomCaches they can go ahead now.
  delete msg;
  if(move_atoms==1)
    bcastAtomsToAtomCache();
  else
    UatomsCacheProxy[thisInstance.proxyOffset].atomsDone();
}
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void AtomsCompute::accept_PIMD_CM(AtomXYZMsg *msg){
  //==============================================================================

  for(int atomnum=0; atomnum<natm; atomnum++){
    atoms[atomnum].xcm = msg->x[atomnum];
    atoms[atomnum].ycm = msg->y[atomnum];
    atoms[atomnum].zcm = msg->z[atomnum];
  }//endfor

  delete msg;
  atomsCMrecv=true;

  if(atomsPIMDXrecv){
#ifdef _CP_DEBUG_ATMS_
    CkPrintf("{%d}[%d] AtomsCompute::accept_PIMD_CM contributing to atomsDone atomsPIMDXrecv is %d iteration %d\n ", thisInstance.proxyOffset, thisIndex, atomsPIMDXrecv, *iteration);     
#endif
    int i=0;
    CkCallback cb(CkIndex_AtomsCompute::atomsDone(NULL), CkArrayIndex1D(0), UatomsComputeProxy[thisInstance.proxyOffset]);
    contribute(sizeof(int),&i,CkReduction::sum_int,cb);
    atomsCMrecv=atomsPIMDXrecv=false;
  }else{
#ifdef _CP_DEBUG_ATMS_
    CkPrintf("{%d}[%d] AtomsCompute::accept_PIMD_CM warning! atomsPIMDXrecv is %d iteration %d\n ", thisInstance.proxyOffset, thisIndex, atomsPIMDXrecv, *iteration);     
#endif
  }//endif

  //==============================================================================
}//end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void AtomsCompute::send_PIMD_u(){
  //==============================================================================
  //  CkPrintf("{%d}[%d] AtomsCompute::send_PIMD_u iteration %d\n ", thisInstance.proxyOffset, CkMyPe(), *iteration);     

  for(int atomnum=0;atomnum<natm;atomnum++){
    UPIBeadAtomsProxy[thisInstance.proxyOffset][atomnum].accept_PIMD_u(
        atoms[atomnum].xu,atoms[atomnum].yu,atoms[atomnum].zu, PIBeadIndex);
  }//endfor

}//end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
// is broadcast to us
void AtomsCompute::accept_PIMD_Fu(double _fxu, double _fyu, double _fzu, int atomI){
  //==============================================================================

  atoms[atomI].fxu    =_fxu; // FastAtoms not used for integration or output
  atoms[atomI].fyu    =_fyu;
  atoms[atomI].fzu    =_fzu;

  acceptCountfu++;
#ifdef _CP_DEBUG_ATMS_
  CkPrintf("{%d}[%d] AtomsCompute::accept_PIMD_fu (%d of %d) iteration %d %d %.5g %.5g %.5g\n", 
      thisInstance.proxyOffset, thisIndex,acceptCountfu, natm, *iteration, atomI, _fxu, _fyu, _fzu);     
#endif
  if(acceptCountfu==natm){
    //  CkPrintf("{%d}[%d] AtomsCompute::accept_PIMD_fu done calling integrator iteration %d\n ", 
    //              thisInstance.proxyOffset, CkMyPe(), *iteration);     
    integrateAtoms();
    acceptCountfu=0;
  }//endif

  //==============================================================================
}//end routine
//==============================================================================

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
// is broadcast to us
void AtomsCompute::accept_PIMD_Fu_and_u(double _fxu, double _fyu, double _fzu, 
    double _xu, double _yu, double _zu, int atomI){
  //==============================================================================

  atoms[atomI].fxu    =_fxu; // FastAtoms not used for integration or output
  atoms[atomI].fyu    =_fyu;
  atoms[atomI].fzu    =_fzu;

  atoms[atomI].xu    =_xu; // FastAtoms not used for integration or output
  atoms[atomI].yu    =_yu;
  atoms[atomI].zu    =_zu;

  acceptCountfu++;
#ifdef _CP_DEBUG_ATMS_
  CkPrintf("{%d}[%d] AtomsCompute::accept_PIMD_fu (%d of %d) iteration %d %d %.5g %.5g %.5g\n", 
      thisInstance.proxyOffset, thisIndex,acceptCountfu, natm, *iteration, atomI, _fxu, _fyu, _fzu);     
#endif
  if(acceptCountfu==natm){
    //  CkPrintf("{%d}[%d] AtomsCompute::accept_PIMD_fu done calling integrator iteration %d\n ", 
    //              thisInstance.proxyOffset, CkMyPe(), *iteration);     
    integrateAtoms();
    acceptCountfu=0;
  }//endif

  //==============================================================================
}//end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void AtomsCompute::accept_PIMD_x(double _x, double _y, double _z, int atomI){
  //==============================================================================
  //  CkPrintf("{%d}[%d] AtomsCompute::accept_PIMD_x iteration %d\n ", thisInstance.proxyOffset, thisIndex, *iteration);     
  //==============================================================================

  atoms[atomI].x=_x;
  atoms[atomI].y=_y;
  atoms[atomI].z=_z;

  fastAtoms.x[atomI]=_x;
  fastAtoms.y[atomI]=_y;
  fastAtoms.z[atomI]=_z;

  acceptCountX++;

  if(acceptCountX==natm){
    acceptCountX=0;
    atomsPIMDXrecv=true;
    if(atomsCMrecv){
      int i=0;
#ifdef _CP_DEBUG_ATMS_
      CkPrintf("{%d}[%d] AtomsCompute::accept_PIMD_x contributing to atomsDone atomsCMrecv is %d iteration %d\n ", thisInstance.proxyOffset, thisIndex, atomsCMrecv, *iteration);
#endif
      CkCallback cb(CkIndex_AtomsCompute::atomsDone(NULL), CkArrayIndex1D(0), UatomsComputeProxy[thisInstance.proxyOffset]);     

      contribute(sizeof(int),&i,CkReduction::sum_int,cb);
      atomsCMrecv=atomsPIMDXrecv=false;
    }else{
#ifdef _CP_DEBUG_ATMS_
      CkPrintf("{%d}[%d] AtomsCompute::accept_PIMD_x warning! atomsCMrecv is %d iteration %d\n",  thisInstance.proxyOffset, thisIndex, atomsCMrecv, *iteration);     
#endif
    }//endif
  }//endif : all data has arrived

  //==============================================================================
}//end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
// done during initialization in 1st iteration
//==============================================================================
void AtomsCompute::send_PIMD_x(){
  //==============================================================================

  for(int atomnum=0;atomnum<natm;atomnum++){
    UPIBeadAtomsProxy[thisInstance.proxyOffset][atomnum].accept_PIMD_x(atoms[atomnum].x,atoms[atomnum].y,atoms[atomnum].z, 
        PIBeadIndex);
  }//endfor

  //==============================================================================
}//endroutine
//==============================================================================

void AtomsCompute::acceptNewTemperature(double temp)
{
  // Hey GLENN do something with your new temperature here
  // when you're done
  int i=1;
  contribute(sizeof(int), &i, CkReduction::sum_int, 
      CkCallback(CkIndex_InstanceController::atomsDoneNewTemp(NULL),CkArrayIndex1D(thisInstance.proxyOffset),instControllerProxy), thisInstance.proxyOffset);
  //==============================================================================
}//end routine
//==============================================================================

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
// done during initialization in 1st iteration
//==============================================================================
void AtomsCompute::accept_PIMD_u(double _xu, double _yu, double _zu, int atomI){
  //==============================================================================
#ifdef _CP_DEBUG_ATMS_
  CkPrintf("{%d}[%d] AtomsCompute::accept_PIMD_u iteration %d\n", thisInstance.proxyOffset, thisIndex, *iteration);     
#endif
  atoms[atomI].xu =_xu;
  atoms[atomI].yu =_yu;
  atoms[atomI].zu =_zu;

  acceptCountu++;
  if(acceptCountu==natm){
    integrateAtoms();
    acceptCountu=0;
  }//endif

  //==============================================================================
}//end routine
//==============================================================================
/*@}*/
#include "Atoms.def.h"

