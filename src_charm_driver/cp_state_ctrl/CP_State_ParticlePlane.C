//=========================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=========================================================================
/** \file CP_State_ParticlePlane.C
 * Life-cycle of a CP_State_ParticlePlane:
 *
 * The particle-plane array is a shadow of the g-space planes. This means
 * that particlePlaneProxy(s, p) exists on the same processor as
 * gSpacePlaneProxy(s, p). Also,  both chare arrays have the same number of
 * objects.
 *
 * The life-cycle of this object involves computation of a portion of
 * the Z-matrix, reduction of the Z matrix over g-space and computing
 * forces using the particle data.  ComputeZ is triggered by the
 * arrival of the structure factor for each atom group at the local
 * sfcache.
 *
 * The computation of the Z-matrix is done in the method computeZ().
 * Immediately after this, the matrix is reduced using the reduceZ() method.
 * Through the reduceZ() method, a portion of the reduced matrix is 
 * available at the zeroth plane of each state.
 * 
 * The reduced matrix is sent from the particle-plane 0 to all the particle
 * planes in the state using the getForces() method, which computes
 * forces and adds these forces to the forces compputed using electron
 * density in the g-space planes. 
 *
 * Important information:
 * 1. The dimensionality of the Z-matrix depends upon how may orbital levels
 * are considered (l = 0,1...). The current code assumes that only l = 0
 * is allowed. This has to be changed.
 * 
 * 2. The special point gx=gy=gz=0 is not dealt with.
 */ 
//=========================================================================

#include "charm++.h"
#include "util.h"
#include "groups.h"
#include "cpaimd.h"
#include "../../include/debug_flags.h"
#include "sim_subroutines.h"
#include "CP_State_Plane.h"
#include "StructFactorCache.h"
#include "ckmulticast.h"
#include "../../src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cpnonlocal.h"
#include "../../src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cplocal.h"

//=========================================================================

extern CProxy_CP_State_GSpacePlane gSpacePlaneProxy;
extern CProxy_AtomsGrp atomsGrpProxy;
extern CProxy_CPcharmParaInfoGrp scProxy;
extern CProxy_CP_State_ParticlePlane particlePlaneProxy;
extern CProxy_StructFactCache sfCacheProxy;
extern CProxy_EnergyGroup egroupProxy; //energy group proxy
extern int nstates;
extern int nchareG;
extern Config config;
extern CkGroupID mCastGrpId;

//=========================================================================


//============================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void printEnl(void *param, void *msg){
  CkReductionMsg *m=(CkReductionMsg *)msg;
  double d = ((double *)m->getData())[0];
  CkPrintf("ENL         = %5.8lf\n", d);
  delete m;

  gSpacePlaneProxy(0,0).computeEnergies(ENERGY_ENL, d);  
}
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
CP_State_ParticlePlane::CP_State_ParticlePlane(int x, int y, int z, 
                                               int _gSpacePPC, int numSfGrps_in,
                                               int natm_nl_in, int natm_nl_grp_max_in)
//============================================================================
  {//begin routine
//============================================================================

  sizeX                = x;
  sizeY                = y;
  sizeZ                = z;
  gSpacePlanesPerChare = _gSpacePPC;
  numSfGrps            = numSfGrps_in;
  natm_nl              = natm_nl_in;
  natm_nl_grp_max      = natm_nl_grp_max_in;
  count                = new int[numSfGrps];
  bzero(count,numSfGrps*sizeof(int));
  haveSFAtmGrp         = new int[numSfGrps];
  memset(haveSFAtmGrp,-1,numSfGrps*sizeof(int));
  doneEnl              = 0;
  doneForces           = 0;
  enl                  = 0.0;
  energy_count         = 0;
  totalEnergy          = 0.0;
  reductionPlaneNum    = calcReductionPlaneNum(thisIndex.x);
  contribute(sizeof(int), &sizeX, CkReduction::sum_int);
  setMigratable(false);
  usesAtSync           = CmiFalse;

  zmatrixSum    = NULL;
  zmatrix       = NULL;
  zmatrix_fx    = NULL; zmatrix_fy    = NULL; zmatrix_fz    = NULL;
  zmatrixSum_fx = NULL; zmatrixSum_fy = NULL; zmatrixSum_fz = NULL;

//============================================================================
//  Specifies the end of the particle plane iteration
//  It is set to 1 in reduceZ.
//  After this only GSpacePlane will begin the next iteration.
//  register with the SFCache

  StructFactCache *sfcache = sfCacheProxy.ckLocalBranch();
  for(int i=0;i<numSfGrps;i++){
    sfcache->registerPP(thisIndex.x, thisIndex.y,i);
  }//endfor

//============================================================================
// Create section proxy for ENL reduction ParticlePlane (any state#, reductionPlaneNum)

  if(thisIndex.x==0 && thisIndex.y==0){

      CkArrayIndexMax *elems = new CkArrayIndexMax[nstates];

      CkArrayIndex2D idx(0, reductionPlaneNum);  // plane# = this plane#
      for (int j = 0; j < nstates; j++) {
	idx.index[0] = j;
	idx.index[1] = calcReductionPlaneNum(j);
	elems[j] = idx;
      }//endfor

      particlePlaneENLProxy = 
	CProxySection_CP_State_ParticlePlane::ckNew(particlePlaneProxy.ckGetArrayID(),
						    elems, nstates);
      CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch();
      particlePlaneENLProxy.ckDelegate(mcastGrp);
      EnlCookieMsg *emsg= new EnlCookieMsg;
      mcastGrp->setSection(particlePlaneENLProxy);
      particlePlaneENLProxy.setEnlCookie(emsg);

  }//endif

//---------------------------------------------------------------------------
   }//end routine
//============================================================================


//============================================================================
/// this function is called from CP_State_GSpacePlane::initGSpace()
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_ParticlePlane::initKVectors(GStateSlab *gss){
    gss->setKVectors(&gSpaceNumPoints, &k_x, &k_y, &k_z);
    myForces = new complex[gSpaceNumPoints]; // forces on coefs
    memset(myForces, 0, gSpaceNumPoints * sizeof(complex));
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
CP_State_ParticlePlane::~CP_State_ParticlePlane(){

	delete [] myForces;
	delete [] k_x;
	delete [] k_y;
	delete [] k_z;
	//delete [] zmatrixSum;

        delete [] zmatrix; //FOO-BAR should i delete this?            
        delete [] zmatrix_fx; 
        delete [] zmatrix_fy; 
        delete [] zmatrix_fz; 
	delete [] count;
	delete [] haveSFAtmGrp;
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_ParticlePlane::pup(PUP::er &p){
	ArrayElement2D::pup(p);
	p|gSpaceNumPoints;
	p|doneGettingForces;
	p|numSfGrps;
        p|natm_nl;
	p|natm_nl_grp_max;
	int zsize=numSfGrps*natm_nl_grp_max;
	if (p.isUnpacking()) {
		k_x = new int[gSpaceNumPoints];
		k_y = new int[gSpaceNumPoints];
		k_z = new int[gSpaceNumPoints];
		myForces = new complex[gSpaceNumPoints]; // forces on coefs
                memset(myForces, 0, gSpaceNumPoints * sizeof(complex));
		// these 
		zmatrixSum = new complex[zsize];
                zmatrixSum_fx = new complex[zsize];
                zmatrixSum_fy = new complex[zsize];
                zmatrixSum_fz = new complex[zsize];
		zmatrix = new complex[zsize];
                zmatrix_fx = new complex[zsize];
                zmatrix_fy = new complex[zsize];
                zmatrix_fz = new complex[zsize];
		count= new int[numSfGrps];
		haveSFAtmGrp= new int[numSfGrps];
	}//endif
	p(k_x, gSpaceNumPoints);
	p(k_y, gSpaceNumPoints);
	p(k_z, gSpaceNumPoints);
	p(count,numSfGrps);
	p(haveSFAtmGrp,numSfGrps);
	p((char*)zmatrixSum,zsize *sizeof(complex));
	p((char*)zmatrixSum_fx,zsize *sizeof(complex));
	p((char*)zmatrixSum_fy,zsize *sizeof(complex));
	p((char*)zmatrixSum_fz,zsize *sizeof(complex));
	p((char*)zmatrix,zsize *sizeof(complex));
	p((char*)zmatrix_fx,zsize *sizeof(complex));
	p((char*)zmatrix_fy,zsize *sizeof(complex));
	p((char*)zmatrix_fz,zsize *sizeof(complex));
	p|reductionPlaneNum;
	// Dont pack zmatrixSum
	p|sizeX;
	p|sizeY;
	p|sizeZ;
	p|gSpacePlanesPerChare;
	p|zsize;
	p|energy_count;
	p|totalEnergy;
	p|enl;
	p|doneEnl;
	p|doneForces;
        p|enl_total;
	p|enlCookie;
	p|particlePlaneENLProxy;
	// "gspace" is not pup'ed since it is always assigned to
}
//============================================================================

//============================================================================
/* computeZ is triggered by the arrival of the structure factor for
   each atom group in the local sfcache. */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_ParticlePlane::computeZ(PPDummyMsg *m){
//============================================================================    

    doneGettingForces = false;    //if you are here, you ain't done.
                                  //It is also reset in sendLambda 
                                  //which is when the results of the previous
                                  //iteration are no longer needed.
    CP_State_GSpacePlane *gsp = 
        	        gSpacePlaneProxy(thisIndex.x, thisIndex.y).ckLocal();
    GStateSlab *gss = &(gsp->gs);
    int atmIndex    = m->atmGrp;
    int sfindex     = m->sfindex;
    delete m;
   

    // you can't add to the energy group, until the previous step is done filling it.
    // since computeZ leads to energies, thats no good. 
    // Since gps is bound to pp, this syncs pp with energy group correctly.
    if(gsp->acceptedPsi){
      if(!gsp->doneNewIter){ // false after send lambda, true after accept energy/atoms
         CkPrintf("Flow of Control Warning in computeZ. \n");
         PPDummyMsg *pmsg = new (8*sizeof(int)) PPDummyMsg;
  	 pmsg->atmGrp     = atmIndex;
         pmsg->sfindex    = sfindex;
         CkSetQueueing(pmsg, CK_QUEUEING_IFIFO);
	 *(int*)CkPriorityPtr(pmsg) = config.sfpriority+config.numSfGrps; 
         //lower than sf and sfcache
         particlePlaneProxy(thisIndex.x, thisIndex.y).computeZ(pmsg);
      }//endif
    }//endif

//============================================================================    
// If you got here and you have the latest psi 

    if(gsp->acceptedPsi){ //need latest Psi to computeZ
      StructFactCache *sfcache = sfCacheProxy.ckLocalBranch();

      complex *structureFactor;
      complex *structureFactor_fx;
      complex *structureFactor_fy;
      complex *structureFactor_fz;
      sfcache->getStructFact(thisIndex.y, atmIndex, &structureFactor, 
			     &structureFactor_fx, &structureFactor_fy, &structureFactor_fz);
      zsize = 0;
      AtomsGrp *ag = atomsGrpProxy.ckLocalBranch(); // find me the local copy
      zsize = natm_nl_grp_max*numSfGrps;

      if(zmatrix==NULL){
	zmatrix    = new complex[zsize];
	zmatrix_fx = new complex[zsize];
	zmatrix_fy = new complex[zsize];
	zmatrix_fz = new complex[zsize];
	memset(zmatrix, 0, sizeof(complex)*zsize);
	memset(zmatrix_fx, 0, sizeof(complex)*zsize);
	memset(zmatrix_fy, 0, sizeof(complex)*zsize);
	memset(zmatrix_fz, 0, sizeof(complex)*zsize);
      }//endif

      int mydoublePack = config.doublePack;
      int zoffset=natm_nl_grp_max * atmIndex;
      CPNONLOCAL::CP_enl_matrix_calc(gSpaceNumPoints,gss->packedPlaneData, k_x,k_y,k_z, 
	     structureFactor,structureFactor_fx,structureFactor_fy,
	     structureFactor_fz,&zmatrix[zoffset],&zmatrix_fx[zoffset],&zmatrix_fy[zoffset],
	     &zmatrix_fz[zoffset],thisIndex.x,mydoublePack,numSfGrps,atmIndex);

      // reduce zmatrices over the planes of each state
      thisProxy(thisIndex.x, reductionPlaneNum).reduceZ(natm_nl_grp_max, atmIndex, 
	&zmatrix[zoffset],&zmatrix_fx[zoffset],&zmatrix_fy[zoffset],&zmatrix_fz[zoffset]);

      /*  overwritten each time so redundant
	  bzero(&zmatrix[zoffset], natm_nl_grp_max * sizeof(complex));
	  bzero(&zmatrix_fx[zoffset], natm_nl_grp_max * sizeof(complex));
	  bzero(&zmatrix_fy[zoffset], natm_nl_grp_max * sizeof(complex));
	  bzero(&zmatrix_fz[zoffset], natm_nl_grp_max * sizeof(complex));
      */

      haveSFAtmGrp[atmIndex]=-1; //this one is done
//============================================================================    
// If you got here and you don't have the latest psi (impossible now)
    }else{
      // we beat Psi here.  Set flag so Psi can kick this off
	haveSFAtmGrp[atmIndex]=sfindex;
    }//endif

//----------------------------------------------------------------------------
  }//end routine
//============================================================================


//============================================================================
/// this function achieves the reduction over gspace for a particular state
//============================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_ParticlePlane::reduceZ(int size, int atmIndex, complex *zmatrix_,
       complex *zmatrix_fx_,complex *zmatrix_fy_,complex *zmatrix_fz_){
//============================================================================
// This method is invoked only CP_State_ParticlePlane(*,reductionPlaneNum)
//  CkPrintf("PP [%d %d] reduceZ\n",thisIndex.x, thisIndex.y);
  count[atmIndex]++;
  int i,j;
  int zsize = natm_nl_grp_max * numSfGrps;
  int zoffset= natm_nl_grp_max * atmIndex;
  if (zmatrixSum==NULL) {
      zmatrixSum = new complex[zsize];
      memset(zmatrixSum, 0, zsize * sizeof(complex));
  }//endif

  if (zmatrixSum_fx==NULL) {
      zmatrixSum_fx = new complex[zsize];
      memset(zmatrixSum_fx, 0, zsize * sizeof(complex));
  }//endif

  if (zmatrixSum_fy==NULL) {
      zmatrixSum_fy = new complex[zsize];
      memset(zmatrixSum_fy, 0, zsize * sizeof(complex));
  }//endif

  if (zmatrixSum_fz==NULL) {
      zmatrixSum_fz = new complex[zsize];
      memset(zmatrixSum_fz, 0, zsize * sizeof(complex));
  }//endif
  
  for (i = 0,j=zoffset; i < size; i++,j++){
    zmatrixSum[j]    += zmatrix_[i];
    zmatrixSum_fx[j] += zmatrix_fx_[i];
    zmatrixSum_fy[j] += zmatrix_fy_[i];
    zmatrixSum_fz[j] += zmatrix_fz_[i];
  }//endif

  if ( ((count[atmIndex] == (sizeX/gSpacePlanesPerChare)/2 - 1)&&(!config.doublePack)) || 
       ((count[atmIndex] == nchareG)&&(config.doublePack)) ) { 
    count[atmIndex] = 0;
    if(doneEnl==0){enl=0.0;}
    doneEnl++;

    // Now we have contributions from all the other gspace planes
    // in this state. Calculate the energies and atoms forces.
    // Normalize ZmatrixSum

    AtomsGrp *ag = atomsGrpProxy.ckLocalBranch();
    Atom *atoms  = ag->atoms;
    int mydoublePack = config.doublePack;
    double myenl=0.0;
    CPNONLOCAL::CP_enl_atm_forc_calc(numSfGrps,atmIndex,atoms,&zmatrixSum[zoffset],
	     &zmatrixSum_fx[zoffset],&zmatrixSum_fy[zoffset],&zmatrixSum_fz[zoffset],
             &myenl,mydoublePack);
    enl+=myenl;

#ifdef _CP_DEBUG_NLMAT_
    for (i = 0; i < zsize; i++){
      CkPrintf("Non-local matrix[%d][%d]={%g, %g}\n",i,thisIndex.x,
	       zmatrix_[i].re,zmatrix_[i].im );}
#endif

    // send the reduced zmatrix to the particleplanes of this state
    // that are non-zero and compute the non-local forces
    for (i = 0; i < nchareG; i ++){
      thisProxy(thisIndex.x, i).getForces(size,atmIndex, &zmatrixSum[zoffset]);
    }//endfor
    //    CkPrintf("PP [%d %d] doneENL %d numSfGrps %d\n",thisIndex.x, thisIndex.y, doneEnl, numSfGrps);
    // done with Z stuff send out our ENL
    if(thisIndex.y==reductionPlaneNum && doneEnl==numSfGrps){
      CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch();
      CkCallback cb=CkCallback(printEnl, NULL);
      mcastGrp->contribute(sizeof(double),(void*) &enl, 
			   CkReduction::sum_double, enlCookie, cb);
      //      CkPrintf("PP [%d %d] contributed ENL %g\n",thisIndex.x, thisIndex.y,enl);
      doneEnl=0;
      enl=0.0;
      bzero(zmatrixSum, zsize * sizeof(complex));
      bzero(zmatrixSum_fx, zsize * sizeof(complex));
      bzero(zmatrixSum_fy, zsize * sizeof(complex));
      bzero(zmatrixSum_fz, zsize * sizeof(complex));
    }//endif

  }/*endif*/

//---------------------------------------------------------------------------
}//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_ParticlePlane::sumEnergies(double energy_) {
        
        enl_total   += energy_;
	totalEnergy += energy_;
	energy_count++;

	if (energy_count == nstates) {
            energy_count = 0;
            totalEnergy = 0;
            enl_total = 0;
	}//endif
}
//============================================================================


void ParticleBarrierClient(void *param, void *msg) {
  CkReductionMsg *m = (CkReductionMsg *)msg;
  //delete m;
  
  printf("After particle plane barrier\n");
  
  gSpacePlaneProxy.startFFT(m);
}


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_ParticlePlane::getForces(int zmatSize, int atmIndex, 
                                       complex *zmatrixSum_loc)
//============================================================================
{//begin routine
//============================================================================
//  CkPrintf("PP [%d %d] getforces \n",thisIndex.x,thisIndex.y);
  doneForces++;
  CP_State_GSpacePlane *gsp = 
       gSpacePlaneProxy(thisIndex.x, thisIndex.y).ckLocal();
  GStateSlab *gss = &(gsp->gs);
  AtomsGrp *ag = atomsGrpProxy.ckLocalBranch();
  int state_ind = gss->istate_ind;                

  // compute the non-local forces using the qc_subroutine
  int mydoublePack = config.doublePack;
  StructFactCache *sfcache = sfCacheProxy.ckLocalBranch();
  complex *structureFactor,*structureFactor_fx,*structureFactor_fy,*structureFactor_fz;
  sfcache->getStructFact(thisIndex.y, atmIndex, &structureFactor, &structureFactor_fx, 
                         &structureFactor_fy, &structureFactor_fz);

  if (gSpaceNumPoints != 0){
    CPNONLOCAL::CP_enl_force_calc(zmatrixSum_loc,
				  gSpaceNumPoints, k_x, k_y, k_z, structureFactor,
				  myForces, state_ind, mydoublePack, numSfGrps, atmIndex);
  }//endif

  // forces have been computed, send them to the corresponding
  // g-space plane
  if(doneForces==numSfGrps){
      doneGettingForces = true;
      doneForces=0;

#if !CP_PARTICLE_BARRIER
      if (gsp->doneDoingIFFT) {
	gsp->combineForcesGetEke(); //calls resume in G very yucky
      }//endif
#else
      int dummy=0;
      contribute(sizeof(int),&dummy,CkReduction::sum_int,
		 CkCallback(ParticleBarrierClient, NULL));
#endif
      
  }//endif : doneforces
  
//----------------------------------------------------------------------------
}//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_ParticlePlane::setEnlCookie(EnlCookieMsg *m){
  CkGetSectionInfo(enlCookie,m);
  //  CkPrintf("PP [%d %d] get enl cookie\n",thisIndex.x, thisIndex.y); 
  //  delete m; 
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_ParticlePlane::ResumeFromSync(){
  if(thisIndex.x==0 && thisIndex.y==0){
      CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch();       
      mcastGrp->resetSection(particlePlaneENLProxy);
  }//endif
}// end routine
//============================================================================

int CP_State_ParticlePlane::calcReductionPlaneNum(int state)
{
  
  int nstatemax=nstates-1;
  int ncharemax=nchareG-1;
  int planeNum= state %ncharemax;
  if(planeNum<0)
    {
      CkPrintf(" PP [%d %d] calc nstatemax %d ncharemax %d state %d planenum %d\n",thisIndex.x, thisIndex.y,nstatemax, ncharemax, state, planeNum); 
      CkExit();
    }
  return planeNum;
}
