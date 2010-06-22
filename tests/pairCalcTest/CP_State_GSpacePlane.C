#include "CP_State_GSpacePlane.h"
#include "pairCalcTest.h"
#include <string>

using namespace std;

//extern MyConfig config;
extern int cp_min_opt;
extern Config config;
extern CProxy_pairCalcTestMain mainProxy;
extern void fastAdd (double*, double*, int);
extern string inputDir;

CP_State_GSpacePlane::CP_State_GSpacePlane()
//CP_State_GSpacePlane::CP_State_GSpacePlane(Config config_in, int cp_min_opt_in)
{
	//cp_min_opt = cp_min_opt_in;

	//memcpy(&config, &config_in, sizeof(Config));

//	printf("CP_State_GSpacePlane: Sizeof(Config) = %d\n",sizeof(Config));

	int nstates = config.nstates;
	int s_grain = config.sGrainSize;
	numPoints = 0;

//	printf("HAPPY----------------------HAPPY1, sgrain=%d and cp_min_opt=%d\n",s_grain,cp_min_opt);

	int ourgrain = thisIndex.x/s_grain*s_grain;
	if(nstates == s_grain){
		AllPsiExpected=1;
	}else{
		//    if(ourgrain<lastgrain){ // corner has no extras
		if(ourgrain<(nstates-s_grain)){ // corner has no extras
			AllPsiExpected=2;
		}else{
			AllPsiExpected=1;
		}//endif
	}//endif
	AllPsiExpected*=config.numChunksSym;

//	printf("HAPPY----------------------HAPPY2\n");

	int numGrains=nstates/s_grain;
	if(config.gSpaceSum){ // no reductions its all coming direct
		AllPsiExpected=numGrains*config.numChunksSym;
	}//endif
	//---------------------------------------------
	// Asymm PC accounting
	if(cp_min_opt==0){    //expect non diagonal column results
		if(nstates == s_grain){
			AllLambdaExpected=1;
		}else {
			AllLambdaExpected=2;
		}//endif
	}else{  //no column reduction results in minimization
		AllLambdaExpected=1;
	}//endif

//	printf("HAPPY----------------------HAPPY3\n");

	if(config.gSpaceSum){ // no reductions its all coming direct
		if(cp_min_opt==0 && numGrains>1)
		{ AllLambdaExpected=(2*numGrains-1)*config.numChunksAsym;}
		else
		{ AllLambdaExpected=numGrains*config.numChunksAsym*AllLambdaExpected;}
	}//endif
	else
	{
		AllLambdaExpected*=config.numChunksAsym;
	}

//	printf("HAPPY----------------------HAPPY4\n");

	countPsi        = 0;
	countLambda     = 0;
	countVPsi       = 0;

	countVPsiO= new int[config.numChunksSym];
	for(int i=0;i<config.numChunksSym;i++)
		countVPsiO[i]=0;

	countPsiO = new int[config.numChunksSym];
	for(int i=0;i<config.numChunksSym;i++)
		countPsiO[i]=0;

	countLambdaO= new int[config.numChunksAsym];
	for(int i=0;i<config.numChunksAsym;i++)
		countLambdaO[i]=0;

	//dumpMatrixDouble("lambdaBf",(double *)force, 1,
    //gnumPoints*2,thisIndex.y,thisIndex.x,thisIndex.x,0,false);

	string path;

	//Load lambda
	path.assign(inputDir)+="lambdaBf";
	packedForceData = loadMatrixDouble(path.c_str(),1,
			thisIndex.y,thisIndex.x,thisIndex.x,0,false);

	//Load Psi
	path.assign(inputDir)+="psiBf";
	packedPlaneData = loadMatrixDouble(path.c_str(),1,
			thisIndex.y,thisIndex.x,thisIndex.x,0,false);

	//Load Psip
	path.assign(inputDir)+="psiBfp";
	packedPlaneDatap = loadMatrixDouble(path.c_str(),1,
				thisIndex.y,thisIndex.x,thisIndex.x,0,false);

	complex *force = packedForceData;
	complex *psi = packedPlaneData;
	complex *psip = packedPlaneDatap;

	if(cp_min_opt==0){
		packedPlaneDataScr   = new complex[numPoints];
		packedPlaneDataTemp  = new complex[numPoints];
		memset(packedPlaneDataScr, 0, sizeof(complex)*numPoints);
	}//endif

	//CkPrintf("Gspace Chare %d x %d re-dumping data %f\n",thisIndex.x,thisIndex.y,((double *)force)[0]);

//	dumpMatrixDouble("dumps2/lambdaBf",(double *)force, 1,
//	                     numPoints*2,thisIndex.y,thisIndex.x,thisIndex.x,0,false);
//
//	dumpMatrixDouble("dumps2/psiBfp",(double *)psi, 1, numPoints*2,
//	                     thisIndex.y,thisIndex.x,thisIndex.x,0,false);

	//CkPrintf("Gspace Chare %d x %d finsihed loading data\n",thisIndex.x,thisIndex.y);


	//Read in gspace configuration data
	char fmt[1000];
	char filename[1000];
	path.assign(inputDir)+="gspaceconfig";
	strncpy(fmt,path.c_str(),999);
	strncat(fmt,"_%d_%d.out",999);
	sprintf(filename,fmt, thisIndex.x, thisIndex.y);
	FILE *loutfile = fopen(filename, "r");
	if(loutfile!=NULL)
	{
		fscanf(loutfile,"%d\n%d\n%d\n",&ihave_kx0,&kx0_strt,&kx0_end);

		fscanf(loutfile,"%d\n",&size_k_x);
		k_x = new int[size_k_x];
		for(int i = 0; i<size_k_x; i++)
			fscanf(loutfile,"%d\n",&k_x[i]);

		fclose(loutfile);
	}
	else
	{
		strncat(filename," not found!",999);
		CkAbort(filename);
	}

	//CkPrintf("ihave_kx0 = %d, kx0_strt=%d, kx0_end=%d, numPoints=%d\n",ihave_kx0,kx0_strt,kx0_end,numPoints);

}

CP_State_GSpacePlane::~CP_State_GSpacePlane()
{
	delete countVPsiO;
	delete countPsiO;
	delete countLambdaO;
	delete k_x;
	if(cp_min_opt==0){
	 delete packedPlaneDataScr;
	 delete	packedPlaneDataTemp;
	}
}

void CP_State_GSpacePlane::acceptPairCalcAIDs(pcSetupMsg *msg)
{
	symmPCmgr  = cp::gspace::PCCommManager(thisIndex, msg->symmCfg, msg->symmIDs);
	asymmPCmgr = cp::gspace::PCCommManager(thisIndex, msg->asymmCfg, msg->asymmIDs);
	myOrtho    = CProxy_Ortho(msg->orthoAID);
	//CkPrintf("TESTFRAMEWORK -- acceptPairCalcAIDs finished\n");
}

void CP_State_GSpacePlane::makePCproxies(){
	lambdaproxy=new CProxySection_PairCalculator[config.numChunksAsym];
	lambdaproxyother=new CProxySection_PairCalculator[config.numChunksAsym];
	psiproxy=new CProxySection_PairCalculator[config.numChunksSym];
	psiproxyother=new CProxySection_PairCalculator[config.numChunksSym];
	//need one proxy per chunk
	if(!config.gSpaceSum){
		for(int chunk=0;chunk<config.numChunksAsym;chunk++){
			lambdaproxy[chunk]=asymmPCmgr.makeOneResultSection_asym(chunk);
			if(AllLambdaExpected/config.numChunksAsym == 2)//additional col. red. in dynamics
				lambdaproxyother[chunk]=asymmPCmgr.makeOneResultSection_asym_column(chunk);
		}//endfor chunk
		for(int chunk=0; chunk < config.numChunksSym ;chunk++){
			psiproxy[chunk]=symmPCmgr.makeOneResultSection_sym1(chunk);
			if(AllPsiExpected / config.numChunksSym > 1)
				psiproxyother[chunk]=symmPCmgr.makeOneResultSection_sym2(chunk);
		}//endfor chunk
	}//endif not gspacesum
	CkPrintf("TESTFRAMEWORK -- makePCproxies finished\n");
}//end routine

void CP_State_GSpacePlane::sendPsi() {
	//==============================================================================
	// Prepare the data : If cp dynamics is going (cp_min_opt == 0), save the non-orthogonal puppies.

	complex *psi   = packedPlaneDatap;
	//int numPoints  = gs.numPoints;

//	if(cp_min_opt==0){
//		int ncoef     = numPoints;
//		complex *scr  = packedPlaneDataScr; //save non-orthog psi
//		CmiMemcpy(scr,psi,sizeof(complex)*ncoef);
//	}//endif
//
//	if(ihave_kx0==1){
//		double rad2i = 1.0/sqrt(2.0);
//		for(int i=kx0_strt; i<kx0_end; i++){psi[i] *= rad2i;}
//	}

	symmPCmgr.sendLeftData(numPoints, psi, false);
	/// Symm loop PC chares in the top left [*,0,0,*] will not receive any right matrix data. Hence, if you're in such a PC's block, dont send right
	if(thisIndex.x >= symmPCmgr.pcCfg.grainSize)
		symmPCmgr.sendRightData(numPoints, psi, false);

}// end routine

void  CP_State_GSpacePlane::sendPsiV() {
	int ncoef     = numPoints;
	complex *data = packedVelData;
	complex *scr  = packedPlaneDataScr;  // replace no-ortho psi
	complex *psi  = packedPlaneDatap;     // by orthonormal psi when norb rotating
	CmiMemcpy(scr,psi,sizeof(complex)*ncoef);

	if(ihave_kx0==1){
		double rad2i = 1.0/sqrt(2.0);
		for(int i=kx0_strt; i<kx0_end; i++){data[i] *= rad2i;}
	}//endif

	//int numPoints = gs.numPoints;
	symmPCmgr.sendLeftData(numPoints,data,true);
	/// Symm loop PC chares in the top left [*,0,0,*] will not receive any right matrix data. Hence, if you're in such a PC's block, dont send right
	if(thisIndex.x >= symmPCmgr.pcCfg.grainSize)
		symmPCmgr.sendRightData(numPoints,data,true);
}// end routine

void  CP_State_GSpacePlane::sendLambda() {
	complex *psi   = packedPlaneData;
	complex *force = packedForceData;

//#ifndef _CP_DEBUG_ORTHO_OFF_
//  if(ihave_kx0==1 && cp_min_opt==0){
//    double rad2i = 1.0/sqrt(2.0);
//    double rad2  = sqrt(2.0);
//    for(int i=kx0_strt; i<kx0_end; i++){psi[i]   *= rad2i;}
//    for(int i=kx0_strt; i<kx0_end; i++){force[i] *= rad2;}
//  }//endif
//#endif

	//int numPoints   = gs.numPoints;
	int toSend = numPoints;
	asymmPCmgr.sendLeftData(toSend,psi,false);
	CmiNetworkProgress();
	asymmPCmgr.sendRightData(toSend,force,false);
}// end routine

void CP_State_GSpacePlane::acceptNewPsi(CkReductionMsg *msg){
	//=============================================================================
	// (0) Nan Check and output

	int N           = msg->getSize()/sizeof(complex);
	complex *data   = (complex *)msg->getData();
	int offset      = msg->getUserFlag();  if(offset<0){offset=0;}
	complex *psi    = packedPlaneDatap;
	int chunksize   = numPoints/config.numChunksSym;
	int chunkoffset = offset*chunksize; // how far into the points this contribution lies

	//=============================================================================
	// (I) Unpack the contribution to newpsi (orthonormal psi)

	int idest=chunkoffset;

	if(countPsiO[offset]<1){
		CmiMemcpy(&(psi[idest]), &(data[0]), N*sizeof(complex)); //slower?
		//for(int i=0; i<N; i++,idest++){psi[idest] = data[i];}
	}else{
		//    for(int i=0; i<N; i++,idest++){psi[idest] += data[i];}
		fastAdd((double *) &psi[idest],(double *)data,N*2);
	}//endif

	delete msg;

	//=============================================================================
	// (II) If you have got it all : Rescale it, produce some output

	countPsi++;//psi arrives in as many as 2 *numblock reductions
	countPsiO[offset]++;//psi arrives in as many as 2
	if(countPsi==AllPsiExpected){
		thisProxy(thisIndex.x,thisIndex.y).doNewPsi();
		//CkPrintf("%d x %d Received NewPsi\n", thisIndex.x, thisIndex.y);
	}//endif

}//end routine

void CP_State_GSpacePlane::acceptNewPsi(partialResultMsg *msg){

	int N           = msg->N;
	complex *data   = msg->result;
	int offset      = msg->myoffset;  if(offset<0){offset=0;}
	complex *psi    = packedPlaneDatap;
	int chunksize   = numPoints/config.numChunksSym;
	int chunkoffset = offset*chunksize; // how far into the points this contribution lies

	//=============================================================================
	// (I) Unpack the contribution to newpsi (orthonormal psi)

	int idest = chunkoffset;
	if(countPsiO[offset]<1){
		CmiMemcpy(&(psi[idest]), &(data[0]), N*sizeof(complex)); //slower?
		//    for(int i=0; i<N; i++,idest++){psi[idest] = data[i];}
	}else{
		//for(int i=0; i<N; i++,idest++){psi[idest] += data[i];}
		fastAdd((double *) &psi[idest],(double *)data,N*2);
	}//endif

	delete msg;

	//==============================================================================
	// When you are done, continue

	countPsi++;         //psi arrives in as many as 2 *numblock * numgrain reductions
	countPsiO[offset]++;//psi arrives in as many as 2 * numgrain
	//
	if(countPsi==AllPsiExpected){
		thisProxy(thisIndex.x,thisIndex.y).doNewPsi();
		//CkPrintf("%d x %d Received NewPsi\n", thisIndex.x, thisIndex.y);
	}//endif
}//end routine

void CP_State_GSpacePlane::doNewPsi()
{

	complex *psi  = packedPlaneDatap;

	if(ihave_kx0==1){
		double rad2 = sqrt(2.0);
		for(int i=kx0_strt; i<kx0_end; i++){psi[i] *= rad2;}
	}//endif

	countPsi = 0;
	bzero(countPsiO,config.numChunksSym*sizeof(int));

	double testvalue=0.00000001;

	string path;
	path.assign(inputDir)+="psiAf";
	complex * savedpsiAf  = loadMatrixDouble(path.c_str(),1,
			thisIndex.y,thisIndex.x,thisIndex.x,0,false);

	for(int i=0;i<numPoints;i++){
		if(fabs(psi[i].re-savedpsiAf[i].re)>testvalue){
			fprintf(stderr, "GSP [%d,%d] %d element psi  %.10g not %.10g\n",
					thisIndex.x, thisIndex.y,i, psi[i].re, savedpsiAf[i].re);
		}//endif
		CkAssert(fabs(psi[i].re-savedpsiAf[i].re)<testvalue);
		CkAssert(fabs(psi[i].im-savedpsiAf[i].im)<testvalue);
	}//endfor

	mainProxy.finishPsi();
}

void CP_State_GSpacePlane::acceptNewPsiV(CkReductionMsg *msg){
	//=============================================================================
	// (I) Local pointers

	int N           = msg->getSize()/sizeof(complex);
	complex *data   = (complex *)msg->getData();
	int offset      = msg->getUserFlag();  if(offset<0){offset=0;}

	complex *vpsi   = packedVelData;
	int chunksize   = numPoints/config.numChunksSym;
	int chunkoffset = offset*chunksize;; // how far into the points this contribution lies

	//=============================================================================
	// (II) Unpack the contribution to newpsi (orthonormal psi)

	int idest=chunkoffset;
	if(countVPsiO[offset]<1){
		CmiMemcpy(&(vpsi[idest]), &(data[0]), N*sizeof(complex));//slower?
		//for(int i=0; i<N; i++,idest++){vpsi[idest] = data[i];}
	}else{
		//    for(int i=0; i<N; i++,idest++){vpsi[idest] += data[i];}
		fastAdd((double *) &vpsi[idest],(double *)data,N*2);
	}//endif

	delete msg;

	//=============================================================================
	// (III) When all has arrive, onward to victory

	countVPsi++;         //psi arrives in as many as 2 reductions
	countVPsiO[offset]++;//psi arrives in as many as 2 reductions

	if(countVPsi==AllPsiExpected){
		//thisProxy(thisIndex.x,thisIndex.y).doNewPsiV();
		countVPsi = 0;
		for(int i=0;i<config.numChunksSym;i++){countVPsiO[i]=0;}
		//CkPrintf("%d x %d Received NewPsiV\n", thisIndex.x, thisIndex.y);
	}//endif

	//----------------------------------------------------------------------------
}//end routine

void CP_State_GSpacePlane::acceptNewPsiV(partialResultMsg *msg){
	//=============================================================================
	// (I) Local pointers

	int N           = msg->N;
	complex *data   = msg->result;
	int offset      = msg->myoffset;  if(offset<0){offset=0;}

	complex *vpsi   = packedVelData;
	int chunksize   = numPoints/config.numChunksSym;
	int chunkoffset = offset*chunksize;; // how far into the points this contribution lies

	//=============================================================================
	// (I) Unpack the contribution to newpsi (orthonormal psi)

	int idest=chunkoffset;
	if(countVPsiO[offset]<1){
		CmiMemcpy(&(vpsi[idest]), &(data[0]), N*sizeof(complex));//slower?
		//for(int i=0; i<N; i++,idest++){vpsi[idest] = data[i];}
	}else{
		//    for(int i=0; i<N; i++,idest++){vpsi[idest] += data[i];}
		fastAdd((double *) &vpsi[idest],(double *)data,N*2);
	}//endif

	delete msg;

	//=============================================================================
	// (II) Continue

	countVPsi++;//psi arrives in as many as 2 reductions
	countVPsiO[offset]++;//psi arrives in as many as 2 reductions

	if(countVPsi==AllPsiExpected){
//		thisProxy(thisIndex.x,thisIndex.y).doNewPsiV();
		countVPsi = 0;
		for(int i=0;i<config.numChunksSym;i++){countVPsiO[i]=0;}
		//CkPrintf("%d x %d Received NewPsiV\n", thisIndex.x, thisIndex.y);
	}//endif
}//end routine

void CP_State_GSpacePlane::acceptLambda(CkReductionMsg *msg) {
//==============================================================================

//  int cp_min_opt    = scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt;
//  eesCache *eesData = UeesCacheProxy[thisInstance.proxyOffset].ckLocalBranch ();
//  int *k_x          = eesData->GspData[iplane_ind].ka;

  int       N    = msg->getSize()/sizeof(complex);
  complex *data  = (complex *)msg->getData();
  int offset     = msg->getUserFlag();  if(offset<0){offset=0;}
  //  CkPrintf("[%d %d] accepts lambda %d \n", thisIndex.x, thisIndex.y,offset);
  complex *force = packedForceData;
  int chunksize  = numPoints/config.numChunksAsym;
  int chunkoffset=offset*chunksize; // how far into the points this contribution lies

  countLambda++;
  //---------------------------------------------------
  // B) Double Pack
  if(config.doublePack==1){
   if(cp_min_opt==1){
#ifdef CMK_BLUEGENEL
#pragma unroll(10)
#endif
     for(int i=0,idest=chunkoffset; i<N; i++,idest++){
       double wght  = (k_x[idest]==0 ? 0.5 : 1);
       CkAssert(idest<=size_k_x);
       force[idest].re -= wght*data[i].re;
       force[idest].im -= wght*data[i].im;
     }//endfor
   }else{
     if(countLambdaO[offset]<1){
#ifdef CMK_BLUEGENEL
#pragma unroll(10)
#endif
       for(int i=0,idest=chunkoffset; i<N; i++,idest++){force[idest]  = data[i]*(-1.0);}
     }else{
#ifdef CMK_BLUEGENEL
#pragma unroll(10)
#endif
       for(int i=0,idest=chunkoffset; i<N; i++,idest++){force[idest]  += data[i]*(-1.0);}
     }// coutlambda
   }//endif : cpmin

  }//endif : double pack

  //---------------------------------------------------
  // C) Double Pack
  if(config.doublePack==0){
#ifdef CMK_BLUEGENEL
#pragma unroll(10)
#endif
    for(int i=0,idest=chunkoffset; i<N; i++,idest){
       force[idest].re -= 0.5*data[i].re;
       force[idest].im -= 0.5*data[i].im;
    }//endfor
  }//endif
  //---------------------------------------------------
  // D) Cleanup
  delete msg;

  countLambdaO[offset]++;
  if(countLambda==AllLambdaExpected){
	  thisProxy(thisIndex.x,thisIndex.y).doLambda();
	  //CkPrintf("%d x %d Received Lambda\n", thisIndex.x, thisIndex.y);
  }//endif
   }//end routine
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::acceptLambda(partialResultMsg *msg) {
//==============================================================================
// 0) unpack the message : pop out variables from groups

  complex *data     = msg->result;
  int       N       = msg->N;
  int offset        = msg->myoffset;
  if(offset<0){offset=0;}

//  int cp_min_opt    = scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt;
//  eesCache *eesData = UeesCacheProxy[thisInstance.proxyOffset].ckLocalBranch ();
//  int *k_x          = eesData->GspData[iplane_ind].ka;

  complex *force    = packedForceData;
  int chunksize     = numPoints/config.numChunksAsym;
  int chunkoffset   = offset*chunksize; // how far into the points this contribution lies

  countLambda++;

  //B) Double Pack
   if(config.doublePack==1){
    if(cp_min_opt==1){
 #ifdef CMK_BLUEGENEL
 #pragma unroll(10)
 #endif
      for(int i=0,idest=chunkoffset; i<N; i++,idest++){
        double wght  = (k_x[idest]==0 ? 0.5 : 1);
        force[idest].re -= wght*data[i].re;
        force[idest].im -= wght*data[i].im;
      }//endfor
    }else{
      if(countLambdaO[offset]<1){
 #ifdef CMK_BLUEGENEL
 #pragma unroll(10)
 #endif
         for(int i=0,idest=chunkoffset; i<N; i++,idest++){force[idest]  = data[i]*(-1.0);}
      }else{
 #ifdef CMK_BLUEGENEL
 #pragma unroll(10)
 #endif
         for(int i=0,idest=chunkoffset; i<N; i++,idest++){force[idest]  += data[i]*(-1.0);}
      }//endif : off set thingy
    }//endif : cp_min_on
   }//endif : double pack


  //----------------------------------------------------------
  //C) Single pack

   if(config.doublePack==0){
 #ifdef CMK_BLUEGENEL
 #pragma unroll(10)
 #endif
     for(int i=0,idest=chunkoffset; i<N; i++,idest){
        force[idest].re -= 0.5*data[i].re;
        force[idest].im -= 0.5*data[i].im;
     }//endfor

   }//endif : single pack

  //----------------------------------------------------------
  //D) Clean up

  delete msg;

    countLambdaO[offset]++;
    if(countLambda==AllLambdaExpected){
    	thisProxy(thisIndex.x,thisIndex.y).doLambda();
    	//CkPrintf("%d x %d Received Lambda\n", thisIndex.x, thisIndex.y);
    }//endif
}//end routine

void CP_State_GSpacePlane::doLambda() {

	complex *force    = packedForceData;

	if(cp_min_opt==0){
		// dynamics scale it out
		if(ihave_kx0==1){
			double rad2i = 1.0/sqrt(2.0);
			for(int i=kx0_strt; i<kx0_end; i++){force[i] *= rad2i;}
		}//endif
	}//endif

	countLambda=0;
	bzero(countLambdaO,config.numChunksAsym*sizeof(int));

	double testvalue=0.00000001;
	string path;
	path.assign(inputDir)+="lambdaAf";
	complex* savedlambdaAf = loadMatrixDouble(path.c_str(),
			1,thisIndex.y,thisIndex.x,thisIndex.x,0,false);

	for(int i=0;i<numPoints;i++){
		CkAssert(fabs(force[i].re-savedlambdaAf[i].re)<testvalue);
		CkAssert(fabs(force[i].im-savedlambdaAf[i].im)<testvalue);
	}//endfor

	mainProxy.finishLambda();
}

void CP_State_GSpacePlane::pup(PUP::er &p) {
	ArrayElement2D::pup(p);
//	UgSpaceDriverProxy[thisInstance.proxyOffset](thisIndex.x,thisIndex.y).pup(p);
//	p|istate_ind;
//	p|iplane_ind;
//	p|registrationFlag;
//	p|initialized;
//	p|istart_typ_cp;
//	p|iteration;
//	p|nrotation;
//	p|exitFlag;
//	p|cleanExitCalled;
//	p|finishedCpIntegrate;
//	p|iRecvRedPsi;
//	p|iSentRedPsi;
//	p|iRecvRedPsiV;
//	p|iSentRedPsiV;
//	p|ireset_cg;
//	p|numReset_cg;
	p|countPsi;
	p|countVPsi;
	p|countLambda;
//	p|countIFFT;
//	p|countFileOut;
//	p|ecount;
//	p|countRedPsi;
//	p|numRecvRedPsi;
	p|AllPsiExpected;
	p|AllLambdaExpected;
//	p|numRDMAlinksSymm;
//	p|numRDMAlinksAsymm;
//	p|doneDoingIFFT;
//	p|doneNewIter;
//	p|acceptedPsi;
//	p|acceptedVPsi;
//	p|acceptedLambda;
//	p|itemp; // 2 temporary variables for debugging in scope
//	p|jtemp;
//
//	p|ehart_total;
//	p|enl_total;
//	p|eke_total;
//	p|fictEke_total;
//	p|fmagPsi_total0;
//	p|fmagPsi_total;
//	p|fmagPsi_total_old;
//	p|fovlap;
//	p|fovlap_old;
//	p|egga_total;
//	p|eexc_total;
//	p|eext_total;
//	p|ewd_total;
//	p|total_energy;
//	p|cpuTimeNow;
	p|numPoints;
	p|ihave_kx0;
	p|kx0_strt;
	p|kx0_end;
//	p|real_proxy;
//	gs.pup(p);
//	p|gSpaceNumPoints;
	if (p.isUnpacking()) {
		lambdaproxy=new CProxySection_PairCalculator[config.numChunksAsym];
		lambdaproxyother=new CProxySection_PairCalculator[config.numChunksAsym];
		psiproxy=new CProxySection_PairCalculator[config.numChunksSym];
		psiproxyother=new CProxySection_PairCalculator[config.numChunksSym];
		countPsiO= new int[config.numChunksSym];
		countVPsiO= new int[config.numChunksSym];
		countLambdaO= new int[config.numChunksAsym];
		k_x = new int[numPoints];
	}//endif
	PUParray(p,lambdaproxy,config.numChunksAsym);
	PUParray(p,lambdaproxyother,config.numChunksAsym);
	PUParray(p,psiproxy,config.numChunksSym);
	PUParray(p,psiproxyother,config.numChunksSym);


//	if (p.isUnpacking()) {
//		packedPlaneDatap     = (complex *)fftw_malloc(numPoints*sizeof(complex));
//		packedForceData     = (complex *)fftw_malloc(numFull*sizeof(complex));
//		packedVelData       = (complex *)fftw_malloc(numPoints*sizeof(complex));
//		if(cp_min_opt==0){
//			packedPlaneDataScr  = (complex *)fftw_malloc(numPoints*sizeof(complex));
//		}//endif
//#ifdef  _CP_DEBUG_UPDATE_OFF_
//		if(cp_min_opt==1){
//			packedPlaneDataTemp = (complex *)fftw_malloc(numPoints*sizeof(complex));
//		}//endif
//#endif
//	}//endif
//
	p((char *) packedPlaneData, numPoints*sizeof(complex));
	p((char *) packedPlaneDatap, numPoints*sizeof(complex));
//	p((char *) packedForceData, numFull*sizeof(complex));
	p((char *) packedVelData, numPoints*sizeof(complex));   //cg under min
//	if(cp_min_opt==0){
//		p((char *) packedPlaneDataScr, numPoints*sizeof(complex));
//	}//endif

	PUParray(p,countPsiO,config.numChunksSym);
	PUParray(p,countVPsiO,config.numChunksSym);
	PUParray(p,countLambdaO,config.numChunksAsym);
	//-------------------------------------------------------
}// end routine : pup

complex * CP_State_GSpacePlane::loadMatrixDouble(const char *infilename, int xdim, int w,int x,int y, int z, bool symmetric)
{
	char fmt[1000];
	char filename[1000];
	int ydim, thisNumPoints;
	complex *cmatrix;
	strncpy(fmt,infilename,999);
	strncat(fmt,"_%d_%d_%d_%d_%d.out",999);
	sprintf(filename,fmt, w, x, y, z, symmetric);
	FILE *loutfile = fopen(filename, "r");
	if(loutfile!=NULL)
	{
		int junk1,junk2;

		//read in ydim and calculate numPoints
		fscanf(loutfile,"%d\n",&ydim);
		thisNumPoints = ydim/2;

		//CkPrintf("Gspace Chare %d x %d has numPoints = %d\n",x,w,thisNumPoints);

		if(numPoints == 0)
		{
			CkAssert(thisNumPoints>0);
			numPoints = thisNumPoints;
		}
		else
			CkAssert(thisNumPoints == numPoints);

		//Allocate space for input
		cmatrix = new complex[numPoints];
		double *matrix = (double *)cmatrix;

		for(int i=0;i<xdim;i++)
			for(int j=0;j<ydim;j++)
			{
				fscanf(loutfile,"%d %d %lf\n",&junk1,&junk2,&(matrix[i*ydim+j]));
				//CkPrintf("Found: %lf\n",((double *)cmatrix)[i*ydim+j]);
			}
		fclose(loutfile);
	}
	else
	{
		strncat(filename," not found!",999);
		CkAbort(filename);
	}

	return cmatrix;
}

void CP_State_GSpacePlane::dumpMatrixDouble(const char *infilename, double *matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric)
{
  char fmt[1000];
  char filename[1000];
  strncpy(fmt,infilename,999);
  strncat(fmt,"_%d_%d_%d_%d_%d.out",999);
  sprintf(filename,fmt, w, x, y, z, symmetric);
  FILE *loutfile = fopen(filename, "w");
  //CkPrintf("File opened!\n");
//#ifdef PAIRCALC_TEST_DUMP
  fprintf(loutfile,"%d\n",ydim);
//#endif
  for(int i=0;i<xdim;i++)
    for(int j=0;j<ydim;j++)
      fprintf(loutfile,"%d %d %.12g\n",x+i,y+j,matrix[i*ydim+j]);
  fclose(loutfile);
}

#include "gStatePlane.def.h"
#include "startupMessages.def.h"
