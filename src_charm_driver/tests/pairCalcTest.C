/*
 * pairCalcTest.C
 *
 *  Created on: Feb 11, 2010
 *      Author: arya
 */

//#include <cmath>
//#include <unistd.h>
//#include <iostream>
//#include <sstream>
//#include <string>
#include "charm++.h"
#include "ckarray.h"
//#include "utility/util.h"
#include "pairCalcTest.h"

using namespace std;

//-----------------------------------------------//
//NECESSARY GLOBALS
//TODO: Add these to the CI file as readonly globals

//MyConfig config;
Config config;
//thisConfig tconfig;

int TimeKeeperID=0;
vector <string> TimeKeeperNames;
CProxy_TimeKeeper TimeKeeperProxy;
CProxy_CPcharmParaInfoGrp scProxy;
CProxy_pairCalcTestMain mainProxy;

CPcharmParaInfo *sim;

TopoManager *topoMgr=NULL;

CkReduction::reducerType sumFastDoubleType=CkReduction::addReducer(sumFastDouble);

int numPes;
bool fakeTorus;

//cp_min_opt == 0 is dynamics, 1 is minimization
int cp_min_opt = 0;
//-----------------------------------------------//

pairCalcTestMain::pairCalcTestMain(CkArgMsg *msg)
{
	printf("TESTFRAMEWORK -- Starting pairCalcTest - NOTE THAT MINIMIZATION MODE RUNS (cp_min_opt == 1) DO NOT CORRECTLY WEIGHT LAMBDA\n");

	mainProxy = thisProxy;

	psiResponses = 0;
	lambdaResponses = 0;

//	printf("PairCalcTest: Sizeof(Config) = %d\n",sizeof(Config));

	TimeKeeperProxy = CProxy_TimeKeeper::ckNew();
	instControllerProxy = CProxy_InstanceController::ckNew();
	printf("TESTFRAMEWORK -- Initialized TimeKeeper and Instance Controller Proxies\n");

	setNumPes();
	printf("TESTFRAMEWORK -- setNumPes() finished\n");
	fillScProxy();
	printf("TESTFRAMEWORK -- fillScProxy() finished\n");
	fillConfig(msg->argv[1]);
	printf("TESTFRAMEWORK -- fillConfig() finished with sGrainSize = %d\n", config.sGrainSize);

	printf("TESTFRAMEWORK -- Initializing CP_State_GSPacePlane Array of dimensions [%d x %d]\n",config.nstates,config.nchareG);
	CkArrayOptions gSpaceOpts(config.nstates,config.nchareG);
	gSpaceProxy = CProxy_CP_State_GSpacePlane::ckNew(gSpaceOpts);
	printf("TESTFRAMEWORK -- Initialized CP_State_GSPacePlane Array of dimensions [%d x %d]\n",config.nstates,config.nchareG);

	buildMap();
	calcBoxSize();

	CkCallback pcHandleCB(CkIndex_CP_State_GSpacePlane::acceptPairCalcAIDs(0), gSpaceProxy);

	createPcConfigs();
	createOrthoConfigs();

	cp::startup::PCCreationManager pcCreator(cfgSymmPC, cfgAsymmPC, orthoCfg);
	pcCreator.build(pcHandleCB, boxSize, peList4PCmapping, &map);

}

void pairCalcTestMain::startTest()
{
	printf("TESTFRAMEWORK -- Done Setup.\n Initiating Psi Testing...\n");
	gSpaceProxy.sendPsi();
}

void pairCalcTestMain::finishPsi()
{
	psiResponses++;
	if(psiResponses >= config.nstates*config.nchareG)
		startTest2();
}

void pairCalcTestMain::startTest2()
{
	printf("TESTFRAMEWORK -- Initiating Lambda Testing...\n");
	gSpaceProxy.sendLambda();
}

void pairCalcTestMain::finishLambda()
{
	lambdaResponses++;
	if(lambdaResponses >= config.nstates*config.nchareG)
		done();
}

void pairCalcTestMain::done()
{
	printf("TESTFRAMEWORK -- All Done\n");
	CkExit();
}

void pairCalcTestMain::fillScProxy()
{
	sim  = new CPcharmParaInfo();

//	int istart_typ_cp   = gensimopts->istart_cp;
//	int gen_wave=0;
//	   if(istart_typ_cp ==0){gen_wave=1;}

	sim->gen_wave = 1;
	sim->tol_norb = 0;
	sim->cp_min_opt = cp_min_opt;
	scProxy  = CProxy_CPcharmParaInfoGrp::ckNew(*sim);
}

void pairCalcTestMain::fillConfig(char* input_name)
{
	int nstates_in,nplanes_in,maxIter_in;
	FILE *loutfile = fopen("dumps/configData.out", "r");
	if(loutfile!=NULL)
	{
		fscanf(loutfile,"%d\n%d\n%d\n",&nstates_in,&nplanes_in,&maxIter_in);
		fclose(loutfile);
	}
	else
	{
		CkAbort("configData.out file not found\n");
	}

	config.readConfig(input_name,nstates_in,nplanes_in,maxIter_in,numPes);

//	config.numPes = CkNumPes();
//
//	int UberImax = 1;
//	int UberJmax = 1;
//	int UberKmax = 1;
//	int numInstances = UberImax * UberJmax * UberKmax; // numIntegrals * numKpoints * numTempers;
//	config.numPesPerInstance = config.numPes / numInstances;
//
//	config.nstates = nstates_in;
//	config.maxIter = maxIter_in;
//
//	//HELP
//	//config.nchareG = NULL;
//
//	config.UberImax = UberImax;
//	config.UberJmax = UberJmax;
//	config.UberKmax = UberKmax;
//
//	config.numInstances = UberImax * UberJmax * UberKmax;
//
//	config.psipriority = 400000000;
//	config.rsfftpriority = 2000000;
//	config.gsfftpriority = 1000000;
//
//	config.conserveMemory = 0;
//
//	//for the following variables, 1 is on, 0 is off
//	config.doublePack = 1;
//	config.usePairDirectSend = 1;
//	config.PCCollectTiles = 0;
//	config.PCdelayBWSend = 0;
//	config.PCstreamBWout = 1;
//	config.useOrthoDirect = 0;
//	config.useOrthoHelpers = 0;
//	config.useOrthoSection = 1;
//	config.useOrthoSectionRed = 1;
//
//	config.invsqr_tolerance = 1e-15;
//	config.invsqr_max_iter = 1000;
//
//	config.PCSpanFactor = 2;
//
//	//HELP
//	config.sGrainSize = nstates_in;
//
//	config.gemmSplitFWk = 16;
//	if(gemmSplitFWk>sGrainSize) gemmSplitFWk=sGrainSize;
//	if(gemmSplitFWk%2!=0) gemmSplitFWk-=1;
//
//	config.gemmSplitFWm = 16;
//	if(gemmSplitFWm>sGrainSize) gemmSplitFWm=sGrainSize;
//	if(gemmSplitFWm%2!=0) gemmSplitFWm-=1;
//
//	config.gemmSplitBW = 16;
//	if(gemmSplitFWm<sGrainSize) gemmSplitFWm=sGrainSize;
//	if(gemmSplitFWm%2!=0) gemmSplitFWm+=1;
//
//	config.gemmSplitOrtho = 32;
//	if(gemmSplitOrtho>orthoGrainSize) gemmSplitOrtho=orthoGrainSize;
//	if(gemmSplitOrtho%2!=0) gemmSplitOrtho-=1;
//
//
//	config.Finale();

}

void pairCalcTestMain::setNumPes()
{
	numPes=CkNumPes();

	fakeTorus = config.fakeTorus>0;
	CkPrintf("for numInstances %d numPes %d numPesPerInstance is %d \n",config.numInstances, numPes, config.numPesPerInstance);
	if(fakeTorus)
		numPes=config.torusDimNX * config.torusDimNY * config.torusDimNZ * config.torusDimNT;
}

void pairCalcTestMain::buildMap()
{
	map.buildMap(config.nstates,config.nchareG);
}

void pairCalcTestMain::calcBoxSize()
{
	int nchareG = config.nchareG;

	CkPrintf("Initializing TopoManager\n");
	if(config.fakeTorus) {
		topoMgr = new TopoManager(config.torusDimNX, config.torusDimNY,
				config.torusDimNZ, config.torusDimNT);
	}
	else {
		topoMgr = new TopoManager();
	}
	CkPrintf("         Torus %d x %d x %d nodes %d x %d x %d VN %d DimNT %d .........\n",
			topoMgr->getDimX(), topoMgr->getDimY(), topoMgr->getDimZ(),
			topoMgr->getDimNX(), topoMgr->getDimNY(), topoMgr->getDimNZ(),
			topoMgr->hasMultipleProcsPerNode(), topoMgr->getDimNT());
	if(config.torusMap==1) {
		//PRINT_LINE_STAR;
		CkPrintf("*****************************************************************\n");
		CkPrintf("         Topology Sensitive Mapping being done for RSMap, GSMap, ....\n");
		CkPrintf("            ......., PairCalc, RhoR, RhoG and RhoGHart .........\n\n");
		//PRINT_LINE_STAR;
		CkPrintf("*****************************************************************\n");
	}
	CkPrintf("Initializing PeList\n");

	PeList *gfoo=NULL;
	PeList *rfoo=NULL;

	if(!config.loadMapFiles && config.useCuboidMap)
	{
		if( config.numPesPerInstance % config.nchareG != 0)
		{
			CkPrintf("To use CuboidMap nchareG %d should be chosen as a factor of numprocs %d / numInstances %d = %d\n",config.nchareG, config.numPes, config.numInstances, config.numPesPerInstance);
			CkExit();
		}
		int procsPerPlane = config.numPesPerInstance / nchareG;
		int bx, by, bz;
		if(config.torusMap == 1) {
			boxSize = procsPerPlane;
			int order;

			// correction to accomodate multiple instances
			int dimNX, dimNY, dimNZ, dimNT;
			int longDim = -1, maxD;

			dimNX = topoMgr->getDimNX();
			dimNY = topoMgr->getDimNY();
			dimNZ = topoMgr->getDimNZ();
			dimNT = topoMgr->getDimNT();
			int x, y, z, x1 = dimNX, y1 = dimNY, z1 = dimNZ;

			maxD = dimNX;
			longDim = 1;

			if (dimNY > maxD) { maxD = dimNY; longDim = 2; }
			if (dimNZ > maxD) { maxD = dimNZ; longDim = 3; }

			for(int i=0; i<config.numInstances; i++) {
				switch(longDim) {
				case 1:
					x = i*(maxD/config.numInstances); y = 0; z = 0;
					x1 = dimNX / config.numInstances;
					break;
				case 2:
					x = 0; y = i*(maxD/config.numInstances); z = 0;
					y1 = dimNY / config.numInstances;
					break;
				case 3:
					x = 0; y = 0; z = i*(maxD/config.numInstances);
					z1 = dimNZ / config.numInstances;
					break;
				}
				// printf("Origin %d %d %d Size %d %d %d\n", x, y, z, x1, y1, z1);
			}

			CkPrintf("         Box per Instance: nodes %d x %d x %d .........\n", x1, y1, z1);
			if(findCuboid(bx, by, bz, order, x1, y1, z1, topoMgr->getDimNT(), boxSize, topoMgr->hasMultipleProcsPerNode()))
			{
				CkPrintf("Using %d, %d, %d dimensions for box %d mapping order %d\n", bx, by, bz, boxSize, order);
				gfoo = new PeList(bx, by, bz, order, x1, y1, z1, dimNT);	// heap it
				peList4PCmapping = PeListFactory(bx,by,bz,order,x1,y1,z1,dimNT);
			}
			else
			{
				CkPrintf("no box for %d\n", boxSize);
				config.useCuboidMap = 0;
				gfoo = new PeList;				// heap it
			}
		}
		else
			gfoo = new PeList;				// heap it
	}
	else
		gfoo = new PeList;				// heap it
	if(!config.loadMapFiles && config.useCuboidMapRS)
		rfoo = new PeList(*gfoo);
	else
		rfoo = new PeList;				// heap it

}

bool pairCalcTestMain::findCuboid(int &x, int &y, int &z, int &order, int maxX, int maxY, int maxZ, int maxT, int volume, int vn){
	//============================================================================
	int maxD=maxX;
	int minD=maxX;
	int minX=maxX;
	int minY=maxY;
	int minZ=maxZ;
	//  vn=0;
	if(config.forceMappingAxis>-1)
	{
		order=config.forceMappingAxis;
		if (config.forceMappingAxis==0)
			minD=minX;
		if (config.forceMappingAxis==1)
			minD=maxY;
		if (config.forceMappingAxis==2)
			minD=maxZ;
	}
	else
	{
		order=0;
		maxD = (maxY>maxD) ? maxY : maxD;
		maxD = (maxZ>maxD) ? maxZ : maxD;
		if(vn)
		{  // using Y as the prism axis seems to suck
			maxD = (maxY>maxD) ? maxY : maxD;
			maxD = (maxZ>maxD) ? maxZ : maxD;
			minD = (maxY<minD) ? maxY : minD;
			minD = (maxZ<maxD) ? maxZ : minD;

		}
		else
		{
			maxD = (maxY>maxD) ? maxY : maxD;
			maxD = (maxZ>maxD) ? maxZ : maxD;
			minD = (maxY<minD) ? maxY : minD;
			minD = (maxZ<minD) ? maxZ : minD;

		}
	}
	// CkPrintf("minD %d maxD %d\n",minD, maxD);
	if(config.useCuboidMapRS)
	{
		CkPrintf("Using long prisms for useCuboidMapRS\n");
	}

	// We were reducing the volume by maxT and then finding the dimensions of the
	// box in terms of the no. of nodes and not processors

	// That worked ok on BG/L, but BG/P maxT of 4 distorts the long axis logic.
	// in TXYZ you have a 32x8x16 for 1 Rack or 32x8x32 for 2 Rack
	// in XYZT you have 8x8x64 for 1 Rack or 8x8x128 for 2 Rack

	// The latter case makes clear the bizarre and unreal distortion created by
	// pretending that maxT is just a multiplier along one dimension when
	// considering actual box size.  What was merely odd on BG/L is now
	// seriously peculiar.  Also default TXYZ mappings and default Order 0
	// mappings to use the X as long axis for prisms fail to work right because
	// we aren't really spanning the X*maxT dimension with these long prisms.

	int redVol = volume / maxT;
	double cubert= cbrt((double) redVol);
	int cubetrunc= (int) cubert;
	x=y=z=cubetrunc;
	if(cubetrunc>minD)
		cubetrunc=minD;
	if(cubetrunc>maxY)
		cubetrunc=maxY;
	if(cubetrunc>maxZ)
		cubetrunc=maxZ;
	if(redVol==x*y*z && !config.useCuboidMapRS)
		return true;
	bool switchSet=false;
	CkAssert(redVol>0);
	switch (redVol) // for the common values we just pick cuboids we like
	{
	case 1:
		x=1; y=1; z=1; switchSet=true; break;
	case 2:
		x=2; y=1; z=1; switchSet=true; break;
	case 3:
		x=3; y=1; z=1; switchSet=true; break;
	case 4:
		x=2; y=2; z=1; switchSet=true; break;
	case 5:
		x=5; y=1; z=1; switchSet=true; break;
	case 6:
		x=3; y=2; z=1; switchSet=true; break;
	case 7:
		x=7; y=1; z=1; switchSet=true; break;
	case 8:
		if(config.useCuboidMapRS)
		{
			if(minD>=2)

			{ x=2; y=2; z=2; switchSet=true; break;}
			/* no evidence that these other schemes help
	  if(minD==4)
	  { x=4; y=2; z=1; switchSet=true; break;}
	  if(minD>=8)
	  { x=8; y=1; z=1; switchSet=true; break;}
			 */
		}
	case 9:
		x=3; y=3; z=1; switchSet=true; break;
	case 10:
		x=5; y=2; z=1; switchSet=true; break;
	case 12:
		x=2; y=3; z=2; switchSet=true; break;
	case 14:
		x=7; y=2; z=1; switchSet=true; break;
	case 15:
		x=5; y=3; z=1; switchSet=true; break;
	case 16:
		if(config.useCuboidMapRS)
		{
			if(minD>=8)
			{ x=8; y=2; z=1; switchSet=true; break;}
		}
		x=4; y=2; z=2; switchSet=true; break;
	case 18:
		x=3; y=3; z=2; switchSet=true; break;
	case 20:
		x=5; y=2; z=2; switchSet=true; break;
	case 21:
		x=7; y=3; z=1; switchSet=true; break;
	case 24:
		x=4; y=3; z=2; switchSet=true; break;
	case 25:
		x=5; y=5; z=1; switchSet=true; break;
	case 27:
		x=3; y=3; z=3; switchSet=true; break;
	case 28:
		x=7; y=2; z=2; switchSet=true; break;
	case 30:
		x=5; y=2; z=2; switchSet=true; break;
	case 32:
		if(config.useCuboidMapRS)
		{
			if(minD==8)
			{ x=8; y=2; z=2; switchSet=true; break;}
			if(minD==16)
			{ x=16; y=2; z=1; switchSet=true; break;}
			if(minD>=32)
			{ x=32; y=1; z=1; switchSet=true; break;}

		}
		x=4; y=2; z=4; switchSet=true; break;
	case 35:
		x=7; y=5; z=1; switchSet=true; break;
	case 36:
		x=4; y=3; z=3; switchSet=true; break;
	case 40:
		x=5; y=4; z=2; switchSet=true; break;
	case 42:
		x=7; y=3; z=2; switchSet=true; break;
	case 43:
		x=7; y=3; z=2; switchSet=true; break;
	case 45:
		x=5; y=3; y=3; switchSet=true; break;
	case 48:
		x=4; y=3; z=4; switchSet=true; break;
	case 50:
		x=5; y=5; z=2; switchSet=true; break;
	case 54:
		x=6; y=3; z=3; switchSet=true; break;
	case 56:
		x=7; y=4; z=2; switchSet=true; break;
	case 60:
		x=5; y= 4; z=3; switchSet=true; break;
	case 64:
		if(config.useCuboidMapRS)
		{
			if(minD==8)
			{ x=8; y=4; z=2; switchSet=true; break;}
			if(minD==16)
			{ x=16; y=2; z=2; switchSet=true; break;}
			if(minD==32)
			{ x=32; y=2; z=1; switchSet=true; break;}
			if(minD>=64)
			{ x=64; y=1; z=1; switchSet=true; break;}
		}
		x=4; y=4; z=4; switchSet=true; break;
	case 128:
		if(config.useCuboidMapRS)
		{
			if(minD==8)
			{ x=8; y=4; z=4; switchSet=true; break;}
			if(minD==16)
			{x=16; y=4; z=2; switchSet=true; break;}
			if(minD==32)
			{  x=32; y=2; z=2; switchSet=true; break;	}
			if(minD==64)
			{  x=64; y=2; z=1; switchSet=true; break;	}
			if(minD>=128)
			{  x=128; y=1; z=1; switchSet=true; break;	}
		}
		x=8; y=4; z=4; switchSet=true; break;
	case 256:
		if(config.useCuboidMapRS)
		{
			if(minD==8)
			{ x=8; y=8; z=4; switchSet=true; break;}
			if(minD==16)
			{ x=16; y=4; z=4; switchSet=true; break;}
			if(minD==32)
			{ x=32; y=4; z=2; switchSet=true; break;}
			if(minD==64)
			{ x=64; y=2; z=2; switchSet=true; break;}
			if(minD>=128)
			{ x=64; y=2; z=1; switchSet=true; break;}
		}
		x=8; y=8; z=4; switchSet=true; break;
	case 512:
		if(config.useCuboidMapRS)
		{
			if(minD==8)
			{ x=8; y=8; z=8; switchSet=true; break;}
			if(minD==16)
			{ x=16; y=4; z=8; switchSet=true; break;}
			if(minD==32)
			{ x=32; y=4; z=4; switchSet=true; break;}
			if(minD==64)
			{ x=64; y=4; z=2; switchSet=true; break;}
			if(minD==128)
			{ x=128; y=2; z=2; switchSet=true; break;}
			if(minD>=256)
			{ x=256; y=2; z=1; switchSet=true; break;}
		}
		x=8; y=8; z=8; switchSet=true; break;

	default:
		break;
	}
	if(switchSet && config.forceMappingAxis>-1)
	{

		switch(config.forceMappingAxis)
		{
		case 0:
			return true; 	  //no change
		case 1:
		{
			order=1; //YXZ
			int swap=x;
			x=y;
			y=swap;
			return true;
		}
		case 2:
		{
			order=2; //ZXY
			int swap=x;
			x=z;
			z=swap;
			return true;
		}
		}
	}
	else if (switchSet)
	{
		// now correct the x,y,z to put long prism axis along the
		// smallest torus dimension which will fit.
		if(x==maxX)
			return true;
		if(x==maxY)
		{ // change to Y
			order=1; //YXZ
			int swap=x;
			x=y;
			y=swap;
			return true;
		}
		if(x==maxZ)
		{ // change to Z
			order=2; //ZXY
			int swap=x;
			x=z;
			z=swap;
			return true;
		}
		// if we're here then we don't have a spanning prism
		// just pick the smallest which will fit.
		if(x<maxX)
			return true;
		if(x<maxY)
		{ // change to Y
			order=1; //YXZ
			int swap=x;
			x=y;
			y=swap;
			return true;
		}
		if(x<maxZ)
		{ // change to Z
			order=2; //ZXY
			int swap=x;
			x=z;
			z=swap;
			return true;
		}
	}
	// its something weird so try a best fit
	int start=cubetrunc-1;
	if(config.useCuboidMapRS)
		x = (redVol>=maxX) ? maxX : cubetrunc;
	else
		x=cubetrunc;
	for(; x<=maxX;x++)
	{
		for(y=start; y<=maxY;y++)
		{
			for(z=start; z<=maxZ;z++)
			{
				if(redVol==x*y*z)
					return true;
			}
		}
	}
	return false;
}

void pairCalcTestMain::createPcConfigs()
{
	int doublePack = config.doublePack;
	int nstates = config.nstates;

	// Stuff it with the actual configurations
	cfgSymmPC.isDynamics         = (sim->cp_min_opt==1)? false: true;

	cfgSymmPC.numPlanes          = config.nchareG;
	cfgSymmPC.numStates          = nstates;
	cfgSymmPC.grainSize          = config.sGrainSize;
	cfgSymmPC.orthoGrainSize     = config.orthoGrainSize;

	cfgSymmPC.conserveMemory     = config.conserveMemory;
	cfgSymmPC.isLBon             = config.lbpaircalc;

	cfgSymmPC.areBWTilesCollected= config.PCCollectTiles;
	cfgSymmPC.isBWstreaming      = config.PCstreamBWout;
	cfgSymmPC.isBWbarriered      = config.useBWBarrier;
	cfgSymmPC.shouldDelayBWsend  = config.PCdelayBWSend;
	cfgSymmPC.isInputMulticast   = !config.usePairDirectSend;
	cfgSymmPC.isOutputReduced    = !config.gSpaceSum;
	cfgSymmPC.inputSpanningTreeFactor = config.PCSpanFactor;

	cfgSymmPC.instanceIndex      = 0; //thisInstance.proxyOffset;

	cfgSymmPC.gemmSplitFWk       = config.gemmSplitFWk;
	cfgSymmPC.gemmSplitFWm       = config.gemmSplitFWm;
	cfgSymmPC.gemmSplitBW        = config.gemmSplitBW;

	cfgAsymmPC = cfgSymmPC;

	// Configurations specific to the symmetric PC instance
	cfgSymmPC.isSymmetric        = true;
	cfgSymmPC.arePhantomsOn      = config.phantomSym;
	cfgSymmPC.numChunks          = config.numChunksSym;
	cfgSymmPC.isDoublePackOn     = doublePack;
	cfgSymmPC.inputMsgPriority   = config.psipriority;
	cfgSymmPC.resultMsgPriority  = config.gsfftpriority;

	// Configurations specific to the asymmetric PC instance
	cfgAsymmPC.isSymmetric        = false;
	cfgAsymmPC.arePhantomsOn      = false;
	cfgAsymmPC.numChunks          = config.numChunksAsym;
	cfgAsymmPC.isDoublePackOn     = 0;
	cfgAsymmPC.inputMsgPriority   = config.lambdapriority;
	cfgAsymmPC.resultMsgPriority  = config.lambdapriority+2;

	// Configure the GSpace entry methods that the PCs will callback
	if(cfgSymmPC.isOutputReduced)
	{
		cfgSymmPC.gSpaceEP        = CkIndex_CP_State_GSpacePlane::__idx_acceptNewPsi_CkReductionMsg;
		cfgSymmPC.PsiVEP          = CkIndex_CP_State_GSpacePlane::__idx_acceptNewPsiV_CkReductionMsg;
	}
	else
	{
		cfgSymmPC.gSpaceEP        = CkIndex_CP_State_GSpacePlane::__idx_acceptNewPsi_partialResultMsg;
		cfgSymmPC.PsiVEP          = CkIndex_CP_State_GSpacePlane::__idx_acceptNewPsiV_partialResultMsg;
	}

	if(cfgAsymmPC.isOutputReduced)
	{
		cfgAsymmPC.gSpaceEP       = CkIndex_CP_State_GSpacePlane::__idx_acceptLambda_CkReductionMsg;
		cfgAsymmPC.PsiVEP         = 0;
	}
	else
	{
		cfgAsymmPC.gSpaceEP       = CkIndex_CP_State_GSpacePlane::__idx_acceptLambda_partialResultMsg;
		cfgAsymmPC.PsiVEP         = 0;
	}

#ifdef _CP_SUBSTEP_TIMING_
	//symmetric AKA Psi
	cfgSymmPC.forwardTimerID      = keeperRegister("Sym Forward");
	cfgSymmPC.backwardTimerID     = keeperRegister("Sym Backward");
	cfgSymmPC.beginTimerCB        = CkCallback(CkIndex_TimeKeeper::collectStart(NULL),0,TimeKeeperProxy);
	cfgSymmPC.endTimerCB          = CkCallback(CkIndex_TimeKeeper::collectEnd(NULL),0,TimeKeeperProxy);
	//asymmetric AKA Lambda AKA Gamma
	cfgAsymmPC.forwardTimerID     = keeperRegister("Asym Forward");
	cfgAsymmPC.backwardTimerID    = keeperRegister("Asym Backward");
	cfgAsymmPC.beginTimerCB       = CkCallback(CkIndex_TimeKeeper::collectStart(NULL),0,TimeKeeperProxy);
	cfgAsymmPC.endTimerCB         = CkCallback(CkIndex_TimeKeeper::collectEnd(NULL),0,TimeKeeperProxy);
#endif

	cfgSymmPC.uponSetupCompletion= CkCallback(CkIndex_InstanceController::doneInit(NULL),0,instControllerProxy.ckGetArrayID());
	cfgAsymmPC.uponSetupCompletion= CkCallback(CkIndex_InstanceController::doneInit(NULL),0,instControllerProxy.ckGetArrayID());

	cfgSymmPC.gSpaceAID    = gSpaceProxy.ckGetArrayID();
	cfgAsymmPC.gSpaceAID   = gSpaceProxy.ckGetArrayID();
}

void pairCalcTestMain::createOrthoConfigs()
{
//	orthoCfg.numStates     = config.nstates;
//	orthoCfg.grainSize     = config.orthoGrainSize;
//	orthoCfg.instanceIndex = 0; //thisInstance.getPO();
//	orthoCfg.uponToleranceFailure = CkCallback(CkCallback::ignore);


	orthoCfg.isDynamics    = (cp_min_opt==1)? false: true;
	orthoCfg.isGenWave     = (sim->gen_wave==1)? true: false;
	orthoCfg.numStates     = config.nstates;
	orthoCfg.grainSize     = config.orthoGrainSize;
	orthoCfg.instanceIndex = 0;
	orthoCfg.maxTolerance  = sim->tol_norb;
	orthoCfg.uponToleranceFailure = CkCallback(CkCallback::ignore);
}

#include "CPcharmParaInfo.def.h"
#include "pairCalcTest.def.h"
//#include "ortho.def.h"
#include "timeKeeper.def.h"
