/*
 * pairCalcTest.h
 *
 *  Created on: Dec 14, 2009
 *      Author: anshu
 */

#ifndef PAIRCALCTEST_H_
#define PAIRCALCTEST_H_

//#include "main/CPcharmParaInfoGrp.h"

#define _CP_SUBSTEP_TIMING_

#include "configure.h"
//#include "this_configure.h"
//#include "MyGStateSlab.h"
#include "main/pcCreationManager.h"
#include "paircalc/pcConfig.h"
#include "orthog_ctrl/orthoConfig.h"
#include "CP_State_GSpacePlane.h"
#include "load_balance/IntMap.h"
#include "main/TimeKeeper.h"
#include "CPcharmParaInfoGrp.h"
#include "CPcharmParaInfo.decl.h"
#include "timeKeeper.decl.h"
#include "instanceController.decl.h"
#include "gStatePlane.decl.h"
#include "pairCalcTest.decl.h"

#include <iostream>
#include <vector>
#include <string>

class pairCalcTestMain : public CBase_pairCalcTestMain {
 public:
	pairCalcTestMain(CkMigrateMessage *m) {};
	pairCalcTestMain(CkArgMsg *);
    ~pairCalcTestMain() {};
    void startTest();
    void startTest2();
    void finishPsi();
    void finishLambda();
    void done();
 private:
    int boxSize;
    PeListFactory peList4PCmapping;
    cp::paircalc::pcConfig cfgSymmPC;
    cp::paircalc::pcConfig cfgAsymmPC;
    cp::ortho::orthoConfig orthoCfg;
    MapType2 map;
    CProxy_InstanceController instControllerProxy;
    CProxy_CP_State_GSpacePlane gSpaceProxy;

    int psiResponses;
    int lambdaResponses;

    void fillScProxy();
    void fillConfig(char* input_name, int nstates_in, int nplanes_in, int maxIter_in);
    void setNumPes();
    void buildMap();
    void calcBoxSize();
    bool findCuboid(int &x, int &y, int &z, int &order, int maxX, int maxY, int maxZ, int maxT, int volume, int vn);
    void createPcConfigs();
    void createOrthoConfigs();
};

//Duplicated from cpaimd.C
inline CkReductionMsg *sumFastDouble(int nMsg, CkReductionMsg **msgs){

  int size0=msgs[0]->getSize();
  int size=size0/sizeof(double);

  double *inmatrix;
  //  int progcount=0;

  double *ret=(double *)msgs[0]->getData();

  for(int i=1; i<nMsg;i++)
    {
#ifdef CMK_BLUEGENEL
      fastAdd(ret, (double *) msgs[i]->getData(), size);
#else
      inmatrix=(double *) msgs[i]->getData();
	for(int d=0;d<size;d++)
	  ret[d]+=inmatrix[d];
#endif
    }
  return CkReductionMsg::buildNew(size0,ret);
}

#endif /* PAIRCALCTEST_H_ */
