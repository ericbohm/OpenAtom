#include "debug_flags.h"
#include "structureFactorCache.decl.h"
#include "StructureFactorMessages.h"
 
#ifndef _StructFactorCache_h_
#define _StructFactorCache_h_

class PlaneAtom {
 public:
  int plane;
  int atom;
  CkVec <CkArrayIndex2D> particles;
  PlaneAtom(){}
  bool operator==(const PlaneAtom& obj) const {
    return (obj.plane==plane && obj.atom ==atom);
  }
  PlaneAtom(int _plane, int _atom): plane(_plane), atom(_atom) {}
  void pup(PUP::er &p) {
    p|plane;
    p|atom;
    p|particles;
  }
      
};

class PlaneCount {
 public: 
    int plane;
    int count;
    int state;
    int updated;
    PlaneCount() {}
    PlaneCount(int p, int c, int _state) : plane(p), count(c), state(_state) { updated = 0;}
};
PUPbytes(PlaneCount);

#include "fft_slab_ctrl/fftCacheSlab.h"

class StructFactCache : public Group {
 public:
    void printCountStructFact();
    int numSfGrps,natm_nl,natm_nl_grp_max, totalsize;
    const UberCollection thisInstance;
    StructFactCache(int numSfGrps_in,int natm_nl_in,int natm_nl_grp_max_in,
		    UberCollection _thisInstance) :thisInstance(_thisInstance)
      {
	numSfGrps       = numSfGrps_in;
	natm_nl         = natm_nl_in;
	natm_nl_grp_max = natm_nl_grp_max_in;
      }
    int getStructFact(int planeIndex, int atmGrp, complex **sf, complex **sf_x, complex **sf_y, complex **sf_z);
    void getStructFactIdx(int sfindex, int atmGrp, complex **sf, complex **sf_x, complex **sf_y, complex **sf_z);
    void acceptStructFact(StructFactorMsg *msg);
    void removeAll();
    void setZero(int); 
    int existStructFactGrp(int planeIndex, int atmGrp);
    int registerPP(int state, int plane, int atmGrp);
    int existsPP(int plane, int atmGrp);

    void pup(PUP::er &p) {
       p|numSfGrps;  p|natm_nl;  p|natm_nl_grp_max, p|totalsize;
    }//routine

    CkVec<GSlabInfo> gSpaceSlabs;
    CkVec<int*> structFactorAtmGrps; 
    CkVec<complex*> structFactorList; 
    CkVec<complex*> structFactorfxList;
    CkVec<complex*> structFactorfyList;
    CkVec<complex*> structFactorfzList;   
    CkVec<int> structFactorSize;
    CkVec<PlaneCount> planeCountList; 
    CkVec <PlaneAtom> ppList;
 private: 
    int existStructFact(int planeIndex);
    int incCountStructFact(int planeIndex);
    int decCountStructFact(int planeIndex);


};

#endif //_StructFactorCache_h_
