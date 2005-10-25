/** \file StructFactorCache.h
 *
 */
 
#ifndef _StructFactorCache_h_
#define _StructFactorCache_h_

class PPDummyMsg: public CMessage_PPDummyMsg {
 public:
  int atmGrp;
  int sfindex;
};


class StructFactorMsg: public CkMcastBaseMsg, public CMessage_StructFactorMsg {
 public:
	int datalen;
	int atmGrpIndex;
	int gsSize;
	int planeIndex; 
	complex *structFactor; 
	complex *structFactor_fx, *structFactor_fy, *structFactor_fz;
};

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

class StructFactCache : public Group {
 public:
    int numSfGrps,natm_nl,natm_nl_grp_max, totalsize;

    StructFactCache(int numSfGrps_in,int natm_nl_in,int natm_nl_grp_max_in) 
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
 private:
    CkVec<GSlabInfo> gSpaceSlabs;
    CkVec<int*> structFactorAtmGrps; 
    CkVec<complex*> structFactorList; 
    CkVec<complex*> structFactorfxList;
    CkVec<complex*> structFactorfyList;
    CkVec<complex*> structFactorfzList;   
    CkVec<int> structFactorSize;
    CkVec<PlaneCount> planeCountList; 
    CkVec <PlaneAtom> ppList;
 
    int existStructFact(int planeIndex);
    int incCountStructFact(int planeIndex);
    int decCountStructFact(int planeIndex);
 public:
    void printCountStructFact();
};

#endif //_StructFactorCache_h_
