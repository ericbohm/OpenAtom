//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \file eesCache.h
 *
 *
 */
//============================================================================

#ifndef _eesCache_h_
#define _eesCache_h_

#include "eesDataClass.h"
#include "fft_types.h"
#include <vector>

//============================================================================
// Group Container class : Only allowed chare data classes have data.
//============================================================================

struct KthSpinStatePlaneTuple {
 public:
  int kpoint, spin, state, plane;
  KthSpinStatePlaneTuple() {}
  KthSpinStatePlaneTuple(int _kpoint, int _spin, int _state, int _plane) : kpoint(_kpoint), spin(_spin), state(_state), plane(_plane) {}
  inline int getKpoint(){return kpoint;}
  inline int getSpin(){return spin;}
  inline int getState(){return state;}
  inline int getPlane(){return plane;}
  inline int setKpoint(int _kpoint){kpoint=_kpoint;return kpoint;}
  inline int setSpin(int _spin){spin=_spin;return spin;}
  inline int setState(int _state){state=_state;return state;}
  inline int setPlane(int _plane){plane=_plane;return plane;}
};

class eesCache: public Group {
  public:
    // Variables
    const UberCollection thisInstance;
    int itimeRPP;
    int itimeRHart;
    int rpp_on;

    int nstates;
    int nchareGSPProcT;           // Total number of GSP on proc
    CkVec <KthSpinStatePlaneTuple> gspKptSpinStatePlaneVec;             // state index of gsp chare on proc
    int nMallSize;

    int nchareGPP;              // Chare Plane array index sizes
    int nchareRPP;
    int nchareGSP;

    int nchareGPPProc;
    int nchareRPPProc;          // # of plane chares of each type on this proc
    int nchareGSPProc;
    int nchareRhoGProc;
    int nkpoint;

    // planes chares on proc : again ignore state index
    int *allowedGppChares;      // lth: nchareGPP   : 1/0 if chare is on/off proc
    int *allowedRppChares;      // lth: nchareRPP   : 1/0 if chare is on/off proc
    int *allowedGspChares; // lth: nchareRHart : 1/0 if chare is on/off proc

    // planes chares on proc : again ignore state index
    int *indGppChares;          // lth: nchareGPPProc  : ind[5]=27 means the 5th chare
    int *indRppChares;          // lth: nchareRPPProc  : on proc is plane index 27
    int *indGspChares;          // lth: nchareGspProc: 

    // plane data on proc : again ignore state index
    GPPDATA      **GppData;      // lth: nchareGPP  : only allowed guys have data
    RPPDATA      **RppData;      // lth: nchareRPP  : only allowed guys have data
    GSPDATA      **GspData;      // lth: nchareGsp  : only allowed guys have data
    std::vector<RHOGHARTDATA> RhoGHartData; // lth: nchareGHart: only allowed guys have data
    std::vector<RHORHARTDATA> RhoRHartData; // lth: nchareRHart: only allowed guys have data

    // Functions
    eesCache(){};
    ~eesCache(){};
    eesCache(int, int, int, int, UberCollection thisInstance);

    // Report into the cache
    void registerCacheGPP  (int ,int, int *,int *,int*);
    void registerCacheRPP  (int );
    int registerCacheGHart(int, std::vector< gridPoint > * myPoints);
    int registerCacheRHart(int, int, int, int, int, int );
    void registerCacheGSP(int, int, int,int);
    void registerCacheRHOG(int );

    // Ask the cache if it is warm
    void queryCacheRPP  (int ,int ,int );
    void queryCacheRHart(int ,int ,int );

};
//============================================================================


//============================================================================

#endif // eesCache not yet defined
