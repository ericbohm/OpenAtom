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

#include "../../include/eesDataClass.h"

//============================================================================
// Group Container class : Only allowed chare data classes have data.
//============================================================================
class eesCache: public Group {
 public:
 public:
  // Variables
   int itimeRPP;
   int itimeRHart;

   int nstates;
   int nchareGSPProcT;           // Total number of GSP on proc
   int *gspStateInd;             // state index of gsp chare on proc
   int *gspPlaneInd;             // state index of gsp chare on proc
   int nMallSize;

   int nchareGPP;              // Chare Plane array index sizes
   int nchareRPP;
   int nchareGHart;
   int nchareRHart;
   int nchareGSP;
   int nchareRhoG;

   int nchareGPPProc;
   int nchareRPPProc;          // # of plane chares of each type on this proc
   int nchareGHartProc;        // two different states with the same plane #
   int nchareRHartProc;        // assigned to this proc counts once here.
   int nchareGSPProc;
   int nchareRhoGProc;

   // planes chares on proc : again ignore state index
   int *allowedGppChares;      // lth: nchareGPP   : 1/0 if chare is on/off proc
   int *allowedRppChares;      // lth: nchareRPP   : 1/0 if chare is on/off proc
   int *allowedRhoGHartChares; // lth: nchareGHart : 1/0 if chare is on/off proc
   int *allowedRhoRHartChares; // lth: nchareRHart : 1/0 if chare is on/off proc
   int *allowedGspChares; // lth: nchareRHart : 1/0 if chare is on/off proc
   int *allowedRhoGChares; // lth: nchareRHart : 1/0 if chare is on/off proc

   // planes chares on proc : again ignore state index
   int *indGppChares;          // lth: nchareGPPProc  : ind[5]=27 means the 5th chare
   int *indRppChares;          // lth: nchareRPPProc  : on proc is plane index 27
   int *indRhoGHartChares;     // lth: nchareGHartProc: 
   int *indRhoRHartChares;     // lth: nchareRHartProc: 
   int *indGspChares;          // lth: nchareGspProc: 
   int *indRhoGChares;         // lth: nchareRhoGProc: 

   // plane data on proc : again ignore state index
   GPPDATA      *GppData;      // lth: nchareGPP  : only allowed guys have data
   RPPDATA      *RppData;      // lth: nchareRPP  : only allowed guys have data
   RHOGHARTDATA *RhoGHartData; // lth: nchareGHart: only allowed guys have data
   RHORHARTDATA *RhoRHartData; // lth: nchareRHart: only allowed guys have data
   GSPDATA      *GspData;      // lth: nchareGsp  : only allowed guys have data
   RHOGDATA     *RhoGData;     // lth: nchareRhoG : only allowed guys have data

  // Functions
   eesCache(){};
  ~eesCache(){};
   eesCache(int , int , int ,int ,int,int);

   // Report into the cache
   void registerCacheGPP  (int ,int, int *,int *,int*);
   void registerCacheRPP  (int );
   void registerCacheGHart(int ,int, int *,int *,int*);
   void registerCacheRHart(int );
   void registerCacheGSP(int,int);
   void registerCacheRHOG(int );

   // Ask the cache if it is warm
   void queryCacheRPP  (int ,int ,int );
   void queryCacheRHart(int ,int ,int );

};
//============================================================================


//============================================================================

#endif // eesCache not yet defined
