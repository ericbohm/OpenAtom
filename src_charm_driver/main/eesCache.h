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
  // Variables
   int itimeRPP;
   int itimeRHart;

   int nchareGPP;
   int nchareRPP;              // Chare array sizes
   int nchareGHart;
   int nchareRHart;

   int nchareGPPProc;
   int nchareRPPProc;          // # of chares of each type on this proc
   int nchareGHartProc;
   int nchareRHartProc;

   int *allowedGppChares;      // lth: nchareGPP   : 1/0 if chare is on/off proc
   int *allowedRppChares;      // lth: nchareRPP   : 1/0 if chare is on/off proc
   int *allowedRhoGHartChares; // lth: nchareGHart : 1/0 if chare is on/off proc
   int *allowedRhoRHartChares; // lth: nchareRHart : 1/0 if chare is on/off proc

   int *indGppChares;          // lth: nchareGPPProc  : ind[5]=27 means the 5th 
   int *indRppChares;          // lth: nchareRPPProc  : allowed chare is index 27
   int *indRhoGHartChares;     // lth: nchareGHartProc: 
   int *indRhoRHartChares;     // lth: nchareRHartProc: 

   GPPDATA      *GppData;      // lth: nchareGPP  : only allowed guys have data
   RPPDATA      *RppData;      // lth: nchareRPP  : only allowed guys have data
   RHOGHARTDATA *RhoGHartData; // lth: nchareGHart: only allowed guys have data
   RHORHARTDATA *RhoRHartData; // lth: nchareRHart: only allowed guys have data

  // Functions
   eesCache(){};
  ~eesCache(){};
   eesCache(int , int , int ,int );

   void registerCacheGPP  (int ,int, int *,int *,int*);
   void registerCacheRPP  (int );
   void registerCacheGHart(int ,int, int *,int *,int*);
   void registerCacheRHart(int );

   void queryCacheRPP  (int ,int ,int );
   void queryCacheRHart(int ,int ,int );

};
//============================================================================

#endif // eesCache not yet defined
