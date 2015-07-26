/** \file MapFile.h
 *  Author: Abhinav S Bhatele
 *  Date Created: December 28th, 2006
 *
 *  This class is used for dumping maps to a file during a program run
 *  and also for loading them for use in initial mapping.
 */

#ifndef _MAPFILE_H_
#define _MAPFILE_H_

#include "../../include/debug_flags.h"
#include "charm++.h"
#include "load_balance/MapTable.h"

class MapFile
{
  private:
    char* mapName;	// name of the map
    int numDim;		// number of dimensions in this map
    int* sizeDim;  	// array of size numDim, size of each dimension

    int numProcs;	// number of processors
    char* mapOrder;	// is the mapping TXYZ or XYZT
    int Xmax;
    int Ymax;
    int Zmax;
    int Tmax;
    int stride;         // stride for non dense indices NOTE evil
    // hardcoding assuming 4D uses stride in inner
    // 2 indices.  There should be a stride for
    // each index.

  public:
    MapFile(const char* name, int numpes);
    MapFile(const char* name, int num, int* size, int numpes, const char *order, int x, int y, int z, int t, int stride=1);
    MapFile();		// default constructor
    ~MapFile();		// destructor

    void setSize(int num, int* size);
    void setAttributes(int num, int* size, char *order, int x, int y, int z, int t, int stride);
    void dumpMap(MapType1 *map, char c) { }
    void dumpMap(MapType2 *map, char c);
    void dumpMap(MapType3 *map, char c);
    void dumpMap(MapType4 *map, char c);

    void dumpMapCoords(MapType1 *map, char c) { }
    void dumpMapCoords(MapType2 *map, char c);
    void dumpMapCoords(MapType3 *map, char c);
    void dumpMapCoords(MapType4 *map, char c);

    int loadMap(const char *filename, MapType1 *map) {
      CkAbort("Loading 1D maps not supported yet.\n");
      return 0;
    }
    int loadMap(const char *filename, MapType2 *map);
    int loadMap(const char *filename, MapType3 *map);
    int loadMap(const char *filename, MapType4 *map);

};

#endif

