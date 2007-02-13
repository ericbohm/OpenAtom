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
#include "cpaimd.h"
#include "MapTable.h"

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

  public:
    MapFile(char* name, int numpes); 
    MapFile(char* name, int num, int* size, int numpes, char *order, int x, int y, int z, int t); 
    MapFile();		// default constructor
    ~MapFile();		// destructor
     
    void setSize(int num, int* size); 
    void setAttributes(int num, int* size, char *order, int x, int y, int z, int t);
    void dumpMap(MapType2 *map);
    void dumpMap(MapType3 *map);
    void dumpMap(MapType4 *map);
    int loadMap(char *filename, MapType2 *map);
    int loadMap(char *filename, MapType3 *map);
    int loadMap(char *filename, MapType4 *map);

};

#endif

