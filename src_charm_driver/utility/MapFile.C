/** \file MapFile.C
 *
 */

#include "TopoManager.h"
#include "MapFile.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <assert.h>
extern TopoManager *topoMgr;

MapFile::~MapFile()
{
  if (mapName !=NULL)
    free(mapName);
  if(sizeDim!=NULL)
    free(sizeDim);
  if(mapOrder!=NULL)
    free(mapOrder);
}
MapFile::MapFile(const char* name, int numpes)
{
  if(topoMgr==NULL) 
    topoMgr= new TopoManager(); // just to keep dump coord from flaking
  mapName = (char *)malloc(sizeof(char)*15);
  strcpy(mapName, name);
  numProcs = numpes;
}


MapFile::MapFile(const char* name, int num, int* size, int numpes, const char *order, int x, int y, int z, int t, int _stride)
{
  if(topoMgr==NULL) 
    topoMgr= new TopoManager(); // just to keep dump coord from flaking
  mapName = (char *)malloc(sizeof(char)*15);
  strcpy(mapName, name);
  numDim = num;
  sizeDim = (int *)malloc(sizeof(int)*numDim);
  for(int i=0; i<numDim; i++)
    sizeDim[i] = size[i];
  numProcs = numpes;
  mapOrder = (char *)malloc(sizeof(char)*5);
  strcpy(mapOrder, (const char*)order);
  Xmax = x;
  Ymax = y;
  Zmax = z;
  Tmax = t;
  stride=_stride;
}


void MapFile::setSize(int num, int* size)
{
  numDim = num;
  if(sizeDim==NULL)
    sizeDim = (int *)malloc(sizeof(int)*numDim);
  for(int i=0; i<numDim; i++)
    sizeDim[i] = size[i];
}


void MapFile::setAttributes(int num, int* size, char *order, int x, int y, int z, int t, int _stride)
{
  numDim = num;
  if(sizeDim==NULL)
    sizeDim = (int *)malloc(sizeof(int)*numDim);
  for(int i=0; i<numDim; i++)
    sizeDim[i] = size[i];
  if(mapOrder==NULL)
    mapOrder = (char *)malloc(sizeof(char)*5);
  strcpy(mapOrder, (const char*)order);
  Xmax = x;
  Ymax = y;
  Zmax = z;
  Tmax = t;
  stride=_stride;
}


void MapFile::dumpMapCoords(MapType2 *map, char c)
{
  char name[100];
  sprintf(name, "%s_inst_%d", mapName, (int)c);
  FILE *fp = fopen(name, "w");
  fprintf(fp, "%s %d ", mapName, numDim);
  for(int i=0; i<numDim; i++)
    fprintf(fp, "%d ", sizeDim[i]);
  fprintf(fp, "\n%d %d %d %d %d\n", numProcs, topoMgr->getDimNX(), topoMgr->getDimNY(), topoMgr->getDimNZ(), topoMgr->getDimNT());
  for(int i=0; i<sizeDim[0]; i++)
    for(int j=0; j<sizeDim[1]; j++)
    {
      int proc=map->get(i, j);
      int X, Y, Z, T;
      topoMgr->rankToCoordinates(proc, X, Y, Z, T);
      fprintf(fp, "%d %d %d %d %d %d %d\n", i, j, proc,X,Y,Z,T);
    }
  fclose(fp);
}


void MapFile::dumpMapCoords(MapType4 *map, char c)
{
  char name[100];
  sprintf(name, "%s_inst_%d", mapName, (int)c);
  FILE *fp = fopen(name, "w");
  fprintf(fp, "%s %d ", mapName, numDim);
  for(int i=0; i<numDim; i++)
    fprintf(fp, "%d ", sizeDim[i]);
  fprintf(fp, "\n%d %d %d %d %d\n", numProcs, topoMgr->getDimNX(), topoMgr->getDimNY(), topoMgr->getDimNZ(), topoMgr->getDimNT());
  for(int i=0; i<sizeDim[0]; i++)
    for(int j=0; j<sizeDim[1]*stride; j+=stride)
      for(int k=0; k<sizeDim[2]*stride; k+=stride)
        for(int l=0; l<sizeDim[3]; l++)
        {
          int proc=map->get(i, j, k, l);
          int X, Y, Z, T;
          topoMgr->rankToCoordinates(proc, X, Y, Z, T);

          fprintf(fp, "%d %d %d %d %d %d %d %d %d\n", i, j, k, l, proc, X, Y, Z, T);
        }
  fclose(fp);
}


void MapFile::dumpMapCoords(MapType3 *map, char c)
{
  char name[100];
  sprintf(name, "%s_inst_%d", mapName, (int)c);
  FILE *fp = fopen(name, "w");
  fprintf(fp, "%s %d ", mapName, numDim);
  for(int i=0; i<numDim; i++)
    fprintf(fp, "%d ", sizeDim[i]);
  fprintf(fp, "\n%d %d %d %d %d\n", numProcs, topoMgr->getDimNX(), topoMgr->getDimNY(), topoMgr->getDimNZ(), topoMgr->getDimNT());
  for(int i=0; i<sizeDim[0]; i++)
    for(int j=0; j<sizeDim[1]; j++)
      for(int k=0; k<sizeDim[2]; k++)
      {
        int proc=map->get(i, j, k);
        int X, Y, Z, T;
        topoMgr->rankToCoordinates(proc, X, Y, Z, T);

        fprintf(fp, "%d %d %d %d %d %d %d %d\n", i, j, k,  proc, X, Y ,Z, T);
      }
  fclose(fp);
}


void MapFile::dumpMap(MapType2 *map, char c)
{
  char name[100];
  sprintf(name, "%s_inst_%d", mapName, (int)c);
  FILE *fp = fopen(name, "w");
  fprintf(fp, "%s %d ", mapName, numDim);
  for(int i=0; i<numDim; i++)
    fprintf(fp, "%d ", sizeDim[i]);
  fprintf(fp, "\n%d %d %d %d %d\n", numProcs, topoMgr->getDimNX(), topoMgr->getDimNY(), topoMgr->getDimNZ(), topoMgr->getDimNT());
  for(int i=0; i<sizeDim[0]; i++)
    for(int j=0; j<sizeDim[1]; j++)
      fprintf(fp, "%d %d %d\n", i, j, map->get(i, j));
  fclose(fp);
}


void MapFile::dumpMap(MapType4 *map, char c)
{
  char name[100];
  sprintf(name, "%s_inst_%d", mapName, (int)c);
  FILE *fp = fopen(name, "w");
  fprintf(fp, "%s %d ", mapName, numDim);
  for(int i=0; i<numDim; i++)
    fprintf(fp, "%d ", sizeDim[i]);
  fprintf(fp, "\n%d %d %d %d %d\n", numProcs, topoMgr->getDimNX(), topoMgr->getDimNY(), topoMgr->getDimNZ(), topoMgr->getDimNT());
  for(int i=0; i<sizeDim[0]; i++)
    for(int j=0; j<sizeDim[1]*stride; j+=stride)
      for(int k=0; k<sizeDim[2]*stride; k+=stride)
        for(int l=0; l<sizeDim[3]; l++)
          fprintf(fp, "%d %d %d %d %d\n", i, j, k, l, map->get(i, j, k, l));
  fclose(fp);
}


void MapFile::dumpMap(MapType3 *map, char c)
{
  char name[100];
  sprintf(name, "%s_inst_%d", mapName, (int)c);
  FILE *fp = fopen(name, "w");
  fprintf(fp, "%s %d ", mapName, numDim);
  for(int i=0; i<numDim; i++)
    fprintf(fp, "%d ", sizeDim[i]);
  fprintf(fp, "\n%d %d %d %d %d\n", numProcs, topoMgr->getDimNX(), topoMgr->getDimNY(), topoMgr->getDimNZ(), topoMgr->getDimNT());
  for(int i=0; i<sizeDim[0]; i++)
    for(int j=0; j<sizeDim[1]; j++)
      for(int k=0; k<sizeDim[2]; k++)
        fprintf(fp, "%d %d %d %d\n", i, j, k,  map->get(i, j, k));
  fclose(fp);
}


int MapFile::loadMap(const char *filename, MapType2 *map)
{
  int x, y, pe;
  FILE *fp = fopen(filename, "r");
  if(fp==NULL)
    return 0;
  assert(fscanf(fp, "%s%d%d%d%d", mapName, &numDim, &sizeDim[0], &sizeDim[1], &numProcs)); 
  for(int i=0; i<sizeDim[0]; i++)
    for(int j=0; j<sizeDim[1]; j++)
    {
      assert(fscanf(fp, "%d%d%d\n", &x, &y, &pe));
      map->set(x, y, pe);
    }
  fclose(fp);
  CkPrintf("%s loaded from file ----\n", filename);
  return 1;
}

int MapFile::loadMap(const char *filename, MapType3 *map)
{
  int x, y, z, w, pe;
  FILE *fp = fopen(filename, "r");
  if(fp==NULL)
    return 0;
  assert(fscanf(fp, "%s%d%d%d%d%d%d", mapName, &numDim, &sizeDim[0], &sizeDim[1], &sizeDim[2], &sizeDim[3], &numProcs)); 
  for(int i=0; i<sizeDim[0]; i++)
    for(int j=0; j<sizeDim[1]; j++)
      for(int k=0; k<sizeDim[2]; k++)
      {	
        assert(fscanf(fp, "%d%d%d%d", &x, &y, &z, &pe));
        map->set(x, y, z, pe);
      }
  fclose(fp);
  CkPrintf("%s loaded from file ----\n", filename);
  return 1;
}


int MapFile::loadMap(const char *filename, MapType4 *map)
{
  int x, y, z, w, pe;
  FILE *fp = fopen(filename, "r");
  if(fp==NULL)
    return 0;
  assert(fscanf(fp, "%s%d%d%d%d%d%d", mapName, &numDim, &sizeDim[0], &sizeDim[1], &sizeDim[2], &sizeDim[3], &numProcs)); 
  for(int i=0; i<sizeDim[0]; i++)
    for(int j=0; j<sizeDim[1]*stride; j+=stride)
      for(int k=0; k<sizeDim[2]*stride; k+=stride)
        for(int l=0; l<sizeDim[3]; l++)
        {	
          assert(fscanf(fp, "%d%d%d%d%d", &x, &y, &z, &w, &pe));
          map->set(x, y, z, w, pe);
        }
  fclose(fp);
  CkPrintf("%s loaded from file ----\n", filename);
  return 1;
}

