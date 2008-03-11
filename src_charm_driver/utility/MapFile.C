/** \file MapFile.C
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "MapFile.h"

MapFile::~MapFile()
{
  if (mapName !=NULL)
    free(mapName);
  if(sizeDim!=NULL)
    free(sizeDim);
  if(mapOrder!=NULL)
    free(mapOrder);
}
MapFile::MapFile(char* name, int numpes)
{
  mapName = (char *)malloc(sizeof(char)*15);
  strcpy(mapName, (const char*)name);
  numProcs = numpes;
}


MapFile::MapFile(char* name, int num, int* size, int numpes, char *order, int x, int y, int z, int t, int _stride)
{
  mapName = (char *)malloc(sizeof(char)*15);
  strcpy(mapName, (const char*)name);
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


void MapFile::dumpMap(MapType2 *map)
{
  FILE *fp = fopen(mapName, "w");
  fprintf(fp, "%s %d ", mapName, numDim);
  for(int i=0; i<numDim; i++)
    fprintf(fp, "%d ", sizeDim[i]);
  fprintf(fp, "\n%d \n", numProcs);
  for(int i=0; i<sizeDim[0]; i++)
    for(int j=0; j<sizeDim[1]; j++)
      fprintf(fp, "%d %d %d\n", i, j, map->get(i, j));
  fclose(fp);
}


void MapFile::dumpMap(MapType4 *map)
{
  FILE *fp = fopen(mapName, "w");
  fprintf(fp, "%s %d ", mapName, numDim);
  for(int i=0; i<numDim; i++)
    fprintf(fp, "%d ", sizeDim[i]);
  fprintf(fp, "\n%d \n", numProcs);
  for(int i=0; i<sizeDim[0]; i++)
    for(int j=0; j<sizeDim[1]*stride; j+=stride)
      for(int k=0; k<sizeDim[2]*stride; k+=stride)
	for(int l=0; l<sizeDim[3]; l++)
	  fprintf(fp, "%d %d %d %d %d\n", i, j, k, l, map->get(i, j, k, l));
  fclose(fp);
}


void MapFile::dumpMap(MapType3 *map)
{
  FILE *fp = fopen(mapName, "w");
  fprintf(fp, "%s %d ", mapName, numDim);
  for(int i=0; i<numDim; i++)
    fprintf(fp, "%d ", sizeDim[i]);
  fprintf(fp, "\n%d \n", numProcs);
  for(int i=0; i<sizeDim[0]; i++)
    for(int j=0; j<sizeDim[1]; j++)
      for(int k=0; k<sizeDim[2]; k++)
	fprintf(fp, "%d %d %d %d\n", i, j, k,  map->get(i, j, k));
  fclose(fp);
}


int MapFile::loadMap(char *filename, MapType2 *map)
{
  int x, y, pe;
  FILE *fp = fopen(filename, "r");
  if(fp==NULL)
    return 0;
  fscanf(fp, "%s%d%d%d%d", mapName, &numDim, &sizeDim[0], &sizeDim[1], &numProcs); 
  for(int i=0; i<sizeDim[0]; i++)
    for(int j=0; j<sizeDim[1]; j++)
    {
      fscanf(fp, "%d%d%d", &x, &y, &pe);
#ifdef USE_INT_MAP
      map->set(x, y, pe);
#else
      map->put(intdual(x, y))=destpe;
#endif
    }
  fclose(fp);
  CkPrintf("%s loaded from file ----\n", filename);
  return 1;
}

int MapFile::loadMap(char *filename, MapType3 *map)
{
  int x, y, z, w, pe;
  FILE *fp = fopen(filename, "r");
  if(fp==NULL)
    return 0;
  fscanf(fp, "%s%d%d%d%d%d", mapName, &numDim, &sizeDim[0], &sizeDim[1], &sizeDim[2],  &numProcs); 
  for(int i=0; i<sizeDim[0]; i++)
    for(int j=0; j<sizeDim[1]; j++)
      for(int k=0; k<sizeDim[2]; k++)
	{	
	  fscanf(fp, "%d%d%d%d", &x, &y, &z, &pe);
#ifdef USE_INT_MAP
          map->set(x, y, z, pe);
#else
	  CkArrayIndex3D idx3d(x, y, z);
	  CmiMemcpy(intidx, idx3d.index, 3*sizeof(int));
	  map->put(inttriple(intidx[0], intidx[1], intidx[2]))=destpe;
#endif
	}
  fclose(fp);
  CkPrintf("%s loaded from file ----\n", filename);
  return 1;
}


int MapFile::loadMap(char *filename, MapType4 *map)
{
  int x, y, z, w, pe;
  FILE *fp = fopen(filename, "r");
  if(fp==NULL)
    return 0;
  fscanf(fp, "%s%d%d%d%d%d%d", mapName, &numDim, &sizeDim[0], &sizeDim[1], &sizeDim[2], &sizeDim[3], &numProcs); 
  for(int i=0; i<sizeDim[0]; i++)
    for(int j=0; j<sizeDim[1]*stride; j+=stride)
      for(int k=0; k<sizeDim[2]*stride; k+=stride)
	for(int l=0; l<sizeDim[3]; l++)
	{	
	  fscanf(fp, "%d%d%d%d%d", &x, &y, &z, &w, &pe);
#ifdef USE_INT_MAP
          map->set(x, y, z, w, pe);
#else
	  CkArrayIndex4D idx4d(x, y, z, w);
	  CmiMemcpy(intidx, idx4d.index, 2*sizeof(int));
	  map->put(intdual(intidx[0], intidx[1]))=destpe;
#endif
	}
  fclose(fp);
  CkPrintf("%s loaded from file ----\n", filename);
  return 1;
}

