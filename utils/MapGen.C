/*****************************************************************************
 * $Source$
 * $Author$
 * $Date$
 * $Revision$
 *****************************************************************************/

/** \file MapGen.C
 *  Author: Abhinav S Bhatele
 *  Date Created: July 10th, 2008
 *
 */

#include <iostream.h>
#include <fstream.h>

#define CUBE_SIZE       2
#define DISTANCE        10
#define OBJS_PER_ROW	10
#define DIST_OBJS	(float)CUBE_SIZE/11.0

int main(int argc, char *argv[]) {
  char *fn, mapName[30];
  int dims, numStates, numPlanes, numPes, numNodes, dimX,dimY, dimZ, dimT;
  int opt;

  if(argc != 3) {
    cout << "Enter OPTION (1 or 2) and the name of the map file!!!\n";
    return 0;
  }

  opt = atoi(argv[1]);
  fn = argv[2];
  ifstream rd (fn, ios::in);

  // read the map file first
  // header
  rd >> mapName >> dims >> numStates >> numPlanes;
  rd >> numPes >> dimX >> dimY >> dimZ >> dimT; 

  cout << mapName << " " << dims << " " << numStates << " " << numPlanes << " " << numPes << "\n";

  int map[numStates][numPlanes][5];
  int state, plane, pe, x, y, z, t;

  for(int i=0; i<numStates; i++)
    for(int j=0; j<numPlanes; j++) {
      rd >> state >> plane >> pe >> x >> y >> z >> t;
      map[state][plane][0] = pe;
      map[state][plane][1] = x;
      map[state][plane][2] = y;
      map[state][plane][3] = z;
      map[state][plane][4] = t;
  
      // cout << state << " " << plane << " " << pe << " " << map[state][plane][1] << " " << map[state][plane][2] << " " << map[state][plane][3] << " " << map[state][plane][4] << "\n";  
    }

  numNodes = numPes/dimT;
  int objsPerPe[numNodes];

  for(int i=0; i<numPes; i++)
    objsPerPe[i] = 0;

  ofstream wr;

  if(opt == 1)
    wr.open("GSMap.vtk", ios::out);
  else
    wr.open("RSMap1.vtk", ios::out);

  // write the header
  wr << "# vtk DataFile Version 2.0\n";
  if(opt == 1)
    wr << "GSMap\n";
  else
    wr << "RSMap\n";
  wr << "ASCII\n";
  wr << "DATASET POLYDATA\n";

  float ox, oy, oz;
  float fx, fy, fz;
  int numPoints;
  
  // points
  numPoints = numStates * numPlanes; 
  wr << "\nPOINTS " << numPoints << " float\n";

  for(int i=0; i<numStates; i++) {
    for(int j=0; j<numPlanes; j++) {
      pe = map[i][j][1] + map[i][j][2]*dimX + map[i][j][3]*dimX*dimY;

      if(opt == 1) {
	ox = map[i][j][1] * DISTANCE;
	oy = map[i][j][2] * DISTANCE + DIST_OBJS;
	oz = map[i][j][3] * DISTANCE + DIST_OBJS;

	fx = ox;
	fy = oy + ((int)(objsPerPe[pe] / OBJS_PER_ROW) * DIST_OBJS);
	fz = oz + ((int)(objsPerPe[pe] % OBJS_PER_ROW) * DIST_OBJS);
      } else {
	ox = map[i][j][1] * DISTANCE + DIST_OBJS;
	oy = map[i][j][2] * DISTANCE + DIST_OBJS;
	oz = map[i][j][3] * DISTANCE;

	fx = ox + ((int)(objsPerPe[pe] % OBJS_PER_ROW) * DIST_OBJS);
	fy = oy + ((int)(objsPerPe[pe] / OBJS_PER_ROW) * DIST_OBJS);
	fz = oz;
      }

      objsPerPe[pe]++;
      // cout << map[i][j][1] << " " << map[i][j][2] << " " << map[i][j][3] << " " << DISTANCE << "\n";
      wr << fx << " " << fy << " " << fz << "    ";
    }
    wr << "\n";
  }

  // cells
  wr << "\nVERTICES " << numPoints << " " << numPoints*2 << " \n";
  for(int i=0; i<numPoints; i++)
    wr << "1 " << i << " ";
  wr << "\n";

  // point data
  wr << "\nPOINT_DATA " << numPoints << "\n";
  wr << "SCALARS thisIndex.x Int32 1\n";
  wr << "LOOKUP_TABLE default\n";

  for(int i=0; i<numStates; i++) {
    for(int j=0; j<numPlanes; j++) {
      wr << j << " ";
    }
    wr << "\n";
  }

  wr.close();
}
