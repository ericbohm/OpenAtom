/*****************************************************************************
 * $Source$
 * $Author$
 * $Date$
 * $Revision$
 *****************************************************************************/

/** \file TorusGen.C
 *  Author: Abhinav S Bhatele
 *  Date Created: July 10th, 2008
 *
 */

#include <iostream.h>
#include <fstream.h>

#define CUBE_SIZE	2
#define DISTANCE	10

int main(int argc, char *argv[]) {
  int dimX, dimY, dimZ, numPoints, pe;

  if(argc != 4) {
    cout << "Dimensions of the torus required!!!\n";
    return 0;
  }
    
  dimX = atoi(argv[1]);
  dimY = atoi(argv[2]);
  dimZ = atoi(argv[3]);

  ofstream wr("Torus-Mesh.vtk", ios::out);

  // write the header
  wr << "# vtk DataFile Version 2.0\n";
  wr << "Torus\n";
  wr << "ASCII\n";
  wr << "DATASET POLYDATA\n";

  // points
  numPoints = dimX * dimY * dimZ;
  wr << "\nPOINTS " << numPoints*8 << " float\n";

  for(int k=0; k < DISTANCE*dimZ; k=k+DISTANCE)
    for(int j=0; j < DISTANCE*dimY; j=j+DISTANCE) 
      for(int i=0; i < DISTANCE*dimX; i=i+DISTANCE) {
	// this whole print is for a single processor
	wr << i		  << " " << j		<< " " << k	      << "    ";
	wr << i+CUBE_SIZE << " " << j		<< " " << k	      << "    ";
	wr << i   	  << " " << j+CUBE_SIZE << " " << k	      << "    ";
	wr << i+CUBE_SIZE << " " << j+CUBE_SIZE << " " << k	      << "    ";
	wr << i   	  << " " << j		<< " " << k+CUBE_SIZE << "    ";
	wr << i+CUBE_SIZE << " " << j		<< " " << k+CUBE_SIZE << "    ";
	wr << i		  << " " << j+CUBE_SIZE	<< " " << k+CUBE_SIZE << "    ";
	wr << i+CUBE_SIZE << " " << j+CUBE_SIZE	<< " " << k+CUBE_SIZE << "\n";
      }

  // polygons
  wr << "\nPOLYGONS " << numPoints*6 << " " << numPoints*30 << "\n";
  for(int k=0; k < dimZ; k=k+1)
    for(int j=0; j < dimY; j=j+1) 
      for(int i=0; i < dimX; i=i+1) {
	// this whole print is for a single processor
	
	pe = (i + j*dimX + k*dimX*dimY) * 8;
	wr << "4 " << pe+0 << " " << pe+1 << " " << pe+3 << " " << pe+2 << "    ";
	wr << "4 " << pe+4 << " " << pe+5 << " " << pe+7 << " " << pe+6 << "    ";
	wr << "4 " << pe+0 << " " << pe+1 << " " << pe+5 << " " << pe+4 << "    ";
	wr << "4 " << pe+2 << " " << pe+3 << " " << pe+7 << " " << pe+6 << "    ";
	wr << "4 " << pe+1 << " " << pe+3 << " " << pe+7 << " " << pe+5 << "    ";
	wr << "4 " << pe+0 << " " << pe+2 << " " << pe+6 << " " << pe+4 << "\n";
      }
 
  wr.close();
  
}
