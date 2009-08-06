/** \file FindProcessor.h
 *  Author: Abhinav S Bhatele
 *  Date Created: November 30th, 2005
 *  Finds the nearest processor from a given processor in different topologies.
 */

#ifndef _FINDPROCESSOR_H
#define _FINDPROCESSOR_H

#include <cmath>
#include <iostream>
/*class inttriple {
 private:
  int x, y,z;
 public:
    inttriple(){x=y=z=0;}

    inttriple(int _x,int _y, int _z) : x(_x), y(_y), z(_z){}
  void pup(PUP::er &p)
      {
	  p|x;
	  p|y;
	  p|z;
      }
    inline int getx(){return x;};
    inline int gety(){return y;};
    inline int getz(){return z;};
    inline CkHashCode hash() const {
	return (CkHashCode)(x+y+z);
    }
    static CkHashCode staticHash(const void *k,size_t){
	return ((inttriple *)k)->hash();
    }
    inline int compare(inttriple &t) const{
	return (t.getx() == x && t.gety() == y && t.getz() ==z);
    }
    static int staticCompare(const void *a,const void *b,size_t){
	return ((inttriple *)a)->compare((*(inttriple *)b));
    }
   
};
*/
class FindProcessor 
{
	public:
		int option;		// input about which topology is it
		int nop;		// total no. of processors
		int nopX, nopY, nopZ;	// processors in X, Y and Z dimension respectively 
		int start[3];		// starting processor
		int next[3];		// next processor w.r.t start
		int w;			// the fourth dimension for bluegene in virtual node mode
		int count;		// no. of processors output till now 
		int negXL, negYL, negZL; //necessary for FindNextInIter
		int posXL, posYL, posZL;
		int distance;
		//		CkVec <inttriple> *avoidList; // list of processors <x,y,z> to avoid
		

	public:
		FindProcessor();			// default constructor
		FindProcessor(int a[]);			// set start[] to a[] 
		void init();  // setup torus constants
		void findNext(int a[]);			// for infinite no. of processors in 3D
		int findNextInMIter(int a[]);		// helper function for 3D mesh
		int findNextInMesh(int a[]);		// for 3D mesh
		int findNextIter(int a[]);		// helper funtion for 3D torus
		int findNextInTorus(int a[]);		// for 3D torus in CO mode
		int findNextInTorusV(int t, int a[]);	// for 3D torus in VN mode
		int compare(int n, int a, int b);	// helper function
		void printSome(int n);			
		void printing(int a, int b, int c);
		void printing_sp(int a, int b, int c);
		void printing(int w, int a, int b, int c);
		//int main(int argc, char*argv[]);
};

#endif
