/** \file FindProcessor.h
 * Finds the nearest processor from a given processor in different topologies.
 *
 */
#ifndef _FINDPROCESSOR_H
#define _FINDPROCESSOR_H


#include <math.h>
#include <iostream.h>

class FindProcessor 
{
	public:
		int option;
		int nop;
		int nopX, nopY, nopZ;
		int start[3];
		int next[3];
		int w;
		int count;
	public:
		FindProcessor();
		FindProcessor(int a[]);
		int main(int argc, char*argv[]);
		void findNext(int a[]);
		void findNextInBluegene(int a[]);
		int findNextInTorus(int a[]);
		int findNextInTorusV(int t, int a[]);
		void printSome(int n);
		void printing(int a, int b, int c);
		void printing_sp(int a, int b, int c);
		void printing(int w, int a, int b, int c);
		int compare(int n, int a, int b);
};

#endif
