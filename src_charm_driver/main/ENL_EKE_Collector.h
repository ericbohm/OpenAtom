/*
** ENL_EKE_Collector.h
** 
** Made by Eric Bohm
** Login   <bohm@alacrity>
** 
** Started on  Apr 14 10:35:54 2010 Eric Bohm
** Last update Apr 14 10:35:54 2010 Eric Bohm
* 
* ENL_EKE_Collector sums the ENL and EKE energies across
* k-points
* spin orbitals
* path-integrals
* note: each tempering would have its own ENL and EKE
* so the number of elements for this array is the number of tempers

* this object is so simple that making it a chare array borders on
* overkill. It receives messages, sums them, and prints.  So its only
* output is a side effect.
*
*/

#ifndef   	ENL_EKE_Collector_H_
# define   	ENL_EKE_Collector_H_



#include "ENL_EKE_Collector.decl.h"

class ENL_EKE_Collector : public CBase_ENL_EKE_Collector
{
 public:
	ENL_EKE_Collector(CkMigrateMessage *m) {}
        ENL_EKE_Collector(int _numEnergyInputs);
	void acceptENL(double _enl);
	void acceptEKE(double _eke);
	~ENL_EKE_Collector(){}
 private:
	int enlIteration, ekeIteration;
	FILE *temperScreenFile;
	void printENL();
	void printEKE();
	int energyExpected;
	int countEKE;
	int countENL;
	double EKE;
	double ENL;
};



#endif 	    /* !ENL_EKE_Collector_H_ */
