/*
** PIBeadAtoms.h
** 
** Made by Eric Bohm
** Login   <bohm@alacrity>
** 
** Started on  Mon Mar  1 10:35:54 2010 Eric Bohm
** Last update Mon Mar  1 10:35:54 2010 Eric Bohm
* 
* PIBeadAtoms is the helper array which accepts the transpose of bead-wise
* atomic positions and forces into atom-wise form.  Atom-wise
* operations are then performed and the result is returned to the
* usual bead-wise atom structures.

* Notionally this is the F_x -> F_u -> F_x process where F_u is
* constructed in PIBeadAtoms

* PIBeadAtoms is triggered during Atom force integration when path
* integrals (UberIMax > 0) are in use.  

* PIBeadAtoms is initialized with atom positions during startup

* PIBeadAtoms receives and sends force contributions each step

* PIBeadAtoms doesn't really have persistent state.  Its purpose is as
* a temporary space for the F_u operations.  As such, it is entirely
* data driven and could theoretically be created and destroyed every
* step, or even computed in the atoms group using some arbitrary
* partioning scheme.  In practice we don't do that because the
* initialization isn't free of cost and it is easier to debug things
* when they have an explicit decomposition.
*
* PIBeadAtoms has no internal synchronization.
*/

#ifndef   	PIBEAD_H_
# define   	PIBEAD_H_



#include "PIBeadAtoms.decl.h"
#include "uber/Uber.h"

class PIBeadAtoms : public CBase_PIBeadAtoms
{
 public:
	PIBeadAtoms(CkMigrateMessage *m) {}
        PIBeadAtoms(UberCollection _thisInstance, int _numBeads);
	void accept_PIMD_Fx(AtomXYZMsg *msg);
	void accept_PIMD_x(double _x, double _y, double _z, int PIBeadIndex);
	void accept_PIMD_u(double _xu, double _yu, double _zu, int PIBeadIndex);
	~PIBeadAtoms(){}
 private:
	const UberCollection thisInstance;
	int numBeads;
	void compute_PIMD_Fu(){}
	void compute_PIMD_u(){}
	void compute_PIMD_x(){}
	// each of these arrays is of length numBeads
	double *x,*y,*z;
	double *xu,*yu,*zu;
	double *fx,*fy,*fz;
	double *fxu,*fyu,*fzu;
	int acceptCount_Fx;
	int acceptCount_u;
	int acceptCount_x;
};




#endif 	    /* !PIBEAD_H_ */
