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
* partioning scheme.  In practice we don't do that because setup isn't
* free of cost and it is easier to debug things when they have an
* explicit decomposition.

*/

#ifndef   	PIBEAD_H_
# define   	PIBEAD_H_

class PIBeadAtoms : public CBase_PIBeadAtoms
{
 public:
	PIBeadAtoms(CkMigrateMessage *m) {}
	PIBeadAtoms(UberCollection) {}
	void acceptForces(){}
	void acceptCoords(){}
	~PIBeadAtoms(){}
 private:
	void compute_Fu(){}
	void integrate_Fu(){}
	void send_force(){}
};




#endif 	    /* !PIBEAD_H_ */
