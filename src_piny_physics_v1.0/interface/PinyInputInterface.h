#ifndef _PINYINPUTINTERFACE_H_
#define _PINYINPUTINTERFACE_H_

#include "BaseInputInterface.h"
#include "leanMD.h"
#include "../main/AtomPosInit.h"

class PinyInputInterface: public BaseInputInterface {
  public:
    PinyInputInterface ();
    ~PinyInputInterface ();

    void calc_max_brnch_root_bond_len (double& max_branch_root_bond_len);

    void setRootID (CkVec<SingleAtom>& atom_list);

    void getConstInfo (double** cut_atm_typ,
        int&     natm_typ,
        int**    iatm_atm_typ,
        int&     nAtoms,
        double** hmat,
        double** hmati,
        double&  box_vol);

    void getAtomPos (int** serial, 
        double** x, double** y, double** z, 
        double** vx, double** vy, double** vz);

    void getAtomInfo (double** mass, 
        double** charge, 
        double** epsilon,
        double** sigma);

    void getIntra (int&      nBonds,
        int&      nAngles,
        int&      nTorsions,
        int&      nonfo,
        Bond**    bonds,
        Angle**   angles,
        Torsion** torsions,
        OneFour** onfo);            

  private:
    int*     _serial;
    double*  _epsilon;
    double*  _sigma;
    Bond*    _bonds;
    Angle*   _angles;
    Torsion* _torsions;
    OneFour* _onfo;

    double* _x;
    double* _y;
    double* _z;
    double* _vx;
    double* _vy;
    double* _vz;

    double* _cut_atm_typ;
};
#endif
