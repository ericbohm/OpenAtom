#include "PinyInputInterface.h"

#define DEBUG //CmiPrintf

void PinyInputInterface::calc_max_brnch_root_bond_len (double& max_branch_root_bond_len) {
  DEBUG ("--> PinyInputInterface::calc_max_brnch_root_bond_len\n");
  // find Max (branch_root_bond_length)
  double max = 0;
  int i;

  if (SimParam.DoConstraints) {
    int     nBonds          = readonly_mdintra.mdbond.npow;
    double* eq_pow          = readonly_mdintra.mdbond.eq_pow;
    int*    j1_pow          = readonly_mdintra.mdbond.j1_pow;
    int*    j2_pow          = readonly_mdintra.mdbond.j2_pow;
    int*    jtyp_pow        = readonly_mdintra.mdbond.jtyp_pow;

    int     nAtoms          = readonly_mdatoms.mdclatoms_info.natm_tot;
    double* mass            = readonly_mdatoms.mdclatoms_info.mass;

    int*    hydrogenFlag    = new int [nAtoms+1];

    double  hydrogenCutMass = SimParam.HydrogenCutMass;


    // check if an atom is hydrogen atom or not
    for (i=1; i<=nAtoms; i++) {
      if (mass [i] <= hydrogenCutMass) {
        hydrogenFlag [i] = 1;
      } else {
        hydrogenFlag [i] = 0;
      }
    }

    for (i=1; i<=nBonds; i++) {
      if ((1 == hydrogenFlag [j1_pow [i]]) || 
          (1 == hydrogenFlag [j2_pow [i]])) {
        if (max < eq_pow [jtyp_pow [i]]) max = eq_pow [jtyp_pow [i]];
      }
    }

    delete [] hydrogenFlag;
  }

  /*
     for (i=1; i<=_rBond.ntyp_21; i++) {
     if (max < eq_21 [1][i]) max = eq [1][i];
     }

     for (i=1; i<=_rBond.ntyp_33; i++) {
     if (max < eq_33 [1][i]) max = eq_33 [1][i];
     if (max < eq_33 [2][i]) max = eq_33 [2][i];
     }

     for (i=1; i<=_rBond.ntyp_46; i++) {
     if (max < eq_46 [1][i]) max = eq_46 [1][i];
     if (max < eq_46 [2][i]) max = eq_46 [2][i];
     if (max < eq_46 [3][i]) max = eq_46 [3][i];
     }
   */  
  max_branch_root_bond_len = max;

  DEBUG ("<-- PinyInputInterface::calc_max_brnch_root_bond_len\n");
}

void PinyInputInterface::setRootID (CkVec<SingleAtom>& atom_list) {

  // assuming that only first atom in rigid bond data-structure is root
  // and all other atoms are branches following calc is done
  DEBUG ("--> PinyInputInterface::setRootID\n");
  int i;
  int  num_21 = readonly_mdintra.mdgrp_bond_con.num_21;
  int* j1_21  = readonly_mdintra.mdgrp_bond_con.j1_21;
  int* j2_21  = readonly_mdintra.mdgrp_bond_con.j2_21;
  int  num_33 = readonly_mdintra.mdgrp_bond_con.num_33;
  int* j1_33  = readonly_mdintra.mdgrp_bond_con.j1_33;
  int* j2_33  = readonly_mdintra.mdgrp_bond_con.j2_33;
  int* j3_33  = readonly_mdintra.mdgrp_bond_con.j3_33;
  int  num_46 = readonly_mdintra.mdgrp_bond_con.num_46;
  int* j1_46  = readonly_mdintra.mdgrp_bond_con.j1_46;
  int* j2_46  = readonly_mdintra.mdgrp_bond_con.j2_46;
  int* j3_46  = readonly_mdintra.mdgrp_bond_con.j3_46;
  int* j4_46  = readonly_mdintra.mdgrp_bond_con.j4_46;

  for (i=1; i<=num_21; i++) {
    // check if it was not assigned a root already
    CkAssert (0 == atom_list [j1_21 [i] - 1].rootID);
    CkAssert (0 == atom_list [j2_21 [i] - 1].rootID);

    CkAssert (j1_21 [i] == atom_list [j1_21 [i] - 1].id);
    CkAssert (j2_21 [i] == atom_list [j2_21 [i] - 1].id);

    // -ve value of rootID indicates that it is root otherwise branch
    atom_list [j1_21 [i] - 1].rootID = -i;
    atom_list [j2_21 [i] - 1].rootID = j1_21 [i];
  }

  for (i=1; i<=num_33; i++) {
    CkAssert (0 == atom_list [j1_33 [i] - 1].rootID);
    CkAssert (0 == atom_list [j2_33 [i] - 1].rootID);
    CkAssert (0 == atom_list [j3_33 [i] - 1].rootID);

    CkAssert (j1_33 [i] == atom_list [j1_33 [i] - 1].id);
    CkAssert (j2_33 [i] == atom_list [j2_33 [i] - 1].id);
    CkAssert (j3_33 [i] == atom_list [j3_33 [i] - 1].id);

    // -ve value of rootID indicates that it is root otherwise branch
    atom_list [j1_33 [i] - 1].rootID = -num_21 - i;
    atom_list [j2_33 [i] - 1].rootID = j1_33 [i];
    atom_list [j3_33 [i] - 1].rootID = j1_33 [i];
  }

  for (i=1; i<=num_46; i++) {
    CkAssert (0 == atom_list [j1_46 [i] - 1].rootID);
    CkAssert (0 == atom_list [j2_46 [i] - 1].rootID);
    CkAssert (0 == atom_list [j3_46 [i] - 1].rootID);
    CkAssert (0 == atom_list [j4_46 [i] - 1].rootID);

    CkAssert (j1_46 [i] == atom_list [j1_46 [i] - 1].id);
    CkAssert (j2_46 [i] == atom_list [j2_46 [i] - 1].id);
    CkAssert (j3_46 [i] == atom_list [j3_46 [i] - 1].id);
    CkAssert (j4_46 [i] == atom_list [j4_46 [i] - 1].id);

    // -ve value of rootID indicates that it is root otherwise branch
    atom_list [j1_46 [i] - 1].rootID = -num_21 - num_33 - i;
    atom_list [j2_46 [i] - 1].rootID = j1_46 [i];
    atom_list [j3_46 [i] - 1].rootID = j1_46 [i];
    atom_list [j4_46 [i] - 1].rootID = j1_46 [i];
  }
  DEBUG ("<-- PinyInputInterface::setRootID\n");
}

void PinyInputInterface::getConstInfo (double** cut_atm_typ,
    int&     natm_typ,
    int**    iatm_atm_typ,
    int&     nAtoms,
    double** hmat,
    double** hmati,
    double&  box_vol) {
  DEBUG ("--> PinyInputInterface::getConstInfo\n");
  // number of atoms
  nAtoms        = readonly_mdatoms.mdclatoms_info.natm_tot;
  natm_typ      = readonly_mdinter.mdinteract.natm_typ;
#ifdef USE_SPLINE_LJ // use spline lj routine
  _cut_atm_typ  = new double [natm_typ];
  for (int i=0; i<natm_typ; i++) {
    _cut_atm_typ [i] = readonly_mdinter.mdinteract.cutoff_max;
  }
  *cut_atm_typ  = _cut_atm_typ;
#else // use when explicit lj is used
  *cut_atm_typ  = &(readonly_mdinter.mdinteract.cut_atm_typ [1]);
#endif
  *iatm_atm_typ = &(readonly_mdatoms.mdatom_maps.iatm_atm_typ [1]);
  *hmat         = &(readonly_general_data.gencell.hmat [1]);
  *hmati        = &(readonly_general_data.gencell.hmati [1]);
  box_vol       = readonly_general_data.gencell.vol;
  DEBUG ("<-- PinyInputInterface::getConstInfo\n");
}

void PinyInputInterface::getAtomPos (int** serial, 
    double** x, double** y, double** z, 
    double** vx, double** vy, double** vz) {
  DEBUG ("--> PinyInputInterface::getAtomPos\n");
  int nAtoms  = readonly_mdatoms.mdclatoms_info.natm_tot;
  _serial = new int [nAtoms];

  // this object holds atom position and velocity data
  AtomPosInit atomPosInit;

  for (int i=0; i<nAtoms; i++) { 
    _serial [i] = i+1;
  }

  _x = new double [nAtoms];
  _y = new double [nAtoms];
  _z = new double [nAtoms];
  _vx = new double [nAtoms];
  _vy = new double [nAtoms];
  _vz = new double [nAtoms];

  // index '1' is used because piny arrays start from '1' and not '0'
  // accessing zeroth element of piny array is invalid and might cause 
  // segmentation violation
  *serial = _serial;
  CmiMemcpy (_x, &(atomPosInit.mdclatoms_pos [1].x [1]), sizeof (double)*nAtoms);
  CmiMemcpy (_y, &(atomPosInit.mdclatoms_pos [1].y [1]), sizeof (double)*nAtoms);
  CmiMemcpy (_z, &(atomPosInit.mdclatoms_pos [1].z [1]), sizeof (double)*nAtoms);
  CmiMemcpy (_vx, &(atomPosInit.mdclatoms_pos [1].vx [1]), sizeof (double)*nAtoms);
  CmiMemcpy (_vy, &(atomPosInit.mdclatoms_pos [1].vy [1]), sizeof (double)*nAtoms);
  CmiMemcpy (_vz, &(atomPosInit.mdclatoms_pos [1].vz [1]), sizeof (double)*nAtoms);

  *x = _x;
  *y = _y,
    *z = _z;

  *vx = _vx;
  *vy = _vy,
    *vz = _vz;
  DEBUG ("<-- PinyInputInterface::getAtomPos\n");
}

void PinyInputInterface:: getAtomInfo (double** mass, 
    double** charge, 
    double** epsilon,
    double** sigma) {
  DEBUG ("--> PinyInputInterface::getAtomInfo\n");
  int nAtoms  = readonly_mdatoms.mdclatoms_info.natm_tot;

  _epsilon       = new double [nAtoms];                                       
  _sigma         = new double [nAtoms];

  int*    iatm_atm_typ  = readonly_mdatoms.mdatom_maps.iatm_atm_typ;
  double* lj_sigma      = readonly_mdinter.mdinteract.lj_sigma;
  double* lj_epsilon    = readonly_mdinter.mdinteract.lj_epsilon;

#ifdef USE_SPLINE_LJ // use when need to use explicit lj
  memset (_epsilon, 0, sizeof (double)*nAtoms);
  memset (_sigma, 0, sizeof (double)*nAtoms);
#else
  if (NULL != lj_epsilon) {
    for (int i=0; i<nAtoms; i++) {
      _epsilon [i] = lj_epsilon [iatm_atm_typ [i+1]];
      _sigma [i]   = lj_sigma [iatm_atm_typ [i+1]];
    }
  }
#endif

  *mass    = &(readonly_mdatoms.mdclatoms_info.mass [1]);
  *charge  = &(readonly_mdatoms.mdclatoms_info.q [1]);
  *epsilon = _epsilon;
  *sigma   = _sigma;
  DEBUG ("<-- PinyInputInterface::getAtomInfo\n");
}

void PinyInputInterface::getIntra (int&      nBonds,
    int&      nAngles,
    int&      nTorsions,
    int&      nonfo,
    Bond**    bonds,
    Angle**   angles,
    Torsion** torsions,
    OneFour** onfo) {
  DEBUG ("--> PinyInputInterface::getIntra\n");
  nBonds    = readonly_mdintra.mdbond.npow;
  nAngles   = readonly_mdintra.mdbend.npow;
  nTorsions = readonly_mdintra.mdtors.npow;
  nonfo     = readonly_mdintra.mdonfo.num;

  _bonds     = NULL;
  _angles    = NULL;
  _torsions  = NULL;
  _onfo      = NULL;

  if (0 != nBonds) {
    _bonds    = new Bond [nBonds];
  }

  if (0 != nAngles) {
    _angles   = new Angle [nAngles];
  }

  if (0 != nTorsions) {
    _torsions = new Torsion [nTorsions];
  }

  if (0 != nonfo) {
    _onfo     = new OneFour [nonfo];
  }

  int i;

  for (i=0; i<nBonds; i++) {
    _bonds [i].atom1    = readonly_mdintra.mdbond.j1_pow [i+1];
    _bonds [i].atom2    = readonly_mdintra.mdbond.j2_pow [i+1];
    _bonds [i].bondType = readonly_mdintra.mdbond.jtyp_pow [i+1];
  }

  for (i=0; i<nAngles; i++) {
    _angles [i].atom1    = readonly_mdintra.mdbend.j1_pow [i+1];
    _angles [i].atom2    = readonly_mdintra.mdbend.j2_pow [i+1];
    _angles [i].atom3    = readonly_mdintra.mdbend.j3_pow [i+1];
    _angles [i].bendType = readonly_mdintra.mdbend.jtyp_pow [i+1];
  }

  for (i=0; i<nTorsions; i++) {
    _torsions [i].atom1   = readonly_mdintra.mdtors.j1_pow [i+1];
    _torsions [i].atom2   = readonly_mdintra.mdtors.j2_pow [i+1];
    _torsions [i].atom3   = readonly_mdintra.mdtors.j3_pow [i+1];
    _torsions [i].atom4   = readonly_mdintra.mdtors.j4_pow [i+1];
    _torsions [i].torType = readonly_mdintra.mdtors.jtyp_pow [i+1];
  }

  for (i=0; i<nonfo; i++) {
    _onfo [i].atom1    = readonly_mdintra.mdonfo.j1 [i+1];
    _onfo [i].atom2    = readonly_mdintra.mdonfo.j2 [i+1];
    _onfo [i].onfoType = readonly_mdintra.mdonfo.jtyp [i+1];
  }

  *bonds    = _bonds;
  *angles   = _angles;
  *torsions = _torsions;
  *onfo     = _onfo;
  DEBUG ("<-- PinyInputInterface::getIntra\n");
}

PinyInputInterface::~PinyInputInterface () {
  DEBUG ("--> PinyInputInterface::~PinyInputInterface \n");
  if (NULL != _serial) {delete [] _serial;}
  if (NULL != _epsilon) {delete [] _epsilon;}
  if (NULL != _sigma) {delete [] _sigma;}
  if (NULL != _bonds) {delete [] _bonds;}
  if (NULL != _angles) {delete [] _angles;}
  if (NULL != _torsions) {delete [] _torsions;}
  if (NULL != _onfo) {delete [] _onfo;}

  if (NULL != _x) {delete [] _x;}
  if (NULL != _y) {delete [] _y;}
  if (NULL != _z) {delete [] _z;}
  if (NULL != _vx) {delete [] _vx;}
  if (NULL != _vy) {delete [] _vy;}
  if (NULL != _vz) {delete [] _vz;}

  if (NULL != _cut_atm_typ) {delete [] _cut_atm_typ;}
  DEBUG ("<-- PinyInputInterface::~PinyInputInterface \n");
}

PinyInputInterface::PinyInputInterface () : _serial (NULL),
  _epsilon (NULL),
  _sigma (NULL),
  _bonds (NULL),
  _angles (NULL), 
  _torsions (NULL),
  _onfo (NULL),
  _x (NULL),
  _y (NULL),
  _z (NULL),
  _vx (NULL),
  _vy (NULL),
  _vz (NULL),
  _cut_atm_typ(NULL) {
    DEBUG ("--> PinyInputInterface::PinyInputInterface\n");
    DEBUG ("<-- PinyInputInterface::PinyInputInterface\n");
  }
