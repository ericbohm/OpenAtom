/*****************************************************************************
 * $Source$
 * $Author$
 * $Date$
 * $Revision$
 *****************************************************************************/

/** \file energy.h
 *
 */

#ifndef _ENERGY_H_
#define _ENERGY_H_

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/* Structure that stores all energies */

struct EnergyStruct {

  int iteration_gsp;       // step upon which energies below computed
  double enl;             // non local - k and spin
  double eext;            // local external energy
  double eke;             // quantum kinetic energy - k and spin
  double ehart;           // hartree energy
  double egga;            // exchange correlation grad corr - check for spin reduce
  double eexc;            // exchange correlation local - check for spin reduce
  double fictEke;         // fict KE from cp dynamics - k and spin
  double fmagPsi;         // coef force magnitude
  double eewald_recip;    // atm(ion)-atm(ion) recip (computed by psi chares)
  double totalElecEnergy; // sum of electronic energies + ewald_recip 
                          // local total? or summed over k and spin?
  // no fict and no ewald_real

  int iteration_atm;       // step upon which energies below computed.
  double eewald_real;     // Real space ewald.
  double eKinetic_atm;    // classical kinetic energy
  double eKineticNhc_atm; // NHC kinetic energy
  double potNhc_atm;      // NHC pot energy
  double fmag_atm;        // magnitude of atm forces
  double potPIMDChain;
  double totalpotPIMDChain;
};
PUPbytes(EnergyStruct);
//============================================================================
#endif
