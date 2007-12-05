/*****************************************************************************
 * $Source$
 * $Author$
 * $Date$
 * $Revision$
 *****************************************************************************/

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
/** \file moduleinit.C
 *
 */
//==============================================================================


extern void _registercommlib(void);
extern void _registerCkMulticast(void);
extern void _registerDummyLB(void);
extern void _registerRefineLB(void);
extern void _registerRefineCommLB(void);
extern void _registerckPairCalculator(void);
void _registerExternalModules(char **argv) {
  _registercommlib();
  _registerCkMulticast();
  _registerDummyLB();
  _registerRefineLB();
  _registerRefineCommLB();
  _registerckPairCalculator();
}
void _createTraces(char **argv) {
}
