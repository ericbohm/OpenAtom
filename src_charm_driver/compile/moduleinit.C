extern void _registercommlib(void);
extern void _registerfftlib(void);
extern void _registerCkMulticast(void);
extern void _registerCkSparseContiguousReducer(void);
extern void _registerDummyLB(void);
extern void _registerRefineLB(void);
extern void _registerRefineCommLB(void);
extern void _registerckPairCalculator(void);
void _registerExternalModules(char **argv) {
  _registercommlib();
  _registerfftlib();
  _registerCkMulticast();
  _registerCkSparseContiguousReducer();
  _registerDummyLB();
  _registerRefineLB();
  _registerRefineCommLB();
  _registerckPairCalculator();
}
void _createTraces(char **argv) {
}
