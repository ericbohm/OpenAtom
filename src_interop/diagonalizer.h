#ifndef _diagonalizer_h_
#define _diagonalizer_h_
#include "paircalc/pcFwdDeclarations.h"

template <typename T>
struct diagData_t {
  int selfsize;
  int totalsize;

  T* plambda;
  T* rlambda;
  T* ilambda;
  T* selflambda;
  int pelements;
};

#endif
