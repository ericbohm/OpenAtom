#include "matrix_test.h"

MatrixTester::MatrixTester(CkArgMsg* msg) {
  thisProxy.test();
}

TestMatrix::TestMatrix(int mr, int mc, int tr, int tc, CkCallback cb) {
  initialize(mr, mc, tr, tc, cb);
  setInitial();
}

void TestMatrix::setInitial(int type) {
  for (int r = 0; r < num_rows; r++) {
    for (int c = 0; c < num_cols; c++) {
      data[r*num_cols + c].re = (double)rand();
      data[r*num_cols + c].im = (double)rand();
    }
  }
}

#include "matrix_test.def.h"
