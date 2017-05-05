#include "matrix_test.h"

MatrixTester::MatrixTester(CkArgMsg* msg) {
  thisProxy.test();
}

TestMatrix::TestMatrix(int mr, int mc, int tr, int tc, int t, CkCallback cb) : CBase_TestMatrix(mr, mc, tr, tc, cb) {
  setInitial(t);
}

void TestMatrix::setInitial(int type) {
  for (int r = 0; r < num_rows; r++) {
    for (int c = 0; c < num_cols; c++) {
      switch (type) {
        case ZERO:
          break;
        case RANDOM:
          data[r*num_cols + c].re = (double)rand();
          data[r*num_cols + c].im = (double)rand();
          break;
        case LOCAL:
          data[r*num_cols + c].re = r;
          data[r*num_cols + c].im = c;
          break;
        case GLOBAL:
          data[r*num_cols + c].re = start_row + r;
          data[r*num_cols + c].im = start_col + c;
          break;
        default:
          CkAbort("Bad initialization type\n");
          break;
      }
    }
  }
}

#include "matrix_test.def.h"
