#include "pmatrix.h"
#include "gw_bse.h"

PMatrix::PMatrix() {
  num_rows = config.rows_per_chare;
  row_size = config.n_elems;
  start_row = thisIndex * num_rows;

  rows = new double*[num_rows];
  for (int i = 0; i < num_rows; i++) {
    rows[i] = new double[row_size];
  }
}

void PMatrix::receiveRowContribution(RowMessage* msg) {
  unsigned local_index = msg->index - start_row;
  for (int i = 0; i < row_size; i++) {
    rows[local_index][i] += msg->row[i];
  }
}
