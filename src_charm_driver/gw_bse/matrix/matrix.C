#include "matrix.h"

#include <cstring> // for memcpy
#include <fstream>
#include <string>

using std::ifstream;
using std::ios;
using std::min;
using std::memcpy;
using std::ofstream;
using std::string;

void matrixCopy(CProxy_Matrix src, CProxy_Matrix dest, CkCallback cb) {
  dest.startCopy(src, cb);
}

void matrixCompare(CProxy_Matrix src, CProxy_Matrix dest, CkCallback cb) {
  dest.startCompare(src, cb);
}

Matrix::Matrix(int mr, int mc, int tr, int tc, CkCallback cb) {
  initialize(mr, mc, tr, tc, cb);
}

void Matrix::initialize(int mr, int mc, int tr, int tc, CkCallback cb) {
  matrix_rows = mr;
  matrix_cols = mc;
  num_rows = tr;
  num_cols = tc;
  start_row = num_rows * thisIndex.x;
  start_col = num_cols * thisIndex.y;

  data_received = 0;
  total_data = num_rows * num_cols;

  data = new complex[num_rows*num_cols];

  contribute(cb);
}

void Matrix::startCopy(CProxy_Matrix src, CkCallback cb) {
  if (thisIndex.x == 0 && thisIndex.y == 0) {
    src.copy(thisProxy, num_rows, num_cols);
  }
  busy = true;
  dataHandler = &Matrix::unpackMsg;
  thisProxy(thisIndex).waitForData(cb);
}

void Matrix::copy(CProxy_Matrix dest, int dest_rows, int dest_cols) {
  sendTo(dest, dest_rows, dest_cols);
}

void Matrix::startCompare(CProxy_Matrix src, CkCallback cb) {
  if (thisIndex.x == 0 && thisIndex.y == 0) {
    src.compare(thisProxy, num_rows, num_cols);
  }
  busy = true;
  dataHandler = &Matrix::compareMsg;
  thisProxy(thisIndex).waitForData(cb);
}

void Matrix::compare(CProxy_Matrix dest, int dest_rows, int dest_cols) {
  sendTo(dest, dest_rows, dest_cols);
}

void Matrix::read(string prefix, CkCallback cb) {
  if (num_cols != matrix_cols) {
    CkAbort("Read/Write only supported for row decomposition\n");
  }
  ifstream infile;
  string filename;
  for (int r = 0; r < num_rows; r++) {
    filename = prefix + std::to_string(start_row+r);
    infile.open(filename, ios::in);
    for (int c = 0; c < num_cols; c++) {
      infile >> data[r * num_cols + c].re;
      infile >> data[r * num_cols + c].im;
    }
    infile.close();
  }
  contribute(cb);
}

void Matrix::write(string prefix, CkCallback cb) {
  if (num_cols != matrix_cols) {
    CkAbort("Read/Write only supported for row decomposition\n");
  }
  ofstream outfile;
  string filename;
  for (int r = 0; r < num_rows; r++) {
    filename = prefix + std::to_string(start_row+r);
    outfile.open(filename, ios::out);
    for (int c = 0; c < num_cols; c++) {
      outfile << data[r * num_cols + c].re << " ";
      outfile << data[r * num_cols + c].im << " ";
    }
    outfile.close();
  }
  contribute(cb);
}

void Matrix::packMsg(DataMessage* msg) {
  int msg_idx = 0;
  int local_idx = (msg->row - start_row) * num_cols + (msg->col - start_col);
  for (int r = 0; r < msg->num_rows; r++) {
    memcpy(&msg->data[msg_idx],&data[local_idx],sizeof(complex)*msg->num_cols);
    msg_idx += msg->num_cols;
    local_idx += num_cols;
  }
}

void Matrix::unpackMsg(DataMessage* msg) {
  int msg_idx = 0;
  int local_idx = (msg->row - start_row) * num_cols + (msg->col - start_col);
  for (int r = 0; r < msg->num_rows; r++) {
    memcpy(&data[local_idx],&msg->data[msg_idx],sizeof(complex)*msg->num_cols);
    msg_idx += msg->num_cols;
    local_idx += num_cols;
  }
}

void Matrix::compareMsg(DataMessage* msg) {
  int msg_idx = 0;
  int local_idx = (msg->row - start_row) * num_cols + (msg->col - start_col);
  for (int r = 0; r < msg->num_rows; r++) {
    for (int c = 0; c < msg->num_cols; c++) {
      if (data[local_idx+c].re != msg->data[msg_idx+c].re ||
          data[local_idx+c].im != msg->data[msg_idx+c].im) {
        CkAbort("Matrices don't match\n");
      }
    }
    msg_idx += msg->num_cols;
    local_idx += num_cols;
  }
}

void Matrix::sendTo(CProxy_Matrix dest, int dest_rows, int dest_cols) {
  int row_idx = start_row;
  int col_idx = start_col;
  int max_row = row_idx + num_rows;
  int max_col = col_idx + num_cols;
  while (row_idx < max_row) {
    int rows = min(dest_rows - (row_idx % dest_rows), max_row - row_idx);
    while (col_idx < max_col) {
      int cols = min(dest_cols - (col_idx % dest_cols), max_col - col_idx);
      DataMessage* msg = new (rows*cols) DataMessage();
      msg->row = row_idx;
      msg->col = col_idx;
      msg->num_rows = rows;
      msg->num_cols = cols;
      packMsg(msg);

      int chare_row = row_idx / dest_rows;
      int chare_col = col_idx / dest_cols;
      dest(chare_row,chare_col).receiveData(msg);
      col_idx += cols;
    }
    row_idx += rows;
  }
}

#include "matrix.def.h"
