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

MatrixConfig convertTo1D(const MatrixConfig& in, int tile_rows) {
  MatrixConfig out = in;
  out.tile_rows = tile_rows;
  out.tile_cols = out.mat_cols;
  return out;
}

void matrixCopy(CProxy_Matrix src, CProxy_Matrix dest, CkCallback cb) {
  dest.copy(src, cb);
}

void matrixCompare(CProxy_Matrix src, CProxy_Matrix dest, CkCallback cb) {
  dest.compare(src, cb);
}

void matrixRead(CProxy_Matrix dest, std::string prefix, CkCallback cb) {
  dest.read(prefix, cb);
}

void matrixWrite(CProxy_Matrix src, std::string prefix, CkCallback cb) {
  src.write(prefix, cb);
}

void matrixVerify(CProxy_Matrix src, string prefix, CkCallback cb) {
  src.verify(prefix, cb);
}

Matrix::Matrix(MatrixConfig config) : config(config) {
  initialize();
}

Matrix::Matrix(int mr, int mc, int tr, int tc) {
  config.mat_rows = mr; config.mat_cols = mc;
  config.tile_rows = tr; config.tile_cols = tc;
  initialize();
}

void Matrix::initialize() {
  start_row = config.tile_rows * thisIndex.x;
  start_col = config.tile_cols * thisIndex.y;

  data_received = 0;
  total_data = config.tile_rows * config.tile_cols;

  data = new complex[total_data];
}

void Matrix::copy(CProxy_Matrix src, CkCallback cb) {
  if (thisIndex.x == 0 && thisIndex.y == 0) {
    src.sendData(thisProxy, config.tile_rows, config.tile_cols);
  }
  busy = true;
  dataHandler = &Matrix::unpackMsg;
  thisProxy(thisIndex).waitForData(cb);
}

void Matrix::compare(CProxy_Matrix src, CkCallback cb) {
  if (thisIndex.x == 0 && thisIndex.y == 0) {
    src.sendData(thisProxy, config.tile_rows, config.tile_cols);
  }
  busy = true;
  dataHandler = &Matrix::compareMsg;
  thisProxy(thisIndex).waitForData(cb);
}

void Matrix::read(string prefix, CkCallback cb) {
  if (config.chareCols() != 1) {
    CkAbort("Read/Write only supported for row decomposition\n");
  }
  ifstream infile;
  string filename;
  for (int r = 0; r < config.tile_rows; r++) {
    filename = prefix + std::to_string(start_row+r);
    infile.open(filename, ios::in);
    for (int c = 0; c < config.tile_cols; c++) {
      infile >> data[r * config.tile_cols + c].re;
      infile >> data[r * config.tile_cols + c].im;
    }
    infile.close();
  }
  contribute(cb);
}

void Matrix::write(string prefix, CkCallback cb) {
  if (config.chareCols() != 1) {
    CkAbort("Read/Write only supported for row decomposition\n");
  }
  ofstream outfile;
  string filename;
  for (int r = 0; r < config.tile_rows; r++) {
    filename = prefix + std::to_string(start_row+r);
    outfile.open(filename, ios::out);
    for (int c = 0; c < config.tile_cols; c++) {
      outfile << data[r * config.tile_cols + c].re << " ";
      outfile << data[r * config.tile_cols + c].im << " ";
    }
    outfile.close();
  }
  contribute(cb);
}

void Matrix::verify(string prefix, CkCallback cb) {
  if (config.chareCols() != 1) {
    CkAbort("Read/Write only supported for row decomposition\n");
  }
  ifstream infile;
  string filename;
  complex tmp;
  for (int r = 0; r < config.tile_rows; r++) {
    filename = prefix + std::to_string(start_row+r);
    infile.open(filename, ios::in);
    for (int c = 0; c < config.tile_cols; c++) {
      infile >> tmp.re;
      infile >> tmp.im;
      if (data[r * config.tile_cols + c].re != tmp.re ||
          data[r * config.tile_cols + c].im != tmp.im) {
        CkAbort("Matrix verification failed!\n");
      }
    }
    infile.close();
  }
  contribute(cb);
}

void Matrix::sendData(CProxy_Matrix dest, int dest_rows, int dest_cols) {
  int max_row = start_row + config.tile_rows;
  int max_col = start_col + config.tile_cols;
  int row_idx = start_row;
  while (row_idx < max_row) {
    int rows = min(dest_rows - (row_idx % dest_rows), max_row - row_idx);
    int col_idx = start_col;
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

void Matrix::packMsg(DataMessage* msg) {
  int msg_idx = 0;
  int local_idx = (msg->row - start_row) * config.tile_cols + (msg->col - start_col);
  for (int r = 0; r < msg->num_rows; r++) {
    memcpy(&msg->data[msg_idx],&data[local_idx],sizeof(complex)*msg->num_cols);
    msg_idx += msg->num_cols;
    local_idx += config.tile_cols;
  }
}

void Matrix::unpackMsg(DataMessage* msg) {
  int msg_idx = 0;
  int local_idx = (msg->row - start_row) * config.tile_cols + (msg->col - start_col);
  for (int r = 0; r < msg->num_rows; r++) {
    memcpy(&data[local_idx],&msg->data[msg_idx],sizeof(complex)*msg->num_cols);
    msg_idx += msg->num_cols;
    local_idx += config.tile_cols;
  }
}

void Matrix::compareMsg(DataMessage* msg) {
  int msg_idx = 0;
  int local_idx = (msg->row - start_row) * config.tile_cols + (msg->col - start_col);
  for (int r = 0; r < msg->num_rows; r++) {
    for (int c = 0; c < msg->num_cols; c++) {
      if (data[local_idx+c].re != msg->data[msg_idx+c].re ||
          data[local_idx+c].im != msg->data[msg_idx+c].im) {
        CkAbort("Matrices don't match\n");
      }
    }
    msg_idx += msg->num_cols;
    local_idx += config.tile_cols;
  }
}

#include "matrix.def.h"
