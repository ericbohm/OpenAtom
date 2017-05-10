#include "standard_include.h"
#include "allclass_gwbse.h"
#include "eps_matrix.h"
#include "messages.h"
#include "pmatrix.h"
#include "controller.h"
#include "states.h"
#include "fft_routines.h"
#include "fft_controller.h"

#include <cstring> // for memcpy
using std::memcpy;

#define eps_rows 20
#define eps_cols 20
#define NSIZE 4
#define IDX_eps(r,c) ((r)*eps_cols + (c))

EpsMatrix::EpsMatrix() {
  GWBSE* gwbse = GWBSE::get();

  // Set some constants
  L = gwbse->gw_parallel.L;
  nfft = gwbse->gw_parallel.fft_nelems;
  qindex = Q_IDX; // Eventually the controller will set this

  data_received = 0;
  total_time = 0.0;
}

EpsMatrix::EpsMatrix(MatrixConfig config) : CBase_EpsMatrix(config) {
  GWBSE* gwbse = GWBSE::get();

  // Set some constants
  L = gwbse->gw_parallel.L;
  nfft = gwbse->gw_parallel.fft_nelems;
  qindex = Q_IDX; // Eventually the controller will set this

  total_time = 0.0;
  data_received = 0;
}

void EpsMatrix::setI(CLA_Matrix_interface mat, bool clean){
  matrix = mat;
  if (clean) {
    delete [] data;
    initialize();
  }
}

void EpsMatrix::receiveFs(Phase3Message* msg) {
  int n = 0;
  // TODO: memcpy
  for(int i=msg->start_i;i<=msg->end_i;i++)
    for(int j=msg->start_j;j<=msg->end_j;j++)
        data[IDX_eps(i,j)] = msg->data[n++];

  data_received+=n;
  if(data_received == eps_cols*eps_rows) {
    CkCallback cb(CkReductionTarget(Controller, epsilon_created), controller_proxy);
    contribute(cb);
  }
}
 
void EpsMatrix::multiply(double alpha, double beta) {
  matrix.multiply(alpha, beta, data, EpsMatrix::done_cb, (void*) this,
       thisIndex.x, thisIndex.y);
}

void EpsMatrix::add_compl_two() {
  int i = 0;
  complex compl_two(2.0, 0.0);
  if(thisIndex.x==thisIndex.y)
  for(int i=0;i<eps_rows;i++)
    data[IDX_eps(i,i)] += compl_two;

  CkCallback cb(CkReductionTarget(Controller, complement_multiplied), controller_proxy);
  contribute(cb);
}

void inline EpsMatrix::round_done(void) {
  CmiMemoryCheck();
#if 0
  if(thisIndex.x==0 && thisIndex.y==0){
    std::stringstream ss;
    for(int i = 0; i<10;i++)
      ss << data[i] << "\n";
    CkPrintf("\nmultiplied data = %s\n", ss.str().c_str());
  }
#endif
  CkCallback cb(CkReductionTarget(Controller, m_multiplied), controller_proxy);
  contribute(cb);
}

void EpsMatrix::scalar_multiply(double alpha) {
  for(int i=0;i<config.tile_rows;i++)
    for(int j=0;j<config.tile_cols;j++)
      data[IDX_eps(i,j)] = alpha*data[IDX_eps(i,j)]; 

  CkCallback cb(CkReductionTarget(Controller, scalar_multiplied), controller_proxy);
  contribute(cb);
}

void EpsMatrix::screenedExchange() {
  complex contribution(0.0,0.0);
  const int N = 2; // TODO: Give an actual value
  FVectorCache* f_cache = fvector_cache_proxy.ckLocalBranch();
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      for (int l = 0; l < L; l++) {
        complex* fi = f_cache->getFVec(i, l, start_row, config.tile_rows);
        complex* fj = f_cache->getFVec(j, l, start_col, config.tile_cols);
        for (int r = 0; r < config.tile_rows; r++) {
          for (int c = 0; c < config.tile_cols; c++) {
            contribution += fi[r]*data[IDX_eps(r,c)]*fj[c];
          }
        }
      }
    }
  }
//  CkCallback cb(CkReductionTarget(Controller, screenedExchangeComplete), controller_proxy);
//  contribute(sizeof(complex), &contribution, CkReduction::sum_double, cb);
}

void EpsMatrix::bareExchange() {
  complex contribution(0.0,0.0);
  FVectorCache* f_cache = fvector_cache_proxy.ckLocalBranch();
  int ndata = nfft[0]*nfft[1]*nfft[2];
  int block_size = ndata/(config.chareRows());
  if(ndata % config.chareRows() > 0) {
    block_size += 2;
  }

  if(thisIndex.x == thisIndex.y) {
    for (int i = 0; i < NSIZE; i++) {//ib = 5 to 8 actually, map to a number from 0
      for (int l = 0; l < L; l++) {
        complex* f = f_cache->getFVec(i, l, thisIndex.x, block_size);//start_row, num_rows);
        int end = block_size;
        if(thisIndex.x == 6) {
          end = ndata%block_size;
        }
        for(int ii=0; ii < end; ii++) {
          complex tmp = f[ii]*f[ii];
          contribution += sqrt(tmp.getMagSqr());
        }
      }
    }
  }
  CkCallback cb(CkReductionTarget(Controller, bareExchangeComplete), controller_proxy);
  contribute(sizeof(complex), &contribution, CkReduction::sum_double, cb);
}

void EpsMatrix::findAlpha() {
  if (config.chareCols() != 1) {
    CkAbort("findAlpha() only implemented for 1D decompositions\n");
  }
  double R = 0;
  for(int i = 0; i < config.tile_cols; i++) {
    R += abs(data[i]);
  }

  CkCallback cb(CkReductionTarget(Controller, found_alpha), controller_proxy);
  contribute(sizeof(long double), &R, CkReduction::max_double, cb);
}

void EpsMatrix::convergence_check(CProxy_EpsMatrix cproxy){
    // TODO: memcpy
    std::vector<complex> data_out(total_data);
    for(int i=0;i<total_data;i++)
      data_out[i] = data[i];

    cproxy(thisIndex.x, thisIndex.y).receiveConvCheck(data_out);
}

void EpsMatrix::receiveConvCheck(std::vector<complex> data_in) {
  double Rmax=0;  // the largest element
  double tmp;
  for(int i=0; i<total_data; i++) {
    tmp = abs(data[i] - data_in[i]);
    if( tmp > Rmax ){ Rmax = tmp; }
  }
  contribute(sizeof(complex), &Rmax, CkReduction::max_double,
      CkCallback(CkReductionTarget(Controller, converge_results), controller_proxy));
}

void EpsMatrix::createTranspose(CProxy_EpsMatrix other, bool todo) {
  std::vector<complex> incoming;
  for(int i=0; i < config.tile_rows; i++) {
    for(int j=0; j < config.tile_cols; j++) {
      if (todo) {
        incoming.push_back(-1*data[IDX_eps(i,j)]);
      } else {
        incoming.push_back(data[IDX_eps(i,j)]);
      }
    }
  }
  if(todo) {
    other(thisIndex.y, thisIndex.x).receiveTranspose(incoming);
  } else {
    other(thisIndex.x, thisIndex.y).receiveTranspose(incoming);
  }
}

void EpsMatrix::receiveTranspose(std::vector<complex> new_data) {
  unsigned n = 0;
  for(int i=0;i<config.tile_rows;i++)
    for(int j=0;j<config.tile_cols;j++)
      data[IDX_eps(i,j)] = new_data[n++];
    
  contribute(CkCallback(CkReductionTarget(Controller, transpose_complete), controller_proxy));  
}

#include "eps_matrix.def.h"
