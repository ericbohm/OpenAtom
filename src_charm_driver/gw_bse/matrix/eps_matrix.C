#include "standard_include.h"
#include "allclass_gwbse.h"
#include "eps_matrix.h"
#include "messages.h"
#include "pmatrix.h"
#include "controller.h"
#include "states.h"
#include "fft_routines.h"
#include "fft_controller.h"

#define eps_rows 20
#define eps_cols 20
#define NSIZE 4
#define IDX_eps(r,c) ((r)*eps_cols + (c))

EpsMatrix2D::EpsMatrix2D(){
  GWBSE* gwbse = GWBSE::get();

  // Set some constants
  L = gwbse->gw_parallel.L;
  nfft = gwbse->gw_parallel.fft_nelems;
  qindex = Q_IDX; // Eventually the controller will set this

  // Grab a local pointer to the fft controller for fft-ing our rows
  // TODO: Is this guaranteed to be safe (is the local branch created for sure)?

  // Figure out what part of the matrix we have and allocate space
  num_rows = eps_rows;
  num_cols = eps_cols;
  num_chares = (matrix_dimension / num_rows) * (matrix_dimension / num_cols);
  start_row = thisIndex.x * num_rows;
  start_col = thisIndex.y * num_cols;

  data = new complex[num_rows * num_cols];
  data_received = 0;
  total_time = 0.0;
}

void EpsMatrix2D::setSize(int matrix_dimension_in){
//  CkPrintf("\nSetting size as %d\n", matrix_dimension_in);
  matrix_dimension = matrix_dimension_in;
  CkCallback *cb = new CkCallback(CkReductionTarget(Controller, set_size), controller_proxy);
  contribute(*cb);
}

EpsMatrix2D::EpsMatrix2D(CLA_Matrix_interface mat){
  GWBSE* gwbse = GWBSE::get();

  // Set some constants
  L = gwbse->gw_parallel.L;
  nfft = gwbse->gw_parallel.fft_nelems;
  qindex = Q_IDX; // Eventually the controller will set this

  num_rows = eps_rows;
  num_cols = eps_cols;
  num_chares = (matrix_dimension / num_rows) * (matrix_dimension / num_cols);
  start_row = thisIndex.x * num_rows;
  start_col = thisIndex.y * num_cols;
   
  data = new complex[num_rows * num_cols];

  for(int i=0;i<num_rows * num_cols;i++){
    data[i].re = 0;
    data[i].im = 0;
  }
  total_time = 0.0;
  data_received = 0;
  
  this->matrix = mat;
}

void EpsMatrix2D::setI(CLA_Matrix_interface mat, bool clean){
  this->matrix = mat;
  if(clean)
    data = new complex[num_rows * num_cols];
}

void EpsMatrix2D::receiveFs(Phase3Message* msg){

  int n = 0;
  for(int i=msg->start_i;i<=msg->end_i;i++)
    for(int j=msg->start_j;j<=msg->end_j;j++)
        data[IDX_eps(i,j)] = msg->data[n++];

  data_received+=n;

    int i = 0;
    if(data_received == eps_cols*eps_rows)
    {
      CkCallback *cb = new CkCallback(CkReductionTarget(Controller, epsilon_created), controller_proxy);
      contribute(sizeof(int), &i, CkReduction::sum_int, *cb);
    }
}
 
void EpsMatrix2D::multiply(double alpha, double beta){
  matrix.multiply(alpha, beta, data, EpsMatrix2D::done_cb, (void*) this,
       thisIndex.x, thisIndex.y);

}

void EpsMatrix2D::add_compl_two(){
  int i = 0;
  complex compl_two(2.0, 0.0);
  if(thisIndex.x==thisIndex.y)
  for(int i=0;i<eps_rows;i++)
    data[IDX_eps(i,i)] += compl_two;

  CkCallback *cb = new CkCallback(CkReductionTarget(Controller, complement_multiplied), controller_proxy);
    contribute(sizeof(int), &i, CkReduction::sum_int,
      *cb
    );
}
 void inline EpsMatrix2D::round_done(void){
    int i=0;
    CmiMemoryCheck();
#if 0
    if(thisIndex.x==0 && thisIndex.y==0){
      std::stringstream ss;
      for(int i = 0; i<10;i++)
        ss << data[i] << "\n";
      CkPrintf("\nmultiplied data = %s\n", ss.str().c_str());
    }
#endif
    CkCallback *cb = new CkCallback(CkReductionTarget(Controller, m_multiplied), controller_proxy);
      contribute(sizeof(int), &i, CkReduction::sum_int,
        *cb
      );
    }


long double EpsMatrix2D::max_fn(int size){

    long double max = -1;

    for(int i =0;i<size;i++){
      double temp = abs(total[i]);
      if(temp < 0)
        temp *= -1;
      if(temp>max)
        max = temp;
    }

    return max;
}

void EpsMatrix2D::scalar_multiply(double alpha){
  for(int i=0;i<num_rows;i++)
    for(int j=0;j<num_cols;j++)
      data[IDX_eps(i,j)] = alpha*data[IDX_eps(i,j)]; 

  CkCallback *cb = new CkCallback(CkReductionTarget(Controller, scalar_multiplied), controller_proxy);
  contribute(*cb);
}

void EpsMatrix2D::sendTo1D(){
  
  for(int i=0;i<num_rows;i++){
    Phase3Message *msg; 
    msg = new(eps_cols)Phase3Message();
    int n = 0;
    for(int j=0;j<num_cols;j++){
      msg->data[n++] =  data[IDX_eps(i,j)];
    }
    msg->start_i = thisIndex.y*eps_cols;
    msg->end_i = (thisIndex.y+1)*eps_cols;
    eps_proxy1D(thisIndex.x*eps_rows+i).receiveData(msg);
  }
  
}

void EpsMatrix2D::screenedExchange() {
  complex contribution(0.0,0.0);
  const int N = 2; // TODO: Give an actual value
  FVectorCache* f_cache = fvector_cache_proxy.ckLocalBranch();
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      for (int l = 0; l < L; l++) {
        complex* fi = f_cache->getFVec(i, l, start_row, num_rows);
        complex* fj = f_cache->getFVec(j, l, start_col, num_cols);
        for (int r = 0; r < num_rows; r++) {
          for (int c = 0; c < num_cols; c++) {
            contribution += fi[r]*data[IDX_eps(r,c)]*fj[c];
          }
        }
      }
    }
  }
//  CkCallback cb(CkReductionTarget(Controller, screenedExchangeComplete), controller_proxy);
//  contribute(sizeof(complex), &contribution, CkReduction::sum_double, cb);
}

void EpsMatrix2D::bareExchange() {
  complex* tile = new complex[num_rows * num_cols];
  complex contribution(0.0,0.0);

  FVectorCache* f_cache = fvector_cache_proxy.ckLocalBranch();
  int ndata = nfft[0]*nfft[1]*nfft[2];
  int block_size = ndata/(matrix_dimension/eps_rows);
  if(ndata%(matrix_dimension/eps_rows) > 0)
    block_size += 2;
  if(thisIndex.x== thisIndex.y)
  for (int i = 0; i < NSIZE; i++) {//ib = 5 to 8 actually, map to a number from 0
      for (int l = 0; l < L; l++) {
      complex* f = f_cache->getFVec(i, l, thisIndex.x, block_size);//start_row, num_rows);
      int end = block_size;
      if(thisIndex.x==6)
        end = ndata%(matrix_dimension/eps_rows);
      for(int ii=0;ii<end;ii++){
        complex tmp = f[ii]*f[ii];
        contribution += sqrt(tmp.getMagSqr());
      }
    }
  }
  CkCallback cb(CkReductionTarget(Controller, bareExchangeComplete), controller_proxy);
  contribute(sizeof(complex), &contribution, CkReduction::sum_double, cb);
}

EpsMatrix1D::EpsMatrix1D(){
  received = 0;
}

void EpsMatrix1D::setSize(int ncols){
  n_1d_cols = ncols;
  data = new complex[n_1d_cols];
}
void EpsMatrix1D::receiveData(Phase3Message *msg){
  int n = 0;
  for(int i=msg->start_i;i<msg->end_i;i++)
    data[i] = msg->data[n++];

  if(++received == n_1d_cols/eps_cols){
    contribute(CkCallback(
        CkReductionTarget(Controller, created1d_complete), controller_proxy)); 
  }
}

void EpsMatrix1D::findAlpha(){
  double R = 0;
  for(int i=0; i<n_1d_cols; i++)
    R += abs(data[i]);

  contribute(sizeof(long double), &R, CkReduction::max_double,
        CkCallback(CkReductionTarget(Controller, found_alpha), controller_proxy)
      );
}

void EpsMatrix2D::findAlpha(){

    for(int i=0;i<num_rows;i++)
        for(int j=0;j<num_cols;j++)
        total[i] += data[IDX_eps(i,j)];

    long double max = 0;
    int chunk_count = matrix_dimension/num_rows;
    if(thisIndex.y==chunk_count-1){
        max = max_fn(num_rows);   //return the result;
    }else
        thisProxy(thisIndex.x, thisIndex.y+1).findAlpha();
 
    contribute(sizeof(long double), &max, CkReduction::max_double,
        CkCallback(CkReductionTarget(Controller, found_alpha), controller_proxy)
      );
}


void EpsMatrix2D::convergence_check(CProxy_EpsMatrix2D cproxy){
    
    std::vector<complex> data_out(num_rows*num_cols);
    for(int i=0;i<num_rows*num_cols;i++)
      data_out[i] = data[i];

    cproxy(thisIndex.x, thisIndex.y).receiveConvCheck(data_out);
}

void EpsMatrix2D::receiveConvCheck(std::vector<complex> data_in){

   double Rmax=0;  // the largest element
   double tmp;
   for(int i=0; i<num_rows*num_cols; i++){
        tmp = abs(data[i] - data_in[i]) ;
        if( tmp > Rmax ){ Rmax = tmp; }
    }
   contribute(sizeof(complex), &Rmax, CkReduction::max_double,
        CkCallback(CkReductionTarget(Controller, converge_results), controller_proxy)
      );
}

void EpsMatrix2D::createTranspose(bool todo){
    std::vector<complex> incoming;
    unsigned n = 0;
    complex* new_data;
    new_data = new complex[num_rows*num_cols];
    for(int i=0;i<num_rows;i++){
        for(int j=0;j<num_cols;j++){
            if(todo)
              new_data[n] = -1*data[IDX_eps(i,j)];//negative values //data[IDX_eps(j,i)];
            else  
              new_data[n] =  data[IDX_eps(i,j)];//data[IDX_eps(j,i)];
            incoming.push_back(new_data[n++]);
        }
    }
    if(todo)
      eps_matrix2D_bproxy(thisIndex.y, thisIndex.x).receiveTranspose(incoming);
    else
      eps_matrix2D_bproxy(thisIndex.x, thisIndex.y).receiveTranspose(incoming);//new_data);
}

void EpsMatrix2D::receiveTranspose(std::vector<complex> new_data){

  unsigned n = 0;
  for(int i=0;i<num_rows;i++)
    for(int j=0;j<num_cols;j++)
      data[IDX_eps(i,j)] = new_data[n++];
    
  contribute(CkCallback(CkReductionTarget(Controller, transpose_complete), controller_proxy));  
}


void EpsMatrix2D::sendTo(CProxy_EpsMatrix2D receiver_proxy){
    std::vector<complex> incoming;
    for(int i=0;i<num_rows;i++){
        for(int j=0;j<num_cols;j++){
            incoming.push_back(data[IDX_eps(i,j)]);
        }
    }
    receiver_proxy(thisIndex.x, thisIndex.y).receiveData(incoming);
}

void EpsMatrix2D::receiveData(std::vector<complex> new_data){

  unsigned n = 0;
  for(int i=0;i<num_rows;i++)
    for(int j=0;j<num_cols;j++)
      data[IDX_eps(i,j)] = new_data[n++];
   
  contribute(CkCallback(CkReductionTarget(Controller, finished_copy), controller_proxy));
}



// TODO: These methods shouldn't be part of PMatrix, and should also just be
// computed once at startup.
void EpsMatrix2D::kqIndex(unsigned ikpt, unsigned& ikq, int* uklapp){
  GWBSE* gwbse = GWBSE::get();

  // temporary space to save k/q/k+q vectors
  double *this_k, *this_q;
  double k_plus_q[3], k_plus_q_orig[3];

  this_k = gwbse->gwbseopts.kvec[ikpt];
  this_q = gwbse->gwbseopts.qvec[qindex];

  for (int i=0; i<3; i++) {
    // calculate k+q vector
    k_plus_q[i] = this_k[i] + this_q[i]; // k+q vector
    k_plus_q_orig[i] = k_plus_q[i]; // save it for Umklapp process
    // if not 0 =< k+q [i] <1, adjust k+q so that k+q[i] is in the Brillouine zone
    if ( k_plus_q[i] >= 1 ) {
      k_plus_q[i] -= 1;
    }
    else if( k_plus_q[i] < 0 ){
      k_plus_q[i] += 1;
    }
  }

  // find k+q vector index
  for (int kk=0; kk < gwbse->gwbseopts.nkpt; kk++) {
    bool match = true;
    this_k = gwbse->gwbseopts.kvec[kk];
    //this_k is now a difference between k and k+q
    for (int i=0; i<3; i++) {
      if (this_k[i] != k_plus_q[i]) {
        match = false;
        break;
      }
    }
    if (match) {
      ikq = kk;
      break;
    }
  }
  // save umklapp scattering information
  for (int i=0; i<3; i++) {
    uklapp[i] = int( k_plus_q_orig[i] - k_plus_q[i] );
  }

}

void EpsMatrix2D::getUmklappFactor(complex* umklapp_factor, int uklpp[3]){

  if (uklpp[0]==0 && uklpp[1]==0 && uklpp[2]==0){
    // do nothing
  }
  else{
    GWBSE *gwbse = GWBSE::get();
    int* nfft;
    nfft = gwbse->gw_parallel.fft_nelems;
    double *a1, *a2, *a3, *b1, *b2, *b3;
    a1 = gwbse->gwbseopts.a1;
    a2 = gwbse->gwbseopts.a2;
    a3 = gwbse->gwbseopts.a3;
    b1 = gwbse->gwbseopts.b1;
    b2 = gwbse->gwbseopts.b2;
    b3 = gwbse->gwbseopts.b3;
    double lattconst = gwbse->gwbseopts.latt;

    double rijk, G0, phase;
    unsigned counter = 0;
    for(int i=0; i<nfft[0]; i++){
      for(int j=0; j<nfft[1]; j++){
        for(int k=0; k<nfft[2]; k++){
          phase = 0;
          for (int l=0; l<3; l++){
            rijk = a1[l]*i/nfft[0] + a2[l]*j/nfft[1] + a3[l]*k/nfft[2];
            G0 = b1[l]*uklpp[0] + b2[l]*uklpp[1] + b3[l]*uklpp[2];
            G0 *= -2*M_PI/lattconst;
            phase += rijk*G0;
          }
          umklapp_factor[counter].re = cos(phase);
          umklapp_factor[counter].im = sin(phase);
          counter += 1;
        }// end k loop
      }// end j loop
    }// end i loop
  }//end if-else statement

}//end function

#include "eps_matrix.def.h"
