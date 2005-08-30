#include "matmul.h"
//#include "ckmulticast.h"


#ifdef FORTRANUNDERSCORE
#define DGEMM dgemm_
#else
#define DGEMM dgemm
#endif

extern "C" {void  DGEMM(char *, char *, int *, int *, int *, double *,
			   double *, int *, double *, int *, double *,
			   double *, int *);}

//#define USER_EVENTS
//#define MULT_EVENT

#ifdef USER_EVENTS
#define MULT_EVENT
#define RECV_EVENTS
#endif

#ifdef RECV_EVENTS
extern int user_event1;
extern int user_event2;
#endif

#ifdef MULT_EVENT
extern int cache_event;
extern int mult_event;
#endif

matmul::matmul(){
	row_neighbors = NULL;
	col_neighbors = NULL;
	dest = NULL;
	row_count = col_count = 0;
	fcb = NULL;
	use_cb = false;
	cb_obj = NULL;
//	abort_calc = false;
}

matmul::matmul(int rows, int chunksize, CkCallback mycb){
	this->rows = rows;
	this->chunksize = chunksize;
	this->chunksqrd = chunksize * chunksize;
	row_neighbors = new double[rows * chunksqrd];
	col_neighbors = new double[rows * chunksqrd];
	row_count = col_count = 0;
	row_group = get_row_group(thisIndex.x, thisIndex.y, rows, thisProxy);
	col_group = get_col_group(thisIndex.x, thisIndex.y, rows, thisProxy);
	dummy = 0;
	fcb = NULL;
	use_cb = false;
	cb_obj = NULL;
//	abort_calc = false;
	mycb.send(NULL);
}

/*
matmul::matmul(int rows, int chunksize, CkCallback mycb, CkGroupID gid){
	this->rows = rows;
	this->chunksize = chunksize;
	this->chunksqrd = chunksize * chunksize;
	row_neighbors = new double[rows * chunksqrd];
	col_neighbors = new double[rows * chunksqrd];
	row_count = col_count = 0;
	row_group = get_row_group(thisIndex.x, thisIndex.y, rows, thisProxy);
	col_group = get_col_group(thisIndex.x, thisIndex.y, rows, thisProxy);
	row_group.ckSectionDelegate(CProxy_CkMulticastMgr(gid).ckLocalBranch());
	col_group.ckSectionDelegate(CProxy_CkMulticastMgr(gid).ckLocalBranch());
	abort_calc = false;
	mycb.send(NULL);
}

matmul::matmul(int rows, int chunksize, CkCallback mycb,
   CProxySection_matmul row_group, CProxySection_matmul col_group){
	this->rows = rows;
	this->chunksize = chunksize;
	this->chunksqrd = chunksize * chunksize;
	this->row_group = row_group;
	this->col_group = col_group;
	row_neighbors = new double[rows * chunksqrd];
	col_neighbors = new double[rows * chunksqrd];
	row_count = col_count = 0;
	abort_calc = false;
	mycb.send(NULL);
}
*/

matmul::matmul(CkMigrateMessage *m){}

matmul::~matmul(){
	delete [] row_neighbors;
	delete [] col_neighbors;
}

void matmul::start_mat_mul(double *A, double *B, double *C, double alpha,
   double beta, CkCallback cb){
/*
	if(abort_calc){
		cb.send(NULL);
		return;
	}
*/
	/* initialize */
	row_count++; col_count++; /* we have our own data */
	this->alpha = alpha;
	this->beta = beta;
	this->dest = C;
	this->cb = cb;
	use_cb = true;

#ifdef MULT_EVENT
	double start_t = CmiWallTimer();
	for(int i = 0; i < chunksqrd; i++){
		dummy += A[i];
		dummy += B[i];
	}
	for(int i = 0; i < chunksqrd * rows; i++){
		dummy += row_neighbors[i];
		dummy += col_neighbors[i];
	}
	traceUserBracketEvent(cache_event, start_t, CmiWallTimer());
#endif

	/* copy my data */
	for(int i = 0; i < chunksize; i++)
		memcpy(&row_neighbors[rows * chunksize * i +
		   chunksize * thisIndex.y], &A[chunksize * i],
		   chunksize * sizeof(double));
	memcpy(col_neighbors + chunksqrd * thisIndex.x, B,
	   sizeof(double) * chunksqrd);

	/* multicast my data */
	mm_matrix_data *msgA = new (chunksqrd, 0)
	   mm_matrix_data(A, thisIndex.y, chunksqrd);
	mm_matrix_data *msgB = new (chunksqrd, 0)
	   mm_matrix_data(B, thisIndex.x, chunksqrd);
	row_group.receive_A(msgA);
	col_group.receive_B(msgB);

	if(row_count == rows && col_count == rows)
		multiply();
}

void matmul::start_mat_mul(double *A, double *B, double *C, double alpha,
   double beta, void (*fptr) (void *), void *obj){
	/* initialize */
	row_count++; col_count++; /* we have our own data */
	this->alpha = alpha;
	this->beta = beta;
	this->dest = C;
	this->cb = cb;
	use_cb = false;
	fcb = fptr;
	cb_obj = obj;

#ifdef MULT_EVENT
	double start_t = CmiWallTimer();
	for(int i = 0; i < chunksqrd; i++){
		dummy += A[i];
		dummy += B[i];
	}
	for(int i = 0; i < chunksqrd * rows; i++){
		dummy += row_neighbors[i];
		dummy += col_neighbors[i];
	}
	traceUserBracketEvent(cache_event, start_t, CmiWallTimer());
#endif

	/* copy my data */
	for(int i = 0; i < chunksize; i++)
		memcpy(&row_neighbors[rows * chunksize * i +
		   chunksize * thisIndex.y], &A[chunksize * i],
		   chunksize * sizeof(double));
	memcpy(col_neighbors + chunksqrd * thisIndex.x, B,
	   sizeof(double) * chunksqrd);

	/* multicast my data */
	mm_matrix_data *msgA = new (chunksqrd, 0)
	   mm_matrix_data(A, thisIndex.y, chunksqrd);
	mm_matrix_data *msgB = new (chunksqrd, 0)
	   mm_matrix_data(B, thisIndex.x, chunksqrd);
	row_group.receive_A(msgA);
	col_group.receive_B(msgB);

	if(row_count == rows && col_count == rows)
		multiply();
}

void matmul::receive_A(mm_matrix_data *msg){
/*
	if(abort_calc)
		return;
*/
#ifdef RECV_EVENTS
	double start_t = CmiWallTimer();
#endif
	/* receive data */
	row_count++;
	for(int i = 0; i < chunksize; i++)
		memcpy(&row_neighbors[rows * chunksize * i +
		   chunksize * msg->from], &msg->data[i * chunksize],
		   chunksize * sizeof(double));
	delete msg;

	/* if we have received all data, multiply and report */
#ifdef RECV_EVENTS
	traceUserBracketEvent(user_event1, start_t, CmiWallTimer());
#endif
	if(row_count == rows && col_count == rows)
		multiply();
}

void matmul::receive_B(mm_matrix_data *msg){
/*
	if(abort_calc)
		return;
*/
#ifdef RECV_EVENTS
	double start_t = CmiWallTimer();
#endif
	/* receive data */
	col_count++;
	memcpy(col_neighbors + msg->from * chunksqrd, msg->data,
	   sizeof(double) * chunksqrd);

	delete msg;

	/* if we have received all data, multiply and report */
#ifdef RECV_EVENTS
	traceUserBracketEvent(user_event2, start_t, CmiWallTimer());
#endif
	if(row_count == rows && col_count == rows)
		multiply();
}

void matmul::multiply(void){
#ifdef MULT_EVENT
	double start_t = CmiWallTimer();
	for(int i = 0; i < chunksqrd; i++)
		dummy += dest[i];
	for(int i = 0; i < chunksqrd * rows; i++){
		dummy += row_neighbors[i];
		dummy += col_neighbors[i];
	}
	traceUserBracketEvent(cache_event, start_t, CmiWallTimer());
	start_t = CmiWallTimer();
#endif
	/* reset counters */
	row_count = col_count = 0;
	/* transpose result matrix (if Beta != 0) */
	if(beta != 0)
		for(int i = 0; i < chunksize; i++)
			for(int j = i + 1; j < chunksize; j++){
				double tmp = dest[i * chunksize + j];
				dest[i * chunksize + j] = dest[j*chunksize + i];
				dest[j * chunksize + i] = tmp;
			}
	

	/* compute */
	int tmp = rows * chunksize;
	char trans = 't';
	DGEMM(&trans, &trans, &chunksize, &chunksize, &tmp, &alpha,
	   row_neighbors, &tmp, col_neighbors, &chunksize, &beta, dest,
	   &chunksize);
	  
	/* transpose result matrix */
	for(int i = 0; i < chunksize; i++)
		for(int j = i + 1; j < chunksize; j++){
			double tmp = dest[i * chunksize + j];
			dest[i * chunksize + j] = dest[j * chunksize + i];
			dest[j * chunksize + i] = tmp;
		}

	/* alert caller */
	if(use_cb)
		cb.send(NULL);
	else
		(*fcb)(cb_obj);
#ifdef MULT_EVENT
	dummy = 0;
	traceUserBracketEvent(mult_event, start_t, CmiWallTimer());
#endif
}

void matmul::pup(PUP::er &p){
	CBase_matmul::pup(p);
	p | rows;
	p | chunksize;
	p | chunksqrd;
	p | row_count;
	p | col_count;
	p | alpha;
	p | beta;
	p | cb;
	p | row_group;
	p | col_group;
	p | dummy;
	p | use_cb;
//	p | abort_calc;
	if(p.isUnpacking()){
		row_neighbors = new double[chunksqrd * rows];
		col_neighbors = new double[chunksqrd * rows];
	}
	p(row_neighbors, chunksqrd * rows);
	p(col_neighbors, chunksqrd * rows);
}

void matmul::reset(double *C){
	this->dest = C;
}

void matmul::reset(double *C, void (*fptr) (void *), void *obj){
	this->dest = C;
	this->fcb = fptr;
	this->cb_obj = obj;
}
/*
void matmul::abort(){
	abort_calc = true;
}
*/

CProxySection_matmul matmul::get_row_group(int xpos, int ypos, int n,
   CProxy_matmul &proxy){
	CkVec<CkArrayIndexMax> elems;
	for(int i = 0, j = 0; i < n; i++)
		if(i != ypos){
			elems.push_back(CkArrayIndex2D(xpos, i));
			j++;
		}
	return CProxySection_matmul::ckNew(proxy, elems.getVec(), elems.size());
}

CProxySection_matmul matmul::get_col_group(int xpos, int ypos, int n,
   CProxy_matmul &proxy){
	CkVec<CkArrayIndexMax> elems;
	for(int i = 0, j = 0; i < n; i++)
		if(i != xpos){
			elems.push_back(CkArrayIndex2D(i, ypos));
			j++;
		}
	return CProxySection_matmul::ckNew(proxy, elems.getVec(), elems.size());
}

mm_matrix_data::mm_matrix_data(){
}

mm_matrix_data::mm_matrix_data(double *data, int from, int size){
	this->from = from;
	this->size = size;
	memcpy(this->data, data, size * sizeof(double));
}

#include "matmul.def.h"
