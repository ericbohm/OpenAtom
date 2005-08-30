#ifndef _MATMUL_H_
#define _MATMUL_H_

#include "matmul.decl.h"

class matmul:public CBase_matmul {
	public:
		matmul();
		matmul(int rows, int chunksize, CkCallback mycb);
//		matmul(int rows, int chunksize, CkCallback mycb, CkGroupID gid);
//		matmul(int rows, int chunksize, CkCallback mycb,
//		   CProxySection_matmul row_group,
//		   CProxySection_matmul col_group);
		matmul(CkMigrateMessage *m);
		~matmul();
		/* callback returns NULL message */
		void start_mat_mul(double *A, double *B, double *C,
		   double alpha, double beta, CkCallback cb);
		void start_mat_mul(double *A, double *B, double *C,
		   double alpha, double beta, void (*fptr) (void *), void *obj);
		void receive_A(mm_matrix_data *msg);
		void receive_B(mm_matrix_data *msg);
//		void abort();
		/* bound class should call reset while unpacking */
		void reset(double *C);
		void reset(double *C, void (*fptr) (void *), void *obj);
		static CProxySection_matmul get_row_group(int xpos, int ypos,
		   int n, CProxy_matmul &proxy);
		static CProxySection_matmul get_col_group(int xpos, int ypos,
		   int n, CProxy_matmul &proxy);
		virtual void pup(PUP::er &p);
	private:
		void multiply(void);
		int rows;
		int chunksize;
		int chunksqrd;
		int row_count;
		int col_count;
		double *row_neighbors;
		double *col_neighbors;
		double alpha, beta;
		double *dest;
		CkCallback cb;
		CProxySection_matmul row_group, col_group;
		double dummy;
		void (*fcb) (void *obj);
		bool use_cb;
		void *cb_obj;
//		bool abort_calc;
};

class mm_matrix_data:public CMessage_mm_matrix_data, public CkMcastBaseMsg{
	public:
		mm_matrix_data();
		mm_matrix_data(double *data, int from, int size);
		int size;
		int from;
		double *data;
};

#endif
