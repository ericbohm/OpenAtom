#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <sys/time.h>
//#include <unistd.h> // for sleep, i guess
#include <mpi.h>
#include "mpi-interoperate.h"
#include "diagonalizer.h"
#include "configure.h"
//#include "preprocessor.h"

using namespace std;

extern "C" {
// BLACS
 void Cblacs_pinfo( int* mypnum, int* nprocs);
 void Cblacs_get( int context, int request, int* value);
 int Cblacs_gridinit( int* context, const char * order, int np_row, int np_col);
 void Cblacs_gridinfo( int context, int*  np_row, int* np_col, int*  my_row, int*  my_col);
 void Cblacs_gridexit( int context);
 void Cblacs_exit( int error_code);

// BLAS
 void dscal_( int *n, double *da, double *dx, int *incx);

// SCALAPACK
 void pdlawrite_( char **filenam, int *m, int *n, double *A, int *ia, int *ja, int *descA, int *irwrit, int *icwrit, double *work);
 void pdelset_( double *A, int *ia, int *ja, int *desca, double *alpha);
 double pdlamch_( int *ictxt, char *cmach);
 int indxg2p_( int *indxglob, int *nb, int *iproc, int *isrcproc, int *nprocs);
 int indxg2l_( int *indxglob, int *nb, int *iproc, int *isrcproc, int *nprocs);
 int numroc_( int *n, int *nb, int *iproc, int *isrcproc, int *nprocs);
 void descinit_( int *desc, int *m, int *n, int *mb, int *nb, int *irsrc, int *icsrc,
                                int *ictxt, int *lld, int *info);
 void pdlaset_( char *uplo, int *m, int *n, double *alpha, double *beta, double *A, int *ia, int *ja, int *descA );
 double pdlange_( char *norm, int *m, int *n, double *A, int *ia, int *ja, int *desca, double *work);
 void pdlacpy_( char *uplo, int *m, int *n, double *a, int *ia, int *ja, int *desca,
                                double *b, int *ib, int *jb, int *descb);
 void pdgesv_( int *n, int *nrhs, double *A, int *ia, int *ja, int *desca, int* ipiv,
                                double *B, int *ib, int *jb, int *descb, int *info);
 void pdgesvd_( char *jobu, char *jobvt, int *m, int *n, double *a, int *ia, int *ja, int *desca,
                                double *s, double *u, int *iu, int *ju, int *descu,
                                double *vt, int *ivt, int *jvt, int *descvt, double *work, int *lwork, int *info);
 void pdgemm_( char *TRANSA, char *TRANSB, int * M, int * N, int * K, double * ALPHA,
                                double * A, int * IA, int * JA, int * DESCA, double * B, int * IB, int * JB, int * DESCB,
                                double * BETA, double * C, int * IC, int * JC, int * DESCC );
 int indxg2p_( int *indxglob, int *nb, int *iproc, int *isrcproc, int *nprocs);


 void pdsyev_( char *jobz, char *uplo, int *n,
                double *a, int *ia, int *ja, int *desca, double *w, double *z, int *iz, int *jz, int *descz,
                double *work, int *lwork, int *info );
}

static int max( int a, int b ){
        if (a>b) return(a); else return(b);
}
static int min( int a, int b ){
        if (a<b) return(a); else return(b);
}

//double *mpiarray;

class Config;
extern diagData_t<double> *diagData;
extern int numOrthosPerDim;
extern int grainSizeOrtho;
extern int totalOrthos;
extern int diagonalization;
extern int numStatesOA;
extern Config config;

void restartcharm();

int main(int argc, char **argv) {
        diagData = NULL;
        int iam, nprocs;
        int myrank_mpi, nprocs_mpi;
        int ictxt, nprow, npcol, myrow, mycol;
        int np, nq, nb, n;
        int mpA, nqA;
        int i, j, k, info, itemp, seed, lwork, min_mn;
        int descA[9], descZ[9];
        double *A, *Z, *work, *W;
        int izero=0,ione=1;
        double mone=(-1.0e0),pone=(1.0e0),dzero=(0.0e0);
        double MPIt1, MPIt2, MPIelapsed, GFLOPS, GFLOPS_per_proc ;
        char jobz, uplo;
        MPI_Init( &argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank_mpi);
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs_mpi);

        MPI_Barrier(MPI_COMM_WORLD);
        CharmLibInit(MPI_COMM_WORLD, argc, argv);
        MPI_Barrier(MPI_COMM_WORLD);

        n = config.nstates; nprow = numOrthosPerDim; npcol = numOrthosPerDim; nb = config.orthoGrainSize; jobz= 'V'; uplo='U';
        //int mpisize = 5;
        //for(int ppp = 0 ; ppp < mpisize ; ppp++) {
        //  printf("mpiarray: <%5.2f>\n", mpiarray[0]);
        //}
        for( i = 1; i < argc; i++ ) {
                if( strcmp( argv[i], "-jobz" ) == 0 ) {
                        if (i+1<argc) {
                                if( strcmp( argv[i+1], "V" ) == 0 ){ jobz = 'V'; i++; }
                                else if( strcmp( argv[i+1], "N" ) == 0 ){ jobz = 'N'; i++; }
                                else if( strcmp( argv[i+1], "A" ) == 0 ){ jobz = 'A'; i++; }
                                else printf(" ** warning: jobu should be set to V, N or A in the command line ** \n");
                        }
                        else
                                printf(" ** warning: jobu should be set to V, N or A in the command line ** \n");
                }
                if( strcmp( argv[i], "-n" ) == 0 ) {
                        n      = atoi(argv[i+1]);
                        i++;
                }
                if( strcmp( argv[i], "-p" ) == 0 ) {
                        nprow  = atoi(argv[i+1]);
                        i++;
                }
                if( strcmp( argv[i], "-q" ) == 0 ) {
                        npcol  = atoi(argv[i+1]);
                        i++;
                }
                if( strcmp( argv[i], "-nb" ) == 0 ) {
                        nb     = atoi(argv[i+1]);
                        i++;
                }
        }
        if (nb>n)
                nb = n;
        if (nprow*npcol>nprocs_mpi){
                if (myrank_mpi==0)
                        printf(" **** ERROR : we do not have enough processes available to make a p-by-q process grid ***\n");
                        printf(" **** Bye-bye                                                                                                                                                            ***\n");
                MPI_Finalize(); exit(1);
        }
        Cblacs_pinfo( &iam, &nprocs ) ;
        Cblacs_get( 0, 0, &ictxt );
        Cblacs_gridinit( &ictxt, "R", nprow, npcol );
        Cblacs_gridinfo( ictxt, &nprow, &npcol, &myrow, &mycol );
        for(int iternum = 1 ; iternum <=4 ; iternum++) {
        printf("Prateek: rank: <%d>, plambda: <%.2f,%.2f,%.2f>\n", myrank_mpi, diagData->plambda[0], diagData->plambda[1], diagData->plambda[2]);
        printf("Prateek: rank: <%d>: pelements: <%d>\n", myrank_mpi, diagData->pelements);
        if ((myrow < nprow)&(mycol < npcol)){
                mpA    = numroc_( &n     , &nb, &myrow, &izero, &nprow );
                nqA    = numroc_( &n     , &nb, &mycol, &izero, &npcol );
                if(mpA!=nb) {
                  exit(1);
                }
                if(nqA!=nb){
                  exit(1);
                }
                min_mn = mpA;
                min_mn = nqA;
                min_mn = n;
                printf("mpa: <%d>, nqa: <%d>, n: <%d>, nb: <%d>, nprow: <%d>, npcol: <%d>, nprocs: <%d>\n", mpA, nqA, n, nb, nprow, npcol, nprocs);
                
                A = (double *)calloc(mpA*nqA,sizeof(double)) ;
                if (A==NULL){ printf("error of memory allocation A on proc %dx%d\n",myrow,mycol); exit(0); }
                Z = (double *)calloc(mpA*nqA,sizeof(double)) ;
                if (Z==NULL){ printf("error of memory allocation VT on proc %dx%d\n",myrow,mycol); exit(0); }
                W = (double *)calloc(min_mn,sizeof(double)) ;
                if (W==NULL){ printf("error of memory allocation S on proc %dx%d\n",myrow,mycol); exit(0); }
                seed = iam*(mpA*nqA*2); srand(seed);
                k = 0;
                if(diagData->pelements != ((nb*nb)+3)){
                  exit(1);
                }
                for (i = 0; i < mpA; i++) {
                        for (j = 0; j < nqA; j++) {
                                //A[k] = ((double) rand()) / ((double) RAND_MAX) - 0.5 ;
                                A[k] = diagData->plambda[3+(i*mpA)+j]; // ((double) rand()) / ((double) RAND_MAX) - 0.5 ;
                                k++;
                        }
                }
                itemp = max( 1, mpA );
                descinit_( descA,  &n, &n, &nb, &nb, &izero, &izero, &ictxt, &itemp, &info );
                descinit_( descZ,  &n, &n, &nb, &nb, &izero, &izero, &ictxt, &itemp, &info );
                work = (double *)calloc(2,sizeof(double)) ;
                if (work==NULL){ printf("error of memory allocation for work on proc %dx%d (1st time)\n",myrow,mycol); exit(0); }
                lwork=-1;
                printf("rank: <%d>, A: <%.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f>, Z: <%.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f>\n", myrank_mpi, A[0], A[1], A[2], A[3], A[64], A[128], A[192], A[256], Z[0], Z[1], Z[2], Z[3], Z[64], Z[128], Z[192], Z[256]);
                pdsyev_( &jobz, &uplo, &n, A, &ione, &ione, descA, W, Z, &ione, &ione, descZ, work, &lwork, &info );
                lwork= (int) work[0];
                free(work);
                work = (double *)calloc(lwork,sizeof(double)) ;
                if (work==NULL){ printf("error of memory allocation work on proc %dx%d\n",myrow,mycol); exit(0); }
                MPIt1 = MPI_Wtime();

                pdsyev_( &jobz, &uplo, &n, A, &ione, &ione, descA, W, Z, &ione, &ione, descZ, work, &lwork, &info );

                diagData->rlambda = new double[(mpA*nqA)+3];
                diagData->rlambda[0] = diagData->plambda[0];
                diagData->rlambda[1] = diagData->plambda[1];
                //double xind = (double) myrank_mpi
                diagData->rlambda[2] = diagData->plambda[2]; //(double) mpA * nqA;
                for (i = 0; i < mpA; i++) {
                        for (j = 0; j < nqA; j++) {
                                //A[k] = ((double) rand()) / ((double) RAND_MAX) - 0.5 ;
                                diagData->rlambda[3+(i*mpA)+j] = Z[(i*mpA)+j];
                        }
                }

                printf("rank: <%d>, A: <%.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f>, Z: <%.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f>\n", myrank_mpi, A[0], A[1], A[2], A[3], A[64], A[128], A[192], A[256], Z[0], Z[1], Z[2], Z[3], Z[64], Z[128], Z[192], Z[256]);
                MPIt2 = MPI_Wtime();
                MPIelapsed=MPIt2-MPIt1;
                free(work);
                if ( iam==0 ){
                        printf("n=%d\t(%d,%d)\t%d\tjobz=%c\t%8.2fs \n",n,nprow,npcol,nb,jobz,MPIelapsed);
                }
                free(W);
                free(Z);
                free(A);
                
        }
        if ( iam==0 ){
                printf("\n");
        }

        MPI_Barrier(MPI_COMM_WORLD);
        restartcharm();
        MPI_Barrier(MPI_COMM_WORLD);
        }

        Cblacs_gridexit( 0 );
        // Cblacs_gridexit(ctxt); - spiros way of doing it
        MPI_Finalize();
        exit(0);
}
