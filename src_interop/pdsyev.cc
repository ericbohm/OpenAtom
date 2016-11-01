#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>
#include <mpi.h>
#include "mpi-interoperate.h"
#include "diagonalizer.h"
#include "configure.h"

using namespace std;

extern "C" {
// BLACS
  void Cblacs_pinfo(int* mypnum, int* nprocs);
  void Cblacs_get(int context, int request, int* value);
  int Cblacs_gridinit(int* context, const char * order, int np_row, int np_col);
  void Cblacs_gridinfo(int context, int* np_row, int* np_col, int* my_row, int* my_col);
  void Cblacs_gridexit(int context);
  void Cblacs_exit(int error_code);

// SCALAPACK
  int numroc_(int *n, int *nb, int *iproc, int *isrcproc, int *nprocs);
  void descinit_(int *desc, int *m, int *n, int *mb, int *nb, int *irsrc, int *icsrc,
                 int *ictxt, int *lld, int *info);
  void pdsyev_(char *jobz, char *uplo, int *n, double *a, int *ia, int *ja,
                 int *desca, double *w, double *z, int *iz, int *jz, int *descz,
                 double *work, int *lwork, int *info);
}

static int max(int a, int b) {
  if (a>b) return(a); else return(b);
}

static int min(int a, int b) {
  if (a<b) return(a); else return(b);
}

class Config;
extern diagData_t<internalType> *diagData;
extern int numOrthosPerDim;
extern int numEOrthosPerDim;
extern int grainSizeOrtho;
extern int totalOrthos;
extern int diagonalization;
extern int numStatesOA;
extern Config config;

void restartcharm();

int main(int argc, char **argv) {
#ifdef CP_DIAGONALIZER_DEBUG
    setbuf(stdout, NULL);
#endif
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
    double MPIt1, MPIt2, MPIelapsed;
    char jobz, uplo;
    MPI_Init( &argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank_mpi);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs_mpi);

    MPI_Barrier(MPI_COMM_WORLD);
    CharmLibInit(MPI_COMM_WORLD, argc, argv);
    MPI_Barrier(MPI_COMM_WORLD);

    n = numStatesOA;
    nprow = numEOrthosPerDim;
    npcol = numEOrthosPerDim;
    nb = grainSizeOrtho;
    jobz= 'V';
    uplo='U';

    if (nprow*npcol>nprocs_mpi) {
        if (myrank_mpi==0) {
          printf(" **** ERROR : we do not have enough processes available to make a p-by-q process grid ***\n");
        }
        printf(" **** Bye-bye ***\n");
        MPI_Finalize();
        exit(1);
    }
    Cblacs_pinfo(&iam, &nprocs);
    Cblacs_get(0, 0, &ictxt);
    Cblacs_gridinit(&ictxt, "R", nprow, npcol);
    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);
    for(int iternum = 1 ; iternum <= 1 ; iternum++) {
        if ((myrow >= 0) && (mycol >= 0) && (myrow < nprow) && (mycol < npcol)) {
            mpA = numroc_(&n, &nb, &myrow, &izero, &nprow);
            nqA = numroc_(&n, &nb, &mycol, &izero, &npcol);
            min_mn = n;
            A = (double *)calloc(mpA*nqA,sizeof(double)) ;
            if (A==NULL){ printf("error of memory allocation A on proc %dx%d\n",myrow,mycol); exit(0); }
            Z = (double *)calloc(mpA*nqA,sizeof(double)) ;
            if (Z==NULL){ printf("error of memory allocation VT on proc %dx%d\n",myrow,mycol); exit(0); }
            W = (double *)calloc(min_mn,sizeof(double)) ;
            if (W==NULL){ printf("error of memory allocation S on proc %dx%d\n",myrow,mycol); exit(0); }
            seed = iam*(mpA*nqA*2); srand(seed);
            k = 0;
            for (i = 0; i < mpA; i++) {
                for (j = 0; j < nqA; j++) {
                    A[k] = diagData->plambda[(i*nqA)+j];
                    k++;
                }
            }
            itemp = max( 1, mpA );
            descinit_( descA,  &n, &n, &nb, &nb, &izero, &izero, &ictxt, &itemp, &info );
            descinit_( descZ,  &n, &n, &nb, &nb, &izero, &izero, &ictxt, &itemp, &info );
            work = (double *)calloc(2,sizeof(double)) ;
            if (work==NULL){ printf("error of memory allocation for work on proc %dx%d (1st time)\n",myrow,mycol); exit(0); }
            lwork=-1;
            pdsyev_( &jobz, &uplo, &n, A, &ione, &ione, descA, W, Z, &ione, &ione, descZ, work, &lwork, &info );
            lwork= (int) work[0];
            free(work);
            work = (double *)calloc(lwork,sizeof(double)) ;
            if (work==NULL){ printf("error of memory allocation work on proc %dx%d\n",myrow,mycol); exit(0); }
            MPIt1 = MPI_Wtime();

            pdsyev_( &jobz, &uplo, &n, A, &ione, &ione, descA, W, Z, &ione, &ione, descZ, work, &lwork, &info );
#ifdef CP_DIAGONALIZER_DEBUG
            printf("step 3 on proc : [%d] : [%d,%d] : info = [%d]\n", myrank_mpi, myrow, mycol, info);
            fflush(stdout);
#endif
            int looptot = mpA * nqA;
            diagData->rlambda = new double[looptot];
            for (i = 0; i < mpA; i++) {
                for (j = 0; j < nqA; j++) {
                    diagData->rlambda[(i*nqA)+j] = Z[(i*nqA)+j];
                }
            }

            MPIt2 = MPI_Wtime();
            MPIelapsed=MPIt2-MPIt1;
            if (iam == 0){
                    printf("n=%d\t(%d,%d)\t%d\tjobz=%c\t%8.2fs \n",n,nprow,npcol,nb,jobz,MPIelapsed);
            }
            free(work);
            free(W);
            free(Z);
            free(A);
        }
        else {
          printf("I <rank: <%d>> did not participate in diagonalization.\n", myrank_mpi);
          fflush(stdout);
        }
        if (iam==0) {
            printf("\n");
        }

        MPI_Barrier(MPI_COMM_WORLD);
        restartcharm();
        MPI_Barrier(MPI_COMM_WORLD);
    }

    if ((myrow >= 0) && (mycol >= 0) && (myrow < nprow) && (mycol < npcol)) {
      Cblacs_gridexit (ictxt);
    }

    MPI_Finalize();
    exit(0);
}
