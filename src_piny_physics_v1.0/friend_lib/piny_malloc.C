#include "standard_include.h"
#include "ckcomplex.h"
#include "../proto_defs/proto_friend_lib_entry.h"

#ifdef DUAL_FFTW
#include <dfftw.h>
#include <drfftw.h>
#else
#include <fftw.h>
#include <rfftw.h>
#endif

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* cmalloc: Careful malloc                                                  */
/*==========================================================================*/
void *cmalloc(size_t len,char *func_name)
{ /* begin routine */
  void *mem_ptr;
  double request;
  
  if(len == 0) return NULL;
  mem_ptr = fftw_malloc(len);
  if(mem_ptr == NULL) {
   request = ((double) len)*1.0e-6;
   PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
   PRINTF("Dude, like you've just requested %g MBytes\n",request);
   PRINTF("of memory -- get real\n\n");
   PRINTF("Error occurred in function %s\n",func_name);
   PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
   FFLUSH(stdout);
   EXIT(1);
  }

  memset(mem_ptr,0,len);

  return mem_ptr;

} /* end routine */

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* cfree: Careful free                                                      */
/*==========================================================================*/

void cfree(void *p,char *func_name)
{/* begin routine */
  if(p==NULL) PRINTF("Freeing unallocated pointer in function %s\n",func_name);
  fftw_free(p);
}/* end routine */


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* crealloc: Careful realloc                                                */
/*==========================================================================*/
void *crealloc(void *ptr,size_t len,char *func_name)
{ /* begin routine */
  void *mem_ptr;
  double request;
  
  if(len==0){
    fftw_free(ptr);
    return NULL;
  }

  mem_ptr = realloc(ptr,len);
  if(mem_ptr == NULL) {
   request = ((double) len)*1.0e-6;
   PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
   PRINTF("Dude, like your reallocating %g MBytes\n",request);
   PRINTF("of memory -- get real (address was %p)\n\n",ptr);
   PRINTF("Error occurred in function %s\n",func_name);
   PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
   FFLUSH(stdout);
   EXIT(1);
  }

  return mem_ptr;

} /* end routine */
/*===========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* cmall_int_mat: Careful malloc a matrix of type int                       */
/*==========================================================================*/
int **cmall_int_mat(long nrl, long nrh, long ncl, long nch,char *func_name)
/* allocate an integer matrix with subscript range m[nrl...nrh][ncl...nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  int **m;

  if(nrow <= 0 || ncol <= 0) return (int **)NULL;
/* allocate pointers to rows */
  m=(int **) fftw_malloc((size_t)((nrow+NR_END)*sizeof(int*)));
  if(!m) 
    {
     PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     PRINTF("allocation failure of row pointers in cmall_mat\n");
     PRINTF("Error occurred in function %s\n",func_name);
     PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     FFLUSH(stdout);
     EXIT(1);
    }
  m += NR_END;
  m -= nrl;

/* allocate rows and set pointers to them */
  m[nrl]=(int *) fftw_malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
  if(!m[nrl]) 
    {
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("allocation failure of rows in cmall_mat\n");
      PRINTF("Error occurred in function %s\n",func_name);
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) 
    m[i]=m[i-1]+ncol;

/* return pointer to array of pointers to rows */
   return m;  
}/* end routine */


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* cmall_mat: Careful malloc a matrix of type double                        */
/*==========================================================================*/
double **cmall_mat(long nrl, long nrh, long ncl, long nch,char *func_name)
/* allocate a double matrix with subscript range m[nrl...nrh][ncl...nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  double **m;

  if(nrow <= 0 || ncol <= 0) return (double **)NULL;
/* allocate pointers to rows */
  m=(double **) fftw_malloc((size_t)((nrow+NR_END)*sizeof(double*)));
  if(!m) 
    {
     PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     PRINTF("allocation failure of row pointers in cmall_mat\n");
     PRINTF("Error occurred in function %s\n",func_name);
     PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     FFLUSH(stdout);
     EXIT(1);
    }
  m += NR_END;
  m -= nrl;

/* allocate rows and set pointers to them */
  m[nrl]=(double *) fftw_malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
  if(!m[nrl]) 
    {
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("allocation failure of rows in cmall_mat\n");
      PRINTF("Error occurred in function %s\n",func_name);
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) 
    m[i]=m[i-1]+ncol;

/* return pointer to array of pointers to rows */
   return m;  
}/* end routine */

/*===========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* cmall_tens3: Careful malloc a rank 3 tensor of type double               */
/*==========================================================================*/
double ***cmall_tens3(long nrl,long nrh,long ncl,long nch,long ndl,long ndh, char *func_name)

/* allocate a double rank3 tensor with subscript range                       */
/* m[nrl...nrh][ncl...nch][ndl...ndh] */
{
  long i,j, nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
  double ***m;

  if(nrow <= 0 || ncol <= 0 || ndep <= 0) return (double ***)NULL;
/* allocate pointers to rows */
  m=(double ***) fftw_malloc((size_t)((nrow+NR_END)*sizeof(double**)));
  if(!m) 
    {
     PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     PRINTF("allocation failure of row pointers cmall_tens3\n");
     PRINTF("Error occurred in function %s\n",func_name);
     PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     FFLUSH(stdout);
     EXIT(1);
    }
  m += NR_END;
  m -= nrl;

/* allocate rows and set pointers to them */
  m[nrl]=(double **) fftw_malloc((size_t)((nrow*ncol+NR_END)*sizeof(double*)));
  if(!m[nrl]) 
    {
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("allocation failure of rows in cmall_tens3\n");
      PRINTF("Error occurred in function %s\n",func_name);
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }
  m[nrl] += NR_END;
  m[nrl] -= ncl;

/* allocate columns and set pointers to them */
  m[nrl][ncl]=(double *) fftw_malloc((size_t)((nrow*ncol*ndep+NR_END)
					 *sizeof(double)));
  if(!m[nrl][ncl]) 
    {
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("allocation failure of columns in cmall_tens3\n");
      PRINTF("Error occurred in function %s\n",func_name);
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }
  m[nrl][ncl] += NR_END;
  m[nrl][ncl] -= ndl;


   for(i=ncl+1;i<=nch;i++) {
     m[nrl][i]=m[nrl][i-1]+ndep;
   }/*endfor*/

   for(i=nrl+1;i<=nrh;i++) {
     m[i]=m[i-1]+ncol;
     m[i][ncl]=m[i-1][ncl]+ncol*ndep;
    for(j=ncl+1;j<=nch;j++) {
      m[i][j]=m[i][j-1]+ndep;
    }/*endfor*/
   }/*endfor*/

/* return pointer to array of pointers to rows */
   return m;  
}/* end routine */

/*===========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* cmall_itens3: Careful malloc a rank 3 tensor of type int               */
/*==========================================================================*/
int ***cmall_itens3(long nrl,long nrh,long ncl,long nch,long ndl,long ndh, char *func_name)

/* allocate a int rank3 tensor with subscript range                       */
/* m[nrl...nrh][ncl...nch][ndl...ndh] */
{
  long i,j, nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
  int ***m;

  if(nrow <= 0 || ncol <= 0 || ndep <= 0) return (int ***)NULL;
/* allocate pointers to rows */
  m=(int ***) fftw_malloc((size_t)((nrow+NR_END)*sizeof(int**)));
  if(!m) 
    {
     PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     PRINTF("allocation failure of row pointers cmall_tens3\n");
     PRINTF("Error occurred in function %s\n",func_name);
     PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     FFLUSH(stdout);
     EXIT(1);
    }
  m += NR_END;
  m -= nrl;

/* allocate rows and set pointers to them */
  m[nrl]=(int **) fftw_malloc((size_t)((nrow*ncol+NR_END)*sizeof(int*)));
  if(!m[nrl]) 
    {
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("allocation failure of rows in cmall_tens3\n");
      PRINTF("Error occurred in function %s\n",func_name);
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }
  m[nrl] += NR_END;
  m[nrl] -= ncl;

/* allocate columns and set pointers to them */
  m[nrl][ncl]=(int *) fftw_malloc((size_t)((nrow*ncol*ndep+NR_END)
					 *sizeof(int)));
  if(!m[nrl][ncl]) 
    {
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("allocation failure of columns in cmall_tens3\n");
      PRINTF("Error occurred in function %s\n",func_name);
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }
  m[nrl][ncl] += NR_END;
  m[nrl][ncl] -= ndl;


   for(i=ncl+1;i<=nch;i++) {
     m[nrl][i]=m[nrl][i-1]+ndep;
   }/*endfor*/

   for(i=nrl+1;i<=nrh;i++) {
     m[i]=m[i-1]+ncol;
     m[i][ncl]=m[i-1][ncl]+ncol*ndep;
    for(j=ncl+1;j<=nch;j++) {
      m[i][j]=m[i][j-1]+ndep;
    }/*endfor*/
   }/*endfor*/

/* return pointer to array of pointers to rows */
   return m;  
}/* end routine */
/*===========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* cmall_tens4: Careful malloc a rank 4 tensor of type double               */
/*==========================================================================*/
double ****cmall_tens4(long nrl,long nrh,long ncl,long nch,
		       long ndl,long ndh,long ntl,long nth,char *func_name)

/* allocate a double rank4 tensor with subscript range                       */
/* m[nrl...nrh][ncl...nch][ndl...ndh][ntl...nth] */
{
  long i,j,k;
  long nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1,ntes=nth-ntl+1;
  double ****m;

  if(nrow <= 0 || ncol <= 0 || ndep <= 0 || ntes <= 0) {
    return (double ****)NULL;
  }

/* allocate pointers to rows */
  m=(double ****) fftw_malloc((size_t)((nrow+NR_END)*sizeof(double***)));
  if(!m) 
    {
     PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     PRINTF("allocation failure of row pointers cmall_tens4\n");
     PRINTF("Error occurred in function %s\n",func_name);
     PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     FFLUSH(stdout);
     EXIT(1);
    }
  m += NR_END;
  m -= nrl;

/* allocate rows and set pointers to them */
  m[nrl]=(double ***) fftw_malloc((size_t)((nrow*ncol+NR_END)*sizeof(double**)));
  if(!m[nrl]) 
    {
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("allocation failure of rows in cmall_tens4\n");
      PRINTF("Error occurred in function %s\n",func_name);
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }
  m[nrl] += NR_END;
  m[nrl] -= ncl;

/* allocate columns and set pointers to them */
  m[nrl][ncl]=(double **) fftw_malloc((size_t)((nrow*ncol*ndep+NR_END)
					 *sizeof(double*)));
  if(!m[nrl][ncl]) 
    {
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("allocation failure of columns in cmall_tens4\n");
      PRINTF("Error occurred in function %s\n",func_name);
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }
  m[nrl][ncl] += NR_END;
  m[nrl][ncl] -= ndl;

/* allocate depths and set pointers to them */
  m[nrl][ncl][ndl]=(double *) fftw_malloc((size_t)((nrow*ncol*ndep*ntes+NR_END)
					       *sizeof(double)));
  if(!m[nrl][ncl][ndl]) 
    {
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("allocation failure of depths in cmall_tens4\n");
      PRINTF("Error occurred in function %s\n",func_name);
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }
  m[nrl][ncl][ndl] += NR_END;
  m[nrl][ncl][ndl] -= ntl;


   for(i=ndl+1;i<=ndh;i++) {
     m[nrl][ncl][i]=m[nrl][ncl][i-1]+ntes;
   }/*endfor*/

   for(i=ncl+1;i<=nch;i++) {
     m[nrl][i]=m[nrl][i-1]+ndep;
     m[nrl][i][ndl]=m[nrl][i-1][ndl]+ndep*ntes;
    for(j=ndl+1;j<=ndh;j++) {
      m[nrl][i][j]=m[nrl][i][j-1]+ntes;
    }/*endfor*/
   }/*endfor*/


   for(i=nrl+1;i<=nrh;i++) {
     m[i]=m[i-1]+ncol;
     m[i][ncl]=m[i-1][ncl]+ncol*ndep;
     m[i][ncl][ndl]=m[i-1][ncl][ndl]+ncol*ndep*ntes;
    for(k=ndl+1;k<=ndh;k++) {
     m[i][ncl][k]=m[i][ncl][k-1]+ndep;
    }/*endfor*/
    for(j=ncl+1;j<=nch;j++) {
      m[i][j]=m[i][j-1]+ndep;
      m[i][j][ndl]=m[i][j-1][ndl]+ncol*ndep;
     for(k=ndl+1;k<=ndh;k++) {
      m[i][j][k]=m[i][j-1][k]+ntes;
      m[i][j][k]=m[i][j][k-1]+ntes;
     }/*endfor*/
    }/*endfor*/
   }/*endfor*/


/* return pointer to array of pointers to rows */
   return m;  
}/* end routine */


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* cmall_tens4: Careful malloc a rank 4 tensor of type int                  */
/*==========================================================================*/
int ****cmall_itens4(long nrl,long nrh,long ncl,long nch,
		     long ndl,long ndh,long ntl,long nth,char *func_name)

/* allocate a int rank4 tensor with subscript range                       */
/* m[nrl...nrh][ncl...nch][ndl...ndh][ntl...nth] */
{
  long i,j,k;
  long nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1,ntes=nth-ntl+1;
  int ****m;

  if(nrow <= 0 || ncol <= 0 || ndep <= 0 || ntes <= 0) {
    return (int ****)NULL;
  }
/* allocate pointers to rows */
  m=(int ****) fftw_malloc((size_t)((nrow+NR_END)*sizeof(int***)));
  if(!m) 
    {
     PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     PRINTF("allocation failure of row pointers cmall_itens4\n");
     PRINTF("Error occurred in function %s\n",func_name);
     PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     FFLUSH(stdout);
     EXIT(1);
    }
  m += NR_END;
  m -= nrl;

/* allocate rows and set pointers to them */
  m[nrl]=(int ***) fftw_malloc((size_t)((nrow*ncol+NR_END)*sizeof(int**)));
  if(!m[nrl]) 
    {
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("allocation failure of rows in cmall_itens4\n");
      PRINTF("Error occurred in function %s\n",func_name);
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }
  m[nrl] += NR_END;
  m[nrl] -= ncl;

/* allocate columns and set pointers to them */
  m[nrl][ncl]=(int **) fftw_malloc((size_t)((nrow*ncol*ndep+NR_END)
					 *sizeof(int*)));
  if(!m[nrl][ncl]) 
    {
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("allocation failure of columns in cmall_itens4\n");
      PRINTF("Error occurred in function %s\n",func_name);
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }
  m[nrl][ncl] += NR_END;
  m[nrl][ncl] -= ndl;

/* allocate depths and set pointers to them */
  m[nrl][ncl][ndl]=(int *) fftw_malloc((size_t)((nrow*ncol*ndep*ntes+NR_END)
					       *sizeof(int)));
  if(!m[nrl][ncl][ndl]) 
    {
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("allocation failure of depths in cmall_itens4\n");
      PRINTF("Error occurred in function %s\n",func_name);
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }
  m[nrl][ncl][ndl] += NR_END;
  m[nrl][ncl][ndl] -= ntl;


   for(i=ndl+1;i<=ndh;i++) {
     m[nrl][ncl][i]=m[nrl][ncl][i-1]+ntes;
   }/*endfor*/

   for(i=ncl+1;i<=nch;i++) {
     m[nrl][i]=m[nrl][i-1]+ndep;
     m[nrl][i][ndl]=m[nrl][i-1][ndl]+ndep*ntes;
    for(j=ndl+1;j<=ndh;j++) {
      m[nrl][i][j]=m[nrl][i][j-1]+ntes;
    }/*endfor*/
   }/*endfor*/


   for(i=nrl+1;i<=nrh;i++) {
     m[i]=m[i-1]+ncol;
     m[i][ncl]=m[i-1][ncl]+ncol*ndep;
     m[i][ncl][ndl]=m[i-1][ncl][ndl]+ncol*ndep*ntes;
    for(k=ndl+1;k<=ndh;k++) {
     m[i][ncl][k]=m[i][ncl][k-1]+ndep;
    }/*endfor*/
    for(j=ncl+1;j<=nch;j++) {
      m[i][j]=m[i][j-1]+ndep;
      m[i][j][ndl]=m[i][j-1][ndl]+ncol*ndep;
     for(k=ndl+1;k<=ndh;k++) {
      m[i][j][k]=m[i][j-1][k]+ntes;
      m[i][j][k]=m[i][j][k-1]+ntes;
     }/*endfor*/
    }/*endfor*/
   }/*endfor*/


/* return pointer to array of pointers to rows */
   return m;  
}/* end routine */


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Reallocate a double matrix                                               */
/*==========================================================================*/

double **creall_mat(double **m_old,long nrl_old,long nrh_old,
                                   long ncl_old,long nch_old,
		                   long nrl_new,long nrh_new,
                                   long ncl_new,long nch_new,char *func_name)
{
  long i,j, nrow=nrh_new-nrl_new+1,ncol=nch_new-ncl_new+1;
  long nrl_t,nrh_t,ncl_t,nch_t;
  double **m_new;

  if(nrow <= 0 || ncol <= 0) return (double **)NULL;
/* allocate pointers to rows */
  m_new=(double **) fftw_malloc((size_t)((nrow+NR_END)*sizeof(double*)));
  if(!m_new) 
    {
     PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@@\n");
     PRINTF("allocation failure of row pointers in creall_mat\n");
     PRINTF("Error occurred in function %s\n",func_name);
     PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@@\n");
     FFLUSH(stdout);
     EXIT(1);
    }
  m_new += NR_END;
  m_new -= nrl_new;

/* allocate rows and set pointers to them */
  m_new[nrl_new]=(double *) fftw_malloc((size_t)((nrow*ncol+NR_END)
                                            *sizeof(double)));
  if(!m_new[nrl_new]) 
    {
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("allocation failure of rows in cmall_mat\n");
      PRINTF("Error occurred in function %s\n",func_name);
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }
  m_new[nrl_new] += NR_END;
  m_new[nrl_new] -= ncl_new;

  for(i=nrl_new+1;i<=nrh_new;i++) 
    m_new[i]=m_new[i-1]+ncol;

/* Copy old matrix onto new matrix */

   nrl_t = MAX(nrl_new,nrl_old); 
   nrh_t = MIN(nrh_new,nrh_old);
   ncl_t = MAX(ncl_new,ncl_old); 
   nch_t = MIN(nch_new,nch_old);

   for(i=nrl_t;i<=nrh_t;i++){
    for(j=ncl_t;j<=nch_t;j++){
     m_new[i][j] = m_old[i][j];
    }
   }

/* Free old pointer */

   cfree_mat(m_old,nrl_old,nrh_old,ncl_old,nch_old);

/* return pointer to array of pointers to rows */
   return m_new;  
  
}/* end routine */

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Reallocate an integer matrix                                             */
/*==========================================================================*/

int **creall_int_mat(int **m_old,long nrl_old,long nrh_old,
                                 long ncl_old,long nch_old,
		                 long nrl_new,long nrh_new,
                                 long ncl_new,long nch_new,char *func_name)
{
  long i,j, nrow=nrh_new-nrl_new+1,ncol=nch_new-ncl_new+1;
  long nrl_t,nrh_t,ncl_t,nch_t;
  int **m_new;

  if(nrow <= 0 || ncol <= 0) return (int **)NULL;
/* allocate pointers to rows */
  m_new=(int **) fftw_malloc((size_t)((nrow+NR_END)*sizeof(int*)));
  if(!m_new) 
    {
     PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@@\n");
     PRINTF("allocation failure of row pointers in creall_mat\n");
     PRINTF("Error occurred in function %s\n",func_name); 
     PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@@\n");
     FFLUSH(stdout);
     EXIT(1);
    }
  m_new += NR_END;
  m_new -= nrl_new;

/* allocate rows and set pointers to them */
  m_new[nrl_new]=(int *) fftw_malloc((size_t)((nrow*ncol+NR_END)
                                            *sizeof(int)));
  if(!m_new[nrl_new]) 
    {
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("allocation failure of rows in cmall_mat\n");
      PRINTF("Error occurred in function %s\n",func_name);
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }
  m_new[nrl_new] += NR_END;
  m_new[nrl_new] -= ncl_new;

  for(i=nrl_new+1;i<=nrh_new;i++) 
    m_new[i]=m_new[i-1]+ncol;

/* Copy old matrix onto new matrix */

   nrl_t = MAX(nrl_new,nrl_old); 
   nrh_t = MIN(nrh_new,nrh_old);
   ncl_t = MAX(ncl_new,ncl_old); 
   nch_t = MIN(nch_new,nch_old);

   for(i=nrl_t;i<=nrh_t;i++){
    for(j=ncl_t;j<=nch_t;j++){
     m_new[i][j] = m_old[i][j];
    }
   }

/* Free old pointer */

   cfree_int_mat(m_old,nrl_old,nrh_old,ncl_old,nch_old);

/* return pointer to array of pointers to rows */
   return m_new;  
  
}/* end routine */

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* cfree_mat: Careful free a double matrix                                  */
/*==========================================================================*/

void cfree_mat(double **m,long nrl, long nrh, long ncl, long nch)
{
  fftw_free((char *) (m[nrl]+ncl-NR_END));
  fftw_free((char *) (m+nrl-NR_END));
}/* end routine */

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* cfree_int_mat: Careful free a integer matrix                             */
/*==========================================================================*/

void cfree_int_mat(int **m,long nrl, long nrh, long ncl, long nch)
{
  fftw_free((char *) (m[nrl]+ncl-NR_END));
  fftw_free((char *) (m+nrl-NR_END));
}/* end routine */

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* cfree_tens3: Careful free a double rank 3 tensor                         */
/*==========================================================================*/

void cfree_tens3(double ***m,long nrl,long nrh,long ncl,long nch,
                             long ndl,long ndh)
{
  fftw_free((char *) (m[nrl][ncl]+ndl-NR_END));
  fftw_free((char *) (m[nrl]+ncl-NR_END));
  fftw_free((char *) (m+nrl-NR_END));
  
}       
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* cfree_tens3: Careful free a int rank 3 tensor                            */
/*==========================================================================*/

void cfree_itens3(int ***m,long nrl,long nrh,long ncl,long nch,
                           long ndl,long ndh)
{
  fftw_free((char *) (m[nrl][ncl]+ndl-NR_END));
  fftw_free((char *) (m[nrl]+ncl-NR_END));
  fftw_free((char *) (m+nrl-NR_END));
  
}       
/*==========================================================================*/
