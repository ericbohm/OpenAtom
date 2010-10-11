#include "matrix2file.h"
#include "charm++.h"//< Just for CkAbort!!

void dumpMatrixDouble(const char *infilename, double *matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric)
{
  char fmt[1000];
  char filename[1000];
  strncpy(fmt,infilename,999);
  strncat(fmt,"_%d_%d_%d_%d_%d.out",999);
  sprintf(filename,fmt, w, x, y, z, symmetric);
  FILE *loutfile = fopen(filename, "w");
#ifdef PAIRCALC_TEST_DUMP
  fprintf(loutfile,"%d\n",ydim);
#endif
  for(int i=0;i<xdim;i++)
    for(int j=0;j<ydim;j++)
      fprintf(loutfile,"%d %d %.12g\n",x+i,y+j,matrix[i*ydim+j]);
  fclose(loutfile);
}



void loadMatrixDouble(const char *infilename, double *matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric)
{
  char fmt[1000];
  char filename[1000];
  strncpy(fmt,infilename,999);
  strncat(fmt,"_%d_%d_%d_%d_%d.out",999);
  sprintf(filename,fmt, w, x, y, z, symmetric);
  FILE *loutfile = fopen(filename, "r");
  if(loutfile!=NULL)
    {
      int junk1,junk2;
      for(int i=0;i<xdim;i++)
	for(int j=0;j<ydim;j++)
	  CkAssert(fscanf(loutfile,"%d %d %lf\n",&junk1,&junk2,&(matrix[i*ydim+j])));
      fclose(loutfile);
    }
  else
    {
      CkAbort(filename);
    }
}

//! NOTE: this uses the evil piny convention
void dumpMatrix2DDouble(const char *infilename, double **matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric)
{
  char fmt[1000];
  char filename[1000];
  strncpy(fmt,infilename,999);
  strncat(fmt,"_%d_%d_%d_%d_%d.out",999);
  sprintf(filename,fmt, w, x, y, z, symmetric);
  FILE *loutfile = fopen(filename, "w");
  for(int i=0;i<xdim;i++)
    for(int j=1;j<=ydim;j++)
      fprintf(loutfile,"%d %d %.12g\n",i,j,matrix[i][j]);
  fclose(loutfile);
}
void loadMatrix2DDouble(const char *infilename, double **matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric)
{
  char fmt[1000];
  char filename[1000];
  strncpy(fmt,infilename,999);
  strncat(fmt,"_%d_%d_%d_%d_%d.out",999);
  sprintf(filename,fmt, w, x, y, z, symmetric);
  FILE *loutfile = fopen(filename, "r");
  if(loutfile!=NULL)
    {
      int junk1,junk2;
      for(int i=0;i<xdim;i++)
	for(int j=1;j<=ydim;j++)
	  CkAssert(fscanf(loutfile,"%d %d %lf\n",&junk1,&junk2,&(matrix[i][j])));
      fclose(loutfile);
    }
  else
    {
      CkAbort(filename);
    }
}


//! NOTE: this uses the evil piny convention
void dumpMatrix2DInt(const char *infilename, int **matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric)
{
  char fmt[1000];
  char filename[1000];
  strncpy(fmt,infilename,999);
  strncat(fmt,"_%d_%d_%d_%d_%d.out",999);
  sprintf(filename,fmt, w, x, y, z, symmetric);
  FILE *loutfile = fopen(filename, "w");
  for(int i=0;i<xdim;i++)
    for(int j=1;j<=ydim;j++)
      fprintf(loutfile,"%d %d %d\n",i,j,matrix[i][j]);
  fclose(loutfile);
}
void loadMatrix2DInt(const char *infilename, int **matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric)
{
  char fmt[1000];
  char filename[1000];
  strncpy(fmt,infilename,999);
  strncat(fmt,"_%d_%d_%d_%d_%d.out",999);
  sprintf(filename,fmt, w, x, y, z, symmetric);
  FILE *loutfile = fopen(filename, "r");
  if(loutfile!=NULL)
    {
      int junk1,junk2;
      for(int i=0;i<xdim;i++)
	for(int j=1;j<=ydim;j++)
	  CkAssert(fscanf(loutfile,"%d %d %d\n",&junk1,&junk2,&(matrix[i][j])));
      fclose(loutfile);
    }
  else
    {
      CkAbort(filename);
    }
}

