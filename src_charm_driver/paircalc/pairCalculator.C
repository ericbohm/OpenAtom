#include "pairCalculator.h"
#include <algorithm>


void dumpMatrixDouble(const char *infilename, double *matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric)
{
  char fmt[1000];
  char filename[1000];
  strncpy(fmt,infilename,999);
  strncat(fmt,"_%d_%d_%d_%d_%d.out",999);
  sprintf(filename,fmt, w, x, y, z, symmetric);
  FILE *loutfile = fopen(filename, "w");
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
	  fscanf(loutfile,"%d %d %lf\n",&junk1,&junk2,&(matrix[i*ydim+j]));
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
	  fscanf(loutfile,"%d %d %lf\n",&junk1,&junk2,&(matrix[i][j]));
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
	  fscanf(loutfile,"%d %d %d\n",&junk1,&junk2,&(matrix[i][j]));
      fclose(loutfile);
    }
  else
    {
      CkAbort(filename);
    }
}



#ifdef ROTATE_LIST
bool reorder_elem_list(CkArrayIndexMax *elems, int numelems, int newstart)
{
  if(newstart>numelems)
    return(false);
  CkArrayIndexMax  swap;
  int where=newstart;
  for(int i=0;i<newstart;i++,where++)
    {
      if(where>=numelems)
	where=0;
      swap=elems[i];
      elems[i]=elems[where];
      elems[where]=swap;
    }
  return(true);
}

#else

bool reorder_elem_list(CkArrayIndexMax *elems, int numelems, int newstart)
{
  //  CkPrintf("reordering list of %d elems\n", numelems);
  std::random_shuffle(elems,elems+numelems);
  return(true);
}

bool reorder_elem_list_4D(CkArrayIndex4D *elems, int numelems, int newstart)
{
  //  CkPrintf("reordering list of %d elems\n", numelems);
  std::random_shuffle(elems,elems+numelems);
  return(true);
}

bool reorder_elem_list_max(CkArrayIndexMax *elems, int numelems, int newstart)
{
  //  CkPrintf("reordering list of %d elems\n", numelems);
  std::random_shuffle(elems,elems+numelems);
  return(true);
}
#endif

