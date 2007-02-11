#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void lst_sort_clean(int , int *, int *);
void sort_commence(int , int *,int *);
void sort_me(int , int *);
void sort_commence_piny(int , int *,int *);

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void lst_sort_clean(int n, int *n_out, int *index)

/*=======================================================================*/
/*            Begin subprogram:                                          */
{/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */

  int m,ir,i,rindex,j;

  if(n<=1){
    (*n_out) = n;
    return;
  }//endif

/*=======================================================================*/
/* I) Sort the list                        */

  m  = n/2+1;
  ir = n;

  for(;;){

/*---------------------------------------------------------------------*/
/*  A)hire rindex */
    if(m>1){ 
      m--;
      rindex = index[(m-1)];
/*--------------------------------------------------------------------*/
/*  B)retire/promote index[1] */
    }else{
      rindex = index[(ir-1)];
      index[(ir-1)]=index[0];
      ir--;
      if(ir==1){
	index[0]=rindex;
	break;
      }/*endif*/
    }/*endif*/
/*---------------------------------------------------------------------*/
/*  C)put rindex in appropriate slot */
    i=m;
    j=2*m;
    while(j<=ir){
      /*    a)compare to rindex to underling */
      if((j<ir) && (index[(j-1)]< index[j])) j++;
      /*    b)demote */
      if(rindex<index[(j-1)]){
	index[(i-1)]=index[(j-1)];
	i=j;
	j=2*j;
      }else{
	/*    c)if no demotations exit while */
	j=ir+1;
      }/*endif*/
    } /*endwhile*/
    /*    d)slot rindex */
    index[(i-1)] = rindex;
  }/*endfor*/

/*==========================================================================*/
/* Clean out repeaters */

  j = 0;
  for(i=1;i<n;i++){
    if(index[i]!=index[j]){
      j++; index[j]=index[i];
    }//endif
  }//endfor

  *n_out=j+1;
/*-----------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  Sort one vector commensorate with another  */
/*==========================================================================*/

void sort_commence(int n, int *index_in,int *jndex_in)

/*=======================================================================*/
/*            Begin subprogram:                                          */
{/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */

  int m,ir,i,j,rindex,rjndex;
  int *index = index_in-1;
  int *jndex = jndex_in-1;

/*=======================================================================*/
/* I) Setup                        */

  m  = n/2+1;
  ir = n;
  if(n==1){return;}

/*=======================================================================*/
/* II) Sort array index keeping jndex commensurrate */

  for(;;){

/*---------------------------------------------------------------------*/
/*  A)hire rindex */
    if(m>1){ 
      m--;
      rindex = index[m];
      rjndex = jndex[m];
/*--------------------------------------------------------------------*/
/*  B)retire/promote index[1] */
    }else{
      rindex = index[ir];
      rjndex = jndex[ir];
      index[ir]=index[1];
      jndex[ir]=jndex[1];
      ir--;
      if(ir==1){
       index[1]=rindex;
       jndex[1]=rjndex;
       break;
      }/*endif*/
    }/*endif*/
/*---------------------------------------------------------------------*/
/*  C)put rindex in appropriate slot */
    i=m;
    j=2*m;
    while(j<=ir){
      /*    a)compare to rindex to underling */
      if((j<ir) && (index[j]< index[(j+1)])) j++;
      /*    b)demote */
      if(rindex<index[j]){
       index[i]=index[j];
       jndex[i]=jndex[j];
       i=j;
       j=2*j;
      }else{
       /*    c)if no demotations exit while */
       j=ir+1;
      }/*endif*/
    } /*endwhile*/
    /*    d)slot rindex */
    index[i] = rindex;
    jndex[i] = rjndex;
  }/*endfor*/

/*-----------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  Sort one vector */
/*==========================================================================*/

void sort_me(int n, int *index_in)

/*=======================================================================*/
/*            Begin subprogram:                                          */
{/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */

  int m,ir,i,j,rindex,rjndex;
  int *index = index_in-1;

/*=======================================================================*/
/* I) Setup                        */

  m  = n/2+1;
  ir = n;
  if(n==1){return;}

/*=======================================================================*/
/* II) Sort array index */

  for(;;){

/*---------------------------------------------------------------------*/
/*  A)hire rindex */
    if(m>1){ 
      m--;
      rindex = index[m];
/*--------------------------------------------------------------------*/
/*  B)retire/promote index[1] */
    }else{
      rindex = index[ir];
      index[ir]=index[1];
      ir--;
      if(ir==1){
       index[1]=rindex;
       break;
      }/*endif*/
    }/*endif*/
/*---------------------------------------------------------------------*/
/*  C)put rindex in appropriate slot */
    i=m;
    j=2*m;
    while(j<=ir){
      /*    a)compare to rindex to underling */
      if((j<ir) && (index[j]< index[(j+1)])) j++;
      /*    b)demote */
      if(rindex<index[j]){
       index[i]=index[j];
       i=j;
       j=2*j;
      }else{
       /*    c)if no demotations exit while */
       j=ir+1;
      }/*endif*/
    } /*endwhile*/
    /*    d)slot rindex */
    index[i] = rindex;
  }/*endfor*/

/*-----------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  Sort one vector commensorate with another  */
/*==========================================================================*/

void sort_commence_piny(int n, int *index,int *jndex)

/*=======================================================================*/
/*            Begin subprogram:                                          */
{/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */

  int m,ir,i,j,rindex,rjndex;

/*=======================================================================*/
/* I) Setup                        */

  m  = n/2+1;
  ir = n;
  if(n==1){return;}

/*=======================================================================*/
/* II) Sort array index keeping jndex commensurrate */

  for(;;){

/*---------------------------------------------------------------------*/
/*  A)hire rindex */
    if(m>1){ 
      m--;
      rindex = index[m];
      rjndex = jndex[m];
/*--------------------------------------------------------------------*/
/*  B)retire/promote index[1] */
    }else{
      rindex = index[ir];
      rjndex = jndex[ir];
      index[ir]=index[1];
      jndex[ir]=jndex[1];
      ir--;
      if(ir==1){
       index[1]=rindex;
       jndex[1]=rjndex;
       break;
      }/*endif*/
    }/*endif*/
/*---------------------------------------------------------------------*/
/*  C)put rindex in appropriate slot */
    i=m;
    j=2*m;
    while(j<=ir){
      /*    a)compare to rindex to underling */
      if((j<ir) && (index[j]< index[(j+1)])) j++;
      /*    b)demote */
      if(rindex<index[j]){
       index[i]=index[j];
       jndex[i]=jndex[j];
       i=j;
       j=2*j;
      }else{
       /*    c)if no demotations exit while */
       j=ir+1;
      }/*endif*/
    } /*endwhile*/
    /*    d)slot rindex */
    index[i] = rindex;
    jndex[i] = rjndex;
  }/*endfor*/

/*-----------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/

