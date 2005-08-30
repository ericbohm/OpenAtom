#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void lst_sort_clean(int , int *, int *);

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
