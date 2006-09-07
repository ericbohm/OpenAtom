//bluegene check
//Authored by Sameer Kumar 2006/8/12
#ifdef CMK_VERSION_BLUEGENE
#include <builtins.h>

//performs the arithmatic a[i] += b[i]
void fastAdd (double *a, double *b, int nelem) {

  //Check alignments

#pragma disjoint (*a, *b)
#pragma disjoint (*b, *a)

  register int i = 0;

  if (((unsigned long)a & 0x0f) == 0 && ((unsigned long)a & 0x0f) == 0
      && nelem > 32) {
    __alignx (16, a);
    __alignx (16, b);

#pragma unroll(16)
    for (i = 0; i < nelem; i ++) {
      a[i] += b[i];
    }
  }
  else {
    for (i = 0; i < nelem; i ++)
      a[i] += b[i];
  }
}
#else 
// this is just to simplify our compilation life
// users of this code should probably do this inline
void fastAdd (double *a, double *b, int nelem) {
#define _UNROLLING_OFF_
#ifdef _UNROLLING_OFF_
    register int i=0;
    for ( ; i < nelem; i ++){a[i] += b[i];}
#else
   int nrem  = (nelem % 5);
   int istrt = (nelem-nrem);
   for(int i=0;i<istrt;i+=5){
     a[i]   += b[i];
     a[i+1] += b[i+1];
     a[i+2] += b[i+2];
     a[i+3] += b[i+3];
     a[i+4] += b[i+4];
   }//endfor
   for(int i=istrt;i<nelem;i++){a[i] += b[i];}
#endif
}
#endif
