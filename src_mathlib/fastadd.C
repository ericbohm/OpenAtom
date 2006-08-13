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
    register int i=0;
    for ( ; i < nelem; i ++)
      a[i] += b[i];
}
#endif
