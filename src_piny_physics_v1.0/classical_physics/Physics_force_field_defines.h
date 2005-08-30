//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                     Lean_MD-PINY_Phyics
//
// Purpose :                                                                         
//      Define simple constants and in-line functions for the physics
//
// Authors :
//      Parallel Programming Laboratory, CS, UIUC
//      Tuckerman Group, Chemistry, New York University
//                                                                          
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================


//==========================================================================
// DEFINE CONSTANTS:  

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define BOLTZ 315777.0
// #define BOLTZ (315773.218) 
#define KCAL 627.50921
#define EV 27.211396
#define PROT_MASS 1822.8885
#define PCONV 3.3989242e-09
#define STENS_CONV 6.423021e-07
#define TIME_CONV 0.0241888

//==========================================================================


//==========================================================================
// DEFINE INLINE FUNCTIONS:  

#ifndef NINT
#ifdef SIMP_NINT
#define NINT(X) ( (int) ((X)>=0.0 ? ((X)+0.5):((X)-0.5)) )
#else
#define MAGIC 6755399441055744.0
#define NINT(X) (((X)+MAGIC)-MAGIC)
#endif
#endif

#define MAX(A,B) (((A)>(B))?(A):(B))
#define MIN(A,B) (((A)<(B))?(A):(B))
#define MAX3(A,B,C) (MAX(MAX((A),(B)),(C)))

//==========================================================================







