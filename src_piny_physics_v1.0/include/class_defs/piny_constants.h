//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                                                                          
//                         PI_MD:                                           
//             The future of simulation technology                          
//             ------------------------------------                         
//                   Module: defines.h                                      
//                                                                          
// File to define simple constants and macros                               
//                                                                          
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================



//==========================================================================
// HARD WIRE SOME ARRAY SIZES   

#define MAX_VALENCE      6
#define MAX_VALENCE1     7
#define MAX_BOND_SITE    4
#define MAX_BOND_SITE1   5
#define NCOEF_GHOST_MAX  5
#define NCOEF_GHOST_MAX1 6
#define MAXWORD          80
#define MAXLINE          100
#define NR_END           0
#define NMEM_MIN         100

//==========================================================================
// DEFINE DATA TYPES   

#ifdef T3E_SCILIB
#define list_int short
#else
#define list_int int
#endif

typedef char NAME[MAXWORD];    // Chr: a name; Lth: MAXWORD           
typedef char LINE[MAXLINE];    // Chr: a line; Lth: MAXLINE           

//==========================================================================
// DEFINE CONSTANTS:  

#ifndef M_PI
#define M_PI       3.14159265358979323846
#endif
#define BOHR       0.529177  // 0.529177724924 
#define BOLTZ      315777.0  // 315773.218 
#define KCAL       627.50921
#define EV         27.211396
#define PROT_MASS  1822.8885
#define PCONV      3.3989242e-09
#define STENS_CONV 6.423021e-07
#define TIME_CONV  0.0241888
#define CP_EMAGIC  0.00050

//==========================================================================
// DEFINE CHARM++ FUNCTIONS

#ifdef CHARM_ON
#define PRINTF CkPrintf
#ifdef PUP_PRINTING_ON
#define PUP_PRINTF CkPrintf
#else
#define PUP_PRINTF
#endif
#define SCANF  CkScanf
#define FFLUSH fflush
#define EXIT(N) {CkExit();}
#endif

#ifdef CHARM_OFF
#define PUP_PRINTF printf
#define PRINTF printf
#define FFLUSH fflush
#define SCANF  scanf
#define EXIT(N) {exit(N);}
#endif

//==========================================================================
// DEFINE FUNCTION MACROS:  


#define PRINT_LINE_STAR {PRINTF("==============================================================================\n");}
#define PRINT_LINE_DASH {PRINTF("------------------------------------------------------------------------------\n");}

#define PUP_PRINT_LINE_STAR {PUP_PRINTF("==============================================================================\n");}
#define PUP_PRINT_LINE_DASH {PUP_PRINTF("------------------------------------------------------------------------------\n");}

#define SKIP_LINE(A) {int ch_; do{ch_=fgetc(A);}while(ch_!='\n'&&ch_!=EOF);}

#ifndef NINT
#ifdef SIMP_NINT
#define NINT(X) ( (int) ((X)>=0.0 ? ((X)+0.5):((X)-0.5)) )
#else
#define MAGIC 6755399441055744.0
#define NINT(X) (((X)+MAGIC)-MAGIC)
#endif
#endif

#ifndef MAX
#define MAX(A,B) (((A)>(B))?(A):(B))
#endif
#ifndef MIN
#define MIN(A,B) (((A)<(B))?(A):(B))
#endif
#define MAX3(A,B,C) (MAX(MAX((A),(B)),(C)))
#define MIN3(A,B,C) (MIN(MIN((A),(B)),(C)))

#define STEP(A) ((A)>0? 1.0:0.0)

#define DEBUG_READ_INT {int iii;PRINTF("Enter an int : ");SCANF("%d",&iii);}
#define DEBUG_WRITE_DBLE(S,X) {PRINTF("%s %g\n",S,X);}
#define DEBUG_WRITE_INT(S,I) {PRINTF("%s %d\n",S,I);}

//-------------------------------------------------------------------------
// Reversible 6th order pade approx to exp(x) : pade(x)*pade(-x)=1         
//           8th order pade would give 10^{-9} -0.5<x<0.5                  
//           |exp(x)-pade(x)| < 10^{-7}  : -0.483<x<0.561                  
//           |exp(x)-pade(x)| < 10^{-8}  : -0.355<x<0.394                  
//           |exp(x)-pade(x)| < 10^{-9}  : -0.259<x<0.279                  
//           |exp(x)-pade(x)| < 10^{-10} : -0.188<x<0.198                  
//           |exp(x)-pade(x)| < 10^{-11} : -0.136<x<0.142                  
//           |exp(x)-pade(x)| < 10^{-12} : -0.987<x<0.101                  
//           |exp(x)-pade(x)| < 10^{-13} : -0.071<x<0.072                  
//           |exp(x)-pade(x)| < 10^{-14} : -0.051<x<0.052                  
//           |exp(x)-pade(x)| < 10^{-15} : -0.035<x<0.036                  
#define C3PA 0.008333333333333333333333333333333
#define FST_EXP(X)((((C3PA*X+0.1)*X+0.5)*X+1.0)/(((-C3PA*X+0.1)*X-0.5)*X+1.0))
//-------------------------------------------------------------------------


//==========================================================================
// DEFINE FORTRAN PROTOCAL:  

#ifdef HP_VECLIB
#ifndef FORTRANUNDERSCORE_OFF
#define FORTRANUNDERSCORE_OFF
#endif
#endif

#ifdef SGI_COMPLIB
#ifndef FORTRANUNDERSCORE
#define FORTRANUNDERSCORE
#endif
#endif

#ifdef IBM_ESSL
#ifndef FORTRANUNDERSCORE_OFF
#define FORTRANUNDERSCORE_OFF
#endif
#endif

#ifdef IBM_NOESSL
#ifndef FORTRANUNDERSCORE_OFF
#define FORTRANUNDERSCORE_OFF
#endif
#endif

#ifdef DEC_ALPHA
#ifndef FORTRANUNDERSCORE
#define FORTRANUNDERSCORE
#endif
#endif

#ifdef SUN_COMPLIB
#ifndef FORTRANUNDERSCORE
#define FORTRANUNDERSCORE
#endif
#endif

#ifdef DARWIN_MAC
#ifndef FORTRANUNDERSCORE_DBLE
#define FORTRANUNDERSCORE_DBLE
#endif
#endif


//==========================================================================



