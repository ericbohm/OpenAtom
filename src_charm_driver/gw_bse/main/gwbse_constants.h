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


//==========================================================================
// DEFINE DATA TYPES   

#define MAXWORD          80
typedef char NAME[MAXWORD];    // Chr: a name; Lth: MAXWORD           

#define MAXLINE          100
typedef char LINE[MAXLINE];    // Chr: a line; Lth: MAXLINE           

//==========================================================================
// DEFINE CONSTANTS:  

#ifndef M_PI
#define M_PI       3.14159265358979323846
#endif
#define BOHR       0.529177  // 0.529177724924 
#define BOLTZ      315777.0  // 315773.218 
#define EV         27.211396

//==========================================================================
// DEFINE CHARM++ FUNCTIONS

#ifdef CHARM_ON
#define PRINTF CkPrintf
#define CKMYPE CkMyPe
#ifdef _DEBUG_PUP_
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
#define CKMYPE() (0) 
#endif
//==========================================================================
// DEFINE FUNCTION MACROS:  


#define PRINT_LINE_STAR {PRINTF("==============================================================================\n");}
#define PRINT_LINE_DASH {PRINTF("------------------------------------------------------------------------------\n");}

#define PUP_PRINT_LINE_STAR {PUP_PRINTF("==============================================================================\n");}
#define PUP_PRINT_LINE_DASH {PUP_PRINTF("------------------------------------------------------------------------------\n");}

#define SKIP_LINE(A) {int ch_; do{ch_=fgetc(A);}while(ch_!='\n'&&ch_!=EOF);}

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



