#ifndef _PARA_GRP_PARSE_
#define _PARA_GRP_PARSE_

#ifndef MAX
#define MAX(A,B) (((A)>(B))?(A):(B))
#endif

#ifndef MIN
#define MIN(A,B) (((A)<(B))?(A):(B))
#endif

#ifndef PRINTF
#define PRINTF printf
//#define PRINTF CkPrintf
#endif

#ifndef EXIT
#define EXIT(N) {exit(N);}
//#define EXIT(N) {CkExit();}
#endif

#include "ckcomplex.h"

class ParaGrpParse{

 //---------------------------------------------------------------------------
 public:

 //---------------------------------------------------------------------------
 //con-destruct:
   ParaGrpParse(){};
  ~ParaGrpParse(){};

//---------------------------------------------------------------------------
// functions
static void get_chareG_line_prms(int , int ,int ,int *,int *,int *,int *,int *);

static void flip_data_set(int , int *, int *, int *, int *,complex *);

//---------------------------------------------------------------------------
   }; //ParaGrpParse
//==========================================================================

#endif
