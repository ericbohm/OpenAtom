//==========================================================================
/** \file para_grp_parse.h
 *
 */
//==========================================================================

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

class ParaGrpParse{

 //---------------------------------------------------------------------------
 public:

 //---------------------------------------------------------------------------
 //con-destruct:
   ParaGrpParse(){};
  ~ParaGrpParse(){};

//---------------------------------------------------------------------------
// functions

static void get_equal_decomp_prms(int , int , int , int *, int *, int *);
static void get_equal_to_plane_info(int, int ,int *,int *, int **, 
                                    int **, int **, int **);
static void get_plane_decomp_prms(int, int ,int ,int *,int *,int *,int *);
static void get_plane_to_equal_info(int, int ,int ,int , int , int *, int **, 
                                    int **, int **, int **);


//---------------------------------------------------------------------------
   }; //ParaGrpParse
//==========================================================================

#endif
