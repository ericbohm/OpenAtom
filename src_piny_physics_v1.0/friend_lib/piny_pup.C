#include "standard_include.h"
#include "ckcomplex.h"
#include "../proto_defs/proto_friend_lib_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* 1d integer arrays                                                        */
/*==========================================================================*/
  void pup1d_int(PUP::er &p,int **vec, int nlen){
   if(nlen > 0){
     if (p.isUnpacking()) {
       *vec = (int *)cmalloc(nlen*sizeof(int),"piny_pup")-1;
     }/*endif*/
    p(&(*vec)[1],nlen);
   }/* endif nlen */
  }/*end routine */
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* 1d double arrays                                                        */
/*==========================================================================*/
  void pup1d_dbl(PUP::er &p,double **vec, int nlen){
   if(nlen > 0){
    if (p.isUnpacking()) {
      *vec = (double *)cmalloc(nlen*sizeof(double),"piny_pup")-1;
    }/*endif*/
    p(&(*vec)[1],nlen);
   }/* endif nlen */
  }/*end routine */
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* 1d complex arrays                                                        */
/*==========================================================================*/
  void pup1d_cpl(PUP::er &p,complex **vec, int nlen){
   if(nlen > 0){
    if (p.isUnpacking()) {
      *vec = (complex *)cmalloc(nlen*sizeof(complex),"piny_pup")-1;
    }/*endif*/
    complex *temp = &(*vec)[1];
    p((char *)temp,nlen*sizeof(complex));
   }/* endif nlen */
  }/*end routine */
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* 1d char arrays                                                           */
/*==========================================================================*/
  void pup1d_char(PUP::er &p,char **vec, int nlen){
   if(nlen > 0){
    if (p.isUnpacking()) {
      *vec = (char *)cmalloc(nlen*sizeof(char),"piny_pup");
    }/*endif*/
    p(*vec,nlen);
   }/* endif */
  }/*end routine */
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* 2d integer arrays                                                        */
/*==========================================================================*/
  void pup2d_int(PUP::er &p,int ***vec, int nlen_1,int nlen_2){

   if(nlen_1 > 0 && nlen_2 > 0){
    int ir,ic;
    int index;
    int  psize    = nlen_1*nlen_2;
    int *pscratch = (int *) cmalloc(psize*sizeof(int),"piny_pup")-1;
    if (p.isUnpacking()) {
      *vec = cmall_int_mat(1,nlen_1,1,nlen_2,"piny_pup");
    }else{
    // packing 2d array into 1d array to be pupped
     index = 1;
     for(ir=1; ir <= nlen_1; ir++){
     for(ic=1; ic <= nlen_2; ic++){
      pscratch[index] = (*vec)[ir][ic];
      ++index;
     }} 
    }// endif
    p(&(pscratch[1]),psize);
    // unpacking 1d array into 2d array
    if (p.isUnpacking()) {
     index = 1;
     for(ir=1; ir <= nlen_1; ir++){
     for(ic=1; ic <= nlen_2; ic++){
      (*vec)[ir][ic] = pscratch[index];
      ++index;
    }} 
    }//endif
     cfree(&(pscratch[1]),"piny_pup");
   }// endif nlen 

  }// end routine 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* 2d double arrays                                                        */
/*==========================================================================*/
  void pup2d_dbl(PUP::er &p,double ***vec, int nlen_1,int nlen_2, char*callerID){
   if(nlen_1 > 0 && nlen_2 > 0){
    int ir,ic;
    int index;
    int  psize       = nlen_1*nlen_2;
    double *pscratch = (double *) cmalloc(psize*sizeof(double),"piny_pup")-1;
    double **vtmp = *vec;
    if (p.isUnpacking()) {
      *vec = cmall_mat(1,nlen_1,1,nlen_2,"piny_pup");
      vtmp = *vec;
    }else{
    // packing 2d array into 1d array to be pupped
     index = 1;
     for(ir=1; ir <= nlen_1;  ir++){
     for(ic=1; ic <= nlen_2; ic++){
      pscratch[index] = vtmp[ir][ic];
      ++index;
     }} 
    }//endif
    p(&(pscratch[1]),psize);
    // unpacking 1d array into 2d array
    if (p.isUnpacking()) {
     index = 1;
     for(ir=1; ir <= nlen_1; ir++){
     for(ic=1; ic <= nlen_2; ic++){
      vtmp[ir][ic] = pscratch[index];
      ++index;
     }} 
    }//endif
    cfree(&(pscratch[1]),"piny_pup");
   }// endif nlen 

  }// end routine 
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* 2d char arrays                                                           */
/*==========================================================================*/
  void pup1d_name(PUP::er &p,NAME **vec, int nlen){
   int i;
   char *scratch;
   int strlen = MAXWORD;

   if(nlen > 0){
     scratch = (char *)cmalloc(MAXWORD*sizeof(char),"piny_pup");
     if (p.isUnpacking()) {
       *vec     = (NAME *)cmalloc(nlen*sizeof(NAME),"piny_pup")-1;
     }/*endif*/
     for(i=1;i<=nlen;i++){
       if (p.isPacking()) {strcpy(scratch,(*vec)[i]);}
       p(scratch,MAXWORD);
       if (p.isUnpacking()) {strcpy((*vec)[i],scratch);}
     }/*endfor*/
     cfree(scratch,"piny_pup");
   }/* endif */

  }/*end routine */
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* 3d double arrays                                                        */
/*==========================================================================*/
  void pup3d_dbl(PUP::er &p,double ****vec, int nlen_1,int nlen_2,int nlen_3){

   if( (nlen_1 > 0) && (nlen_2 > 0) && (nlen_3 > 0) ){
    int ir,ic,id;
    int index;
    int psize = nlen_1*nlen_2*nlen_3;
    double *pscratch = (double *) cmalloc(psize*sizeof(double),"piny_pup")-1;
    if (p.isUnpacking()) {
     (*vec) = (double ***) cmalloc(nlen_1*sizeof(double **),"piny_pup")-1;
      for(ir=1; ir<= nlen_1; ir++){
       (*vec)[ir] = (double **) cmalloc(nlen_2*sizeof(double *),"piny_pup")-1;
       for(ic=1; ic<=nlen_2; ic++){
         (*vec)[ir][ic] = (double *) cmalloc(nlen_3*sizeof(double ),"piny_pup")-1;
       }//endfor
      }//endfor
    }else{
    // packing 3d array into 1d array to be pupped
     index = 1;
     for(ir=1; ir <= nlen_1; ir++){
     for(ic=1; ic <= nlen_2; ic++){
     for(id=1; id <= nlen_3; id++){
      pscratch[index] = (*vec)[ir][ic][id];
      ++index;
     }}} 
    }//endif
    p(&(pscratch[1]),psize);
    // unpacking 1d array into 2d array
    if (p.isUnpacking()) {
     index = 1;
     for(ir=1; ir <= nlen_1; ir++){
     for(ic=1; ic <= nlen_2; ic++){
     for(id=1; id <= nlen_3; id++){
      (*vec)[ir][ic][id] = pscratch[index];
      ++index;
     }}}
    }//endif
    cfree(&(pscratch[1]),"piny_pup");
   }// endif nlen 

  }// end routine 
//==========================================================================
