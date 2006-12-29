/*
** orthoHelper.h
** 
**
**
** Made by Eric Bohm
** Login   <bohm@localhost.localdomain>
** 
** Started on  Mon Sep 11 12:12:53 2006 Eric Bohm
** Last update Mon Sep 11 12:12:53 2006 Eric Bohm
**
** orthoHelper is just a chare array to host the 2nd multiply of the
** S->T method.  It should be mapped adjacent to ortho such that
** ortho(x,y) and orthoHelper(x,y) are on different processors which
** are only one away.  This assumes that ortho itself is mapped with a
** stride greater than 1. Orthohelper should only be used if numProcs
** >= 2 x numOrtho chares, otherwise we expect the parallelism gains
** would be lost to communication overhead.
**
** OrthoHelper is triggered by being given the result of multiply 1 by
** Ortho.  It then computes multiply 2 and returns it to the ortho
** chares who sent it the multiply 1 input.  This is a simple point to
** point communication ortho(x,y) <-> orthoHelper(x,y).   Ortho then
** proceeds as normal through the S->T process.
*/

#ifndef   	ORTHOHELPER_H_
# define   	ORTHOHELPER_H_
#include "CLA_Matrix.h"

extern MapType2 OrthoHelperImaptable;
extern CkHashtableT <intdual, int> OrthoHelpermaptable;

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class OrthoHelperMsg : public CMessage_OrthoHelperMsg {
public:
  double *A;
  double *B;
  int size; 
  double factorA;
  double factorB;
  double factorC;
  void init(int _size, double *_A, double *_B, double _factorA, double _factorB, double _factorC)
    {
      size=_size;
      memcpy(A,_A,size*sizeof(double));
      memcpy(B,_B,size*sizeof(double));
      factorA=_factorA;
      factorB=_factorB;
      factorC=_factorC;
    }

#ifdef CMK_VERSION_BLUEGENE
  // if we use our own allocator we can get 16 byte alignment
  // to please BGL
#define ALIGN16(x)        (int)((~15)&((x)+15))
 static  void *alloc(int msgnum, size_t sz, int *sizes, int pb) {
    int offsets[3];
    offsets[0] = ALIGN16(sz);
    if(sizes==0)
      offsets[1] = offsets[0];
    else
      offsets[1] = offsets[0] + ALIGN16(sizeof(double)*sizes[0]);
    if(sizes==0)
      offsets[2] = offsets[0];
    else
      offsets[2] = offsets[1] + ALIGN16(sizeof(double)*sizes[1]);
    OrthoHelperMsg *newmsg = (OrthoHelperMsg *) CkAllocMsg(msgnum, offsets[2], pb);
    newmsg->A = (double *) ((char *)newmsg + offsets[0]);
    newmsg->B = (double *) ((char *)newmsg + offsets[1]);
    return (void *) newmsg;
  }

#endif
 friend class CMessage_OrthoHelperMsg;
};
//============================================================================


class OrthoHelper : public CBase_OrthoHelper 
{
 public:
  OrthoHelper(){}
  OrthoHelper(int _m, int _n, CLA_Matrix_interface matA2,
	      CLA_Matrix_interface matB2, CLA_Matrix_interface matC2):
    m(_m), n(_n), matA (matA2), matB(matB2), matC(matC2)
    {
      C= new double[m*n];
      A=NULL;
      B=NULL;
      trigger=NULL;
    }
  OrthoHelper(CkMigrateMessage *m){}
  ~OrthoHelper(){
    delete [] C;
  }

  void recvAB(OrthoHelperMsg *msg)
    {  
      trigger=msg;
      A=trigger->A;  // no sense in making extra copies
      B=trigger->B;
      matA.multiply(trigger->factorA, 0, A, OrthoHelper::step_cb, (void*) this,
		     thisIndex.x, thisIndex.y);
      CmiNetworkProgress();
      matB.multiply(trigger->factorB, 0, B, OrthoHelper::step_cb, (void*) this,
		     thisIndex.x, thisIndex.y);
      CmiNetworkProgress();
      matC.multiply(trigger->factorC, 0, C, OrthoHelper::step_cb, (void*) this,
		     thisIndex.x, thisIndex.y);

    }

  void sendMatrix();

  virtual void pup(PUP::er &p){
//    CBase_Ortho::pup(p);
    ArrayElement2D::pup(p);
    p | m;
    p | n;
    if(p.isUnpacking()){
      C = new double[m * n];
    }
    p(C, m * n);
  }
  static inline void step_cb(void *obj){
    ((OrthoHelper *) obj)->sendMatrix();
  }
  int m;
  int n;
  double *A, *B, *C;
  OrthoHelperMsg *trigger;
  CLA_Matrix_interface matA, matB, matC;
};

class OrthoHelperMap : public CkArrayMapTable2 {
  public:
    OrthoHelperMap()
    {
#ifdef USE_INT_MAP
      maptable= &OrthoHelperImaptable;
#else
      maptable= &OrthoHelpermaptable;
#endif
    }

    ~OrthoHelperMap() { }
    
    void pup(PUP::er &p)
    {
      CkArrayMapTable2::pup(p);
#ifdef USE_INT_MAP
      maptable= &OrthoHelperImaptable;
#else
      maptable= &OrthoHelpermaptable;
#endif
    }
    
    inline int procNum(int, const CkArrayIndex &iIndex)
    {
      int *index=(int *) iIndex.data();
#ifdef USE_INT_MAP
      return(maptable->get(index[0],index[1]));
#else
      return(maptable->get(intdual(index[0],index[1])));
#endif
    }
};

#endif 	    /* !ORTHOHELPER_H_ */
