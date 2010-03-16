/*****************************************************************************
 * $Source$
 * $Author$
 * $Date$
 * $Revision$
 *****************************************************************************/

/** \file orthoHelper.h
 *
 *  Made by Eric Bohm
 *  Login   <bohm@localhost.localdomain>
 *  
 *  Started on  Mon Sep 11 12:12:53 2006 Eric Bohm
 *  Last update Mon Sep 11 12:12:53 2006 Eric Bohm
 *
 *  orthoHelper is just a chare array to host the 2nd multiply of the
 *  S->T method.  It should be mapped adjacent to ortho such that
 *  ortho(x,y) and orthoHelper(x,y) are on different processors which
 *  are only one away.  This assumes that ortho itself is mapped with a
 *  stride greater than 1. Orthohelper should only be used if numProcs
 *  >= 2 x numOrtho chares, otherwise we expect the parallelism gains
 *  would be lost to communication overhead.
 *
 *  OrthoHelper is triggered by being given the result of multiply 1 by
 *  Ortho.  It then computes multiply 2 and returns it to the ortho
 *  chares who sent it the multiply 1 input.  This is a simple point to
 *  point communication ortho(x,y) <-> orthoHelper(x,y).   Ortho then
 *  proceeds as normal through the S->T process.
 */

#include "ortho.decl.h"

#ifndef   	ORTHOHELPER_H_
# define   	ORTHOHELPER_H_

#include "main/CLA_Matrix.h"

extern CkVec <MapType2> OrthoHelperImaptable;
extern CkHashtableT <intdual, int> OrthoHelpermaptable;
extern bool fakeTorus;

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
  void init(int insize, double * inA, double * inB, double infactorA, double infactorB, double infactorC)
    {
      size=insize;
      memcpy(A,inA,size*sizeof(double));
      memcpy(B,inB,size*sizeof(double));
      factorA=infactorA;
      factorB=infactorB;
      factorC=infactorC;
    }

#ifdef CMK_BLUEGENEL
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
	      CLA_Matrix_interface matB2, CLA_Matrix_interface matC2, 
	      CkCallback _orthoCB):
    m(_m), n(_n), matA (matA2), matB(matB2), matC(matC2), uponCompletion(_orthoCB)
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
      trigger = msg;
      A = trigger->A;  // no sense in making extra copies
      B = trigger->B;
      matA.multiply(trigger->factorA, 0, A, OrthoHelper::step_cb, (void*) this,
		     thisIndex.x, thisIndex.y);
      CmiNetworkProgress();
      matB.multiply(trigger->factorB, 0, B, OrthoHelper::step_cb, (void*) this,
		     thisIndex.x, thisIndex.y);
      CmiNetworkProgress();
      matC.multiply(trigger->factorC, 0, C, OrthoHelper::step_cb, (void*) this,
		     thisIndex.x, thisIndex.y);
      // DO NOT DELETE MSG we're using that memory in A and B
    }

  void sendMatrix()
    {
        if(trigger!=NULL)
            delete trigger;
        uponCompletion.send(m*n, C);
    }


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
  /// Callback to the owner ortho chare array to be used at the end of my work (step 2)
  CkCallback uponCompletion;
};

class OrthoHelperMap : public CkArrayMapTable2 {
  public:
    OrthoHelperMap(UberCollection _instance)
    {
      thisInstance=_instance;
#ifdef USE_INT_MAP
      maptable= &OrthoHelperImaptable[thisInstance.getPO()];
#else
      maptable= &OrthoHelpermaptable;
#endif
    }

    ~OrthoHelperMap() { }
    
    void pup(PUP::er &p)
    {
      CkArrayMapTable2::pup(p);
#ifdef USE_INT_MAP
      maptable= &OrthoHelperImaptable[thisInstance.getPO()];
#else
      maptable= &OrthoHelpermaptable;
#endif
    }
    
    inline int procNum(int, const CkArrayIndex &iIndex)
    {
      int *index=(int *) iIndex.data();
      int proc;
#ifdef USE_INT_MAP
      proc=maptable->get(index[0],index[1]);
#else
      proc=maptable->get(intdual(index[0],index[1]));
#endif
      if(fakeTorus)
	return(proc%CkNumPes());
      else
	return(proc);

    }
};

#endif 	    /* !ORTHOHELPER_H_ */
