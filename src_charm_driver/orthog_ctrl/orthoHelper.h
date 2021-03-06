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

#include "CLA_Matrix.h"

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class OrthoHelperMsg : public CMessage_OrthoHelperMsg {
  public:
    internalType *A;
    internalType *B;
    int size; 
    double factorA;
    double factorB;
    double factorC;
    void init(int insize, internalType* inA, internalType* inB, double infactorA, double infactorB, double infactorC)
    {
      size=insize;
      memcpy(A,inA,size*sizeof(internalType));
      memcpy(B,inB,size*sizeof(internalType));
      factorA=infactorA;
      factorB=infactorB;
      factorC=infactorC;
    }
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
    C= new internalType[m*n];
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
      matB.multiply(trigger->factorB, 0, B, OrthoHelper::step_cb, (void*) this,
          thisIndex.x, thisIndex.y);
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
        C = new internalType[m * n];
      }
      PUParray(p, C, m * n);
    }
    static inline void step_cb(void *obj){
      ((OrthoHelper *) obj)->sendMatrix();
    }
    int m;
    int n;
    internalType *A, *B, *C;
    OrthoHelperMsg *trigger;
    CLA_Matrix_interface matA, matB, matC;
    /// Callback to the owner ortho chare array to be used at the end of my work (step 2)
    CkCallback uponCompletion;
};

#endif 	    /* !ORTHOHELPER_H_ */
