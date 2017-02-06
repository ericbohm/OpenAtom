/*
 ** PIBeadAtoms.C
 ** 
 ** Made by (Eric Bohm)
 ** Login   <bohm@alacrity>
 ** 
 ** Started on  Tue Mar  2 10:15:33 2010 Eric Bohm
 ** Last update Sun May 12 01:17:25 2002 Speed Blue

 */


#include "AtomsCompute.h" 
#include "PIBeadAtoms.h"

extern CkVec <CProxy_AtomsCompute>                   UatomsComputeProxy;

// NOTE: thisIndex == atomIndex

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
PIBeadAtoms::PIBeadAtoms(UberCollection _thisInstance, int _numBeads, int _natm) :thisInstance(_thisInstance), 
  numBeads(_numBeads), natm (_natm)
                                   //============================================================================
{// begin routine
  //============================================================================
  //  CkPrintf("{%d}[%d] PIBeadAtoms::PIBeadAtoms numbeads %d\n",thisInstance.proxyOffset,thisIndex, numBeads);

  startAtm = (thisIndex * natm)/config.numBeadAtomChares;
  endAtm = ((thisIndex + 1) * natm)/config.numBeadAtomChares;
  numAtm = endAtm - startAtm;

  x= new double[numAtm * numBeads];
  y= new double[numAtm * numBeads];
  z= new double[numAtm * numBeads];
  xu= new double[numAtm * numBeads];
  yu= new double[numAtm * numBeads];
  zu= new double[numAtm * numBeads];
  fx= new double[numAtm * numBeads];
  fy= new double[numAtm * numBeads];
  fz= new double[numAtm * numBeads];
  fxu= new double[numAtm * numBeads];
  fyu= new double[numAtm * numBeads];
  fzu= new double[numAtm * numBeads];

  rat1 = new double[numBeads];
  rat2 = new double[numBeads];
  veig = new double[numBeads];

  //============================================================================
  // Initialize magic ratios

  rat1[0] = 0.0;
  rat2[0] = 0.0;
  veig[0] = 0.0;
  for(int ip=2;ip<=numBeads;ip++){
    int ip1 = ip-1;
    rat1[ip1] = ((double)(ip1))/((double)(ip));
    rat2[ip1] = 1.0/((double)(ip));
    veig[ip1] = ((double)(ip))/((double)(ip1));
  }//endfor

  //============================================================================
  // Initialize communication counters

  acceptCount_Fx=0;
  acceptCount_u=0;
  acceptCount_x=0;

  //============================================================================
}//end routine
//============================================================================



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void PIBeadAtoms::accept_PIMD_Fx(AtomXYZMsg *msg)
  //============================================================================
{// begin routine
  int curAtm = startAtm;
  int offset = msg->index;
  for(int i = 0; i < numAtm; i++) {
    fx[offset] = msg->x[curAtm];
    fy[offset] = msg->y[curAtm];
    fz[offset] = msg->z[curAtm];
    curAtm++;
    offset += numBeads;
  }
  delete msg;
  acceptCount_Fx++;

  if(acceptCount_Fx==numBeads){
    compute_PIMD_Fu();
    acceptCount_Fx=0;
    UberCollection instance=thisInstance;
    for(int bead=0;bead<numBeads; bead++){
      instance.idxU.x=bead;
      int proxyOffset=instance.setPO();
      //	  CkPrintf("{%d}[%d] PIBeadAtoms::accept_PIMD_Fx sending atomsGrp{%d}.accept_PIMD_fu\n",thisInstance.proxyOffset,thisIndex, proxyOffset);
      //	  UatomsGrpProxy[proxyOffset][atomdest].accept_PIMD_fu(fxu[bead],
      //	  fyu[bead], fzu[bead], thisIndex);
      // this atom index has to send the Fu to everyone
      AtomXYZMsg * toSend = new (numAtm, numAtm, numAtm,  8*sizeof(int)) AtomXYZMsg;
      toSend->index  = thisIndex;
      int offset = bead;
      for(int i = 0; i < numAtm; i++) {
        toSend->x[i] = fxu[offset];
        toSend->y[i] = fyu[offset];
        toSend->z[i] = fzu[offset];
        offset += numBeads;
      }
      UatomsComputeProxy[proxyOffset].accept_PIMD_Fu(toSend);
    }//endfor
  }//endif
  //============================================================================
}//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void PIBeadAtoms::accept_PIMD_Fx_and_x(AtomXYZMsg *msg)
  //============================================================================
{// begin routine
  //============================================================================
  int curAtm = startAtm;
  int offset = msg->index;
  for(int i = 0; i < numAtm; i++) {
    fx[offset] = msg->x[curAtm];
    fy[offset] = msg->y[curAtm];
    fz[offset] = msg->z[curAtm];
    x[offset] = msg->x[(curAtm+natm)];
    y[offset] = msg->y[(curAtm+natm)];
    z[offset] = msg->z[(curAtm+natm)];
    offset += numBeads;
    curAtm++;
  }
  delete msg;

  //============================================================================
  //  atomdest[PIBeadIndex]=atomdest;

  acceptCount_Fx++;
#ifdef _CP_DEBUG_ATMS_
  CkPrintf("{%d}[%d] PIBeadAtoms::accept_PIMD_Fx_and_x (%d of %d)\n",thisInstance.proxyOffset,thisIndex, acceptCount_Fx, numBeads);
#endif

  if(acceptCount_Fx==numBeads){
    compute_PIMD_Fu();
    compute_PIMD_u();
#ifdef _CP_DEBUG_ATMS_
    output_PIMD_u("c");
#endif
    acceptCount_Fx=0;
    UberCollection instance=thisInstance;
    for(int bead=0;bead<numBeads; bead++){
      instance.idxU.x=bead;
      int proxyOffset=instance.setPO();
#ifdef _CP_DEBUG_ATMS_
      CkPrintf("{%d}[%d] PIBeadAtoms::accept_PIMD_Fx_and_x sending atomsGrp{%d}.accept_PIMD_Fu_and_u\n",thisInstance.proxyOffset,thisIndex, proxyOffset);
#endif
      // this atom index has to send the Fu to everyone
      AtomXYZMsg * toSend = new (2*numAtm, 2*numAtm, 2*numAtm,  8*sizeof(int)) AtomXYZMsg;
      toSend->index  = thisIndex;
      int offset = bead;
      for(int i = 0; i < numAtm; i++) {
        toSend->x[i] = fxu[offset];
        toSend->y[i] = fyu[offset];
        toSend->z[i] = fzu[offset];
        toSend->x[numAtm + i] = xu[offset];
        toSend->y[numAtm + i] = yu[offset];
        toSend->z[numAtm + i] = zu[offset];
        offset += numBeads;
      }
      UatomsComputeProxy[proxyOffset].accept_PIMD_Fu_and_u(toSend);
    }//endfor
  }//endif
  //============================================================================
}//end routine
//============================================================================



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void PIBeadAtoms::accept_PIMD_u(AtomXYZMsg *msg)
  //============================================================================
{// begin routine
  //============================================================================
  int offset = msg->index;
  for(int i = 0; i < numAtm; i++) {
    xu[offset] = msg->x[i];
    yu[offset] = msg->y[i];
    zu[offset] = msg->z[i];
    offset += numBeads;
  }
#ifdef _CP_DEBUG_ATMS_
  CkPrintf("{%d}[%d] PIBeadAtoms::accept_PIMD_u (%d of %d) index %d\n",thisInstance.proxyOffset,thisIndex, acceptCount_u, numBeads, msg->index );
#endif
  delete msg;
  acceptCount_u++;
  if(acceptCount_u == numBeads){
#ifdef _CP_DEBUG_ATMS_
    output_PIMD_x("p");
#endif
    compute_PIMD_x();
#ifdef _CP_DEBUG_ATMS_
    output_PIMD_u("r");
    output_PIMD_x("c");
#endif
    acceptCount_u=0;
    UberCollection instance=thisInstance;
    for(int bead=0;bead<numBeads; bead++){
      instance.idxU.x=bead;
      int proxyOffset=instance.setPO();
      AtomXYZMsg * toSend = new (numAtm, numAtm, numAtm,  8*sizeof(int)) AtomXYZMsg;
      toSend->index  = thisIndex;
      int offset = bead;
      for(int i = 0; i < numAtm; i++) {
        toSend->x[i] = x[offset];
        toSend->y[i] = y[offset];
        toSend->z[i] = z[offset];
        offset += numBeads;
      }
      UatomsComputeProxy[proxyOffset].accept_PIMD_x(toSend);
    }//endfor
  }//endif

  //============================================================================
}//end routine
//============================================================================



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void PIBeadAtoms::accept_PIMD_x(AtomXYZMsg *msg)
  //============================================================================
{// begin routine
  //============================================================================
  int offset = msg->index;
  for(int i = 0; i < numAtm; i++) {
    x[offset] = msg->x[i];
    y[offset] = msg->y[i];
    z[offset] = msg->z[i];
    offset += numBeads;
  }
  delete msg;

  acceptCount_x++;
#ifdef _CP_DEBUG_ATMS_
  CkPrintf("{%d}[%d] PIBeadAtoms::accept_PIMD_x (%d of %d) \n",thisInstance.proxyOffset,thisIndex, acceptCount_x, numBeads);
#endif
  if(acceptCount_x==numBeads){
    compute_PIMD_u();
#ifdef _CP_DEBUG_ATMS_
    output_PIMD_u("c");
#endif
    acceptCount_x=0;
    UberCollection instance=thisInstance;
    for(int bead=0;bead<numBeads; bead++){
      instance.idxU.x=bead;
      int proxyOffset=instance.setPO();
      AtomXYZMsg * toSend = new (numAtm, numAtm, numAtm,  8*sizeof(int)) AtomXYZMsg;
      toSend->index  = thisIndex;
      int offset = bead;
      for(int i = 0; i < numAtm; i++) {
        toSend->x[i] = xu[offset];
        toSend->y[i] = yu[offset];
        toSend->z[i] = zu[offset];
        offset += numBeads;
      }
      UatomsComputeProxy[proxyOffset].accept_PIMD_u(toSend);
    }//endfor
  }//endif

  //============================================================================
}//end routine
//============================================================================



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void PIBeadAtoms::compute_PIMD_Fu()
  //============================================================================
{// begin routine
  //============================================================================

  if(numBeads<=0){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("The number of beads must be >=0\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
    exit(0);
  }//endif

  //============================================================================

  if(numBeads>1){
    for(int i = 0; i < numAtm; i++) {
      int offset = i * numBeads;
      fxu[offset]  = 0.0;  fyu[offset]  = 0.0;  fzu[offset]  = 0.0;
      for(int ip = 0; ip < numBeads; ip++) {
        fxu[offset] += fx[offset + ip]; 
        fyu[offset] += fy[offset + ip]; 
        fzu[offset] += fz[offset + ip];
      }/*endfor*/
      fxu[offset + 1] = fx[offset + 1];  
      fyu[offset + 1] = fy[offset + 1];  
      fzu[offset + 1] = fz[offset + 1];
      for(int ip = 1; ip < numBeads-1; ip++){
        int ip1 = ip+1;
        fxu[offset + ip1] = fx[offset + ip1] + rat1[ip]*fxu[offset + ip];
        fyu[offset + ip1] = fy[offset + ip1] + rat1[ip]*fyu[offset + ip];
        fzu[offset + ip1] = fz[offset + ip1] + rat1[ip]*fzu[offset + ip];
      }/*endfor*/
    }
  }else{
    for(int i = 0; i < numAtm; i++) {
      fxu[i]  = fx[i];  fyu[i]  = fy[i];  fzu[i] = fz[i];
    }
  }/*endif*/

  //============================================================================
}//end routine
//============================================================================




//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void PIBeadAtoms::compute_PIMD_u()
  //============================================================================
{// begin routine
  //============================================================================

  if(numBeads<=0){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("The number of beads must be >=0\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
    exit(0);
  }//endif

  //============================================================================

  if(numBeads>1){
    for(int i = 0; i < numAtm; i++) {
      int offset = i * numBeads;
      xu[offset]   = x[offset];  
      yu[offset]   = y[offset];  
      zu[offset]   = z[offset];
      for(int ip = 1; ip < numBeads - 1; ip++){
        int ip1 = ip+1;
        double xstar = rat1[ip] * x[offset + ip1] + rat2[ip]*x[offset];
        double ystar = rat1[ip] * y[offset + ip1] + rat2[ip]*y[offset];
        double zstar = rat1[ip] * z[offset + ip1] + rat2[ip]*z[offset];
        xu[offset + ip] =  x[offset + ip] - xstar;
        yu[offset + ip] =  y[offset + ip] - ystar;
        zu[offset + ip] =  z[offset + ip] - zstar;
      }/*endfor:ip*/
      int ip = numBeads-1;
      xu[offset + ip] =  x[offset + ip] - x[offset];
      yu[offset + ip] =  y[offset + ip] - y[offset];
      zu[offset + ip] =  z[offset + ip] - z[offset];
    }
  }else{
    for(int i = 0; i < numAtm; i++) {
      xu[i]  = x[i];  yu[i]  = y[i];  zu[i] = z[i];
    }
  }//endif

  //============================================================================
}//end routine
//============================================================================




//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void PIBeadAtoms::compute_PIMD_x()
  //============================================================================
{// begin routine
  //============================================================================

  if(numBeads<=0){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("The number of beads must be >=0\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
    exit(0);
  }//endif

  //============================================================================


  if(numBeads>1){
    for(int i = 0; i < numAtm; i++) {
      int offset = i * numBeads;
      x[offset] = xu[offset];  y[offset] = yu[offset];  z[offset] = zu[offset];
      int ip = numBeads-1;
      x[offset + ip] = xu[offset + ip] + x[offset]; 
      y[offset + ip] = yu[offset + ip] + y[offset]; 
      z[offset + ip] = zu[offset + ip] + z[offset];
      for(int ip=numBeads-2;ip>=1;ip--){
        int ip1 = ip+1;
        double xadd = rat1[ip] * x[offset + ip1] + xu[offset] * rat2[ip];
        double yadd = rat1[ip] * y[offset + ip1] + yu[offset] * rat2[ip];
        double zadd = rat1[ip] * z[offset + ip1] + zu[offset] * rat2[ip];
        x[offset + ip] = xu[offset + ip] + xadd;
        y[offset + ip] = yu[offset + ip] + yadd;
        z[offset + ip] = zu[offset + ip] + zadd;
      }/*endfor*/
    }
  }else{
    for(int i = 0; i < numAtm; i++) {
      int offset = i * numBeads; // technically correct, if silly as numBeads=0
      x[offset]  = xu[offset];  y[offset]  = yu[offset];  z[offset] = zu[offset];
    }
  }//endif

  //============================================================================
}//end routine
//============================================================================



//============================================================================
// useful for debugging
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void PIBeadAtoms::output_PIMD_x(char *s)
  //============================================================================
{// begin routine
  //============================================================================

  CkPrintf("\n=================================\n");
  for(int atm = 0; atm < numAtm; atm++) {
    for(int i =0; i< numBeads;i++) {
      CkPrintf("%s x[%d,%d]: %g %g %g\n", s, startAtm + atm, i, x[atm * numBeads + i], 
        y[atm * numBeads + i], z[atm * numBeads + i]);
    }//endfor
  }
  CkPrintf("=================================\n\n");

  //============================================================================
}//end routine
//============================================================================


//============================================================================
// useful for debugging
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void PIBeadAtoms::output_PIMD_u( char *s)
  //============================================================================
{// begin routine
  //============================================================================

  CkPrintf("\n=================================\n");
  for(int atm = 0; atm < numAtm; atm++) {
    for(int i =0;i<numBeads;i++){
      CkPrintf("%s u[%d,%d]: %g %g %g\n", s,startAtm + atm, i, xu[atm * numBeads + i],       
        yu[atm * numBeads + i], zu[atm * numBeads + i]);
    }//endfor
  }
  CkPrintf("=================================\n\n");

  //============================================================================
}//end routine
//============================================================================


//============================================================================
// useful for debugging
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void PIBeadAtoms::zero_PIMD_u()
  //============================================================================
{// begin routine
  //============================================================================

  for(int atm = 0; atm < numAtm; atm++) {
    for(int i =0;i<numBeads;i++){
      xu[atm * numBeads + i] = 0.0; 
      yu[atm * numBeads + i] = 0.0; 
      zu[atm * numBeads + i] = 0.0;;
    }//endfor
  }

  //============================================================================
}//end routine
//============================================================================

//============================================================================
// useful for debugging
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void PIBeadAtoms::zero_PIMD_fu()
  //============================================================================
{// begin routine
  //============================================================================

  for(int atm = 0; atm < numAtm; atm++) {
    for(int i =0;i<numBeads;i++){
      fxu[atm * numBeads + i] = 0.0; 
      fyu[atm * numBeads + i] = 0.0; 
      fzu[atm * numBeads + i] = 0.0;;
    }//endfor
  }

  //============================================================================
}//end routine
//============================================================================


//============================================================================
// useful for debugging
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void PIBeadAtoms::zero_PIMD_x()
  //============================================================================
{// begin routine
  //============================================================================

  for(int atm = 0; atm < numAtm; atm++) {
    for(int i = 0; i < numBeads; i++){
      x[atm * numBeads + i] = 0.0; 
      y[atm * numBeads + i] = 0.0; 
      z[atm * numBeads + i] = 0.0;;
    }//endfor
  }

  //============================================================================
}//end routine
//============================================================================


//============================================================================
// useful for debugging
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void PIBeadAtoms::energy_PIMD_x()
  //============================================================================
{// begin routine
  //============================================================================

  int ip  = 0;
  int ip1 = numBeads-1;
  double ep;
  for(int atm = 0; atm < numAtm; atm++) {
    int offset = atm * numBeads;
    ep = ( (x[offset + ip]-x[offset + ip1])*(x[offset + ip]-x[offset + ip1])
        +(y[offset + ip]-y[offset + ip1])*(y[offset + ip]-y[offset + ip1])
        +(z[offset + ip]-z[offset + ip1])*(z[offset + ip]-z[offset + ip1]) 
        );
    for(int ip=1;ip<numBeads;ip++){
      int ip1 = ip-1;
      ep += ( (x[offset + ip]-x[offset + ip1])*(x[offset + ip]-x[offset + ip1])
          +(y[offset + ip]-y[offset + ip1])*(y[offset + ip]-y[offset + ip1])
          +(z[offset + ip]-z[offset + ip1])*(z[offset + ip]-z[offset + ip1]) 
          );
    }//endfor
    CkPrintf("\n=================================\n");
    CkPrintf("x[%d] energy %.10g\n", startAtm + atm, ep);
    CkPrintf("=================================\n");
  }
  //============================================================================
}//end routine
//============================================================================


//============================================================================
// useful for debugging
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void PIBeadAtoms::energy_PIMD_u()
  //============================================================================
{// begin routine
  //============================================================================

  for(int atm = 0; atm < numAtm; atm++) {
    double ep = 0.0;
    int offset = atm * numBeads;
    for(int ip =1;ip<numBeads;ip++){
      double fk = veig[ip];
      ep += fk*(xu[offset + ip]*xu[offset + ip]
          +yu[offset + ip]*yu[offset + ip]
          +zu[offset + ip]*zu[offset + ip]
          );
    }//endfor
    CkPrintf("\n=================================\n");
    CkPrintf("u[%d] energy %.10g\n", startAtm + atm, ep);
    CkPrintf("=================================\n");
  }

  //============================================================================
}//end routine
//============================================================================



//============================================================================
// useful for debugging
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void PIBeadAtoms::checkUforce()
  //============================================================================
{// begin routine
  //============================================================================
#if 0
  printf("\n=================================\n");

  double potpx,potmx;
  double potpy,potmy;
  double potpz,potmz;
  double delta = 1.e-05;
  for(int ip =0;ip<numBeads;ip++){

    xu[ip] += delta;
    compute_PIMD_x();
    modelpot_PIMD_x(&potpx);
    xu[ip] -= 2.0*delta;
    compute_PIMD_x();
    modelpot_PIMD_x(&potmx);
    xu[ip] += delta;

    yu[ip] += delta;
    compute_PIMD_x();
    modelpot_PIMD_x(&potpy);
    yu[ip] -= 2.0*delta;
    compute_PIMD_x();
    modelpot_PIMD_x(&potmy);
    yu[ip] += delta;

    zu[ip] += delta;
    compute_PIMD_x();
    modelpot_PIMD_x(&potpz);
    zu[ip] -= 2.0*delta;
    compute_PIMD_x();
    modelpot_PIMD_x(&potmz);
    zu[ip] += delta;
    printf("fu[%d]: %.10g %.10g : %.10g %.10g : %.10g %.10g\n",
        ip,fxu[ip],0.5*(potmx-potpx)/delta,
        fyu[ip],0.5*(potmy-potpy)/delta,
        fzu[ip],0.5*(potmz-potpz)/delta);
  }//endfor

  printf("=================================\n");
#endif

  //============================================================================
}//end routine
//============================================================================


//============================================================================
// useful for debugging
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void PIBeadAtoms::modelpot_PIMD_x(double *pot_out)
  //============================================================================
{// begin routine
  //============================================================================

  for(int atm = 0; atm < numAtm; atm++) {
    double pot = 0.0;
    double fk = 0.5;
    int offset = atm * numBeads;
    for(int ip =0;ip<numBeads;ip++){
      pot += fk*(x[offset + ip]*x[offset + ip]
          +y[offset + ip]*y[offset + ip]
          +z[offset + ip]*z[offset + ip]
          );
    }//endfor

    pot_out[offset] = pot;
  }

  //============================================================================
}//end routine
//============================================================================

#include "PIBeadAtoms.def.h"
