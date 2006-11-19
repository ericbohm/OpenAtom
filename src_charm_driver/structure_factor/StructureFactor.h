/** \file StructureFactor.h
 *
 */

#ifndef _StructureFactor_h_
#define _StructureFactor_h_
void *fftw_malloc(size_t );
void fftw_free(void *);

/******************************************************************************
 * Three dimensions: atm_grp, plane, dup
 * Calculates structure factor, sends to sfcache group.
 * 
 * 
 *****************************************************************************/


extern Config config;

class SFDummyMsg: public CMessage_SFDummyMsg {
 public:
  int iteration_src;
};


class StructureFactor : public CBase_StructureFactor
{
 public:
  StructureFactor(CkMigrateMessage *m);
  StructureFactor(int _numSfGrps, int _numSfDups, int _natm_nl_grp_max, int _numdest, int *_destinations) : numSfGrps(_numSfGrps), numSfDups(_numSfDups), natm_nl_grp_max(_natm_nl_grp_max), numdest(_numdest)
    {
      k_x=NULL;
      k_y=NULL;
      k_z=NULL;
      structFactor=NULL;   
      structFactor_fx=NULL;
      structFactor_fy=NULL;
      structFactor_fz=NULL;
      destinations =new int[numdest];
      memcpy(destinations,_destinations, numdest *sizeof(int));
      
      gsSize=0;
//      usesAtSync=CmiTrue; not until we come up with a migration scheme
      setMigratable(false);
#ifdef _CP_DEBUG_SF_CALC_
  CkPrintf("[%d %d %d] created SF\n",thisIndex.x, thisIndex.y, thisIndex.z);
#endif
    }

  ~StructureFactor(){
    if(k_x!=NULL)
      {
	fftw_free(k_x);
	fftw_free(k_y);
	fftw_free(k_z);
      }
    if(structFactor!=NULL)
      {
	fftw_free(structFactor);
	fftw_free(structFactor_fx);
	fftw_free(structFactor_fy);
	fftw_free(structFactor_fz);
      }
    if(destinations!=NULL)
      delete [] destinations;
  }

  void acceptDestination(int _numdest, int *_destinations)
    {
#ifdef _CP_DEBUG_SF_CALC_
  CkPrintf("[%d %d %d] SF has %d destinations\n",thisIndex.x, thisIndex.y, thisIndex.z, _numdest);
#endif
      numdest=_numdest;
      destinations =new int[numdest];
      memcpy(destinations,_destinations, numdest *sizeof(int));
    }
  // compute and send
  void computeSF(SFDummyMsg *msg);

  // accept initializing kvector from gspace
  void acceptKVectors(int n, int *_k_x, int *_k_y, int *_k_z)
    {
      gsSize=n;
      k_x= (int *)fftw_malloc(gsSize*sizeof(int));
      k_y= (int *)fftw_malloc(gsSize*sizeof(int));
      k_z= (int *)fftw_malloc(gsSize*sizeof(int));
      memcpy(k_x,_k_x,gsSize*sizeof(int));
      memcpy(k_y,_k_y,gsSize*sizeof(int));
      memcpy(k_z,_k_z,gsSize*sizeof(int));
    }
  void pup(PUP::er &p)
    {
      p|gsSize;
      p|numSfGrps;
      p|numSfDups;
      p|natm_nl_grp_max;
      p|numdest;
      p((char*)structFactor,gsSize*natm_nl_grp_max*sizeof(complex));
      p((char*)structFactor_fx,gsSize*natm_nl_grp_max*sizeof(complex));
      p((char*)structFactor_fy,gsSize*natm_nl_grp_max*sizeof(complex));
      p((char*)structFactor_fz,gsSize*natm_nl_grp_max*sizeof(complex));
      p(k_x,gsSize);
      p(k_y,gsSize);
      p(k_z,gsSize);
      p(destinations, numdest);
    }

 private:	
  int numSfGrps;
  int numSfDups;
  int natm_nl_grp_max;
  int gsSize;
  int *k_x;
  int *k_y;
  int *k_z;
  complex *structFactor;   
  complex *structFactor_fx;
  complex *structFactor_fy;
  complex *structFactor_fz;
  int numdest;
  int *destinations;
};




#endif //_StructureFactor_h_
