#ifndef _PSNONLOCAL_
#define _PSNONLOCAL_

#include "ckcomplex.h"

//==========================================================================
class PSNONLOCAL{
//-------------------------------------------------------------------------
  public:

   int natm_tot;            //Num: Total number of atoms
   int natm;                //Num: Number of non-local atoms <= total # of atoms
   int natm_typ;            //Num: Total num of atom types with ANY non-local channel
   int ees_eext_on;         //Opt: EES option on for external energy
   int ees_on;              //Opt: EES option on for non-local energy
   int n_interp;            //Order of interpolation : same for NL and EEXT EES
   int n_interp2;           //Num: n_interp^2
   int nk1,nk2,nk3;         //Num: kspace shape (-nk1,nk1)  : NL space
                            //                  (-nk2,nk2)
                            //                  (-nk3,nk3)
   int ngrid_a,ngrid_b,ngrid_c; //Num: EES real space mesh  : NL sizes
   int nfft;                //Num: ngrid_a*ngrid_b*ngrid_c

   int nl_iter;             //Num: Total number NL ees iterations
   int ntot_zmat;           //Num: Total zmat size
   int nmax_zmat;           //Num: Total zmat size
   int nl_max;              //Num: maximum occupied l-chan     
   int nl_max1;             //Num: maximum occupied l-chan+1  : spawns a lang channel loop
   int *natm_typ_lang;      // # of atom typs in the lth channel : spawns ityp loop
   int **iatm_typ_lang;     // atom type lookup: iatm_typ_lang[ityp][lang] = ityp_true
   int *natm_lang;          // # of atms a given type: natm_lang[ityp_true]
   int *iatm_str_lang;      // where atm of this type start in nl-atm-list:
   int *map_nl;             // map[i]=j : the ith atom in nl-list is the jth atm in the
                            // list of all atoms.
   int *lang_v;             // l-angular moment of ith iteration (lang)
   int *mang_v;             // m-angular moment of ith iteration (mang)
   int *ityp_v;             // ityp of ith iteration  (ityp)
   int *n_zmat;             // size of zmatrix for this iteration
   int *ioff_zmat;          // offset into the zmatrix for this iteration

   double ecut;             // Num: Cutoff in Hartree : FOR NL
   double fft_size_scale;   // Num: Increase in FFT size tuning EES accuracy : 
                            //      same for NL and EEXT
  double *vnorm_0;


//-------------------------------------------------------------------------
//con-destruct:
   PSNONLOCAL(){
     natm_tot    = 0;
     natm        = 0;
     natm_typ    = 0;
     ees_eext_on = 0;
     ees_on      = 0;
     n_interp    = 0;
     n_interp2   = 0;
     nk1         = 0;
     nk2         = 0;
     nk3         = 0;
     ngrid_a     = 0;
     ngrid_b     = 0;
     ngrid_c     = 0;
     nfft        = 0; 
     nl_max      = 0;
     nl_max1     = 0;
     nl_iter    = 0; 
     ntot_zmat   = 0;
     nmax_zmat   = 0;
     ecut        = 0.0;
     fft_size_scale = 0.0;
   };
  ~PSNONLOCAL(){};

//-------------------------------------------------------------------------
#ifdef PUP_ON
 //pupping
  void pup(PUP::er &p){
    //pupping ints
    p | natm_tot;
    p | natm;
    p | natm_typ;
    p | ees_on;
    p | ees_eext_on;
    p | n_interp;
    p | n_interp2;
    p | nl_max;
    p | nl_max1;
    p | nk1;      
    p | nk2;      
    p | nk3;
    p | ngrid_a;
    p | ngrid_b;
    p | ngrid_c;
    p | nfft; 
    p | nl_max;
    p | nl_max1;
    p | nl_iter; 
    p | ntot_zmat;
    p | nmax_zmat;
    //pupping dbles
    p | ecut;
    p | fft_size_scale;


     // the pups also malloc
    if(natm>0){
      pup1d_int(p,&natm_typ_lang,nl_max1);
      pup1d_int(p,&map_nl,natm);
      pup1d_int(p,&natm_lang,natm_typ);
      pup1d_int(p,&iatm_str_lang,natm_typ);
      pup1d_int(p,&lang_v,nl_iter);
      pup1d_int(p,&mang_v,nl_iter);
      pup1d_int(p,&ityp_v,nl_iter);
      pup1d_int(p,&n_zmat,nl_iter);
      pup1d_int(p,&ioff_zmat,nl_iter);
      pup2d_int(p,&iatm_typ_lang,natm_typ,nl_max1);
      pup1d_dbl(p,&vnorm_0,natm);
    }//endif


  }//end pupping
#endif
//-------------------------------------------------------------------------
}; //PSNONLOCAL
//==========================================================================

class PSSCRATCH
{
 public:
  // SF scratch for non-ees non-local method
  double *x,*y,*z;
  complex *ei_inc;
  complex *ti_inc;
  double *vtemp;
  int *index_atm;

  // Scratch for ees method : non-local AND eext scratch are here!!
  double *aj;
  double *rn;
  double *rn1;
  int *index_a;
  int *index_b;
  int *igrid_at;
  int *igrid_bt;
  int *iatemp;
  int *ibtemp;
  int *ictemp;
  double *frac_a;
  double *frac_b;
  double *frac_c;
  int **igrid_a;
  int **igrid_b;
  int **igrid_c;
  double **mn_a;
  double **mn_b;
  double **mn_c;
  double **ua;
  double **ub;
  double **uc;
  double **dmn_a;
  double **dmn_b;
  double **dmn_c;
  PSNONLOCAL *psnonlocal; // local ptr to the readonly stuff

  PSSCRATCH(PSNONLOCAL *_psnonlocal, CPATOM_MAPS *maps);
  PSSCRATCH(){};
  ~PSSCRATCH(){}; 
  // Hey Glenn, you never implemented a proper destructor for PSNONLOCAL.
  /* These are all local scratches, therefore we don't define a pup
     for them because actual usage is entirely ephemeral.  The only
     reason we don't malloc/free them on the fly is to avoid the
     overhead of doing those operations over and over again.  If you
     find yourself thinking we need to migrate these then you should
     seek counseling for extreme pupperphilia.
  */

}; //PSSCRATCH

#ifdef PUP_ON
PUPmarshall(PSNONLOCAL);
#endif

#endif
