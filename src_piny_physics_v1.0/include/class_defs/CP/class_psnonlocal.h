#ifndef _PSNONLOCAL_
#define _PSNONLOCAL_

#include "ckcomplex.h"

//==========================================================================
class PSNONLOCAL{
//-------------------------------------------------------------------------
  public:
   int natm;                //Num: Number of non-local atoms 
   int natm_typ;            //Num: Total num of atom types   
   int nl_max;              //Num: maximum occupied l-chan   
   int nl_max1;             //Num: maximum occupied l-chan+1
   int *map_nl;             // map[i]=j jth atm is ith non-local atom 
   int *natm_lang;          // # atms of this type: iatm_typ_lang     
   int *iatm_str_lang;      // where atm is in list: iatm_typ_lang    
   int *natm_typ_lang;      // # of non-local atm typs in lth channel 
   int **iatm_typ_lang;     // atom type lookup: [ityp][lang]
   double *x,*y,*z;
   double *vnorm_0;
   complex *ei_inc;
   complex *ti_inc;

//-------------------------------------------------------------------------
//con-destruct:
   PSNONLOCAL(){
     natm     = 0;
     natm_typ = 0;
     nl_max   = 0;
     nl_max1  = 0;
   };
  ~PSNONLOCAL(){};

//-------------------------------------------------------------------------
#ifdef PUP_ON
 //pupping
  void pup(PUP::er &p){
    //pupping ints
    p | natm;
    p | natm_typ;
    p | nl_max;
    p | nl_max1;
    if(natm>0){
      pup1d_int(p,&natm_typ_lang,nl_max1);
      pup1d_int(p,&map_nl,natm);
      pup1d_int(p,&natm_lang,natm_typ);
      pup1d_int(p,&iatm_str_lang,natm_typ);
      pup2d_int(p,&iatm_typ_lang,natm_typ,nl_max1);
      pup1d_dbl(p,&x,natm);
      pup1d_dbl(p,&y,natm);
      pup1d_dbl(p,&z,natm);
      pup1d_dbl(p,&vnorm_0,natm);
      pup1d_cpl(p,&ei_inc,natm);
      pup1d_cpl(p,&ti_inc,natm);
    }//endif
  }//end pupping
#endif
//-------------------------------------------------------------------------
}; //PSNONLOCAL
//==========================================================================

#ifdef PUP_ON
PUPmarshall(PSNONLOCAL);
#endif

#endif
