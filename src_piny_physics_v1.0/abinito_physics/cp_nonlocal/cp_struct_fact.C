#include "standard_include.h"
#include "ckcomplex.h"

#include "../../../include/Atoms.h"
#include "../class_defs/allclass_gen.h"
#include "../class_defs/allclass_cp.h"
#include "../class_defs/allclass_mdintegrate.h"
#include "../class_defs/allclass_mdatoms.h"
#include "../class_defs/allclass_mdinter.h"
#include "../class_defs/allclass_mdintra.h"

#include "../class_defs/CP_OPERATIONS/class_cplocal.h"
#include "../class_defs/CP_OPERATIONS/class_cpnonlocal.h"

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//
// this function gets the following as input:
//
// a portion of gspace with gSpaceSize coefficients
//
// k_x, k_y, k_z are used to provide the g_x, g_y, g_z values for the
// coefficients
//
// atoms is an array of Atom classes. The atom class contains x, y, z 
// coordinates and the charge q.
//   
// The output is into the array zmatrix
// This should contain the z matrix reduced over the given portion of g-space
// The reduction could be done outside this function, but then the first thing
// after the function call would be a summation over the portion of g-space.
// Right now we assume l=0 and m=1. So zsize will be equal to natoms.
//
//============================================================================

void CPNONLOCAL::CP_calc_Struct_Fact(int gSpaceSize, 
                           int *k_x, int *k_y, int *k_z, 
                           complex *StructFact,
                           complex *StructFact_fx,
                           complex *StructFact_fy,
                           complex *StructFact_fz,
                           FastAtoms *atoms,
                           int mydoublePack, int numSfGrps, int indexSfGrp)

//============================================================================
  {// Begin function
//============================================================================

  MDINTEGRATE  *mdintegrate  = MDINTEGRATE::get();
  MDATOMS      *mdatoms      = MDATOMS::get();
  MDINTER      *mdinter      = MDINTER::get();
  MDINTRA      *mdintra      = MDINTRA::get();
  GENERAL_DATA *general_data = GENERAL_DATA::get();
  CP           *cp           = CP::get();

#include "../class_defs/allclass_strip_cp.h"
#include "../class_defs/allclass_strip_mdatoms.h"
#include "../class_defs/allclass_strip_gen.h"

  /*------------------*/
  /* Atom information */

  int natm_nl        = cppseudo->nonlocal.natm;
  int natm_typ       = cppseudo->nonlocal.natm_typ;
  int *map_nl        = cppseudo->nonlocal.map_nl;
  int *natm_lang     = cppseudo->nonlocal.natm_lang;
  int *iatm_str_lang = cppseudo->nonlocal.iatm_str_lang;
  int *natm_typ_lang = cppseudo->nonlocal.natm_typ_lang;
  int **iatm_typ_lang= cppseudo->nonlocal.iatm_typ_lang;

  double *x          = cppseudo->nonlocal.x;
  double *y          = cppseudo->nonlocal.y;
  double *z          = cppseudo->nonlocal.z;
  complex *ei_inc    = cppseudo->nonlocal.ei_inc;
  complex *s_now     = cppseudo->nonlocal.ti_inc;

  /*--------------------------------*/
  /* Cell and pressure information  */

  double *hmati     = gencell->hmati;  
  double vol        = gencell->vol;

  /*----------------------------------------------*/
  /* Pseudo-potential and Reduced Periodicity info*/

  int nsplin_g      = cppseudo->nsplin_g;
  double gmin_spl   = cppseudo->gmin_spl;
  double dg_spl     = cppseudo->dg_spl;
  double *vps0      = cppseudo->vps0;
  double *vps1      = cppseudo->vps1;
  double *vps2      = cppseudo->vps2;
  double *vps3      = cppseudo->vps3;
  double *vpsnorm   = cppseudo->vpsnorm;
  double *gzvps0    = cppseudo->gzvps0;
  int n_ang_max     = cppseudo->n_ang_max;

  double wght,wght_now;
// Local variable declarations   
  double tpi = 2.0*M_PI;
  double fpi = 4.0*M_PI;

  double *xatm      = atoms->x;
  double *yatm      = atoms->y;
  double *zatm      = atoms->z;

#ifdef GJM_DEBUG_SIZE
   PRINTF(" %d : coefs in CP_calc_Struct_Fact: %d\n",gSpaceSize,mydoublePack);
#endif

//============================================================================
// Check sizes, find your piece of the action

   if(indexSfGrp>=numSfGrps ||  indexSfGrp<0 || numSfGrps>natm_nl){
     PRINTF("Incorrect SF index %d %d %d\n",indexSfGrp,numSfGrps,natm_nl);
     EXIT(1);
   }//endif

   int natm_nl_grp,istrt,iend,ioff;
   get_grp_params(natm_nl,numSfGrps,indexSfGrp,&natm_nl_grp,&istrt,&iend,&ioff);

//============================================================================
// Strip out the non-local atoms and compute some helpful vectors 

   for(int iatm=istrt;iatm<=iend;iatm++){
     int ind = map_nl[iatm]-1;
     int i   = iatm-istrt+1;
     x[i]    = xatm[ind];
     y[i]    = yatm[ind];
     z[i]    = zatm[ind];
   }//endfor

   for(int iatm=1;iatm<=natm_nl_grp;iatm++){
     double arg   = tpi*(hmati[1]*x[iatm]+hmati[4]*y[iatm]+hmati[7]*z[iatm]);
     ei_inc[iatm] = complex(cos(arg),sin(arg));
   } // endfor

   double y00 = 1.0/sqrt(fpi);

   int kx_old = k_x[0];
   int ky_old = k_y[0];
   int kz_old = k_z[0];

   wght = 1.0;  if(mydoublePack==1){wght = 2.0;}
   int igo = 0;

//============================================================================
// Loop over this piece of g-space and all nonlocal atoms; Construct this
// piece of the structure factor (i.e. for the current g-space plane) 

   int ic = 0;
   for(int i=0;i<gSpaceSize;i++){

     if(i % 8 == 0)
       CmiNetworkProgress();

      double gx,gy,gz,g2,g;
      gx = tpi*(k_x[i]*hmati[1] + k_y[i]*hmati[2] + k_z[i]*hmati[3]);
      gy = tpi*(k_x[i]*hmati[4] + k_y[i]*hmati[5] + k_z[i]*hmati[6]);
      gz = tpi*(k_x[i]*hmati[7] + k_y[i]*hmati[8] + k_z[i]*hmati[9]);
      g2 = gx*gx + gy*gy + gz*gz;
      g  = sqrt(g2);

      //--------------------------------------------------------------------
      if(kx_old != k_x[i]-1 || ky_old != k_y[i] || kz_old != k_z[i] || igo==0){
        for(int iatm=1;iatm<=natm_nl_grp;iatm++){
          double arg  = x[iatm]*gx + y[iatm]*gy + z[iatm]*gz;
          double si   = sin(arg);
          double ca   = cos(arg);
          s_now[iatm] = complex(ca,si);
	}//endfor*/
      }//endif
      kx_old = k_x[i];   ky_old = k_y[i];  kz_old = k_z[i];
      igo = 1;

      wght_now = (k_x[i]==0 ? 1.0 : wght);      
      ic++;
      if(k_x[i]==0 && k_y[i]<0){ic--;}
      if(k_x[i]==0 && k_y[i]==0 && k_z[i]<0){ic--;}

      for(int iatm=1;iatm<=natm_nl_grp;iatm++){
        int s_ind  = gSpaceSize*(iatm-1) + i;
        StructFact[s_ind] = s_now[iatm]*wght_now;
        complex tmp = complex(-s_now[iatm].im,s_now[iatm].re);
        tmp = tmp*wght_now;
        StructFact_fx[s_ind]  = tmp*gx;
        StructFact_fy[s_ind]  = tmp*gy;
        StructFact_fz[s_ind]  = tmp*gz;
      }// endfor

      for(int iatm=1;iatm<=natm_nl_grp;iatm++){
        s_now[iatm]=s_now[iatm]*ei_inc[iatm];
      }//endfor

      //---------------------------------------------------------------------
      // hard wire l=0
      int lang = 0; int lang1 = 1;
      int natm_typ = natm_typ_lang[lang1];
      if(natm_typ>0){

        int iatm=0;
        for(int ityp=1;ityp<=natm_typ;ityp++){      // l=0 atm types
          int iatm_typ = iatm_typ_lang[ityp][lang1];// true atom type index
          int natm_now = natm_lang[iatm_typ];       // # atms of this type 
          int iatm_str = iatm_str_lang[iatm_typ];   // where atms strt in list
          int iatm_end = iatm_str+natm_now-1;       // where atms end in list
          double proj_00;
          get_rad_proj(lang,iatm_typ,g,dg_spl,gmin_spl,nsplin_g,n_ang_max,
                       vps0,vps1,vps2,vps3,gzvps0,&proj_00);
          proj_00 *= y00;
          if(iatm_str<=iend && iatm_end>=istrt){
            int jstrt = MAX(iatm_str,istrt);
            int jend  = MIN(iatm_end,iend);
            for(int jatm=jstrt;jatm<=jend;jatm++){        // atms of this type
              iatm++;
              int s_ind  = gSpaceSize*(iatm-1)+i;
              StructFact[s_ind]    = StructFact[s_ind]*proj_00;
              StructFact_fx[s_ind] = StructFact_fx[s_ind]*proj_00;
              StructFact_fy[s_ind] = StructFact_fy[s_ind]*proj_00;
              StructFact_fz[s_ind] = StructFact_fz[s_ind]*proj_00;
  	    }//endfor
	  }//endif
	}//endfor
        if(iatm!=natm_nl_grp){
           PRINTF("Internal error in SF %d %d\n",iatm,natm_nl_grp);
           EXIT(1);
	}//endif

      }//endif : l=0 atom types
     //-----------------------------------------------------------------------
   }// endfor : g-space

//============================================================================

   if(mydoublePack==1){
     int nsize = natm_nl_grp*gSpaceSize;
     for(int i=0;i<nsize;i++){StructFact[i].im    = -StructFact[i].im;}
     for(int i=0;i<nsize;i++){StructFact_fx[i].im = -StructFact_fx[i].im;}
     for(int i=0;i<nsize;i++){StructFact_fy[i].im = -StructFact_fy[i].im;}
     for(int i=0;i<nsize;i++){StructFact_fz[i].im = -StructFact_fz[i].im;}
   }

//============================================================================

#ifdef GJM_DEBUG_SIZE
   PRINTF(" %d : coefs in CP_calc_Struct_Fact: %d\n",ic,mydoublePack);
#endif

//============================================================================
  } // End function
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================

void CPNONLOCAL::get_rad_proj(int lang,int ityp,double g,
                  double dg_spl,double gmin_spl,int nsplin_g,
                  int n_ang_max, double *vps0,double *vps1,double *vps2,
                  double *vps3,double *gzvps0,double *proj_ret)

//============================================================================
  {// Begin function
//============================================================================
 
    double vnow;

    if(g!=0.0){

      int iii   = ((int)((g-gmin_spl)/dg_spl)) + 1;
      iii       = MIN(iii,nsplin_g);
      iii       = MAX(iii,1);
      double h0 = ((double)(iii-1))*dg_spl+gmin_spl;
      double h  = g-h0;

      int ind_now = (ityp-1)*nsplin_g*(n_ang_max+1) + lang*nsplin_g + iii;
      double v0   = vps0[ind_now];
      double v1   = vps1[ind_now];
      double v2   = vps2[ind_now];
      double v3   = vps3[ind_now];
      vnow        = ((v3*h+v2)*h+v1)*h+v0;

    }else{

      vnow       = gzvps0[ityp];

    }//endif

    (*proj_ret) = vnow;

//---------------------------------------------------------------------------
  } // End function
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================

void CPNONLOCAL::get_grp_params(int natm_nl, int numSfGrps, int indexSfGrp,
                                int *n_ret, int *istrt_ret, int *iend_ret,
                                int *ioff_ret)

//============================================================================
{

   int n     = (natm_nl/numSfGrps);
   int m     = (natm_nl % numSfGrps);

   int n_max = n;
   if(m!=0){n_max++;}
   int ioff  = n_max*indexSfGrp;

   int istrt = n*indexSfGrp + 1;
   if(indexSfGrp>=m){istrt += m;}
   if(indexSfGrp<m) {istrt += indexSfGrp;}
   if(indexSfGrp<m) {n++;}
   int iend  = n+istrt-1;
   if(numSfGrps>natm_nl)
     if(m>=indexSfGrp)
       {
	 n=0;
	 istrt=1;
	 iend=0;
       }

   (*n_ret)     = n;
   (*istrt_ret) = istrt;
   (*iend_ret)  = iend;
   (*ioff_ret)  = ioff;
}
//============================================================================
