//========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//========================================================================

#include "standard_include.h"
#include "ckcomplex.h"
#include "../class_defs/typedefs_par.h"

#include "../class_defs/allclass_gen.h"
#include "../class_defs/allclass_cp.h"
#include "../class_defs/allclass_mdatoms.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_handle_entry.h"

#include "../class_defs/CP/class_gen_wave.h"

//========================================================================


//========================================================================
// Simple constructor
//========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//========================================================================
GEN_WAVE::GEN_WAVE(){

   natm_typ_cp = 0;
   nsplin      = 0;
   cp_lda      = 0;
   cp_lsda     = 0;

   dg          = 0.0;
   gmin        = 0.0;
   gmax        = 0.0;

   n_ang           = NULL;  // natm_typ_cp
   iatm_atm_typ_cp = NULL;  // nab_initio
   iatm_state_str  = NULL;  // nab_initio
   iatm_state_end  = NULL;  // nab_initio


   gpsi0           = NULL;     // natm_typ_cp x 3 x nsplin
   gpsi1           = NULL;
   gpsi2           = NULL;
   gpsi3           = NULL;
   gpsi_now        = NULL;    // natm_typ_cp x 3
   gpsi00          = NULL;    // natm_typ_cp

//========================================================================
   }// end routine
//========================================================================


//========================================================================
// Simple destructor
//========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//========================================================================
GEN_WAVE::~GEN_WAVE( ){

  if(natm_typ_cp>0){

    cfree(&(n_ang[1]),"gen_wave_destructor");
    cfree(&(iatm_atm_typ_cp[1]),"gen_wave_destructor");
    cfree(&(iatm_state_str[1]),"gen_wave_destructor");
    cfree(&(iatm_state_end[1]),"gen_wave_destructor");

    cfree(&(gpsi00[1]),"gen_wave_destructor");
    cfree_mat(gpsi_now,1,natm_typ_cp,1,3);

    int i,j;
    for(i=1; i <= natm_typ_cp; i++){
     for(j=1; j<=3; j++){
       cfree(&(gpsi0[i][j][1]),"gen_wave_destructor");
       cfree(&(gpsi1[i][j][1]),"gen_wave_destructor");
       cfree(&(gpsi2[i][j][1]),"gen_wave_destructor");
       cfree(&(gpsi3[i][j][1]),"gen_wave_destructor");
     }//endfor
     cfree(&(gpsi0[i][1]),"gen_wave_destructor");
     cfree(&(gpsi1[i][1]),"gen_wave_destructor");
     cfree(&(gpsi2[i][1]),"gen_wave_destructor");
     cfree(&(gpsi3[i][1]),"gen_wave_destructor");
    }//endfor
    cfree(&(gpsi0[1]),"gen_wave_destructor");
    cfree(&(gpsi1[1]),"gen_wave_destructor");
    cfree(&(gpsi2[1]),"gen_wave_destructor");
    cfree(&(gpsi3[1]),"gen_wave_destructor");

  }//endif

//========================================================================
  } // end routine destructor
//========================================================================


//========================================================================
// Initialize the gen_wave data structures
//========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//========================================================================
void GEN_WAVE::fill_gw_gpsi(CPATOM_MAPS * cpatom_maps,CPCOEFFS_INFO *cpcoeffs_info,
                            int nsplin_in, double gmin_in,double gmax_in, 
			    int *iatm_atm_typ, int natm_typ,NAME *vps_name,
                            int cp_lda_in, int cp_lsda_in)
//========================================================================
  { //begin routine
//========================================================================
//             Local variable declarations                               
  
  int i,j,k,iatm;
  int ind1,ind2,n;
  int iflag;

  nsplin            = nsplin_in;
  gmin              = gmin_in;
  gmax              = gmax_in;
  cp_lda            = cp_lda_in;
  cp_lsda           = cp_lsda_in;
  nab_initio        = cpatom_maps->nab_initio;

  int   *cp_atm_lst = cpatom_maps->cp_atm_lst;
  int   *cp_vlnc_up = cpatom_maps->cp_vlnc_up;
  int   *cp_vlnc_dn = cpatom_maps->cp_vlnc_dn;

  int nstate_up     = cpcoeffs_info->nstate_up;
  int nstate_dn     = cpcoeffs_info->nstate_dn; 

  double *g;

//====================================================================

  PRINTF("\n");
  PRINT_LINE_STAR
  PRINTF("Initializing GenWave\n");
  PRINT_LINE_DASH;PRINTF("\n");

//=======================================================================
// Check Sum of cp_valences is equal to the number of states in set file  

  if(cp_lda!=1 && cp_lsda !=0){
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     PRINTF("GenWave does not support lsda\n");
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     EXIT(1);
  }//endif

  if(nab_initio == 0){
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     PRINTF("There are no ab initio atoms \n");
     PRINTF("Please check whether you should have assigned some atoms \n");
     PRINTF("to be ab initio.  To proceed would be pointless.\n");
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     EXIT(1);
  }//endif

  int  nstate_up_gw = 0;
  int  nstate_dn_gw = 0;
  for(i=1; i<= nab_initio; i++){
   iatm = cp_atm_lst[i];
   nstate_up_gw += cp_vlnc_up[iatm];
   nstate_dn_gw += cp_vlnc_dn[iatm];
  }//endfor

  if( nstate_up_gw != nstate_up ){
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     PRINTF("The number of states up not equal to what is in the set file\n");
     PRINTF("%d here and %d from the set file\n",nstate_up_gw, nstate_up);
     PRINTF("Please check the values of cp_valence_up in the parm files\n");
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     EXIT(1);
  }//endif

  if( nstate_dn_gw != nstate_dn  && cp_lsda==1){
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     PRINTF("The number of states dn not equal to what is in the set file\n");
     PRINTF("%d here and %d from the set file\n",nstate_dn_gw, nstate_dn);
     PRINTF("Please check the values of cp_valence_dn in the parm files\n");
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     EXIT(1);
  }//endif

  if( nstate_dn_gw != nstate_dn  && cp_lda==1){
     PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
     PRINTF("The number of states dn not equal to value in the set file.\n");
     PRINTF("This is OK only if you have some unit occupation numbers. \n");
     PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
     EXIT(1);
  }//endif

  if( nstate_dn_gw > nstate_dn  && cp_lda==1){
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     PRINTF("The number of dn nstates is assumed <= up nstates \n");
     PRINTF("in gen_wave under LDA. Since you already have the warning\n" );
     PRINTF("about occupation numbers and states what you have to do is\n" );
     PRINTF("check the values of cp_valence_(dn/up) in the parm files\n");
     PRINTF("to ensure that the up states get the extra occupation \n");
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     EXIT(1);
  }//endif

//=======================================================================
// check that the value assigned for cp_valence_up is same for 
// for all atoms of the same atm_typ                         

  for(i=1; i <= nab_initio; i++){
    ind1 = cp_atm_lst[i];
    for(j=i+1; j <= nab_initio;j++){
     ind2 = cp_atm_lst[j];
     if( iatm_atm_typ[ind2] == iatm_atm_typ[ind1] &&
         cp_vlnc_up[ind2] != cp_vlnc_up[ind1]){
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     PRINTF("There are two different assigned values for    \n");
     PRINTF("cp_valence_up for the same atom type           \n");
     PRINTF("atom numbers %d and %d \n",ind1,ind2);
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     EXIT(1);
     }//endif

   if( cp_lsda == 1 && iatm_atm_typ[ind2] == iatm_atm_typ[ind1] &&
       cp_vlnc_dn[ind2] != cp_vlnc_dn[ind1]){
       PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       PRINTF("There are two different assigned values for    \n");
       PRINTF("cp_valence_dn for the same atom type           \n");
       PRINTF("atom numbers %d and %d \n",ind1,ind2);
       PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       EXIT(1);
      }//endif
    }//endfor
  }//endfor

//=======================================================================
// Count the number CP atom types

  natm_typ_cp = 1;
  for(i=2; i<= nab_initio; i++){
    iflag = 0;
    ind1 = cp_atm_lst[i];
    for(j=i-1; j >= 1; j--){
      ind2 = cp_atm_lst[j];
      if( iatm_atm_typ[ind2] == iatm_atm_typ[ind1] ) {
        iflag = 0;break;
      }else{
        iflag = 1;
      }//endif
    }//endfor
    natm_typ_cp += iflag;
  }//endfor

//=======================================================================
// malloc 
   
  iatm_atm_typ_cp = (int *) cmalloc(nab_initio*sizeof(int),"gen_wave_init") -1;
  iatm_state_str  = (int *) cmalloc(nab_initio*sizeof(int),"gen_wave_init") -1;
  iatm_state_end  = (int *) cmalloc(nab_initio*sizeof(int),"gen_wave_init") -1;

  NAME *fname_ps = (NAME *) cmalloc(natm_typ_cp*sizeof(NAME),"gen_wave_init")-1;

//=======================================================================
// Assign cp_atm_types 

  int tag = 1;  // used as counter for unique cp atom types 
  iatm_atm_typ_cp[1] = 1;
  strcpy(fname_ps[1],vps_name[iatm_atm_typ[cp_atm_lst[1]]]);

  for(i=2; i<= nab_initio; i++){
    iflag = 0;
    ind1 = cp_atm_lst[i];
    for(j=i-1; j >= 1; j--){
      ind2 = cp_atm_lst[j];
      if( iatm_atm_typ[ind1] == iatm_atm_typ[ind2] ){
        iatm_atm_typ_cp[i] = iatm_atm_typ_cp[j];
        iflag = 0;
        break;
      }else{
        iflag = 1;
      }//endif
    }//endfor
    tag += iflag;
    if(iflag == 1){
      iatm_atm_typ_cp[i] = tag;
      strcpy(fname_ps[tag],vps_name[iatm_atm_typ[ind1]]);
    }//endif
  }//endfor

//=========================================================================
//  Malloc memory for Bessel transform of radial wavefunctions.

  n_ang    = (int *) cmalloc(natm_typ_cp*sizeof(int ),"gen_wave_init") -1;
  gpsi_now = cmall_mat(1,natm_typ_cp,1,3,"gen_wave_init");
  gpsi00   = (double *) cmalloc(natm_typ_cp*sizeof(double),"gen_wave_init") -1;

  gpsi0 = (double ***) cmalloc(natm_typ_cp*sizeof(double **),"gen_wave_init")-1;
  gpsi1 = (double ***) cmalloc(natm_typ_cp*sizeof(double **),"gen_wave_init")-1;
  gpsi2 = (double ***) cmalloc(natm_typ_cp*sizeof(double **),"gen_wave_init")-1;
  gpsi3 = (double ***) cmalloc(natm_typ_cp*sizeof(double **),"gen_wave_init")-1;

  for(i=1; i<= natm_typ_cp; i++){
    gpsi0[i] = (double **) cmalloc(3*sizeof(double *),"gen_wave_init")-1;
    gpsi1[i] = (double **) cmalloc(3*sizeof(double *),"gen_wave_init")-1;
    gpsi2[i] = (double **) cmalloc(3*sizeof(double *),"gen_wave_init")-1;
    gpsi3[i] = (double **) cmalloc(3*sizeof(double *),"gen_wave_init")-1;
    for(j=1; j<=3; j++){
      gpsi0[i][j] = (double *) cmalloc(nsplin*sizeof(double ),"gen_wave_init")-1;
      gpsi1[i][j] = (double *) cmalloc(nsplin*sizeof(double ),"gen_wave_init")-1;
      gpsi2[i][j] = (double *) cmalloc(nsplin*sizeof(double ),"gen_wave_init")-1;
      gpsi3[i][j] = (double *) cmalloc(nsplin*sizeof(double ),"gen_wave_init")-1;
    }//endfor
  }//endfor

//===========================================================================
// assign the occupation numbers 

  double *occ_dn = (double *) cmalloc(nstate_up*sizeof(double),"gen_wave_init") -1;
  double *occ_up = (double *) cmalloc(nstate_up*sizeof(double),"gen_wave_init") -1;

  for(i=1; i<= nstate_up; i++){
    occ_up[i] = 0.0;
    occ_dn[i] = 0.0;
  }//endfor

  for(i=1; i<= nstate_up_gw; i++){
    occ_up[i] = 1.0;
  }//endfor

  for(i=1; i<= nstate_dn_gw; i++){
    occ_dn[i] = 1.0;
  }//endfor

  if(cp_lda==1){
    for(i=1; i<= nstate_up; i++){
      occ_up[i] += occ_dn[i];
    }//endfor
  }//endif cp_lda

//=========================================================================
// bessel transform the radial wave functions and spline the result        
// nsplin rows 3 columns natm_typ dimensions  [d][c][r]                    

  dg = (gmax-gmin)/(double)(nsplin);
  g  = (double *) cmalloc(nsplin*sizeof(double),"splin_btrans") -1;

  for(i=1; i<= nsplin; i++){
    g[i] = dg*(double)(i-1) + gmin;
  }//endfor

  for(i=1; i<= natm_typ_cp; i++){
    splin_btrans(g,gpsi0[i],gpsi1[i],gpsi2[i],gpsi3[i],
                &(gpsi00[i]),&(n_ang[i]),fname_ps[i]);
  }//endfor
 
  cfree(&g[1],"gen_wave");
  cfree(&(fname_ps[1]),"fill_gw_gpsi");

//===========================================================================
// Create map : Which states are assigned to each atom's radial psi.

  iatm_state_str[1] = 1;
  iatm_state_end[1] = iatm_state_str[1] + cp_vlnc_up[1]-1;

  for(i=2; i<= nab_initio; i++){
    iatm_state_str[i] = iatm_state_end[i-1] + 1;
    iatm_state_end[i] = iatm_state_str[i] + cp_vlnc_up[i]-1;
  }//endfor

//===========================================================================

  PRINTF("\n");PRINT_LINE_DASH
  PRINTF("Completed Initializing GEN_WAVE\n");
  PRINT_LINE_STAR;PRINTF("\n");

//==========================================================================
  } //  end routine
//==========================================================================


//========================================================================
// Create a KS state in g-space using radial wavefunctions
//========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//========================================================================
void GEN_WAVE::create_coefs(int *k_x,int *k_y,int *k_z,
		            int gSpaceSize,int state_ind_m1,complex *gspace_coefs,
                            double *xfull, double *yfull, double *zfull)
//========================================================================
   {//begin routine
//========================================================================

  GENERAL_DATA *general_data = GENERAL_DATA::get();
  CP           *cp           = CP::get();

#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_cp.h"

//----------------------------------------------------------------------

  int nab_initio_g = cpatom_maps->nab_initio;
  int *cp_atm_lst  = cpatom_maps->cp_atm_lst;

  double rt_fpi    = cpylm_cons->rt_fpi;
  double pi        = M_PI;
  double tpi       = 2.0*pi;
  double rad2      = sqrt(2.0);

  double psi_r[21],psi_i[21];
  double ylmr[21],ylmi[21];

  if(nab_initio_g!=nab_initio){
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     PRINTF("GenWave :  bad number of ab initio atoms %d %d\n",
               nab_initio_g,nab_initio);
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     EXIT(1);
  }//endif

//----------------------------------------------------------------------
  double *cp_box_center     =  gencell->cp_box_center;
  double *cp_box_center_rel =  gencell->cp_box_center_rel;
  double *hmat              =  gencell->hmat_cp;
  double *hmati             =  gencell->hmati_cp;
  double *hmat_big          =  gencell->hmat;
  double *hmati_big         =  gencell->hmati;
  double  vol               =  gencell->vol_cp;
  double volrt              =  sqrt(vol);

  int state_ind             =  state_ind_m1+1;

//=========================================================================
// I need a random seed

  long seed=103481;
  double xx;
  for(int i=1;i<=state_ind;i++){xx = altRandom(&seed);}
  seed = (long)(1.e6*xx);
  while(seed<10){seed*=10.0;}

//=========================================================================
// Which atom's radial psi forms this KS state

  int my_atom = 0;
  for(int i=1;i<=nab_initio;i++){
    if(state_ind>=iatm_state_str[i] &&
       state_ind<=iatm_state_end[i]){my_atom=i;}
  }//endfor

  if(my_atom==0){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("YIKES GEN_WAVE can't find an atomic psi for the present state\n");
    PRINTF("atom= %d state=%d nab_atom=%d\n",my_atom,state_ind_m1,nab_initio);
    for(int i=1;i<=nab_initio;i++){
      PRINTF("  iatm=%d state_low=%d state_high=%d\n",i,iatm_state_str[i],iatm_state_end[i]);
    }//endfor
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    EXIT(1);
  }//endif

  int ipart = my_atom;
  int is    = state_ind - iatm_state_str[ipart] + 1;
  int ityp  = iatm_atm_typ_cp[ipart];
  int iatm  = cp_atm_lst[ipart];

  double x  = xfull[iatm];
  double y  = yfull[iatm];
  double z  = zfull[iatm];

//=========================================================================
// get the wave functions in g space                                       

  double anorm = 0.0;
  for(int i=0; i < gSpaceSize; i++){

//-------------------------------------------------------------------------
// get g vectors                                                           

    double gx = tpi*(k_x[i]*hmati[1] + k_y[i]*hmati[2] + k_z[i]*hmati[3]);
    double gy = tpi*(k_x[i]*hmati[4] + k_y[i]*hmati[5] + k_z[i]*hmati[6]);
    double gz = tpi*(k_x[i]*hmati[7] + k_y[i]*hmati[8] + k_z[i]*hmati[9]);
    double g2 = gx*gx + gy*gy + gz*gz;
    double g  = sqrt(g2);
  
//-------------------------------------------------------------------------
//  get the spherical bessel tranform of the radial wave functions         
//  at this g space point using the spline and the spherical harmonics

    if(k_x[i]==0 && k_y[i]==0 && k_z[i]==0){   
       gpsi_now[ityp][1] = gpsi00[ityp];
       gpsi_now[ityp][2] = 0.0;
       gpsi_now[ityp][3] = 0.0;
       for(int j=1;j<=8;j++){ylmr[j]=0;ylmi[j]=0;}
       ylmr[1] = rt_fpi;
    }else{
       get_gpsi(g,gpsi0[ityp],gpsi1[ityp],gpsi2[ityp],gpsi3[ityp],
                gpsi_now[ityp],n_ang[ityp]);
       get_ylm(gx,gy,gz,g,ylmr,ylmi,cpylm_cons);
    }//endif

//  S STATE                                                                
    psi_r[1] = ylmr[1]*gpsi_now[ityp][1]/volrt; 
    psi_i[1] = 0.0;

// SPHERICALIZED P BAND                                                   
    double p_a,p_b,p_c;
    if(n_ang[ityp] >= 1){
      int itemp;
      itemp = (int) (3.0*altRandom(&seed));
      itemp = MAX(itemp,0);
      itemp = MIN(itemp,2);
      switch(itemp){
       case 0:
        p_c = -ylmr[2]*gpsi_now[ityp][2]/volrt;
        p_b = -rad2*ylmr[3]*gpsi_now[ityp][2]/volrt;
        p_a = -rad2*ylmi[3]*gpsi_now[ityp][2]/volrt;
       break;
       case 1:
        p_a = -ylmr[2]*gpsi_now[ityp][2]/volrt;
        p_c = -rad2*ylmr[3]*gpsi_now[ityp][2]/volrt;
        p_b = -rad2*ylmi[3]*gpsi_now[ityp][2]/volrt;
       break;
       case 2:
        p_b = -ylmr[2]*gpsi_now[ityp][2]/volrt;
        p_a = -rad2*ylmr[3]*gpsi_now[ityp][2]/volrt;
        p_c = -rad2*ylmi[3]*gpsi_now[ityp][2]/volrt;
       break;
      }//end switch
      psi_r[2] = 0.0;
      psi_i[2] = (p_c + (p_b+p_a)/rad2)/rad2;
      psi_r[3] = 0.0;
      psi_i[3] = (p_c - (p_b+p_a)/rad2)/rad2;
      psi_r[4] = 0.0;
      psi_i[4] = (p_b-p_a)/rad2;
    }//endif
//  D BAND                                                                  
    if(n_ang[ityp] >= 2){
      psi_r[5] = -ylmr[5]*gpsi_now[ityp][3]/volrt;
      psi_i[5] = 0.0;
      psi_r[6] = -rad2*ylmr[6]*gpsi_now[ityp][3]/volrt;
      psi_i[6] = 0.0;
      psi_r[7] = -rad2*ylmi[6]*gpsi_now[ityp][3]/volrt;
      psi_i[7] = 0.0;
      psi_r[8] = -rad2*ylmr[8]*gpsi_now[ityp][3]/volrt;
      psi_i[8] = 0.0;
      psi_r[9] = -rad2*ylmi[8]*gpsi_now[ityp][3]/volrt;
      psi_i[9] = 0.0;
    }//endif
 
//---------------------------------------------------------------------------
//  structure factor                                                         

    double helr,heli;
    if(k_x[i]==0 && k_y[i]==0 && k_z[i]==0){
      helr = 1.00;
      heli = 0.00;
    }else{
      double dx,dy,dz;
      dx  = x;
      dy  = y;
      dz  = z;
#ifdef _QM_MM_
      dx  = x - cp_box_center[1];
      dy  = y - cp_box_center[2];
      dz  = z - cp_box_center[3];

      double asx,asy,asz;
      asx = dx*hmati_big[1]+dy*hmati_big[4]+dz*hmati_big[7];
      asy = dx*hmati_big[2]+dy*hmati_big[5]+dz*hmati_big[8];
      asz = dx*hmati_big[3]+dy*hmati_big[6]+dz*hmati_big[9];

      double sx,sy,sz;
      sx  = asx - NINT(asx);
      sy  = asy - NINT(asy);
      sz  = asz - NINT(asz);
      dx  = sx*hmat_big[1] + sy*hmat_big[4] + sz*hmat_big[7] + cp_box_center_rel[1];
      dy  = sx*hmat_big[2] + sy*hmat_big[5] + sz*hmat_big[8] + cp_box_center_rel[2];
      dz  = sx*hmat_big[3] + sy*hmat_big[6] + sz*hmat_big[9] + cp_box_center_rel[3];
#endif

      double atemp,btemp,ctemp,arg;
      atemp = dx*hmati[1] + dy*hmati[4] + dz*hmati[7];
      btemp = dx*hmati[2] + dy*hmati[5] + dz*hmati[8];
      ctemp = dx*hmati[3] + dy*hmati[6] + dz*hmati[9];

      arg   = (k_x[i]*atemp + k_y[i]*btemp + k_z[i]*ctemp)*tpi;

      helr  = cos(arg);
      heli  = sin(arg);
    }//endif G=0 

//-------------------------------------------------------------------------
//  construct the coeff and the contribution to the norm 

    gspace_coefs[i].re = helr*psi_r[is] - heli*psi_i[is];
    gspace_coefs[i].im = heli*psi_r[is] + helr*psi_i[is];

    double wght = 2.0;
    if(k_x[i]==0)wght=1.0;
    anorm += (gspace_coefs[i].re*gspace_coefs[i].re+
              gspace_coefs[i].im*gspace_coefs[i].im)*wght;
    if(i%30==0)
      CmiNetworkProgress();
  }//endfor:gSpaceSize

//=========================================================================
// Give this guy norm=2

#ifdef _DEBUG_GEN_AVE_
  PRINTF("state %d seed %d ipart %d is %d ityp %d: %g %g %g : %g : %g\n",
          state_ind,seed,ipart,is,ityp,x,y,z,anorm,vol);
#endif

  anorm  = sqrt(2.0/anorm);
  for(int i=0; i < gSpaceSize; i++){
    gspace_coefs[i].re *= anorm;
    gspace_coefs[i].im *= anorm;
  }//endfor

//==========================================================================
  }// end routine
//==========================================================================


//==========================================================================
// Control the spline fit of the radial wavefunction
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void GEN_WAVE::splin_btrans(double *g,double **gpsi0, double **gpsi1,
                            double **gpsi2, double **gpsi3,
         	            double *gpsi00, int *n_ang_out, char *fname_ps)
//==========================================================================
  {//begin routine
//==========================================================================

  double *r,*rphi; // length nr

  int iii;
  int i,ir,iang_now;
  int iang,nr;
  int n_ang_now;
  double rmax,xx,dr;
  int n_ang1;
  FILE *fp_name_ps;

//========================================================================
// Set up the r's 

  fp_name_ps = cfopen(fname_ps,"r");

    fscanf(fp_name_ps,"%d %lg %d ",&nr,&rmax,&n_ang_now);
    readtoendofline(fp_name_ps);
    readtoendofline(fp_name_ps);

    n_ang_out[0] = n_ang_now;

    r    = (double *) cmalloc(nr*sizeof(double),"splin_btrans") -1;
    rphi = (double *) cmalloc(nr*sizeof(double),"splin_btrans") -1;

    dr = rmax/(double)(nr);

    for(ir=1; ir<= nr; ir++){r[ir] = (double)(ir-1)*dr;}
  
//========================================================================
// Set up the gpsi's 
   
    n_ang1 = n_ang_now + 1;

    for(iang=1; iang <= n_ang1; iang++){

      for(ir=1; ir<= nr; ir++){
        fscanf(fp_name_ps,"%lg %lg ",&xx,&(rphi[ir]));
      }//endfor

      iang_now = iang-1;
      bess_trans(rphi,nr,dr,r,gpsi0[iang],g,iang_now,gpsi00);
      fit_spline(gpsi0[iang],gpsi1[iang],gpsi2[iang],gpsi3[iang],g);

    }//endfor

  fclose(fp_name_ps);

//========================================================================
// free locally assigned memory                                          

  cfree(&(r[1]),"spline-b-trans");
  cfree(&(rphi[1]),"spline-b-trans");

//==========================================================================
  }//end routine
//==========================================================================


//==========================================================================
// Spherical bessel transform of the radial wavefunction
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void GEN_WAVE::bess_trans(double *v_rphi,int nr,double dr,double *r,
                          double *fv_rphi,double *g,int iang,double *gzero)
//==========================================================================
  {//begin routine
//==========================================================================


//  local variables                                                         
   double fpidr,rj0,rj1,rj2,rj3,arg,fpi,pi,tpi;
   int ig,ir;
   int iii;

// -------------------------------------------------------------------------
//  slow spherical bessel fourier transform                                 

   pi  = M_PI;
   tpi = 2.0*pi;
   fpi = 4.0*pi;
   fpidr = fpi*dr;

  for(ig=1; ig <= nsplin; ig++){
    fv_rphi[ig] = 0.0;
  }//endfor

  if(iang == 0){
// l=0 or local ,g=0 ----------------------------------------------------- 
   *gzero = 0.0;
  for(ir=1; ir <= nr; ir++){
    *gzero += fpidr*r[ir]*v_rphi[ir];
  }//endfor

// l=0 or local, g ne 0 ----------------------------------------------------
   for(ig=1; ig <= nsplin; ig++){
     for(ir=2; ir <= nr; ir++){
        arg = r[ir]*g[ig];
        rj0 = sin(arg)/arg*r[ir];
        fv_rphi[ig] += fpidr*rj0*v_rphi[ir];
      }//endfor
   }//endfor
 }//endif

   if(iang == 1){
// l=1, g ne 0 --------------------------------------------------------------
   for(ig=1; ig <= nsplin; ig++){
     for(ir=2 ; ir <= nr; ir++){
       arg = r[ir]*g[ig];
       rj1 = (sin(arg)/arg - cos(arg))/arg*r[ir];
       fv_rphi[ig] +=  fpidr*rj1*v_rphi[ir];
     }//endfor
   }//endfor
  }//endif

   if(iang == 2){
//  l=2, g ne 0 ------------------------------------------------------------
   for(ig=1; ig <= nsplin; ig++){
     for(ir=2; ir <= nr; ir++){  
        arg = r[ir]*g[ig];
        rj2 = ((3.0/(arg*arg)-1.0)*sin(arg)-3.0*cos(arg)/arg)/arg*r[ir];
        fv_rphi[ig] += fpidr*rj2*v_rphi[ir];
     }//endfor
   }//endfor
  }//endif

  if(iang == 3){
//  l=3, g ne 0 ---------------------------------------------------------- 
   for(ig=1; ig <= nsplin; ig++){
     for(ir=2; ir <= nr; ir++){
       arg = r[ir]*g[ig];
       rj3 = ((15.0/(arg*arg) - 6.0)*sin(arg)/arg + 
             (1.0 - 15.0/(arg*arg))*cos(arg))/arg*r[ir];
      fv_rphi[ig] =+  fpidr*rj3*v_rphi[ir];
     }//endfor
   }//endfor
  }//endif


//==========================================================================
  }//end routine
//==========================================================================

//==========================================================================
// Do the spline lookup of the radial wavefunction
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
 void GEN_WAVE::get_gpsi(double g,double **gpsi0,double **gpsi1,
                         double **gpsi2,double **gpsi3,double *gpsi_now,
                         int n_ang_now)
//==========================================================================
  {//begin routine
//==========================================================================


   double partem1,partem2,partem3,partem4,h,h0;
   int iii,iang,i;

   if(g<gmin || g>gmax){
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     PRINTF("g value out of range in spline %g %g %g\n",g,gmin,gmax);
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     EXIT(1);
   }//endif

   iii = int((g-gmin)/dg + 1);
   iii = MIN(iii,nsplin);
   iii = MAX(iii,1);

   h0  = (double)(iii-1)*dg+gmin;
   h = g-h0;

   for(iang =1; iang <= (n_ang_now+1); iang++){
     partem1 = gpsi0[iang][iii];
     partem2 = gpsi1[iang][iii];
     partem3 = gpsi2[iang][iii];
     partem4 = gpsi3[iang][iii];

     gpsi_now[iang] = ((partem4*h+partem3)*h+partem2)*h+partem1;
   }//endfor

//==========================================================================
   }//end routine
//==========================================================================


//==========================================================================
// Spline the radial wavefunctions
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void GEN_WAVE::fit_spline(double *c0i,double *c1i,double *c2i,
                          double *c3i,double *xi)
//==========================================================================
   {//begin routine
//==========================================================================
// temporary vectors                                                       
   double *d,*diag; //length nsplin

// temporary scalars                                                       
   double c0,c1,c2,c3,g,divdf1,divdf3,dx;

// temporary integers                                                      
   int iii;
   int m,mm1,mp1,i,ip1,n;

//malloc local memory 
   d    = (double *) cmalloc(nsplin*sizeof(double ),"fit_spline") -1;
   diag = (double *) cmalloc(nsplin*sizeof(double ),"fit_spline") -1;

//==========================================================================
// fit the spline                                                          

//  1st approximate initial and final derivatives                          
   n = nsplin-1;
   c1i[1]   = (c0i[2]-c0i[1])/(xi[2]-xi[1]);
   c1i[n+1] = (c0i[n+1]-c0i[n])/(xi[n+1]-xi[n]);
   c0 = 0.0;
   c1 = 1.0;
   c2 = 2.0;
   c3 = 3.0;
   diag[1] = c1;
   d[1]    = c0;
  for(m=2; m<= n+1; m++){
    mm1 = m-1;
    d[m] = xi[m]-xi[mm1];
    diag[m] = (c0i[m]-c0i[mm1])/d[m];
  }//endfor
  for(m=2; m<=n; m++){
     mp1 = m+1;
     c1i[m] = c3*(d[m]*diag[mp1]+d[mp1]*diag[m]);
     diag[m] = c2*(d[m]+d[mp1]);
  }//endfor
  for(m=2; m<=n; m++){
    mp1 = m+1;
    mm1 = m-1;
    g = -d[mp1]/diag[mm1];
    diag[m] = diag[m]+g*d[mm1];
    c1i[m] = c1i[m]+g*c1i[mm1];
  }//endfor
  for(m=n; m>=2; m--){
    mp1 = m+1;
    c1i[m] = (c1i[m]-d[m]*c1i[mp1])/diag[m];
  }//endfor

// calculate all other coefficients                                       
  for(i=1; i<=n; i++){
    ip1 = i+1;
    dx = xi[ip1]-xi[i];
    divdf1 = (c0i[ip1]-c0i[i])/dx;
    divdf3 = c1i[i]+c1i[ip1]-c2*divdf1;
    c2i[i] = (divdf1-c1i[i]-divdf3)/dx;
    c3i[i] = divdf3/(dx*dx);
  }//endfor

//==========================================================================
// free locally assigned memory 
  cfree(&(d[1]),"fit_spline");
  cfree(&(diag[1]),"fit_spline");

//==========================================================================
  }//end routine
//==========================================================================


//==========================================================================
// Compute some spherical harmonics
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void GEN_WAVE::get_ylm(double xk,double yk,double zk,double g,
                       double *ylmr,double *ylmi,CPYLM_CONS *ylm_cons)
//========================================================================
  {//begin routine
//========================================================================

 // Local pointers 
  double rt_fpi     = ylm_cons->rt_fpi;
  double rt_thrfpi  = ylm_cons->rt_thrfpi;
  double rt_threpi  = ylm_cons->rt_threpi;
  double hrt_fivfpi = ylm_cons->hrt_fivfpi;
  double rt_fiftepi = ylm_cons->rt_fiftepi;
  double hrt_sevfpi = ylm_cons->hrt_sevfpi;
  double hrt_toepi  = ylm_cons->hrt_toepi;
  double hrt_ohffpi = ylm_cons->hrt_ohffpi;
  double hrt_tfepi  = ylm_cons->hrt_tfepi;

//             Local variable declarations                                

  double y00r,y00i,y10r,y10i,y11r,y11i,y20r,y20i,y21r,y21i;
  double y22r,y22i,y30r,y30i,y31r,y31i,y32r,y32i,y33r,y33i;
  double ctheta,stheta,cphi,sphi,c2phi,s2phi,c3phi,s3phi;
  double xydist;

//==========================================================================
// I) Calculate polar angles of the vector g                                

  ctheta = zk/g;
  stheta = sqrt(xk*xk + yk*yk)/g;
  xydist = sqrt(xk*xk + yk*yk);
  if(xydist==0){
    cphi = 1.0;
    sphi = 0.0;
  }else{
    cphi = xk/xydist;
    sphi = yk/xydist;
  }//endif
  c2phi = cphi*cphi - sphi*sphi;
  s2phi = 2.0*cphi*sphi;
  c3phi = cphi*c2phi - sphi*s2phi;
  s3phi = cphi*s2phi + c2phi*sphi;

//==========================================================================
// I.i) l=0                                                                 

  y00r    = rt_fpi;
  y00i    = 0.0;
  ylmr[1] = y00r;
  ylmi[1] = y00i;

//==========================================================================
// I.ii) l=1 (phi derivatives have stheta divided out)                      

  y10r        =  rt_thrfpi*ctheta;
  y10i        = 0.0;

  y11r        =  rt_threpi*stheta*cphi;
  y11i        =  rt_threpi*stheta*sphi;

  ylmr[2] = y10r;
  ylmi[2] = y10i;
  ylmr[3] = y11r;
  ylmi[3] = y11i;
  ylmr[4] = y11r;
  ylmi[4] = -y11i;

//==========================================================================
// I.iii) l=2 (phi derivatives have stheta divided out)                     

  y20r    = hrt_fivfpi*(3.0*ctheta*ctheta - 1.0);
  y20i    = 0.0;

  y21r    = rt_fiftepi*stheta*ctheta*cphi;
  y21i    = rt_fiftepi*stheta*ctheta*sphi;

  y22r    = 0.50*rt_fiftepi*stheta*stheta*c2phi;
  y22i    = 0.50*rt_fiftepi*stheta*stheta*s2phi;

  ylmr[5] = y20r;
  ylmi[5] = y20i;
  ylmr[6] = y21r;
  ylmi[6] = y21i;
  ylmr[7] = y21r;
  ylmi[7] = -y21i;
  ylmr[8] = y22r;
  ylmi[8] = y22i;
  ylmr[9] = y22r;
  ylmi[9] = -y22i;

//==========================================================================
// I.iv) l=3 (phi derivatives have stheta divided out)                      

  y30r        = hrt_sevfpi*(5.0*ctheta*ctheta*ctheta - 3.0*ctheta);
  y30i        = 0.0;

  y31r        = hrt_toepi*stheta*(5.0*ctheta*ctheta - 1.0)*cphi;
  y31i        = hrt_toepi*stheta*(5.0*ctheta*ctheta - 1.0)*sphi;

  y32r        = hrt_ohffpi*stheta*stheta*ctheta*c2phi;
  y32i        = hrt_ohffpi*stheta*stheta*ctheta*s2phi;

  y33r        = hrt_tfepi*stheta*stheta*stheta*c3phi;
  y33i        = hrt_tfepi*stheta*stheta*stheta*s3phi;

  ylmr[10] = y30r;
  ylmi[10] = y30i;
  ylmr[11] = y31r;
  ylmi[11] = y31i;
  ylmr[12] = y31r;
  ylmi[12] = -y31i;
  ylmr[13] = y32r;
  ylmi[13] = y32i;
  ylmr[14] = y32r;
  ylmi[14] = -y32i;
  ylmr[15] = y33r;
  ylmi[15] = y33i;
  ylmr[16] = y33r;
  ylmi[16] = -y33i;

//--------------------------------------------------------------------------
  }//end routine
//==========================================================================

