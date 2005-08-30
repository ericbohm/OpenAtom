
class complex;
class AtomsGrp;
class CPATOM_MAPS;
class CPCOEFFS_INFO;
class MDATOM_MAPS;


//==========================================================================
//    Class used to generate initial guess for wavefunction 
//             
//             


#ifndef _GEN_WAVE_
#define _GEN_WAVE_

class GEN_WAVE{
  public:
   GEN_WAVE();
  ~GEN_WAVE();
   void fill_gw_gpsi(CPATOM_MAPS * ,CPCOEFFS_INFO *,
		     int ,double ,double ,
                     MDATOM_MAPS *mdatom_maps,NAME *vps_name);
   void create_coefs(const int *k_x,const int *k_y,const int *k_z,
                     int gSpaceSize,int nstate_in,
                     complex *gspace_coefs);

   int natm_typ_cp,nsplin;
   int *iatm_state_str;
   int *iatm_state_end;
   int *iatm_atm_typ_cp;
   int *n_ang;
   int  *nstate_up_atm,*nstate_dn_atm;
   double dg;
   double ***gpsi0,***gpsi1,***gpsi2,***gpsi3;
   double **gpsi_now;    //(maxatmtyp,maxatmtyp)
   double *gpsi00;       //(maxatmtyp)

//-------------------------------------------------------------------------
#ifdef PUP_ON
 //pupping
  void pup(PUP::er &p){
    //pupping ints
      p | natm_typ_cp;
      p | nsplin;
    //pupping dbles
      p | dg;
    //pupping arrays
      if(natm_typ_cp>0){
        pup1d_int(p,&n_ang,natm_typ_cp);
        pup1d_int(p,&nstate_up_atm,natm_typ_cp);
        pup1d_int(p,&nstate_dn_atm,natm_typ_cp);
        pup1d_dbl(p,&gpsi00,natm_typ_cp);

        pup2d_dbl(p,&gpsi_now,natm_typ_cp,3);

        pup3d_dbl(p,&gpsi0,natm_typ_cp,3,nsplin); 
        pup3d_dbl(p,&gpsi1,natm_typ_cp,3,nsplin); 
        pup3d_dbl(p,&gpsi2,natm_typ_cp,3,nsplin); 
        pup3d_dbl(p,&gpsi3,natm_typ_cp,3,nsplin); 
      }//endif
#ifdef _PARALLEL_DEBUG_        
    if (p.isUnpacking())
     state_class_out ();
#endif     
  } // end pup
#endif
//----------------------------------------------------------------------
  void state_class_out(){
     int i,j,k;
     char fileName [255];
     sprintf (fileName, "%d_gen_wave.state", CkMyPe());
     FILE *fp;  fp = fopen(fileName,"w");

    // ints
     fprintf(fp,"natm_typ_cp %d\n",natm_typ_cp);
     fprintf(fp,"dg %lg\n",dg);
     if(natm_typ_cp>0){
      for(i=1;i<=natm_typ_cp;i++){fprintf(fp,"n_ang[%d] %d\n",i,n_ang[i]);}
      for(i=1;i<=natm_typ_cp;i++){fprintf(fp,"nstate_up_atm[%d] %d\n",
                                               i,nstate_up_atm[i]);}
      for(i=1;i<=natm_typ_cp;i++){fprintf(fp,"nstate_dn_atm[%d] %d\n",
                                               i,nstate_dn_atm[i]);}
      for(i=1;i<=natm_typ_cp;i++){fprintf(fp,"gpsi00[%d] %lg\n",i,gpsi00[i]);}

      for(i=1; i<= natm_typ_cp;i++){
      for(j=1; j<= 3;j++){
	fprintf(fp,"gpsi_now[%d][%d] %lg \n",i,j,gpsi_now[i][j]);
      }}

      for(i=1; i<= natm_typ_cp;i++){
      for(j=1; j<= 3;j++){
      for(k=1; k<= nsplin; k++){
	fprintf(fp,"gpsi0[%d][%d][%d] %lg\n",i,j,k,gpsi0[i][j][k]);
	fprintf(fp,"gpsi1[%d][%d][%d] %lg\n",i,j,k,gpsi1[i][j][k]);
	fprintf(fp,"gpsi2[%d][%d][%d] %lg\n",i,j,k,gpsi2[i][j][k]);
	fprintf(fp,"gpsi3[%d][%d][%d] %lg\n",i,j,k,gpsi3[i][j][k]);
      }}}//endfor

   }//endif
   fclose(fp);
  }// end routine 

//-------------------------------------------------------------------------
  }; //GEN_WAVE
//==========================================================================

#ifdef PUP_ON
PUPmarshall(GEN_WAVE);
#endif

#endif
//==========================================================================
