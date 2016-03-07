//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                                                                          
//                         PI_MD:                                           
//             The future of simulation technology                          
//             ------------------------------------                         
//                     Module: write_gen_header                             
//                                                                          
//                                                                          
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//               Header:                                                    

#include "standard_include.h"
#include "../../../include/Atoms.h"
#include "../../../include/debug_flags.h"
#include "../class_defs/allclass_gen.h"
#include "../class_defs/allclass_cp.h"
#include "../class_defs/allclass_mdintegrate.h"
#include "../class_defs/allclass_mdatoms.h"
#include "../class_defs/allclass_mdinter.h"
#include "../class_defs/allclass_mdintra.h"
#include "../class_defs/ATOM_OPERATIONS/class_atomoutput.h"

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void ATOMOUTPUT::ctrl_piny_output(int itime,int natm,int len_nhc,int pi_beads_in, 
    int myid, Atom *atoms, AtomNHC *atomsNHC,
    int *iwrite_atm_ret,int output_on,
    int itemper, int ibead)
  //==========================================================================
{//begin routine
  //=======================================================================

  GENERAL_DATA *general_data = GENERAL_DATA::get();
#include "../class_defs/allclass_strip_gen.h"
  int low_lim_par      = genfilenames->low_lim_par-1;
  int high_lim_par     = genfilenames->high_lim_par;
  int iwrite_confp     = genfilenames->iwrite_confp;
  int iwrite_par_confp = genfilenames->iwrite_par_confp;
  int iwrite_dump      = genfilenames->iwrite_dump;

  char *dump_dir       = genfilenames->atm_crd_dir_out;
  char *cpname         = genfilenames->cpname;
  char *cpparname      = genfilenames->cpparname;
  char *dname          = genfilenames->dname;
  int ntime            = gentimeinfo->ntime; //correct number of time steps
  int pi_beads         = 1;

  char temp_ext[1000];

  if(pi_beads_in!=1){
    PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("   Input value of pi_ beads be unity %d\n",pi_beads_in);
    PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    EXIT(1);
  }//endif

  //==========================================================================
  // Write the file now.

  int iwrite_atm=0;
  if(itime>0 && output_on==1){

    if( (itime % iwrite_confp)==0 ){
      iwrite_atm++;
      if(myid==0){
        int low = 0; int high = natm;
        sprintf (temp_ext,"%s/Bead.%d_Temper.%d/%s",dump_dir,ibead,itemper,cpname);
        write_atom_output_conf(low,high,pi_beads,atoms,temp_ext);
      }//endif
    }//endif

    if( (itime % iwrite_par_confp)==0 && low_lim_par<high_lim_par){
      iwrite_atm++;
      if(myid==0){
        sprintf (temp_ext,"%s/Bead.%d_Temper.%d/%s",dump_dir,ibead,itemper,cpparname);
        write_atom_output_conf(low_lim_par,high_lim_par,pi_beads,atoms,temp_ext);
      }//endif
    }//endif

    if( (itime % iwrite_dump)==0 || itime==ntime || itime==ntime-1 ){
      iwrite_atm++;
      if(myid==0){
        sprintf (temp_ext,"%s/Bead.%d_Temper.%d/%s",dump_dir,ibead,itemper,dname);
        write_atom_output_dump(natm,len_nhc,pi_beads,itime,atoms,atomsNHC,temp_ext);
      }//endif
    }//endif

  }//endif


  (*iwrite_atm_ret)=iwrite_atm;
  //==========================================================================
}//end routine
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void ATOMOUTPUT::write_atom_output_dump(int natm, int len_nhc,int pi_beads,
    int itime,Atom *atoms,AtomNHC *atomsNHC,
    char *dname)
  //==========================================================================
{//begin routine
  //=======================================================================

  GENERAL_DATA *general_data = GENERAL_DATA::get();
  MDATOMS      *mdatoms      = MDATOMS::get();
  MDINTEGRATE  *mdintegrate  = MDINTEGRATE::get();

#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_mdatoms.h"
#include "../class_defs/allclass_strip_mdintegrate.h"

  int n,i,ip,i1,j;
  int num_nhc       = mdtherm_info->num_nhc;
  int iextended_on  = mdtherm_info->iextended_on;
  double *hmat      = gencell->hmat;
  double *hmat_ewd  = gencell->hmat_ewd;
  double vol        = gencell->vol;
  int *iatm_atm_typ = mdatom_maps->iatm_atm_typ;
  int *iatm_res_typ = mdatom_maps->iatm_res_typ;
  int *iatm_mol_typ = mdatom_maps->iatm_mol_typ;
  int *iatm_mol_num = mdatom_maps->iatm_mol_num;
  int *iatm_res_num = mdatom_maps->iatm_res_num;
  NAME *atm_typ     = mdatom_maps->atm_typ;
  NAME *res_typ     = mdatom_maps->res_typ;
  NAME *mol_typ     = mdatom_maps->mol_typ;

  //==========================================================================

  FILE *fp = cfopen(dname,"o");

  //==========================================================================
  // I)Atm positions                                                      


  fprintf(fp,"natm_tot restart_typ itime pi_beads\n");
  fprintf(fp,"%d restart_all %d %d\n",natm,itime,pi_beads);
  fprintf(fp,"atm pos, atm_typ, mol_typ mol_num\n");
  for(ip=1;ip<=pi_beads;ip++){
    for(i=1,i1=0;i<=natm;i++,i1++){
      fprintf(fp,"%.13g %.13g %.13g %s %s %s %d %d\n",
          atoms[i1].xold,atoms[i1].yold,atoms[i1].zold,
          atm_typ[iatm_atm_typ[i]],res_typ[iatm_res_typ[i]],
          mol_typ[iatm_mol_typ[i]],iatm_mol_num[i],iatm_res_num[i]);
    }//endfor
  }//endfor

  //==========================================================================
  //  II)Cell shape                                                         

  fprintf(fp,"h matrix\n");
  fprintf(fp,"%.13g %.13g %.13g\n",hmat[1],hmat[4],hmat[7]);
  fprintf(fp,"%.13g %.13g %.13g\n",hmat[2],hmat[5],hmat[8]);
  fprintf(fp,"%.13g %.13g %.13g\n",hmat[3],hmat[6],hmat[9]);
  fprintf(fp,"h matrix for Ewald setup\n");
  fprintf(fp,"%.13g %.13g %.13g\n",hmat_ewd[1],hmat_ewd[4],hmat_ewd[7]);
  fprintf(fp,"%.13g %.13g %.13g\n",hmat_ewd[2],hmat_ewd[5],hmat_ewd[8]);
  fprintf(fp,"%.13g %.13g %.13g\n",hmat_ewd[3],hmat_ewd[6],hmat_ewd[9]);
  fprintf(fp,"1/3 log( Vol(t)/Vol(0) )\n");
  fprintf(fp,"0.0\n");

  //==========================================================================
  //  III)Atm and Atm NHC Velocities                                         

  fprintf(fp,"atm vel\n");
  for(ip=1;ip<=pi_beads;ip++){
    for(i=0;i<natm;i++){
      fprintf(fp,"%.13g %.13g %.13g\n",atoms[i].vxold,atoms[i].vyold,atoms[i].vzold);
    }//endfor
  }//endfor

  fprintf(fp,"number of atm nhc, length of nhc\n");
  fprintf(fp,"0 %d\n",len_nhc);
  fprintf(fp,"atm nhc velocities\n");

  //==========================================================================
  // IV)Vol and Vol NHC Velocities                                            

  fprintf(fp,"vol velocities\n");
  fprintf(fp,"0.0 0.0 0.0\n");
  fprintf(fp,"0.0 0.0 0.0\n");
  fprintf(fp,"0.0 0.0 0.0\n");
  fprintf(fp,"log(vol) velocity\n");
  fprintf(fp,"0.0\n");
  fprintf(fp,"vol nhc velocities\n");
  for(i=1;i<=len_nhc;i++){
    fprintf(fp,"0.0\n");
  }//endfor

  //==========================================================================
  // V)Misc                                                               

  fprintf(fp,"dt=%.13g\n",gentimeinfo->dt);
  fprintf(fp,"nfree=%d\n",mdclatoms_info->nfree);
  fprintf(fp,"nve=%d nvt=%d npt_i=%d npt_f=%d nst=%d\n",
      genensopts->nve,genensopts->nvt,genensopts->npt_i,
      genensopts->npt_f,genensopts->nst);
  fprintf(fp,"nbond_free=0 nbend_free=0 ntors_free=0\n");
  fprintf(fp,"t_ext=%.13g,pext=%.13g,stens_ext=%.13g\n",
      genstatepoint->t_ext,genstatepoint->pext,genstatepoint->stens_ext);

  //==========================================================================
  // VI)Close                                                               

  fflush(fp);
  fclose(fp); 

  //==========================================================================
}//end routine
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void ATOMOUTPUT::write_atom_output_conf(int low, int high, int pi_beads, 
    Atom *atoms,char *fname)
  //==========================================================================
{//begin routine
  //==========================================================================

  GENERAL_DATA *general_data = GENERAL_DATA::get();

#include "../class_defs/allclass_strip_gen.h"

  int n,i,ip;
  int iwrite_conf_binary = genfilenames->iwrite_conf_binary;
  double *hmat           = gencell->hmat;

  FILE *fp  = cfopen(fname,"a");

  //=====================================================================
  // I) Write to the atm config file in ascii 

  if(iwrite_conf_binary==0){
    for(ip=1;ip<=pi_beads;ip++){
      for(i=low;i<high;i++){
        fprintf(fp,"%.12g  %.12g  %.12g\n",atoms[i].xold,atoms[i].yold,atoms[i].zold);
      }//endfor
    }//endfor
    for(i=1;i<=9;i+=3){
      fprintf(fp,"%.13g %.13g %.13g\n",hmat[i],hmat[(i+1)],hmat[(i+2)]);
    }//endfor
    fflush(fp);
  }//endif

  //=====================================================================
  // II) Write to the atm config file in binary

  if(iwrite_conf_binary==1){
    n=1;
    for(ip=1;ip<=pi_beads;ip++){
      for(i=low;i<high;i++){ 
        double x=atoms[i].xold;
        double y=atoms[i].yold;
        double z=atoms[i].zold;
        fwrite(&x,sizeof(double),n,fp);
        fwrite(&y,sizeof(double),n,fp);
        fwrite(&z,sizeof(double),n,fp);
      }//endfor
    }//endfor
    for(i=1;i<=9;i+=3){ 
      fwrite(&hmat[i],sizeof(double),n,fp);
      fwrite(&hmat[i+1],sizeof(double),n,fp);
      fwrite(&hmat[i+2],sizeof(double),n,fp);
    }//endfor
  }//endif

  fclose(fp);

  //==========================================================================
}//end routine
//==========================================================================

