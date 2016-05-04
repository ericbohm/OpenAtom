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
void ATOMOUTPUT::initialize_piny_output()
  //==========================================================================
{ //begin routine
  //=======================================================================
  // Piny read onlys 

  CP           *cp           = CP::get();
  GENERAL_DATA *general_data = GENERAL_DATA::get();
  MDATOMS      *mdatoms      = MDATOMS::get();
  MDINTEGRATE  *mdintegrate  = MDINTEGRATE::get();

#include "../class_defs/allclass_strip_cp.h"
#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_mdatoms.h"
#include "../class_defs/allclass_strip_mdintegrate.h"

  NAME file_typ;
  int ibinary          = genfilenames->iwrite_conf_binary;
  int low_lim_par      = genfilenames->low_lim_par;
  int high_lim_par     = genfilenames->high_lim_par;
  int iwrite_confp     = genfilenames->iwrite_confp;
  int iwrite_par_confp = genfilenames->iwrite_par_confp;

  char *dname          = genfilenames->dname;             // restart file (both pos and vel)
  char *cpname         = genfilenames->cpname;            // position trajectory file
  char *cpparname      = genfilenames->cpparname;         // partial position trajectory file
  char *cvname         = genfilenames->cvname;            // velocity trajectory file
  char *dump_dir       = genfilenames->atm_crd_dir_out;   // directory where the files go

  int ntemper         = gensimopts->ntemper;
  int pi_beads        = gensimopts->pi_beads;

  char temp_ext[1000];

  int itemper, ibead;

  //=======================================================================
  // Invoke the header for the appropriate files

  for(ibead=0;ibead<pi_beads;ibead++){
    for(itemper=0;itemper<ntemper;itemper++){

      sprintf (temp_ext,"%s/Bead.%d_Temper.%d/%s",dump_dir,ibead,itemper,dname);
      FILE *fp = cfopen(temp_ext,"w");
      fclose(fp);

      sprintf (temp_ext,"%s/Bead.%d_Temper.%d/%s",dump_dir,ibead,itemper,cpname);
      strcpy(file_typ,"pos_file");
      write_gen_header_cp(ibinary,iwrite_confp,file_typ,temp_ext);

      sprintf (temp_ext,"%s/Bead.%d_Temper.%d/%s",dump_dir,ibead,itemper,cvname);
      strcpy(file_typ,"pos_file");
      write_gen_header_cp(ibinary,iwrite_confp,file_typ,temp_ext);

      if(low_lim_par<=high_lim_par){
        sprintf (temp_ext,"%s/Bead.%d_Temper.%d/%s",dump_dir,ibead,itemper,cpparname);
        strcpy(file_typ,"par_file");
        write_gen_header_cp(ibinary,iwrite_par_confp,file_typ,temp_ext);
      }//endif

    }}//endfor

  //==========================================================================
}//end routine
//=======================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void ATOMOUTPUT::write_gen_header(int ibinary,int iwrite_freq,
    NAME file_typ, char *fname)
  //==========================================================================
{ //begin routine
  //=======================================================================
  // Piny read onlys and local variables

  CP           *cp           = CP::get();
  GENERAL_DATA *general_data = GENERAL_DATA::get();
  MDATOMS      *mdatoms      = MDATOMS::get();
  MDINTEGRATE  *mdintegrate  = MDINTEGRATE::get();

#include "../class_defs/allclass_strip_cp.h"
#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_mdatoms.h"
#include "../class_defs/allclass_strip_mdintegrate.h"

  int iwrite_freq_now = iwrite_freq;
  int ibinary_now=ibinary;
  int i,n,ifound;
  int low,high;
  int iwarn=0;
  NAME name_scr;
  strcpy(name_scr,file_typ);

  //=======================================================================
  // I) Formatted file                                                       

  if(ibinary==0){
    FILE *fp = cfopen(fname,"w");
    fprintf(fp,
        "%d leanCP %s %d %.10g %d %d %d %d %d %d %d %d %d %d %d %d %.10g %.10g %.10g %d\n",
        ibinary,
        file_typ,
        gentimeinfo->ntime,
        gentimeinfo->dt,
        iwrite_freq_now,
        mdclatoms_info->pi_beads,
        mdatom_maps->nmol_typ,
        mdatom_maps->nres_typ,
        mdatom_maps->natm_typ,
        mdclatoms_info->natm_tot,
        mdclatoms_info->nfree,
        genensopts->nve,
        genensopts->nvt,
        genensopts->npt_i,
        genensopts->npt_f,
        genensopts->nst,
        genstatepoint->t_ext,
        genstatepoint->pext,
        genstatepoint->stens_ext,
        gencell->iperd);
    if(strcasecmp(file_typ,"ins_file")==0){
      fprintf(fp,"%d %d %d %d\n",
          mdtherm_info->num_nhc,mdtherm_info->len_nhc,
          mdtherm_info_bead->num_nhc,mdtherm_info_bead->len_nhc);
    }//endif

    if(strcasecmp(file_typ,"par_file")==0){
      fprintf(fp,"%d %d \n",
          genfilenames->low_lim_par,
          genfilenames->high_lim_par);
    }//endif

    if(strcasecmp(file_typ,"vel_file")==0||
        strcasecmp(file_typ,"pos_file")==0||
        strcasecmp(file_typ,"cen_file")==0||
        strcasecmp(file_typ,"par_file")==0){
      for(i=1;i<=mdatom_maps->nmol_typ;i++){
        fprintf(fp,"%s\n",mdatom_maps->mol_typ[i]);
      }//endfor
      for(i=1;i<=mdatom_maps->nres_typ;i++){
        fprintf(fp,"%s\n",mdatom_maps->res_typ[i]);
      }//endfor
      for(i=1;i<=mdatom_maps->natm_typ;i++){
        fprintf(fp,"%s\n",mdatom_maps->atm_typ[i]);
      }//endfor
      low = 1;high=mdclatoms_info->natm_tot;
      if(strcasecmp(file_typ,"par_file")==0){
        low  = genfilenames->low_lim_par;
        high = genfilenames->high_lim_par;
      }//endif
      for(i=low;i<=high;i++){
        fprintf(fp,"%g %g %d %d %d %d %d\n",
            mdclatoms_info->mass[i],
            mdclatoms_info->q[i],
            mdatom_maps->iatm_mol_typ[i],
            mdatom_maps->iatm_mol_num[i],
            mdatom_maps->iatm_res_typ[i],
            mdatom_maps->iatm_res_num[i],
            mdatom_maps->iatm_atm_typ[i]);
      }//endfor
    }//endif
    fflush(fp);
    fclose(fp);
  }//endif

  //=======================================================================
  // II) Unformatted file                                                    

  if(ibinary==1){
    FILE *fp = cfopen(fname,"w");
    n=1;
    fwrite(&ibinary_now,sizeof(int),n,fp);
    fwrite(&name_scr,sizeof(NAME),n,fp); 
    fwrite(&name_scr,sizeof(NAME),n,fp); 
    fwrite(&(gentimeinfo->ntime),sizeof(int),n,fp);
    fwrite(&(gentimeinfo->dt),sizeof(double),n,fp);
    fwrite(&iwrite_freq_now,sizeof(int),n,fp);
    fwrite(&(mdclatoms_info->pi_beads),sizeof(int),n,fp);
    fwrite(&(mdatom_maps->nmol_typ),sizeof(int),n,fp);
    fwrite(&(mdatom_maps->nres_typ),sizeof(int),n,fp);
    fwrite(&(mdatom_maps->natm_typ),sizeof(int),n,fp);
    fwrite(&(mdclatoms_info->natm_tot),sizeof(int),n,fp);
    fwrite(&(mdclatoms_info->nfree),sizeof(int),n,fp);
    fwrite(&(genensopts->nve),sizeof(int),n,fp);
    fwrite(&(genensopts->nvt),sizeof(int),n,fp);
    fwrite(&(genensopts->npt_i),sizeof(int),n,fp);
    fwrite(&(genensopts->npt_f),sizeof(int),n,fp);
    fwrite(&(genensopts->nst),sizeof(int),n,fp);
    fwrite(&(genstatepoint->t_ext),sizeof(double),n,fp);
    fwrite(&(genstatepoint->pext),sizeof(double),n,fp);
    fwrite(&(genstatepoint->stens_ext),sizeof(double),n,fp);
    fwrite(&(gencell->iperd),sizeof(int),n,fp);

    if(strcasecmp(file_typ,"ins_file")==0){
      fwrite(&(mdtherm_info->num_nhc),sizeof(int),n,fp);
      fwrite(&(mdtherm_info->len_nhc),sizeof(int),n,fp);
      fwrite(&(mdtherm_info_bead->num_nhc),sizeof(int),n,fp);
      fwrite(&(mdtherm_info_bead->len_nhc),sizeof(int),n,fp);
    }//endif

    if(strcasecmp(file_typ,"par_file")==0){
      fwrite(&(genfilenames->low_lim_par),sizeof(int),n,fp);
      fwrite(&(genfilenames->high_lim_par),sizeof(int),n,fp);
    }//endif

    if(strcasecmp(file_typ,"vel_file")==0||
        strcasecmp(file_typ,"pos_file")==0||
        strcasecmp(file_typ,"cen_file")==0||
        strcasecmp(file_typ,"par_file")==0){
      for(i=1;i<=mdatom_maps->nmol_typ;i++){
        strcpy(name_scr,mdatom_maps->mol_typ[i]);
        fwrite(&(name_scr),sizeof(NAME),n,fp);
      }//endfor
      for(i=1;i<=mdatom_maps->nres_typ;i++){
        strcpy(name_scr,mdatom_maps->res_typ[i]);
        fwrite(&(name_scr),sizeof(NAME),n,fp);
      }//endfor
      for(i=1;i<=mdatom_maps->natm_typ;i++){
        strcpy(name_scr,mdatom_maps->atm_typ[i]);
        fwrite(&(name_scr),sizeof(NAME),n,fp);
      }//endfor
      low = 1;high=mdclatoms_info->natm_tot;
      if(strcasecmp(file_typ,"par_file")==0){
        low  = genfilenames->low_lim_par;
        high = genfilenames->high_lim_par;
      }//endif
      n = 1;
      for(i=low;i<=high;i++){
        fwrite(&(mdclatoms_info->mass[i]),sizeof(double),n,fp);
        fwrite(&(mdclatoms_info->q[i]),sizeof(double),n,fp);
        fwrite(&(mdatom_maps->iatm_mol_typ[i]),sizeof(int),n,fp);
        fwrite(&(mdatom_maps->iatm_mol_num[i]),sizeof(int),n,fp);
        fwrite(&(mdatom_maps->iatm_res_typ[i]),sizeof(int),n,fp);
        fwrite(&(mdatom_maps->iatm_res_num[i]),sizeof(int),n,fp);
        fwrite(&(mdatom_maps->iatm_atm_typ[i]),sizeof(int),n,fp);
      }//endfor
    }//endif
    fclose(fp);
  }//endif

  //==========================================================================
}//end routine
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void ATOMOUTPUT::write_gen_header_cp(int ibinary,int iwrite_freq,
    NAME file_typ,char *fname)
  //==========================================================================
{//begin routine
  //=======================================================================

  CP           *cp           = CP::get();
  GENERAL_DATA *general_data = GENERAL_DATA::get();
  MDATOMS      *mdatoms      = MDATOMS::get();
  MDINTEGRATE  *mdintegrate  = MDINTEGRATE::get();
#include "../class_defs/allclass_strip_cp.h"
#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_mdatoms.h"
#include "../class_defs/allclass_strip_mdintegrate.h"

  int iwrite_freq_now = iwrite_freq;
  int ibinary_now     = ibinary;
  int i,n,low,high,ifound;
  NAME name_scr;
  strcpy(name_scr,file_typ);

  //=======================================================================
  // I) Formatted file                                                       

  if(ibinary==0){
    FILE *fp = cfopen(fname,"w");
    fprintf(fp,
        "%d leanCP %s %d %.10g %d %d %d %d %d %d %d %d %d %d %d %d %.10g %.10g %.10g %d\n",
        ibinary,
        file_typ,
        gentimeinfo->ntime,
        gentimeinfo->dt,
        iwrite_freq_now,
        mdclatoms_info->pi_beads,
        mdatom_maps->nmol_typ,
        mdatom_maps->nres_typ,
        mdatom_maps->natm_typ,
        mdclatoms_info->natm_tot,
        mdclatoms_info->nfree,
        genensopts->nve,
        genensopts->nvt,
        genensopts->npt_i,
        genensopts->npt_f,
        genensopts->nst,
        genstatepoint->t_ext,
        genstatepoint->pext,
        genstatepoint->stens_ext,
        gencell->iperd);

    fprintf(fp,"%d %d %d %.10g %.10g %s %s %s %d %d %d %d %d\n",
        cpcoeffs_info->ncoef, 
        cpcoeffs_info->nstate_up, 
        cpcoeffs_info->nstate_dn, 
        cpcoeffs_info->ecut, 
        cppseudo->gga_cut,
        cppseudo->vxc_typ,
        cppseudo->ggax_typ,
        cppseudo->ggac_typ,
        cpopts->cp_lda,
        cpopts->cp_lsda,
        cpopts->cp_sic,
        cpopts->cp_gga,
        cpopts->cp_nonint);

    if(strcasecmp(file_typ,"ins_file")==0){
      fprintf(fp,"%d %d %d %d\n",
          mdtherm_info->num_nhc,mdtherm_info->len_nhc,
          mdtherm_info_bead->num_nhc,mdtherm_info_bead->len_nhc);
    }//endif

    if(strcasecmp(file_typ,"par_file")==0){
      fprintf(fp,"%d %d \n",
          genfilenames->low_lim_par,
          genfilenames->high_lim_par);
    }//endif

    if(strcasecmp(file_typ,"vel_file")==0||
        strcasecmp(file_typ,"pos_file")==0||
        strcasecmp(file_typ,"cen_file")==0||
        strcasecmp(file_typ,"par_file")==0){
      for(i=1;i<=mdatom_maps->nmol_typ;i++){
        fprintf(fp,"%s\n",mdatom_maps->mol_typ[i]);
      }//endfor
      for(i=1;i<=mdatom_maps->nres_typ;i++){
        fprintf(fp,"%s\n",mdatom_maps->res_typ[i]);
      }//endfor
      for(i=1;i<=mdatom_maps->natm_typ;i++){
        fprintf(fp,"%s\n",mdatom_maps->atm_typ[i]);
      }//endfor
      low = 1;high=mdclatoms_info->natm_tot;
      if(strcasecmp(file_typ,"par_file")==0){
        low  = genfilenames->low_lim_par;
        high = genfilenames->high_lim_par;
      }//endif
      for(i=low;i<=high;i++){
        fprintf(fp,"%g %g %d %d %d %d %d\n",
            mdclatoms_info->mass[i],
            mdclatoms_info->q[i],
            mdatom_maps->iatm_mol_typ[i],
            mdatom_maps->iatm_mol_num[i],
            mdatom_maps->iatm_res_typ[i],
            mdatom_maps->iatm_res_num[i],
            mdatom_maps->iatm_atm_typ[i]);
      }//endfor
    }//endif
    fflush(fp);
    fclose(fp);
  }//endif

  //=======================================================================
  // II) Unformatted file                                                    

  if(ibinary==1){
    FILE *fp = cfopen(fname,"w");
    n=1;
    fwrite(&ibinary_now,sizeof(int),n,fp);
    fwrite(&name_scr,sizeof(NAME),n,fp);
    fwrite(&name_scr,sizeof(NAME),n,fp);
    fwrite(&(gentimeinfo->ntime),sizeof(int),n,fp);
    fwrite(&(gentimeinfo->dt),sizeof(double),n,fp);
    fwrite(&iwrite_freq_now,sizeof(int),n,fp);
    fwrite(&(mdclatoms_info->pi_beads),sizeof(int),n,fp);
    fwrite(&(mdatom_maps->nmol_typ),sizeof(int),n,fp);
    fwrite(&(mdatom_maps->nres_typ),sizeof(int),n,fp);
    fwrite(&(mdatom_maps->natm_typ),sizeof(int),n,fp);
    fwrite(&(mdclatoms_info->natm_tot),sizeof(int),n,fp);
    fwrite(&(mdclatoms_info->nfree),sizeof(int),n,fp);
    fwrite(&(genensopts->nve),sizeof(int),n,fp);
    fwrite(&(genensopts->nvt),sizeof(int),n,fp);
    fwrite(&(genensopts->npt_i),sizeof(int),n,fp);
    fwrite(&(genensopts->npt_f),sizeof(int),n,fp);
    fwrite(&(genensopts->nst),sizeof(int),n,fp);
    fwrite(&(genstatepoint->t_ext),sizeof(double),n,fp);
    fwrite(&(genstatepoint->pext),sizeof(double),n,fp);
    fwrite(&(genstatepoint->stens_ext),sizeof(double),n,fp);
    fwrite(&(gencell->iperd),sizeof(int),n,fp);

    fwrite(&(cpcoeffs_info->ncoef),sizeof(int),n,fp);
    fwrite(&(cpcoeffs_info->nstate_up),sizeof(int),n,fp);
    fwrite(&(cpcoeffs_info->nstate_dn),sizeof(int),n,fp);
    fwrite(&(cpcoeffs_info->ecut),sizeof(double),n,fp);
    fwrite(&(cppseudo->gga_cut),sizeof(double),n,fp);
    n = MAXWORD;
    fwrite(&(cppseudo->vxc_typ),sizeof(char),n,fp);
    fwrite(&(cppseudo->ggax_typ),sizeof(char),n,fp);
    fwrite(&(cppseudo->ggac_typ),sizeof(char),n,fp);
    n = 1;
    fwrite(&(cpopts->cp_lda),sizeof(int),n,fp);
    fwrite(&(cpopts->cp_lsda),sizeof(int),n,fp);
    fwrite(&(cpopts->cp_sic),sizeof(int),n,fp);
    fwrite(&(cpopts->cp_gga),sizeof(int),n,fp);
    fwrite(&(cpopts->cp_nonint),sizeof(int),n,fp);

    if(strcasecmp(file_typ,"ins_file")==0){
      fwrite(&(mdtherm_info->num_nhc),sizeof(int),n,fp);
      fwrite(&(mdtherm_info->len_nhc),sizeof(int),n,fp);
      fwrite(&(mdtherm_info_bead->num_nhc),sizeof(int),n,fp);
      fwrite(&(mdtherm_info_bead->len_nhc),sizeof(int),n,fp);
    }//endif

    if(strcasecmp(file_typ,"par_file")==0){
      fwrite(&(genfilenames->low_lim_par),sizeof(int),n,fp);
      fwrite(&(genfilenames->high_lim_par),sizeof(int),n,fp);
    }//endif

    if(strcasecmp(file_typ,"vel_file")==0||
        strcasecmp(file_typ,"pos_file")==0||
        strcasecmp(file_typ,"cen_file")==0||
        strcasecmp(file_typ,"par_file")==0){
      for(i=1;i<=mdatom_maps->nmol_typ;i++){
        strcpy(name_scr,mdatom_maps->mol_typ[i]);
        fwrite(&(name_scr),sizeof(NAME),n,fp);
      }//endfor
      for(i=1;i<=mdatom_maps->nres_typ;i++){
        strcpy(name_scr,mdatom_maps->res_typ[i]);
        fwrite(&(name_scr),sizeof(NAME),n,fp);
      }//endfor
      for(i=1;i<=mdatom_maps->natm_typ;i++){
        strcpy(name_scr,mdatom_maps->atm_typ[i]);
        fwrite(&(name_scr),sizeof(NAME),n,fp);
      }//endfor
      low = 1;high=mdclatoms_info->natm_tot;
      if(strcasecmp(file_typ,"par_file")==0){
        low  = genfilenames->low_lim_par;
        high = genfilenames->high_lim_par;
      }//endif
      for(i=low;i<=high;i++){
        fwrite(&(mdclatoms_info->mass[i]),sizeof(double),n,fp);
        fwrite(&(mdclatoms_info->q[i]),sizeof(double),n,fp);
        fwrite(&(mdatom_maps->iatm_mol_typ[i]),sizeof(int),n,fp);
        fwrite(&(mdatom_maps->iatm_mol_num[i]),sizeof(int),n,fp);
        fwrite(&(mdatom_maps->iatm_res_typ[i]),sizeof(int),n,fp);
        fwrite(&(mdatom_maps->iatm_res_num[i]),sizeof(int),n,fp);
        fwrite(&(mdatom_maps->iatm_atm_typ[i]),sizeof(int),n,fp);
      }//endfor
    }//endif
    fclose(fp);
  }//endif

  //==========================================================================
}//end routine
//==========================================================================


