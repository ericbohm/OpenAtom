//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                                                                          
//                         PI_MD:                                           
//             The future of simulation technology                          
//             ------------------------------------                         
//                     Module: read_coord                                   
//                                                                          
// This subprogram reads atm-atm_NHC vol-vol_NHC input for a MD on a        
// LD-classical potential energy surface (LD-PES)                           
//                                                                          
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

#include "standard_include.h"

#include "../class_defs/allclass_gen.h"
#include "../class_defs/allclass_cp.h"
#include "../class_defs/allclass_mdintegrate.h"
#include "../class_defs/allclass_mdatoms.h"
#include "../class_defs/allclass_mdinter.h"
#include "../class_defs/allclass_mdintra.h"
#include "../class_defs/typedefs_par.h"
#include "../class_defs/ATOMS/class_mdclatoms_pos.h"
#include "../class_defs/INTEGRATE/class_mdtherm_pos.h"

#include "../proto_defs/proto_coords_entry.h"
#include "../proto_defs/proto_coords_local.h"
#include "../proto_defs/proto_handle_entry.h"

#define PATH_INTEGRALS_NOT_IMPLEMENTED


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

void  read_coord(MDINTEGRATE *mdintegrate,MDATOMS *mdatoms,MDINTER *mdinter,
                 MDINTRA *mdintra,GENERAL_DATA *general_data, CP *cp,
                 MDCLATOMS_POS *clatoms_pos,MDTHERM_POS *therm_class,
                 MDTHERM_POS *therm_bead)

//======================================================================
  {   //begin routine 
//======================================================================
//               Local variable declarations                            

  int iii,upper,lower,igloc,ighost,igo;
  int i,j,ip,nmall,ip_now;       
  int natm_tot_now,istart_now,pi_beads_now;
  int imol_num_now,num_nhc_now,len_nhc_now;
  int itime_dump; 

  double vol,h1,h2,h3,deth;
  double *x_tmp,*y_tmp,*z_tmp;
  double *v_nhc_tmp;

  NAME restart_type_now,atm_typ_now,res_typ_now,mol_typ_now;
  NAME restart_type_spec;
  char *line;

  FILE *fp_dnamei;  

#include "../class_defs/allclass_strip_mdatoms.h"
#include "../class_defs/allclass_strip_mdintegrate.h"
#include "../class_defs/allclass_strip_mdinter.h"

#include "../class_defs/allclass_strip_mdintra.h"
#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_cp.h"

 //---------------------------------------------------------------------
 // Local Pointers 

  double *x,*y,*z;        // assigned below
  double *vx,*vy,*vz;
  double **v_nhc,**x_nhc;

  int cp_dual_grid_opt= cpopts->cp_dual_grid_opt;

  char *dnamei        = genfilenames->dnamei;

  int natm_tot        = mdclatoms_info->natm_tot;
  int pi_beads        = mdclatoms_info->pi_beads; 

  int iseed           = mdvel_samp->iseed;
  int iseed2          = mdvel_samp->iseed2;
  double qseed        = mdvel_samp->qseed;

  int class_num_nhc   = mdtherm_info->num_nhc;
  int class_len_nhc   = mdtherm_info->len_nhc;
  int bead_len_nhc    = mdtherm_info_bead->len_nhc;
  int bead_num_nhc    = mdtherm_info_bead->num_nhc;

  NAME *atm_typ       = mdatom_maps->atm_typ;
  int *iatm_atm_typ   = mdatom_maps->iatm_atm_typ;
  NAME *res_typ       = mdatom_maps->res_typ;
  int *iatm_res_typ   = mdatom_maps->iatm_res_typ;
  NAME *mol_typ       = mdatom_maps->mol_typ;
  int *iatm_mol_typ   = mdatom_maps->iatm_mol_typ;
  int *iatm_mol_num   = mdatom_maps->iatm_mol_num;

  int nfreeze         = mdconstrnt->nfreeze;
  int pimd_freez_typ  = mdconstrnt->pimd_freez_typ;
  int *freeze_map     = mdconstrnt->freeze_map;

  int nghost_tot      = mdghost_atoms->nghost_tot;
  int *ighost_map     = mdghost_atoms->ighost_map;

  int initial_spread_opt = gensimopts->initial_spread_opt;

  int nvt             = genensopts->nvt;
  int npt_i           = genensopts->npt_i;
  int npt_f           = genensopts->npt_f;

  double *vgmat       = mdpar_rahman->vgmat;

  double *v_vol_nhc   = mdbaro->v_vol_nhc;
  double *x_vol_nhc   = mdbaro->x_vol_nhc;

  int hmat_int_typ    = gencell->hmat_int_typ;
  int hmat_cons_typ   = gencell->hmat_cons_typ;
  int iperd           = gencell->iperd;
  double *hmat        = gencell->hmat;

  int ip_start        = 1;       // hard coded because PIMD not ready
  int ip_end          = pi_beads;

  int istart          = gensimopts->istart;

  int pimd_on         =  (gensimopts->pimd + gensimopts->cp_pimd 
                        + gensimopts->cp_wave_pimd 
                        + gensimopts->cp_wave_min_pimd  
                        + gensimopts->debug_pimd 
                        + gensimopts->debug_cp_pimd);


//========================================================================
// I)Write to screen:                                                   


  PRINT_LINE_STAR;
  PRINTF("Reading user specified atm coordinate file %s\n",dnamei);
   if(istart==1){printf("using the `initial' restart option\n");}
   if(istart==2){printf("using the `restart_pos' restart option\n");}
   if(istart==3){printf("using the `restart_posvel' restart option\n");}
   if(istart==4){printf("using the `restart_all' restart option\n");}
  PRINT_LINE_DASH;printf("\n");

#ifndef PATH_INTEGRALS_IMPLEMENTED
  if( (pimd_on==1) || (pi_beads > 1) ){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Path integrals not yet implemented in read_coord\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }//endif
#endif

//========================================================================
// II) Open the file and malloc:                                          



  fp_dnamei = cfopen(dnamei,"r");


  line      = (char *)cmalloc(MAXLINE*sizeof(char),"read_coord");
  x_tmp     = (double *)cmalloc(natm_tot*sizeof(double),"read_coord")-1;
  y_tmp     = (double *)cmalloc(natm_tot*sizeof(double),"read_coord")-1;
  z_tmp     = (double *)cmalloc(natm_tot*sizeof(double),"read_coord")-1;

  nmall     = MAX(class_num_nhc,bead_num_nhc);
  if(nmall==0)nmall++;
  v_nhc_tmp = (double *)cmalloc(nmall*sizeof(double),"read_coord")-1;
  nmall     = MAX(class_len_nhc,1);
  v_vol_nhc = (double *)cmalloc(nmall*sizeof(double),"read_coord")-1;

//========================================================================
//  III)Read in header                                                    


//-------------------
// A) Type 1 start : 
    if(istart==1){
      if(fscanf(fp_dnamei,"%d %d %d",&natm_tot_now,&istart_now,
                                   &pi_beads_now)!=3){
        PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("Error reading start type and number of atoms \n");
        PRINTF("in file %s\n",dnamei);
        PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        FFLUSH(stdout);
        EXIT(1);
      }//endif
      readtoendofline(fp_dnamei);
    }//endif

//---------------------
// B) Type 2,3,4 start 
    if(istart>1){
      if(fgetc(fp_dnamei)==EOF){
        PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("Error reading header information \n");
        PRINTF("in file %s\n",dnamei);
        PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        FFLUSH(stdout);
        EXIT(1);
      }//endif
      readtoendofline(fp_dnamei);
      if(fscanf(fp_dnamei,"%d %s %d %d",&natm_tot_now,restart_type_now,
                                        &itime_dump,&pi_beads_now)!=4){
        PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("Error reading start type and number of atoms \n");
        PRINTF("in file %s\n",dnamei);
        PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      }//endif
      readtoendofline(fp_dnamei);

      istart_now = 0;
      if(strcasecmp(restart_type_now,"initial") == 0){        istart_now = 1;}
      if(strcasecmp(restart_type_now,"restart_pos") == 0){    istart_now = 2;}
      if(strcasecmp(restart_type_now,"restart_posvel") == 0){ istart_now = 3;}
      if(strcasecmp(restart_type_now,"restart_all") == 0){    istart_now = 4;}

      if(istart==1){strcpy(restart_type_spec,"initial");}
      if(istart==2){strcpy(restart_type_spec,"restart_pos");}
      if(istart==3){strcpy(restart_type_spec,"restart_posvel");}
      if(istart==4){strcpy(restart_type_spec,"restart_all");}

      if(istart_now==0){
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("Start up %s option in ",restart_type_now);
        PRINTF("user specified coordinate file %s\n",dnamei);
        PRINTF("not supported. Supported genre: \n");
        PRINTF("initial, restart_pos, restart_posvel, restart_all-> \n");     
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        FFLUSH(stdout);
        EXIT(1);
      } // endif 
    } // endif: istart 
  
//--------------------- 
// C) General checks   
    if(istart_now < istart ) {
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("Start up option, %s, in ",restart_type_now);
      PRINTF("user specified coordinate file %s\n",dnamei);
      PRINTF("Incompatible with class setup,%s\n",restart_type_spec);
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    } // endif 
 
    if(natm_tot_now != natm_tot) {
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("Number of particles in\n");
      PRINTF("user specified coordinate file %s \n",dnamei);
      PRINTF("incompatible with class setup\n");
      PRINTF("%d vs %d\n",natm_tot_now,natm_tot);
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    } // endif 
    if(pi_beads_now != pi_beads) {
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("Number of path integral beads in\n");
      PRINTF("user specified coordinate file %s \n",dnamei);
      PRINTF("incompatible with class setup\n");
      PRINTF("%d vs %d\n",pi_beads_now,pi_beads);
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    } // endif 
    if(istart>2 && pi_beads > 1 && initial_spread_opt == 1){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("The spread option is to be used only with \n");
      PRINTF("the options : restart_pos or initial  \n");
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }//endif

//========================================================================
// IV)istart = 1 (initial)                                              


  if(istart == 1) {

//-----------------------------------------------------------------------
// A)Atm positions                                                  
    PRINTF("Reading in coordinates\n");
    upper = pi_beads;
    if(initial_spread_opt == 1){upper = 1;}
    for(ip=1;ip<=upper;ip++){
      for(i=1;i<=natm_tot;i++){
        if(fscanf(fp_dnamei,"%lf %lf %lf",
                        &(x_tmp[i]),
                        &(y_tmp[i]),
                        &(z_tmp[i])) != 3) {
            PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
            PRINTF("Error while reading in the %d atom coordinate\n",i);
            PRINTF("in file \"%s\"\n",dnamei);
            PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
            FFLUSH(stdout);
            EXIT(1);
        }//endif
        readtoendofline(fp_dnamei);
        x_tmp[i] /= BOHR;
        y_tmp[i] /= BOHR;
        z_tmp[i] /= BOHR;
      }//endfor:atoms
      if( (ip>=ip_start) && (ip <=ip_end) ){
         ip_now = ip-ip_start + 1;
         x = clatoms_pos[ip_now].x;
         y = clatoms_pos[ip_now].y;
         z = clatoms_pos[ip_now].z;
         for(i=1;i<=natm_tot;i++){
           x[i] = x_tmp[i];          
           y[i] = y_tmp[i];          
           z[i] = z_tmp[i];          
         }//endfor
      }//endif
    }//endfor : beads

#ifdef PATH_INTEGRALS_IMPLEMENTED
    if(initial_spread_opt == 1 && pi_beads>1){
      spread_coord(mdclatoms_info,clatoms_pos,x_tmp,y_tmp,z_tmp,
                   &iseed,&iseed2,&qseed,mdatommaps);
    }//endif
#endif
 
//------------------------------------------------------------------
//skip over the box
    for(i=0;i<3;i++){
      if(fscanf(fp_dnamei,"%lf %lf %lf",&h1,&h2,&h3) != 3){
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("Error while reading in the %d cell vector \n",i+1);
        PRINTF("in file \"%s\"\n",dnamei);
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        FFLUSH(stdout);
        EXIT(1);
      }//endif
      readtoendofline(fp_dnamei);
    }//endfor

  }//endif: start=1 

//========================================================================
// V)istart = 2 (restart_pos)                                             


  if(istart >= 2) {

//----------------------------------------------------------------------
//     A)Atm positions                                                  
    PRINTF("Reading in coordinates\n");
    if(fgets(line,MAXLINE,fp_dnamei)==NULL){
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("EOF before particle coordinates \n");
        PRINTF("in file \"%s\"\n",dnamei);
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        FFLUSH(stdout);
        EXIT(1);
    }//endif
    upper = pi_beads;
    if(initial_spread_opt == 1){upper = 1;}
    for(ip=1;ip<=upper;ip++){
      for(i=1;i<=natm_tot;i++){
        if(fscanf(fp_dnamei,"%lf %lf %lf %s %s %s %d",
                            &(x_tmp[i]),
                            &(y_tmp[i]),
                            &(z_tmp[i]),
                            atm_typ_now,res_typ_now,mol_typ_now,
                            &imol_num_now) != 7) {
          PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
          PRINTF("Error while reading in the %d atom coordinate\n",i);
          PRINTF("in file \"%s\"\n",dnamei);
          PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
          FFLUSH(stdout);
          EXIT(1);
        }//endif
        readtoendofline(fp_dnamei);
        if(strcasecmp(atm_typ_now,atm_typ[iatm_atm_typ[i]]) != 0){
          PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
          PRINTF("Atom type mismatch for particle %d\n",i);
          PRINTF("in user specified coordinate file %s \n",dnamei);
          PRINTF("File says %s program expects %s\n",atm_typ_now,
                  atm_typ[iatm_atm_typ[i]]);
          PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
          FFLUSH(stdout);
          EXIT(1);
        } // endif 
        if(strcasecmp(res_typ_now,res_typ[iatm_res_typ[i]]) != 0){
          PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
          PRINTF("Residue type mismatch for particle %d\n",i);
          PRINTF("in user specified coordinate file %s \n",dnamei);
          PRINTF("File says %s program expects %s \n",res_typ_now,
                  res_typ[iatm_res_typ[i]]);
          PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
          FFLUSH(stdout);
          EXIT(1);
        } // endif 
        if(strcasecmp(mol_typ_now,mol_typ[iatm_mol_typ[i]]) != 0){
          PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
          PRINTF("Molecule type mismatch for particle %d\n",i);
          PRINTF("in user specified coordinate file %s \n",dnamei);
          PRINTF("File says %s program expects %s \n",mol_typ_now,
                  mol_typ[iatm_mol_typ[i]]);
          PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
          FFLUSH(stdout);
          EXIT(1);
        } // endif 
        if(imol_num_now != iatm_mol_num[i]) {
          PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
          PRINTF("Molecule number mismatch for particle %d\n",i);
          PRINTF("in user specified coordinate file %s \n",dnamei);
          PRINTF("File says %d program expects %d \n",imol_num_now,
                  iatm_mol_num[i]);
          PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
          FFLUSH(stdout);
          EXIT(1);
        } // endif 
      }//endfor:atoms
      if(ip>=ip_start && ip <=ip_end){
        ip_now = ip - ip_start + 1;
        x = clatoms_pos[ip_now].x;
        y = clatoms_pos[ip_now].y;
        z = clatoms_pos[ip_now].z;
        for(i=1;i<=natm_tot;i++){
          x[i] = x_tmp[i];          
          y[i] = y_tmp[i];          
          z[i] = z_tmp[i];          
        }//endfor
      }//endif
    }//endfor:pi_beads

#ifdef PATH_INTEGRALS_IMPLEMENTED
    if(initial_spread_opt == 1&&pi_beads>1){
      spread_coord(mdclatoms_info,clatoms_pos,x_tmp,y_tmp,z_tmp,
                   &iseed,&iseed2,&qseed,mdatommaps);
    }//endif
#endif

//------------------------------------------------------------------
// Skip over the box
    PRINTF("Reading in cell shape information\n");
    if(fgets(line,MAXLINE,fp_dnamei)==NULL){
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("EOF before cell vectors \n");
        PRINTF("in file \"%s\"\n",dnamei);
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        FFLUSH(stdout);
        EXIT(1);
    }//endif
    for(i=0;i<3;i++){
      if(fscanf(fp_dnamei,"%lf %lf %lf",&h1,&h2,&h3) != 3){
          PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
          PRINTF("Error reading in cell vector %d\n",i);
          PRINTF("in file \"%s\"\n",dnamei);
          PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
          FFLUSH(stdout);
          EXIT(1);
      }// endif
      readtoendofline(fp_dnamei);
    }// endfor
    if(cp_dual_grid_opt >= 1){
      readtoendofline(fp_dnamei);
      readtoendofline(fp_dnamei);
      readtoendofline(fp_dnamei);
      readtoendofline(fp_dnamei);
      readtoendofline(fp_dnamei);
      readtoendofline(fp_dnamei);
    }// endif cp_dual_grid_opt
    readtoendofline(fp_dnamei);
    readtoendofline(fp_dnamei);
    readtoendofline(fp_dnamei);
    readtoendofline(fp_dnamei);
    readtoendofline(fp_dnamei);
    readtoendofline(fp_dnamei);

  }//endif: start>=2

//========================================================================
// VI)Istart = 3 (restart_posvel)                                         


  if(istart >= 3) {

//----------------------------------------------------------------------
// A)Atm Velocities                                                 
    PRINTF("Reading in atom velocities\n");
    if(fgets(line,MAXLINE,fp_dnamei)==NULL){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("EOF before particle velocities \n");
      PRINTF("in file \"%s\"\n",dnamei);
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }//endif
    for(ip=1;ip<=pi_beads;ip++){
      for(i=1;i<=natm_tot;i++){
        if(fscanf(fp_dnamei,"%lf %lf %lf",
                               &(x_tmp[i]),
                               &(y_tmp[i]),
                               &(z_tmp[i]) )!=3){
          PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
          PRINTF("Error while reading in the %d atom velocities\n",i);
          PRINTF("in file \"%s\"\n",dnamei);
          PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
          FFLUSH(stdout);
          EXIT(1);
        }//endif
        readtoendofline(fp_dnamei);
      }//endfor:atoms
      if( (ip>=ip_start) && (ip <=ip_end)){
         ip_now = ip-ip_start + 1;
         vx = clatoms_pos[ip_now].vx;
         vy = clatoms_pos[ip_now].vy;
         vz = clatoms_pos[ip_now].vz;
         for(i=1;i<=natm_tot;i++){
           vx[i] = x_tmp[i];          
           vy[i] = y_tmp[i];          
           vz[i] = z_tmp[i];          
         }//endfor
      }//endif
    }//endfor : pi_beads 

//----------------------------------------------------------------------
// B) Zero velocities of ghost atoms if any                            

    if(nghost_tot > 0) {
      for(ip=1;ip<=pi_beads;ip++){
        vx = clatoms_pos[ip].vx;
        vy = clatoms_pos[ip].vy;
        vz = clatoms_pos[ip].vz;
        for(ighost=1;ighost <= nghost_tot;ighost++){
          igloc = ighost_map[ighost];
          vx[igloc] = 0.0;
          vy[igloc] = 0.0;
          vz[igloc] = 0.0;
        }//endfor
      }//endfor : pi_beads
    }//endif

//----------------------------------------------------------------------
// C) Zero velocities of freeze atoms if any                         

    if(nfreeze > 0) {
      igo = 0;
      upper = (ip_start == 1 ? 1 : 0);
      if(pimd_freez_typ == 2){upper=pi_beads;igo=1;}
      if(igo == 1){
        for(ip=1;ip<=upper;ip++){
        vx = clatoms_pos[ip].vx;
        vy = clatoms_pos[ip].vy;
        vz = clatoms_pos[ip].vz;
        for(i=1;i <= nfreeze;i++){
          igloc = freeze_map[i];
          vx[igloc] = 0.0;
          vy[igloc] = 0.0;
          vz[igloc] = 0.0;
        }//endfor
        }//endfor
      }//endif
    }//endif

  }// endif : istart >=3

//========================================================================
//   VII)Istart = 4 (restart_all)                                         


  if(istart == 4) {

//------------------------------------------------------------------
//  A)Atm NHC Velocities                                         

    PRINTF("Reading Nose-Hoover chain (NHC) velocities\n");
    if(fgets(line,MAXLINE,fp_dnamei)==NULL){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("EOF before NHC information \n");
      PRINTF("in file \"%s\"\n",dnamei);
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }//endif
    if(fscanf(fp_dnamei,"%d %d",&num_nhc_now,&len_nhc_now)!=2){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("Error while reading NHC information\n");
      PRINTF("in file \"%s\"\n",dnamei);
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }//endif 
    readtoendofline(fp_dnamei);

    if(num_nhc_now != class_num_nhc) {
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("Mismatched number of Nose-Hoover chains\n");
      PRINTF("%d vs %d\n",class_num_nhc, num_nhc_now);
      PRINTF("in user specified coordinate file %s \n",dnamei);
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }// endif 
    if(len_nhc_now != class_len_nhc) {
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("Mismatched length of Nose-Hoover chains\n");
      PRINTF("%d vs %d\n",class_len_nhc,len_nhc_now);
      PRINTF("in user specified coordinate file %s \n",dnamei);
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    } // endif 
    if(fgets(line,MAXLINE,fp_dnamei)==NULL){
       PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       PRINTF("EOF before NHC velocities \n");
       PRINTF("in file \"%s\"\n",dnamei);
       PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       FFLUSH(stdout);
       EXIT(1);
    }//endif

    v_nhc = therm_class->v_nhc;
    for(j=1;j<=class_len_nhc;j++){
      for(i=1;i<=class_num_nhc;i++){
        if(fscanf(fp_dnamei,"%lf",&(v_nhc_tmp[i]))!=1){
          PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
          PRINTF("Error while reading in the %d %d NHC velocity \n",j,i);
          PRINTF("in file \"%s\"\n",dnamei);
          PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
          FFLUSH(stdout);
          EXIT(1);
        }//endif
        readtoendofline(fp_dnamei);
      }//endfor
    }//endif
    if(ip_start==1){
      for(i=1;i<=class_num_nhc;i++){
        v_nhc[j][i] = v_nhc_tmp[i];
      }//endfor
    }//endfor : length of nhc 

//------------------------------------------------------------------
//  B)Bead NHC Velocities                                         

    if(pi_beads>1){
      PRINTF("Reading Nose-Hoover bead (NHC) velocities\n");
      if(fgets(line,MAXLINE,fp_dnamei)==NULL){
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("EOF before NHC information \n");
        PRINTF("in file \"%s\"\n",dnamei);
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        FFLUSH(stdout);
        EXIT(1);
      }//endif
      fgets(line,MAXLINE,fp_dnamei);
      if(sscanf(line,"%d %d \n",&num_nhc_now,&len_nhc_now)!=2){
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("Error while reading NHC information\n");
        PRINTF("in file \"%s\"\n",dnamei);
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        FFLUSH(stdout);
        EXIT(1);
      }//endif 
      if(num_nhc_now != bead_num_nhc) {
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("Mismatched number of Nose-Hoover chains\n");
        PRINTF("%d vs %d\n",bead_num_nhc, num_nhc_now);
        PRINTF("in user specified coordinate file %s \n",dnamei);
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        FFLUSH(stdout);
        EXIT(1);
      }// endif 
      if(len_nhc_now != bead_len_nhc) {
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("Mismatched length of Nose-Hoover chains\n");
        PRINTF("%d vs %d\n",bead_len_nhc,len_nhc_now);
        PRINTF("in user specified coordinate file %s \n",dnamei);
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        FFLUSH(stdout);
        EXIT(1);
      } // endif 
      if(fgets(line,MAXLINE,fp_dnamei)==NULL){
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("EOF before NHC velocities \n");
        PRINTF("in file \"%s\"\n",dnamei);
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        FFLUSH(stdout);
        EXIT(1);
      }//endif

      for(ip=2;ip<=pi_beads;ip++){
        for(j=1;j<=bead_len_nhc;j++){
          for(i=1;i<=bead_num_nhc;i++){
            if(fgets(line,MAXLINE,fp_dnamei)==NULL){
              PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
              PRINTF("EOF while reading in the %d %d NHC velocity \n",j,i);
              PRINTF("in file \"%s\"\n",dnamei);
              PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
              FFLUSH(stdout);
              EXIT(1);
            }//endif
            if(sscanf(line,"%lf \n",&v_nhc_tmp[i])!=1){
              PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
              PRINTF("Error reading in the %d %d bead NHC velocity \n",j,i);
              PRINTF("in file \"%s\"\n",dnamei);
              PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
              FFLUSH(stdout);
              EXIT(1);
            }//endif
          }//endfor:num_nhc
          if(ip>=ip_start && ip <=ip_end){
            ip_now = ip-ip_start + 1;
            v_nhc  = therm_bead[ip_now].v_nhc;
            for(i=1;i<=bead_num_nhc;i++){
              v_nhc[j][i] = v_nhc_tmp[i];          
            }//endfor
          }//endif : ip in range
        }//endfor:len_nhc
      }//endfor:pi_beads

    }//endif:pi_beads > 1

//------------------------------------------------------------------
//  C)Vol and Vol NHC Velocities                                    

    PRINTF("Reading Cell velocities\n");
    if(fgets(line,MAXLINE,fp_dnamei)==NULL){
       PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       PRINTF("EOF before Cell velocities \n");
       PRINTF("in file \"%s\"\n",dnamei);
       PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       FFLUSH(stdout);
       EXIT(1);
    }//endif
    for(i=0;i<3;i++){
      if(fscanf(fp_dnamei,"%lf %lf %lf",
         &(vgmat[(1+i)]),&(vgmat[(4+i)]),&(vgmat[(7+i)])) != 3){
         PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
         PRINTF("EOF reading in cell vector velocity %d\n",i);
         PRINTF("in file \"%s\"\n",dnamei);
         PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
         FFLUSH(stdout);
         EXIT(1);
      }//endif
      readtoendofline(fp_dnamei);
    }//endfor
    if(hmat_int_typ==0){
      if((vgmat[2]!=vgmat[4])&&
         (vgmat[3]!=vgmat[7])&&
         (vgmat[6]!=vgmat[8])){
         PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
         PRINTF("The Cell velocity matrix must be symmetric\n");
         PRINTF("in file \"%s\"\n",dnamei);
         PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
         FFLUSH(stdout);
         EXIT(1);
      }//endif
    }//endif
    if(hmat_int_typ==1){
      if((vgmat[2]!=0.0)&&
         (vgmat[3]!=0.0)&&
         (vgmat[6]!=0.0)){
         PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
         PRINTF("The Cell velocity matrix must be upper triangular\n");
         PRINTF("in file \"%s\"\n",dnamei);
         PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
         FFLUSH(stdout);
         EXIT(1);
      }//endif
    }//endif
    if(hmat_cons_typ==1){
      if( (vgmat[2] != 0.0) || (vgmat[3] != 0.0) ||
          (vgmat[4] != 0.0) || (vgmat[6] != 0.0) || 
          (vgmat[7] != 0.0) || (vgmat[8] != 0.0) ){
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("The Cell velocity matrix should contain no off-diagonal \n"); 
        PRINTF("coupling with the orthorhombic constraint \n");
        PRINTF("in file \"%s\"\n",dnamei);
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        FFLUSH(stdout);
        EXIT(1);      
      }// endif  
    }// endif  
    if(hmat_cons_typ==2){
      if( (vgmat[3] != 0.0) || (vgmat[6] != 0.0) ||
          (vgmat[7] != 0.0) || (vgmat[8] != 0.0) ){
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("The Cell velocity matrix should contain no c-coupling with\n");
        PRINTF("the monoclinic constraint in file \"%s\"\n",dnamei);
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        FFLUSH(stdout);
        EXIT(1);      
      }// endif  
    }// endif  
    if(iperd==2){
      if((vgmat[7]!=0.0)&&
         (vgmat[8]!=0.0)&&
         (vgmat[3]!=0.0)&&
         (vgmat[6]!=0.0)){
         PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
         PRINTF("The Cell velocity matrix must not contain c-coupling\n");
         PRINTF("in file \"%s\"\n",dnamei);
         PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
         FFLUSH(stdout);
         EXIT(1);
      }//endif
    }//endif
    PRINTF("Reading Volume velocity\n");
    if(fgets(line,MAXLINE,fp_dnamei)==NULL){
       PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       PRINTF("EOF before Volume velocity\n");
       PRINTF("in file \"%s\"\n",dnamei);
       PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       FFLUSH(stdout);
       EXIT(1);
    }//endif
    if(fscanf(fp_dnamei,"%lf",&(mdbaro->v_lnv)) != 1){
       PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       PRINTF("Error reading in volume velocity\n");
       PRINTF("in file \"%s\"\n",dnamei);
       PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       FFLUSH(stdout);
       EXIT(1);
    }//endif
    readtoendofline(fp_dnamei);
    PRINTF("Reading Volume NHC velocities\n");
    if(fgets(line,MAXLINE,fp_dnamei)==NULL){
       PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       PRINTF("EOF before Vol. NHC velocities \n");
       PRINTF("in file \"%s\"\n",dnamei);
       PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       FFLUSH(stdout);
       EXIT(1);
    }//endif
    for(i=1;i<=class_len_nhc;i++){
      if(fscanf(fp_dnamei,"%lf",&(v_vol_nhc[i])) != 1){
         PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
         PRINTF("Error reading in volume NHC velocity %d\n",i);
         PRINTF("in file \"%s\"\n",dnamei);
         PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
         FFLUSH(stdout);
         EXIT(1);
      }//endif
      readtoendofline(fp_dnamei);
    }//endfor

  }//endif : start==4 

//========================================================================
//  VIII) Assign NHC positions                       


  if( ((nvt+npt_i+npt_f)==1)  && (ip_start ==1) ){
     x_nhc = therm_class->x_nhc;
     for(j=1;j<=class_len_nhc;j++){
       for(i=1;i<=class_num_nhc;i++){
          x_nhc[j][i] = 0.0;
       }//endfor
     }//endfor
  }//endif : bead and ensemble 

  if( pi_beads>1 ){
    lower = (ip_start==1 ? 2 : 1);
    for(ip=lower;ip<=pi_beads;ip++){
      x_nhc = therm_bead[ip].x_nhc;
      for(j=1;j<=bead_len_nhc;j++){
        for(i=1;i<=bead_num_nhc;i++){
          x_nhc[j][i] = 0.0;
        }//endfor:num_nhc
      }//endfor:len_nhc
    }//endfor:pi_beads
  }//endif:pi_beads

  if( ((npt_i+npt_f)==1) && (ip_start==1) ){
    for(j=1;j<=class_len_nhc;j++){
      x_vol_nhc[j] = 0.0;
    }//endfor
  }//endif

//========================================================================
//  IX) Calculate the spread                                          

#ifdef PATH_INTEGRALS_IMPLEMENTED
  mdinteract->spread     = 0.0;
  mdinteract->spread_now = 0.0;
  if( (pi_beads>1) || (pimd_on==1) ){
    get_pimd_spread(mdclatoms_info,clatoms_pos,&(mdinteract->spread_now));
    mdinteract->spread = mdinteract->spread_now;
  }//endif
#endif

//========================================================================
//  X) Close file, free arrays                                             
  
  fclose(fp_dnamei); 

  cfree(&x_tmp[1],"read_coord");
  cfree(&y_tmp[1],"read_coord");
  cfree(&z_tmp[1],"read_coord");
  cfree(&v_nhc_tmp[1],"read_coord");
  cfree(&v_vol_nhc[1],"read_coord");
  cfree(line,"read_coord");

//========================================================================
//  XI) Output to the screen                                               


  PRINTF("\n");
  PRINT_LINE_DASH;
  PRINTF("Done reading user specified atm coordinate file %s\n",dnamei);
  PRINT_LINE_STAR;PRINTF("\n");

//--------------------------------------------------------------------------
    }// end routine 
//==========================================================================



