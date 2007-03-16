//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                                                                          
//                         PI_MD:                                           
//             The future of simulation technology                          
//             ------------------------------------                         
//                     Module: parse                                        
//                                                                          
// This subprogram reads in all user inputs                                 
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
#include "../../../include/Atoms.h"
#include "../class_defs/ATOM_OPERATIONS/class_atomoutput.h"

#include "../proto_defs/proto_parse_entry.h"
#include "../proto_defs/proto_parse_local.h"
#include "../proto_defs/proto_sim_params_entry.h"
#include "../proto_defs/proto_mol_params_entry.h"
#include "../proto_defs/proto_intra_params_entry.h"
#include "../proto_defs/proto_cp_ewald_entry.h"
#include "../proto_defs/proto_surf_params_entry.h"
#include "../proto_defs/proto_inter_params_entry.h"
#include "../proto_defs/proto_inter_params_local.h"
#include "../proto_defs/proto_path_init_entry.h"
#include "../proto_defs/proto_coords_entry.h"
#include "../proto_defs/proto_vps_params_entry.h"
#include "../proto_defs/proto_lists_entry.h"
#include "../proto_defs/proto_coords_cp_entry.h"

#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_friend_lib_entry.h"

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//  Parse: note there is a noncommuting order of calls                      
//==========================================================================

void parse(MDINTEGRATE *mdintegrate, MDATOMS *mdatoms, MDINTER *mdinter,
           MDINTRA *mdintra, GENERAL_DATA *general_data,CP *cp, 
           char *input_name)

//========================================================================
//             Begin Routine                                              
         {//Begin subprogram: 
//========================================================================
//             Local variable declarations                                

  CLASS_PARSE      class_parse;
  CP_PARSE         cp_parse;
  FILENAME_PARSE   filename_parse;
  FREE_PARSE       free_parse;
  SPLINE_PARSE     spline_parse;
  NULL_INTER_PARSE null_inter_parse;

//------------------------------------------------------------------------

  int iii,nstate2,nstate,ncoef;
  int iextended_on,ewald_on;
  int cp_on,pimd_on,cp_md,isurf_on;
  int nchrg,iperd,pi_beads,pme_on;
  int cp_dual_grid_opt_on;
  int hmat_cons_typ;

//------------------------------------------------------------------------
// local pointers to classes in container classes stored in 

#include "../class_defs/allclass_strip_mdatoms.h"
#include "../class_defs/allclass_strip_mdintegrate.h"
#include "../class_defs/allclass_strip_mdinter.h"
#include "../class_defs/allclass_strip_mdintra.h"
#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_cp.h"

//========================================================================
//   I) Set the sim parameters: Done first to get input file names       
//               (interface/sim_params/control_sim_params.c)              

 
//   printf("Before sim_params \n");   DEBUG_READ_INT

   filename_parse.input_name  = (char *)cmalloc(MAXWORD*sizeof(char),"parse");
   strcpy(filename_parse.input_name,input_name);

   control_sim_params(mdintegrate,mdatoms,mdinter,general_data,mdintra,
                      cp,&class_parse,&cp_parse,&filename_parse);

//========================================================================
//   II) Read in atom and CP parameters: Done second to get info         
//                                        needed for set_intra;           
//                                        some atom mallocing             
//               (interface/mol_params/control_mol_params.c)              

//   printf("Before mol_params \n");   DEBUG_READ_INT

   control_mol_params(mdintegrate,mdatoms,mdinter,general_data,mdintra,
                       cp,&class_parse,&cp_parse,&free_parse,&filename_parse);

//========================================================================
//   III) Assign Flags                                                     


  pimd_on = gensimopts->pimd         + gensimopts->cp_pimd 
          + gensimopts->cp_wave_pimd + gensimopts->cp_wave_min_pimd 
          + gensimopts->debug_pimd   + gensimopts->debug_cp_pimd;

  cp_on = gensimopts->cp_min   + gensimopts->cp_wave_min
         +gensimopts->cp       + gensimopts->cp_wave
         +gensimopts->cp_pimd  + gensimopts->cp_wave_pimd
         +gensimopts->debug_cp + gensimopts->debug_cp_pimd
         +gensimopts->cp_wave_min_pimd;

  cp_md = gensimopts->cp       + gensimopts->cp_wave
         +gensimopts->debug_cp + gensimopts->cp_pimd
         +gensimopts->cp_wave_pimd;

  iextended_on = genensopts->nvt+genensopts->npt_i+genensopts->npt_f;

  cp_dual_grid_opt_on = cpopts->cp_dual_grid_opt;
  isurf_on            = mdsurface->isurf_on;
  iperd               = gencell->iperd;
  pi_beads            = mdclatoms_info->pi_beads;
  pme_on              = mdpart_mesh->pme_on;

//========================================================================
//  IV) Molecule connectivity data: Done before setting     
//                                  therms;                
//                                  intramol mallocing     
//               (interface/intra_params/control_intra_params.c)          

//   printf("Before intra_params \n");   DEBUG_READ_INT

   control_intra_params(tot_memory,mdclatoms_info,mdclatoms_pimd,
                        cpatom_maps,mdatom_maps,
                        mdintra,&filename_parse,&free_parse,
                        &class_parse,&null_inter_parse,
                        gensimopts,isurf_on);


  assign_classes_int_data(mdclatoms_info,mdatom_maps,mdclatoms_pimd,
                          mdbrnch_root,mdinteract,mdverlist,mdexcl,
                          mdnbr_list,mdconstrnt,mdghost_atoms,
                          cpatom_maps,cppseudo,mdtherm_info);


   nchrg = mdclatoms_info->nchrg;

//========================================================================
//    V) Read in hmat. Do before set_cp_ewald                           
//                (interface/coords/read_coord.c)                         

//  printf("Before read_hmat \n");  DEBUG_READ_INT

  read_hmat(mdintegrate,mdatoms,mdinter,
            general_data,&filename_parse,class_parse.istart,
            cp_dual_grid_opt_on,&(cpewald->dbox_rat),
            &(cpewald->box_rat));

  hmat_cons_typ = gencell->hmat_cons_typ;

//========================================================================
// VI) Set up the ewald/cp: Done before setting intermol PE             
//                            CP/Ewald mallocing                          
//                (interface/cp_ewald/control_set_cp_ewald                

//  printf("Before cpewald \n");  DEBUG_READ_INT

  ewald_on = 0;
  if((nchrg > 0 && iperd > 0) || cp_on==1){
     ewald_on = 1;
     control_set_cp_ewald(gensimopts,gencell,cpcoeffs_info,cpopts,genewald,
                          cpewald,&cp_parse,
                          &(cppseudo->gmin_true),&(cppseudo->gmin_spl),
                          &(cppseudo->gmax_spl),(class_parse.kmax_ewd),
                          (class_parse.kmax_res),tot_memory,
                          (gentimeinfo->int_res_ter),mdpart_mesh,mdecor,
                          (cpopts->cp_lsda),(genminopts->cp_min_diis),
                          cp_dual_grid_opt_on,&(cppseudo->nonlocal),cppseudo); 

     if(cp_on==1){
       cpopts->te_ext         /= 4.1e6; // empirical scaling factor for h2o 32-70
	 //                   (double) ( (cpcoeffs_info->nstate_up
	 //                              *cpcoeffs_info->ncoef) );
       cpvel_samp->div_scal = (double) (2*(cpcoeffs_info->nstate_up)
                                              *cpcoeffs_info->ncoef);
       if(cpopts->cp_lsda==1){
         cpvel_samp->div_scal += (double) (2*(cpcoeffs_info->nstate_dn)
                                                 *cpcoeffs_info->ncoef);
       }//endif
       PRINTF("\nYour fictious CP temperature is: %gK : different scaling then PINY\n\n",
                                             cpopts->te_ext);
     }//endif : cp_on

  }//endif : charged periodic systems and/or CP
  genewald->ewald_on = ewald_on;


//========================================================================
//  VII) Setup intermolecular potential stuff: interspline mallocing      
//                (interface/inter_params/control_inter_params.c)         

//  printf("Before inter_params \n");  //DEBUG_READ_INT

  control_inter_params(mdinteract,&spline_parse,
                       &filename_parse,genewald->alp_ewd,
                       nchrg,mdinter->mdenergy_ctrl.lj_coul_explicit_opt,
                       mdclatoms_info->natm_tot,
                       mdatom_maps->atm_typ,
                       mdatom_maps->iatm_atm_typ,iperd,
                       class_parse.ishift_pot,tot_memory,
                       gentimeinfo->int_res_ter);

//========================================================================
//   VIII) Setup the surface potential if needed                          

  if(isurf_on == 1){
     control_surf_params(mdsurface, &(filename_parse),
                         mdatom_maps->natm_typ, mdatom_maps->atm_typ, 
                         tot_memory); 
  }// endif  


//========================================================================
//    IX) Setup pseudopotential stuff: pseudospline mallocing            
//                (interface/interparams/control_vps_params.c)            

//  printf("Before vps_params \n");  //DEBUG_READ_INT

  if(cp_on==1){

    set_ylm_cons(&(cp->cpylm_cons));
    make_cp_atom_list(cpatom_maps, cppseudo, mdclatoms_info->natm_tot,
                      mdatom_maps->natm_typ,mdclatoms_info->q,
                      mdatom_maps->iatm_atm_typ);
    control_vps_params(cppseudo,gencell,&filename_parse,
                       &spline_parse,mdatom_maps->iatm_atm_typ,
                       mdatom_maps->natm_typ,
                       mdatom_maps->atm_typ,tot_memory,
                       mdclatoms_info->natm_tot,cpatom_maps->nab_initio,
                       cpopts->cp_ptens_calc,cp_dual_grid_opt_on,
                       cp_parse.cp_ecut);

   //----------------------------------------------------------------------
   // These are obsolete and can be nuked eventually

      ncoef   = cpcoeffs_info->ncoef;
      nstate  = cpcoeffs_info->nstate_up;
      nstate2 = nstate*nstate;
      cpcoeffs_info->occ_up      =  (double *)
                     cmalloc(nstate*sizeof(double),"parse:hack")-1;
      cpcoeffs_info->occ_dn      = (double *)
                     cmalloc(nstate*sizeof(double),"parse:hack")-1;
  }//endif : cp_on

//========================================================================
//    X) Set particle exclusions and ewald corrections                 
//                (interface/lists/set_exclude.c)                         

//  printf("Before exclude \n");  DEBUG_READ_INT


  set_exclude(mdclatoms_info,mdghost_atoms,mdintra,
              mdexcl,&null_inter_parse,cp->cpatom_maps.cp_atm_flag,
              iperd,tot_memory,genewald->alp_ewd);

//========================================================================
//   XI) Set thermostats: Done before reading the coordinates;           
//                        atm NHC mallocing                              
//                (interface/coords/set_atm_NHC.c                        

//  printf("Before atm_NHC \n");  DEBUG_READ_INT

  if(iextended_on==1){
      set_atm_NHC(genensopts,genstatepoint,gensimopts,mdclatoms_info,
                  mdclatoms_pimd,mdghost_atoms,mdatom_maps,
                  mdconstrnt,mdtherm_info_bead,mdtherm_info,
                  mdbaro,mdpar_rahman,&class_parse,tot_memory,pimd_on,
                  hmat_cons_typ);
  }//endif

//========================================================================
//  XII) Spline the ewald corrections

//  printf("Before ecor \n");  DEBUG_READ_INT

  if((nchrg > 0 && iperd > 0)){

    splin_ecor(mdecor,genewald,pi_beads,tot_memory);

  }else{
    genewald->self_erf = 1.0;
  }//endif

//========================================================================
// XIII) Set thermostats: Done before reading the coeffs                 
//                        CP NHC mallocing                                
//                (interface/coords/set_coef_NHC.c)                       

//  printf("Before coef_nhc \n");  DEBUG_READ_INT

  if(cp_on==1){

     if(cp_md==1){

       set_coef_NHC(cpcoeffs_info,cpopts,cptherm_info,&cp_parse,tot_memory);

     } else {

       cptherm_info->num_c_nhc      = 0;
       cptherm_info->massiv_flag    = 0;
 
     }//endif : cp_md_on

  }//endif : cp_on

//========================================================================
//   XIV) Initialize path integral transformations:

//  printf("Before path ini \n");  DEBUG_READ_INT

  if(pi_beads>1||pimd_on==1){
    path_integral_init(mdclatoms_info,mdclatoms_pimd,mdghost_atoms, 
                       gensimopts,mdatom_maps,mdconstrnt);
  }//endif


//========================================================================
//   XV) malloc neigbor list memory                                     
//                (interface/lists/mall_make_lists.c)                     

// Right now we only have coloumbs law
// This can be fixed
#ifdef DONE
  get_cut_skin(mdinteract->cutskin,
               mdinteract->cutskin_root,
               mdinteract->cutoff,
               mdinteract->cutskin_res,
               mdinteract->cutskin_root_res,
               mdinteract->cutoff_res,
               mdinteract->skin,
               mdinteract->spread,
               mdinteract->brnch_root_skin,
               mdinteract->nter_typ);
#endif

//========================================================================
// Initialize the output configuration/trajectory files

   ATOMOUTPUT::initialize_piny_output();
  
//========================================================================
// XVI) Flush the buffers                                                 

  FFLUSH(stdout);
  FFLUSH(stderr);

//------------------------------------------------------------------------
   }//end routine 
//==========================================================================




//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

void assign_classes_int_data(MDCLATOMS_INFO *mdclatoms_info,
                             MDATOM_MAPS *mdatom_maps,
                             MDCLATOMS_PIMD *mdclatoms_pimd,
                             MDBRNCH_ROOT *mdbrnch_root,
                             MDINTERACT *mdinteract,
                             MDVERLIST *mdverlist,
                             MDEXCL *mdexcl,
                             MDNBR_LIST *mdnbr_list,
                             MDCONSTRNT *mdconstrnt,
                             MDGHOST_ATOMS *mdghost_atoms,
                             CPATOM_MAPS *cpatom_maps,
                             CPPSEUDO *cppseudo,
                             MDTHERM_INFO *mdtherm_info)

//========================================================================
      {//Begin subprogram: 
//========================================================================

  int natm_tot = mdclatoms_info->natm_tot;

  cpatom_maps->natm_tot    = natm_tot;
  mdatom_maps->natm_tot    = natm_tot;
  mdclatoms_pimd->natm_tot = natm_tot;
  mdtherm_info->natm_tot   = natm_tot;

  mdnbr_list->natm_tot     = natm_tot;
  mdverlist->natm_tot      = natm_tot;
  mdbrnch_root->natm_tot   = natm_tot;
  mdexcl->natm_tot         = natm_tot; 

  mdghost_atoms->natm_tot  = natm_tot;
  mdconstrnt->natm_tot     = natm_tot;

  //---------------------------------------------
  int pi_beads = mdclatoms_info->pi_beads;

  mdclatoms_pimd->pi_beads = pi_beads;

  //---------------------------------------------
  int natm_typ = mdatom_maps->natm_typ;

  mdinteract->natm_typ = natm_typ;
  cppseudo->natm_typ   = natm_typ;

  //---------------------------------------------

  int nmol_typ = mdatom_maps->nmol_typ;

  mdclatoms_info->nmol_typ = nmol_typ;
  mdconstrnt->nmol_typ     = nmol_typ;

//------------------------------------------------------------------------
   }//end routine 
//==========================================================================


