#include "standard_include.h"
#include "../class_defs/allclass_gen.h"
#include "../class_defs/allclass_cp.h"
#include "../class_defs/allclass_mdintegrate.h"
#include "../class_defs/allclass_mdatoms.h"
#include "../class_defs/allclass_mdinter.h"
#include "../class_defs/allclass_mdintra.h"
#include "../class_defs/typedefs_par.h"
#include "../proto_defs/proto_parse_local.h"


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

void control_zero_cntr(MDINTEGRATE *mdintegrate, MDATOMS *mdatoms, 
                       MDINTER *mdinter,MDINTRA *bonded,
                       GENERAL_DATA *general_data, CP *cp,
                       CLASS_PARSE *class_parse,
                       NULL_INTER_PARSE *null_inter_parse)

//========================================================================
//             Begin Routine                                              
         {//Begin subprogram: 
//========================================================================

 general_data->tot_memory = 0.0;

 // PRINTF("Before zero sys\n"); DEBUG_READ_INT
 zero_sys(mdintegrate,mdatoms,mdinter,general_data);

 // PRINTF("Before zero bnd\n"); DEBUG_READ_INT
 zero_bnd(bonded);

 // PRINTF("Before zero cp\n"); DEBUG_READ_INT
 zero_cp(cp);

 // PRINTF("Before zero par\n"); DEBUG_READ_INT
 zero_par(class_parse,null_inter_parse);

//------------------------------------------------------------------------
   }//end routine 
//==========================================================================



//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

 void zero_bnd(MDINTRA *bonded)

//==========================================================================
{/*begin routine*/
//==========================================================================

int i;

/* zero element to ghost_atoms.ghost_atoms_info  */
    bonded->mdghost_atoms.nghost_tot    = 0;        
    bonded->mdghost_atoms.natm_comp_max = 0;     
    bonded->mdghost_atoms.nghost_old    = 0;
    bonded->mdghost_atoms.ncomp_old     = 0;  

/* zero the int elements of bond_info  */

    bonded->mdbond.npow               = 0;                   
    bonded->mdbond.ntyp_pow           = 0;               
    bonded->mdbond.ncon               = 0;                   
    bonded->mdbond.ntyp_con           = 0;               


/* zero the int elements of group_bond_con_info  */
    bonded->mdgrp_bond_con.max_iter          = 0;                
    bonded->mdgrp_bond_con.num_33            = 0;                 
    bonded->mdgrp_bond_con.num_21            = 0;                  
    bonded->mdgrp_bond_con.num_43            = 0;                  
    bonded->mdgrp_bond_con.num_23            = 0;                  
    bonded->mdgrp_bond_con.num_46            = 0;                  
    bonded->mdgrp_bond_con.ntyp_33           = 0;                 
    bonded->mdgrp_bond_con.ntyp_46           = 0;                 
    bonded->mdgrp_bond_con.ntyp_23           = 0;
    bonded->mdgrp_bond_con.ntyp_21           = 0;
    bonded->mdgrp_bond_con.ntyp_43           = 0; 



/* zero elements of bond_free_info  */
  bonded->mdbond_free.num     = 0;                  
  bonded->mdbond_free.j1      = 0;
  bonded->mdbond_free.j2      = 0;                 
  bonded->mdbond_free.npow    = 0;                  
  bonded->mdbond_free.nhist   = 0;                 



/*======================================================================*/
/* zero elements of  bend   */
    bonded->mdbend.npow                  = 0;                 
    bonded->mdbend.ntyp_pow              = 0;             
    bonded->mdbend.ncon                  = 0;                 
    bonded->mdbend.ntyp_con              = 0;             


/* zero elements of bend_free_info  */

    bonded->mdbend_free.num     = 0;        
    bonded->mdbend_free.j1      = 0;
    bonded->mdbend_free.j2      = 0;
    bonded->mdbend_free.j3      = 0;   
    bonded->mdbend_free.npow    = 0;       
    bonded->mdbend_free.nhist   = 0;      

/* zero elements of bend_bnd_info  */

    bonded->mdbend_bnd.num                 = 0;  
    bonded->mdbend_bnd.ntyp                = 0; 



/* zero tors_info */
    bonded->mdtors.npow               = 0;        
    bonded->mdtors.ntyp_pow           = 0;    
    bonded->mdtors.nimpr              = 0;       
    bonded->mdtors.ncon               = 0;       
    bonded->mdtors.ntyp_con           = 0;   

/* zero tors_free_info  */
    bonded->mdtors_free.num    = 0;
    for(i=1;i<=2;i++){
     bonded->mdtors_free.j1[i] = 0;
     bonded->mdtors_free.j2[i] = 0;
     bonded->mdtors_free.j3[i] = 0;
     bonded->mdtors_free.j4[i] = 0; 
     bonded->mdtors_free.eq[i] = 0.0; 
    }

    bonded->mdtors_free.fk     = 0;        
    bonded->mdtors_free.del    = 1;        
    bonded->mdtors_free.npow   = 0;        
    bonded->mdtors_free.nhist  = 0;       

/* zero onfo_info  */
    bonded->mdonfo.num            = 0;                
    bonded->mdonfo.ntyp           = 0;               


/* zero ecor_info  */
    bonded->mdecor.num         = 0;                    
    bonded->mdecor.alp_ewd     = 0.0;                    


}/*endroutine*/
//==========================================================================



//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

  void zero_sys(MDINTEGRATE *mdintegrate, MDATOMS *mdatoms, MDINTER *mdinter,
                GENERAL_DATA *general_data)

//==========================================================================
{/*begin routine*/


//------------------------------------------------------------------------
// local pointers to classes in container classes stored in 
#include "../class_defs/allclass_strip_mdatoms.h"
#include "../class_defs/allclass_strip_mdintegrate.h"
#include "../class_defs/allclass_strip_mdinter.h"

/* zero cell->cell_info  */
  general_data->gencell.iperd       = 0;     
  general_data->gencell.intra_perds = 0;

/* zero general_data->genstat_avg  */

  general_data->genstat_avg.vintert = 0;
  general_data->genstat_avg.aivintert = 0;
  general_data->genstat_avg.avintert = 0;
  general_data->genstat_avg.vintrat = 0;
  general_data->genstat_avg.aivintrat = 0;
  general_data->genstat_avg.avintrat = 0;
  general_data->genstat_avg.vbondt = 0;
  general_data->genstat_avg.vbendt = 0;
  general_data->genstat_avg.vtorst = 0;
  general_data->genstat_avg.vonfot = 0;
  general_data->genstat_avg.vbend_bndt = 0;
  general_data->genstat_avg.vbend_bnd_bond = 0;
  general_data->genstat_avg.vbend_bnd_bend = 0;
  general_data->genstat_avg.vreal = 0;
  general_data->genstat_avg.vrecip = 0;
  general_data->genstat_avg.vvdw = 0;
  general_data->genstat_avg.vcoul = 0;
  general_data->genstat_avg.vlong = 0;
  general_data->genstat_avg.vbond_free = 0;
  general_data->genstat_avg.vbend_free = 0;
  general_data->genstat_avg.vtors_free = 0;
  general_data->genstat_avg.kinet = 0;
  general_data->genstat_avg.aikinet = 0;
  general_data->genstat_avg.akinet = 0;
  general_data->genstat_avg.kinet_v = 0;
  general_data->genstat_avg.aikinet_v = 0;
  general_data->genstat_avg.akinet_v = 0;
  general_data->genstat_avg.vol = 0;
  general_data->genstat_avg.aivol = 0;
  general_data->genstat_avg.avol = 0;
  general_data->genstat_avg.kinet_nhc = 0;
  general_data->genstat_avg.aikinet_nhc = 0;
  general_data->genstat_avg.akinet_nhc = 0;
  general_data->genstat_avg.kinet_nhc_bead = 0;
  general_data->genstat_avg.aikinet_nhc_bead = 0;
  general_data->genstat_avg.akinet_nhc_bead = 0;
  general_data->genstat_avg.vpot_v = 0;
  general_data->genstat_avg.vpotnhc = 0;
  general_data->genstat_avg.aiter_shake = 0;
  general_data->genstat_avg.aiter_ratl = 0;
  general_data->genstat_avg.iter_23 = 0;
  general_data->genstat_avg.iter_33 = 0;
  general_data->genstat_avg.iter_46 = 0;
  general_data->genstat_avg.iter_43 = 0;
  general_data->genstat_avg.iter_21 = 0;
  general_data->genstat_avg.aiter_23 = 0;
  general_data->genstat_avg.aiter_33 = 0;
  general_data->genstat_avg.aiter_46 = 0;
  general_data->genstat_avg.aiter_21 = 0;
  general_data->genstat_avg.aiter_43 = 0;
  general_data->genstat_avg.iter_23r = 0;
  general_data->genstat_avg.iter_33r = 0;
  general_data->genstat_avg.iter_46r = 0;
  general_data->genstat_avg.iter_43r = 0;
  general_data->genstat_avg.iter_21r = 0;
  general_data->genstat_avg.aiter_23r = 0;
  general_data->genstat_avg.aiter_33r = 0;
  general_data->genstat_avg.aiter_46r = 0;
  general_data->genstat_avg.aiter_21r = 0;
  general_data->genstat_avg.aiter_43r = 0;
  general_data->genstat_avg.acella = 0;
  general_data->genstat_avg.acellb = 0;
  general_data->genstat_avg.acellc = 0;
  general_data->genstat_avg.aicella = 0;
  general_data->genstat_avg.aicellb = 0;
  general_data->genstat_avg.aicellc = 0;
  general_data->genstat_avg.acellab = 0;
  general_data->genstat_avg.acellbc = 0;
  general_data->genstat_avg.acellac = 0;
  general_data->genstat_avg.aicellab = 0;
  general_data->genstat_avg.aicellbc = 0;
  general_data->genstat_avg.aicellac = 0;
  general_data->genstat_avg.apress = 0;
  general_data->genstat_avg.aipress = 0;
  general_data->genstat_avg.press_inter = 0;
  general_data->genstat_avg.press_intra = 0;
  general_data->genstat_avg.apress_inter = 0;
  general_data->genstat_avg.aipress_inter = 0;
  general_data->genstat_avg.apress_intra = 0;
  general_data->genstat_avg.aipress_intra = 0;
  general_data->genstat_avg.press_kin = 0;
  general_data->genstat_avg.apress_kin = 0;
  general_data->genstat_avg.aipress_kin = 0;
  general_data->genstat_avg.econv0 = 0;
  general_data->genstat_avg.econv = 0;
  general_data->genstat_avg.cpu1 = 0;
  general_data->genstat_avg.cpu2 = 0;
  general_data->genstat_avg.acpu = 0;
  general_data->genstat_avg.cpu_now = 0;
  general_data->genstat_avg.updates = 0;
  general_data->genstat_avg.updates_w = 0;
  general_data->genstat_avg.kinet_cp = 0;
  general_data->genstat_avg.aikinet_cp = 0;
  general_data->genstat_avg.akinet_cp = 0;
  general_data->genstat_avg.kinet_nhc_cp = 0;
  general_data->genstat_avg.aikinet_nhc_cp = 0;
  general_data->genstat_avg.akinet_nhc_cp = 0;
  general_data->genstat_avg.vpotnhc_cp = 0;
  general_data->genstat_avg.cp_ehart = 0;
  general_data->genstat_avg.aicp_ehart = 0;
  general_data->genstat_avg.acp_ehart = 0;
  general_data->genstat_avg.cp_eext = 0;
  general_data->genstat_avg.aicp_eext = 0;
  general_data->genstat_avg.acp_eext = 0;
  general_data->genstat_avg.cp_exc = 0;
  general_data->genstat_avg.cp_muxc = 0;
  general_data->genstat_avg.aicp_exc = 0;
  general_data->genstat_avg.acp_exc = 0;
  general_data->genstat_avg.cp_eke = 0;
  general_data->genstat_avg.aicp_eke = 0;
  general_data->genstat_avg.acp_eke = 0;
  general_data->genstat_avg.cp_enl = 0;
  general_data->genstat_avg.aicp_enl = 0;
  general_data->genstat_avg.acp_enl = 0;
  general_data->genstat_avg.aiter_shake_cp = 0;
  general_data->genstat_avg.aiter_ratl_cp = 0;
  general_data->genstat_avg.maxfc = 0;
  general_data->genstat_avg.maxf = 0;
  general_data->genstat_avg.pi_ke_prim = 0;
  general_data->genstat_avg.pi_ke_vir = 0;
  general_data->genstat_avg.api_ke_prim = 0;
  general_data->genstat_avg.api_ke_vir = 0;
  general_data->genstat_avg.aipi_ke_prim = 0;
  general_data->genstat_avg.aipi_ke_vir = 0;
  general_data->genstat_avg.kin_harm = 0;
  general_data->genstat_avg.akin_harm = 0;
  general_data->genstat_avg.aikin_harm = 0;
  general_data->genstat_avg.iter_shake = 0;


/* zero clatoms_info */
  mdclatoms_info->natm_tot           = 0;        
  mdclatoms_info->pi_beads           = 0;        
  mdclatoms_info->nfree              = 0;              
  mdclatoms_info->nchrg              = 0;              
  mdclatoms_info->pi_beads_full_ter  = 0;  
  mdclatoms_info->pi_beads_res_ter   = 0;   
  mdclatoms_info->pi_beads_full_tra  = 0;  
  mdclatoms_info->pi_beads_res_tra   = 0;   




/* zero elements of atommaps.atommaps_info  */
    mdatom_maps->nmol_typ      = 0;            
    mdatom_maps->nres_typ      = 0;            
    mdatom_maps->natm_typ      = 0;            
    mdatom_maps->nfreeze       = 0;             
    mdatom_maps->nres_max      = 0;    
    mdatom_maps->nres_typ_max  = 0;
    mdatom_maps->nres_tot      = 0;
    mdatom_maps->nres_sum      = 0;

/* zero therm_info */
  mdtherm_info->num_nhc   = 0;    
  mdtherm_info->len_nhc   = 0;    
  mdtherm_info->nres_nhc  = 0;   
  mdtherm_info->nyosh_nhc = 0;  

  mdtherm_info_bead->num_nhc   = 0;    
  mdtherm_info_bead->len_nhc   = 0;    
  mdtherm_info_bead->nres_nhc  = 0;   
  mdtherm_info_bead->nyosh_nhc = 0;  

/* zero nbr_list.nbr_list_info */
    mdnbr_list->nolst            = 0;
    mdnbr_list->iver             = 0;                
    mdnbr_list->mdverlist.iver_init        = 0;           
    mdnbr_list->mdverlist.iver_fill        = 0;  
    mdnbr_list->mdverlist.iver_count       = 0;
    mdnbr_list->mdverlist.nolst_ver_update = 0;
    mdnbr_list->mdverlist.lnk_ver_update   = 0; 
    mdnbr_list->mdverlist.nver_lst         = 0;           
    mdnbr_list->mdverlist.nver_lst_res     = 0;       
    mdnbr_list->mdverlist.jver_pad         = 0;           
    mdnbr_list->mdverlist.nmem_min_lst     = 0;       
    mdnbr_list->mdlnklist.lnk_excl_safe    = 0;      
    mdnbr_list->mdlnklist.lnk_vol_safe     = 0;       
    mdnbr_list->ilnk             = 0;               
    mdnbr_list->mdlnklist.ilnk_init        = 0;          
    mdnbr_list->mdlnklist.lnk_for_odd      = 0;        
    mdnbr_list->mdlnklist.ncell_div_avg    = 0;     
    mdnbr_list->mdlnklist.ncell_a          = 0;
    mdnbr_list->mdlnklist.ncell_b          = 0;
    mdnbr_list->mdlnklist.ncell_c          = 0;
    mdnbr_list->mdlnklist.natm_cell_max    = 0;          
    mdnbr_list->mdlnklist.lnk_list_dim     = 0;           
    mdnbr_list->mdlnklist.nshft_lnk        = 0;              
    mdnbr_list->mdlnklist.nshft_lnk_res    = 0;          

/* zero interact.interact_info  */
    mdinteract->nter_typ        = 0;   
    mdinteract->nsplin          = 0;     
    mdinteract->nsplin_tot      = 0; 
    mdinteract->iswit_vdw       = 0;  
    mdinteract->dielectric_opt  = 0; 
  
/* zero ewald.ewald_info  */
    general_data->genewald.nktot          = 0; 
    general_data->genewald.nktot_res      = 0; 

/* zero part_mesh.part_mesh_info  */
    mdpart_mesh->pme_on           = 0;                   
    mdpart_mesh->kmax_pme         = 0;                 
    mdpart_mesh->n_interp         = 0;                 
    mdpart_mesh->nktot_pme        = 0;                
    mdpart_mesh->ngrid_a          = 0;
    mdpart_mesh->ngrid_b          = 0;
    mdpart_mesh->ngrid_c          = 0;
    mdpart_mesh->pme_res_on       = 0;
    mdpart_mesh->kmax_pme_res     = 0;                 
    mdpart_mesh->n_interp_res     = 0;                  
    mdpart_mesh->nktot_pme_res    = 0;            
    mdpart_mesh->ngrid_a_res      = 0;
    mdpart_mesh->ngrid_b_res      = 0;
    mdpart_mesh->ngrid_c_res      = 0;
    mdpart_mesh->nlen_pme         = 0;                 

}/*end routine*/
//==========================================================================



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

 void zero_cp(CP *cp)

/*==========================================================================*/
{/*begin routine*/
/*==========================================================================*/

/* zero elements of  cpcoeffs_info  */
	
  cp->cpcoeffs_info.ncoef                = 0;           
  cp->cpcoeffs_info.ncoef_l              = 0;         
  cp->cpcoeffs_info.pi_beads             = 0;        
  cp->cpcoeffs_info.nstate_up            = 0;
  cp->cpcoeffs_info.nstate_dn            = 0;
  cp->cpcoeffs_info.icmass_unif          = 0;         
  cp->cpcoeffs_info.ks_rot_on            = 0;           
  cp->cpcoeffs_info.n_ks_rot             = 0;               

/* zero cptherm  */
  cp->cptherm_info.num_c_nhc       = 0;          
  cp->cptherm_info.len_c_nhc       = 0;          
  cp->cptherm_info.nres_c_nhc      = 0;         
  cp->cptherm_info.nyosh_c_nhc     = 0;        

/* zero cpewald->cpewald_info */
  cp->cpewald.nktot_sm         = 0;   


/* zero pseudo->pseudo_info  */
  cp->cppseudo.n_ang_max     = 0;       
  cp->cppseudo.nsplin_g      = 0;        
  cp->cppseudo.nsplin_g_tot  = 0;    
  cp->cppseudo.num_nl_lst    = 0;      
  cp->cppseudo.nl_cut_on     = 0;       


}/*endroutine*/
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

 void zero_par(CLASS_PARSE *class_parse,NULL_INTER_PARSE *null_inter_parse)

/*==========================================================================*/
   {/*begin routine*/
/*==========================================================================*/
/* zero class_parse->class_parse_info */

  class_parse->ivx_smpl           = 0;
  class_parse->ivnhc_smpl         = 0; 
  class_parse->kmax_ewd           = 0;
  class_parse->kmax_res           = 0;   
  class_parse->istart             = 0;              
  class_parse->ishift_pot         = 0;          
  class_parse->initial_spread_opt = 0;  
  class_parse->zero_com_vel       = 0;        



/* zero resbond_parse->resbond_parse_info  */
  /*  resbond_parse->ionfo         = 0;               
  resbond_parse->iconv         = 0;               
  resbond_parse->nresidue      = 0;
  resbond_parse->nresidue_max  = 0; 
  resbond_parse->nres_bond     = 0;
  resbond_parse->nres_bond_max = 0;*/

/* zero null_inter_parse->null_inter_parse_info */
  null_inter_parse->nbond_nul = 0;
  null_inter_parse->nbend_nul = 0;  
  null_inter_parse->ntors_nul = 0;
  null_inter_parse->nonfo_nul = 0;  


/* zero build_intra->build_intra_info */
  /*      build_intra->build_intra_info.natm_1res_pure_now  = 0;
      build_intra->natm_tot_max        = 0;
      build_intra->nfreeze_max         = 0;
      build_intra->nfreeze_now         = 0;
      build_intra->natm_1res_now       = 0;
      build_intra->natm_1res_max       = 0;
      build_intra->natmind_1res_now    = 0;
      build_intra->natmind_1res_max    = 0;
      build_intra->nghost_now          = 0;
      build_intra->nghost_tot_max      = 0;      
      build_intra->nbond_pow_max       = 0;       
      build_intra->nbond_con_max       = 0; 
      build_intra->nbond_nul_max       = 0;
      build_intra->nbond_typ_pow_max   = 0; 
      build_intra->nbond_typ_con_max   = 0;
      build_intra->nbend_pow_max       = 0;       
      build_intra->nbend_con_max       = 0; 
      build_intra->nbend_nul_max       = 0;
      build_intra->nbend_typ_pow_max   = 0; 
      build_intra->nbend_typ_con_max   = 0;
      build_intra->ntors_pow_max       = 0;       
      build_intra->ntors_con_max       = 0; 
      build_intra->ntors_nul_max       = 0;
      build_intra->ntors_typ_pow_max   = 0; 
      build_intra->ntors_typ_con_max   = 0;
      build_intra->nonfo_max           = 0;           
      build_intra->nonfo_nul_max       = 0;
      build_intra->nonfo_typ_max       = 0;
      build_intra->nbend_bnd_max       = 0;       
      build_intra->nbend_bnd_typ_max   = 0; 
      build_intra->ngrp_33_max         = 0;      
      build_intra->ngrp_21_max         = 0; 
      build_intra->ngrp_43_max         = 0; 
      build_intra->ngrp_23_max         = 0; 
      build_intra->ngrp_46_max         = 0;
      build_intra->ngrp_typ_21_max     = 0; 
      build_intra->ngrp_typ_43_max     = 0; 
      build_intra->ngrp_typ_33_max     = 0; 
      build_intra->ngrp_typ_23_max     = 0;
      build_intra->ngrp_typ_46_max     = 0;
      build_intra->natm_typ_max        = 0;
      build_intra->nres_typ_max        = 0;

      */
}/*end routine*/












