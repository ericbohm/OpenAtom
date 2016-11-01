#include "../include/class_defs/AtomPosInit.h"
#include "../include/class_defs/vx_smpl.h"

AtomPosInit::AtomPosInit () {

  PRINTF ("AtomPosInit object is being created\n");

  MDINTEGRATE  *mdintegrate  = MDINTEGRATE::get();
  MDATOMS      *mdatoms      = MDATOMS::get();
  MDINTER      *mdinter      = MDINTER::get();
  MDINTRA      *mdintra      = MDINTRA::get();
  GENERAL_DATA *general_data = GENERAL_DATA::get();
  CP           *cp           = CP::get();

  pi_beads      = (mdatoms->mdclatoms_info.pi_beads);
  natm_tot      = (mdatoms->mdclatoms_info.natm_tot);
  iextended_on  = (mdintegrate->mdtherm_info.iextended_on);

  /*
     If pi_beads is > 1, you have pi_beads independent classical systems.
     That is, if you have 256 water molecules and pi_beads=25,
     you have 256 x 25 particle positions. Now, the 25 copies don't interact
     via the intermolecule, intramolecular potentials, so far as leanMD
     is concerned you will have spawn  25 independent calculations.  The
     interactions among the copies is very simple

     PHI_PI_BEADS       = sum_ij f_i (pos[j].x[i] -pos
     [j+1].x[i])^2
     pos[pi_beads+1]. x  = pos[1].x

     Anyway, with some luck, we will get to this before too long. The
     parallelism is interesting 
   */
  mdclatoms_pos = (MDCLATOMS_POS *)
    cmalloc(pi_beads*sizeof(MDCLATOMS_POS),"main.C")-1;

  for(ip=1;ip<=pi_beads;ip++){
    mdclatoms_pos[ip].natm_tot = natm_tot;
    mdclatoms_pos[ip].x = (double *)cmalloc(natm_tot*sizeof(double),"main.C")-1;
    mdclatoms_pos[ip].y = (double *)cmalloc(natm_tot*sizeof(double),"main.C")-1;
    mdclatoms_pos[ip].z = (double *)cmalloc(natm_tot*sizeof(double),"main.C")-1;
    mdclatoms_pos[ip].vx =(double *)cmalloc(natm_tot*sizeof(double),"main.C")-1;
    mdclatoms_pos[ip].vy =(double *)cmalloc(natm_tot*sizeof(double),"main.C")-1;
    mdclatoms_pos[ip].vz =(double *)cmalloc(natm_tot*sizeof(double),"main.C")-1;
  }// endfor

  therm_bead    = (MDTHERM_POS *)
    cmalloc(pi_beads*sizeof(MDTHERM_POS),"main.C")-1;

  if(iextended_on==1){
    num_nhc = (mdintegrate->mdtherm_info.num_nhc);
    len_nhc = (mdintegrate->mdtherm_info.len_nhc);
    therm_class.num_nhc = num_nhc;
    therm_class.len_nhc = len_nhc;
    therm_class.x_nhc   = cmall_mat(1,len_nhc,1,num_nhc,"main.C");
    therm_class.v_nhc   = cmall_mat(1,len_nhc,1,num_nhc,"main.C");
  }// endif
  CkAbort("do we get here?\n");
  read_coord(mdintegrate,mdatoms,mdinter,mdintra,general_data,cp,
      mdclatoms_pos,&therm_class,therm_bead);

  for(ip=1;ip<=pi_beads;ip++){
    vx_smpl vsampler (mdclatoms_pos[ip].vx,
        mdclatoms_pos[ip].vy,
        mdclatoms_pos[ip].vz,
        mdatoms->mdclatoms_info.mass,
        mdatoms->mdclatoms_info.text_atm,
        mdatoms->mdclatoms_info.natm_tot,
        &(general_data->genensopts),
        &(mdintegrate->mdvel_samp),
        &(mdintra->mdghost_atoms),
        &(mdintra->mdconstrnt)); 
  }

}

AtomPosInit::~AtomPosInit () {

  PRINTF ("AtomPosInit object is being destroyed\n");

  if(iextended_on==1){
    cfree_mat(therm_class.x_nhc,1,len_nhc,1,num_nhc);
    cfree_mat(therm_class.v_nhc,1,len_nhc,1,num_nhc);
  }// endif

  cfree (&(therm_bead[1]),"~AtomPosInit()");

  for(ip=1;ip<=pi_beads;ip++){
    cfree (&(mdclatoms_pos[ip].x[1]),"~AtomPosInit()");
    cfree (&(mdclatoms_pos[ip].y[1]),"~AtomPosInit()");
    cfree (&(mdclatoms_pos[ip].z[1]),"~AtomPosInit()");
    cfree (&(mdclatoms_pos[ip].vx[1]),"~AtomPosInit()");
    cfree (&(mdclatoms_pos[ip].vy[1]),"~AtomPosInit()");
    cfree (&(mdclatoms_pos[ip].vz[1]),"~AtomPosInit()");
  }// endfor

  cfree (&(mdclatoms_pos[1]),"~AtomPosInit()");
}
