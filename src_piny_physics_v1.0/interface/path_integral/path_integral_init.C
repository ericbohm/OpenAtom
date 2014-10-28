/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: transform                                    */
/*                                                                          */
/* This subprogram transforms between cartesian and normal mode coords      */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/



#include "standard_include.h"
#include "../class_defs/allclass_gen.h"
#include "../class_defs/allclass_mdintegrate.h"
#include "../class_defs/allclass_mdatoms.h"
#include "../class_defs/allclass_mdinter.h"
#include "../class_defs/allclass_mdintra.h"
#include "../proto_defs/proto_path_init_entry.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_friend_lib_entry.h"


#define DEBUG_OFF


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void path_integral_init(MDCLATOMS_INFO *clatoms_info,
    MDCLATOMS_PIMD *clatoms_pimd,
    MDGHOST_ATOMS *ghost_atoms, GENSIMOPTS *simopts,
    MDATOM_MAPS *atommaps,MDCONSTRNT *mdconstrnt)

  /*==========================================================================*/
  /*            Begin Routine */
{/*begin routine */
  /*==========================================================================*/
  /*           Local variable declarations                                */

  int ip,i,igloc,ighost;
  double dpi_beads;
  double beta,tau;

  /*--------------------------------------------------------------------------*/
  /*  Assign local pointers                                                   */

  int pi_md_typ        = simopts->pi_md_typ;

  int *ighost_map      = ghost_atoms->ighost_map;
  int nghost           = ghost_atoms->nghost_tot;

  int nfreeze          = mdconstrnt->nfreeze;
  int *freeze_map      = mdconstrnt->freeze_map;

  int pimd_freez_typ   = clatoms_pimd->pimd_freez_typ;

  int natm_tot         = clatoms_info->natm_tot;
  int pi_beads         = clatoms_info->pi_beads;
  double *text_atm     = clatoms_info->text_atm;
  double *mass         = clatoms_info->mass;

  double *veig,*prekf;     /* assigned after malloc*/
  double *rat1,*rat2;

  /*==========================================================================*/
  /* Malloc and assign copies of simple variables */

  clatoms_pimd->natm_tot = natm_tot;
  clatoms_pimd->pi_beads = pi_beads;
  clatoms_pimd->path_eig = (double *)cmalloc(pi_beads*sizeof(double),"path_integral_init")-1;
  clatoms_pimd->prekf    = (double *)cmalloc(natm_tot*sizeof(double),"path_integral_init")-1;
#ifdef TRANSFORMATIONS_IMPLEMENTED
  clatoms_pimd->x_trans  = (double *)cmalloc(pi_beads*sizeof(double),"path_integral_init")-1;
  clatoms_pimd->y_trans  = (double *)cmalloc(pi_beads*sizeof(double),"path_integral_init")-1;
  clatoms_pimd->z_trans  = (double *)cmalloc(pi_beads*sizeof(double),"path_integral_init")-1;
#endif

  dpi_beads = (double)pi_beads;
  prekf         = clatoms_pimd->prekf;
  veig          = clatoms_pimd->path_eig;


  /*==========================================================================*/
  /* II) Set pre kinetic energy constants */

  for(i=1;i<=natm_tot;i++){
    beta = BOLTZ/(text_atm[i]);
    tau = BOLTZ/(text_atm[i]*dpi_beads);
    prekf[i]  = mass[i]/(2.0*tau*beta);
  }/*endfor*/

  for(ighost=1;ighost <= nghost;ighost++){
    igloc = ighost_map[ighost];
    prekf[igloc] = 0.0;
  }/*endfor*/

  if(pimd_freez_typ==2){
    for(i=1;i <= nfreeze;i++){
      igloc = freeze_map[i];
      prekf[igloc] = 0.0;      
    }/*endif*/
  }/*endif*/

  /*==========================================================================*/
  /* III) Set centroid eigenvalues and use them to define the masses          */

  if(pi_md_typ == 2){
    PRINTF("Sorry no centroids, yet. Technical support busy charm++ing\n");
    EXIT(1);
  }

  /*==========================================================================*/
  /* III) Set staging eigenvalues and use them to define the masses          */

  if(pi_md_typ == 1){
    clatoms_pimd->rat1_stag = 
      (double *)cmalloc(clatoms_info->pi_beads*sizeof(double),"path_integral_init")-1;
    clatoms_pimd->rat2_stag = 
      (double *)cmalloc(clatoms_info->pi_beads*sizeof(double),"path_integral_init")-1;
    rat1 = clatoms_pimd->rat1_stag;
    rat2 = clatoms_pimd->rat2_stag;
    for(ip=2;ip<=pi_beads;ip++){
      rat1[ip] = ((double)(ip-1))/((double)(ip));
      rat2[ip] = 1.0/((double)(ip));
    }/*endfor*/
    veig[1]  = 1.0;   /* funny value for mass assignment */
    for(ip=2;ip<=pi_beads;ip++){
      veig[ip] = ((double)(ip))/((double)(ip-1));
    }/*endfor*/

    veig[1]  = 0.0;      /*True eigenvalue is 0*/

  }/*endif*/

  /*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/




