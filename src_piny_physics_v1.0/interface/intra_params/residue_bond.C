#include "standard_include.h"
#include "../class_defs/allclass_gen.h"
#include "../class_defs/allclass_mdintegrate.h"
#include "../class_defs/allclass_mdatoms.h"
#include "../class_defs/allclass_mdinter.h"
#include "../class_defs/allclass_mdintra.h"
#include "../class_defs/allclass_cp.h"
#include "../class_defs/typedefs_par.h"
#include "../proto_defs/proto_intra_params_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"

#define DEBUG_OFF
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*   residue_bonding routine:                                               */
/*==========================================================================*/

void residue_bond(MDCLATOMS_INFO *clatoms_info,
    MDATOM_MAPS *atommaps,MDINTRA *mdintra,
    RESBOND_PARSE *resbond_parse,
    BUILD_INTRA *build_intra,int jmol_typ,
    int mol_hydrog_con_opt)


  /*=======================================================================*/
  /*              Begin Routine                                          */
{ /*begin routine */
  /*=======================================================================*/
  /*              Local Variables                                          */
  int i,j,k,iii,ind[5],map[5];
  int itype,itype1,itype2,itype3,itype4;
  int ires_off,iresidue;
  int nbond_pow_new,ibond;
  int ipow,icon;
  int nbond_con_new;
  int nbend_bnd_new,ibend_bnd;
  int ntors_pow_new,itors;
  int nonfo_new,ionfo;
  int ind1,ind2,icon_opt;
  int mass_now1,mass_now2;

#include "../class_defs/allclass_strip_mdintra.h"

  /*=======================================================================*/
  /* I) Form the res-bonds */
  /*-----------------------------------------------------------------------*/
  /*  i) Count them up */

  nbond_pow_new = mdbond->npow;
  nbond_con_new = mdbond->ncon;
  for(i=1;i<=resbond_parse->nres_bond;i++){
    ind1 = resbond_parse->resbond_prm[i].iatm_res1_site;
    ind2 = resbond_parse->resbond_prm[i].iatm_res2_site;
    mass_now1 =
      (int)NINT(clatoms_info->mass[ind1]);
    mass_now2 =
      (int)NINT(clatoms_info->mass[ind2]);
    icon_opt = resbond_parse->resbond_prm[i].opt;
    /* check to see if this is any atom-H bond                      */
    if((mass_now1<=2) || (mass_now2<=2)){
      if(mol_hydrog_con_opt==1){icon_opt=1;}
    }/*endif*/
    /* check to see if this is a polar atom-H bond (polar=N,O,S)    */
    if((mass_now1<=2) && ((mass_now2==14) || (mass_now2==16)
          || (mass_now2==32))){
      if(mol_hydrog_con_opt==2){icon_opt=1;}
    }/*endif*/
    if((mass_now2<=2) && ((mass_now1==14) || (mass_now1==16)
          || (mass_now1==32))){
      if(mol_hydrog_con_opt==2){icon_opt=1;}
    }/*endif*/
    if(icon_opt==0){nbond_pow_new++;}
    if(icon_opt==1){
      nbond_con_new++;
      atommaps->nfree_1mol_jmol_typ[jmol_typ] -= 1;
      mdconstrnt->icons_jmol_typ[jmol_typ] = 1;

      ires_off = atommaps->jres_jmol_typ_strt[jmol_typ]-1;
      iresidue = resbond_parse->resbond_prm[i].res1_index;
      mdconstrnt->icons_jres_jmol_typ[(ires_off+iresidue)] = 2;
      iresidue = resbond_parse->resbond_prm[i].res2_index;
      mdconstrnt->icons_jres_jmol_typ[(ires_off+iresidue)] = 2;
    } /*endif*/
  }/*endfor*/

  /*-----------------------------------------------------------------------*/
  /*  ii) Malloc them */

  if(nbond_pow_new > build_intra->nbond_pow_max){
    build_intra->nbond_pow_max += MAX((nbond_pow_new-mdbond->npow),
        NMEM_MIN);
    mdbond->j1_pow = (int *)crealloc(&(mdbond->j1_pow)[1],
        build_intra->nbond_pow_max*sizeof(int),"residue_bond")-1;
    mdbond->j2_pow = (int *)crealloc(&(mdbond->j2_pow)[1],
        build_intra->nbond_pow_max*sizeof(int),"residue_bond")-1;
    mdbond->jtyp_pow=(int *)crealloc(&(mdbond->jtyp_pow)[1],
        build_intra->nbond_pow_max*sizeof(int),"residue_bond")-1;
  }/*endif*/
  if(nbond_con_new > build_intra->nbond_con_max){
    build_intra->nbond_con_max += MAX((nbond_con_new-mdbond->ncon),
        NMEM_MIN);
    mdbond->j1_con = (int *)crealloc(&(mdbond->j1_con)[1],
        build_intra->nbond_con_max*sizeof(int),"residue_bond")-1;
    mdbond->j2_con = (int *)crealloc(&(mdbond->j2_con)[1],
        build_intra->nbond_con_max*sizeof(int),"residue_bond")-1;
    mdbond->jtyp_con=(int *)crealloc(&(mdbond->jtyp_con)[1],
        build_intra->nbond_con_max*sizeof(int),"residue_bond")-1;
  }/*endif*/

  /*-----------------------------------------------------------------------*/
  /*  iii) Fill the pow list and check the type */
  /*       Done this way so label can be transfered to type               */

  ipow = 0;
  for(i=1;i<=resbond_parse->nres_bond;i++){
    ind1 = resbond_parse->resbond_prm[i].iatm_res1_site;
    ind2 = resbond_parse->resbond_prm[i].iatm_res2_site;
    mass_now1 =
      (int)NINT(clatoms_info->mass[ind1]);
    mass_now2 =
      (int)NINT(clatoms_info->mass[ind2]);
    icon_opt = resbond_parse->resbond_prm[i].opt;
    /* check to see if this is any atom-H bond                      */
    if((mass_now1<=2) || (mass_now2<=2)){
      if(mol_hydrog_con_opt==1){icon_opt=1;}
    }/*endif*/
    /* check to see if this is a polar atom-H bond (polar=N,O,S)    */
    if((mass_now1<=2) && ((mass_now2==14) || (mass_now2==16)
          || (mass_now2==32))){
      if(mol_hydrog_con_opt==2){icon_opt=1;}
    }/*endif*/
    if((mass_now2<=2) && ((mass_now1==14) || (mass_now1==16)
          || (mass_now1==32))){
      if(mol_hydrog_con_opt==2){icon_opt=1;}
    }/*endif*/
    if(icon_opt==0){
      ipow++;
      ibond = ipow + mdbond->npow;
      mdbond->j1_pow[ibond]=resbond_parse->resbond_prm[i].iatm_res1_site;
      mdbond->j2_pow[ibond]=resbond_parse->resbond_prm[i].iatm_res2_site;
      itype1 = atommaps->iatm_atm_typ[(mdbond->j1_pow[ibond])];
      itype2 = atommaps->iatm_atm_typ[(mdbond->j2_pow[ibond])];
      strcpy(build_intra->cbond_typ_now->atm1,atommaps->atm_typ[itype1]);
      strcpy(build_intra->cbond_typ_now->atm2,atommaps->atm_typ[itype2]);
      strcpy(build_intra->cbond_typ_now->label,
          resbond_parse->resbond_prm[i].label);
      itype = mdbond->ntyp_pow+1;
      for(j=1;j<=(mdbond->ntyp_pow);j++){
        if((strcasecmp(build_intra->cbond_typ_pow[j].atm1,
                build_intra->cbond_typ_now->atm1)==0)
            &&(strcasecmp(build_intra->cbond_typ_pow[j].atm2,
                build_intra->cbond_typ_now->atm2)==0)
            &&(strcasecmp(build_intra->cbond_typ_pow[j].label,
                build_intra->cbond_typ_now->label)==0)){itype=j;}
        if((strcasecmp(build_intra->cbond_typ_pow[j].atm1,
                build_intra->cbond_typ_now->atm2)==0)
            &&(strcasecmp(build_intra->cbond_typ_pow[j].atm2,
                build_intra->cbond_typ_now->atm1)==0)
            &&(strcasecmp(build_intra->cbond_typ_pow[j].label,
                build_intra->cbond_typ_now->label)==0)){itype=j;}
      }/*endfor*/
      if(itype>build_intra->nbond_typ_pow_max){
        build_intra->nbond_typ_pow_max += NMEM_MIN;
        build_intra->cbond_typ_pow      = (CBOND *) 
          crealloc(&(build_intra->cbond_typ_pow)[1],
              build_intra->nbond_typ_pow_max*sizeof(CBOND),"residue_bond")-1;
      }/*endif*/
      if(itype==mdbond->ntyp_pow+1){
        mdbond->ntyp_pow++;
        strcpy(build_intra->cbond_typ_pow[itype].atm1,
            build_intra->cbond_typ_now->atm1);
        strcpy(build_intra->cbond_typ_pow[itype].atm2,
            build_intra->cbond_typ_now->atm2);
        strcpy(build_intra->cbond_typ_pow[itype].label,
            build_intra->cbond_typ_now->label);
      }/*endif*/
      mdbond->jtyp_pow[ibond] = itype;
    }/*endif:pow bond*/
  }/*endfor*/
  mdbond->npow = nbond_pow_new; 

  /*-----------------------------------------------------------------------*/
  /*  iii) Fill the con list and check the type */
  /*       Done this way so label can be transfered to type               */

  icon = 0;
  for(i=1;i<=resbond_parse->nres_bond;i++){
    ind1 = resbond_parse->resbond_prm[i].iatm_res1_site;
    ind2 = resbond_parse->resbond_prm[i].iatm_res2_site;
    mass_now1 =
      (int)NINT(clatoms_info->mass[ind1]);
    mass_now2 =
      (int)NINT(clatoms_info->mass[ind2]);
    icon_opt = resbond_parse->resbond_prm[i].opt;
    /* check to see if this is any atom-H bond                      */
    if((mass_now1<=2) || (mass_now2<=2)){
      if(mol_hydrog_con_opt==1){icon_opt=1;}
    }/*endif*/
    /* check to see if this is a polar atom-H bond (polar=N,O,S)    */
    if((mass_now1<=2) && ((mass_now2==14) || (mass_now2==16)
          || (mass_now2==32))){
      if(mol_hydrog_con_opt==2){icon_opt=1;}
    }/*endif*/
    if((mass_now2<=2) && ((mass_now1==14) || (mass_now1==16)
          || (mass_now1==32))){
      if(mol_hydrog_con_opt==2){icon_opt=1;}
    }/*endif*/
    if(icon_opt==1){
      icon++;
      ibond = icon + mdbond->ncon;
      mdbond->j1_con[ibond]=resbond_parse->resbond_prm[i].iatm_res1_site;
      mdbond->j2_con[ibond]=resbond_parse->resbond_prm[i].iatm_res2_site;
      ibond = i;
      itype1 = atommaps->iatm_atm_typ[(mdbond->j1_con[ibond])];
      itype2 = atommaps->iatm_atm_typ[(mdbond->j2_con[ibond])];
      strcpy(build_intra->cbond_typ_now->atm1,atommaps->atm_typ[itype1]);
      strcpy(build_intra->cbond_typ_now->atm2,atommaps->atm_typ[itype2]);
      strcpy(build_intra->cbond_typ_now->label,
          resbond_parse->resbond_prm[i].label);
      itype = mdbond->ntyp_con+1;
      for(j=1;j<=mdbond->ntyp_con;j++){
        if((strcasecmp(build_intra->cbond_typ_con[j].atm1,
                build_intra->cbond_typ_now->atm1)==0)
            &&(strcasecmp(build_intra->cbond_typ_con[j].atm2,
                build_intra->cbond_typ_now->atm2)==0)
            &&(strcasecmp(build_intra->cbond_typ_con[j].label,
                build_intra->cbond_typ_now->label)==0)){itype=j;}
        if((strcasecmp(build_intra->cbond_typ_con[j].atm1,
                build_intra->cbond_typ_now->atm2)==0)
            &&(strcasecmp(build_intra->cbond_typ_con[j].atm2,
                build_intra->cbond_typ_now->atm1)==0)
            &&(strcasecmp(build_intra->cbond_typ_con[j].label,
                build_intra->cbond_typ_now->label)==0)){itype=j;}
      }/*endfor*/
      if(itype>build_intra->nbond_typ_con_max){
        build_intra->nbond_typ_con_max += NMEM_MIN;
        build_intra->cbond_typ_con      = (CBOND *) 
          crealloc(&(build_intra->cbond_typ_con)[1],
              build_intra->nbond_typ_con_max*sizeof(CBOND),"residue_bond")-1;
      }/*endif*/
      if(itype==(mdbond->ntyp_con)+1){
        mdbond->ntyp_con++;
        strcpy(build_intra->cbond_typ_con[itype].atm1,
            build_intra->cbond_typ_now->atm1);
        strcpy(build_intra->cbond_typ_con[itype].atm2,
            build_intra->cbond_typ_now->atm2);
        strcpy(build_intra->cbond_typ_con[itype].label,
            build_intra->cbond_typ_now->label);
      }/*endif*/
      mdbond->jtyp_con[ibond] = itype;
    }/*endif:con bond*/
  }/*endfor*/
  mdbond->ncon = nbond_con_new; 

  /*=======================================================================*/
  /* II) Form the res-bend_bnds */ 
  /*-----------------------------------------------------------------------*/
  /* i) Count the bend_bnds */

  nbend_bnd_new = mdbend_bnd->num;
  for(i=1;i<=resbond_parse->nres_bond;i++){
    for(j=1;j<=resbond_parse->resbond_prm[i].natm_res1_1back;j++){
      nbend_bnd_new++;
    }/*endfor*/
    for(j=1;j<=resbond_parse->resbond_prm[i].natm_res2_1back;j++){
      nbend_bnd_new++;
    } /*endfor*/
  }/*endfor*/

  /*-----------------------------------------------------------------------*/
  /* ii) Malloc */

  if(nbend_bnd_new > build_intra->nbend_bnd_max){
    build_intra->nbend_bnd_max += 
      MAX((nbend_bnd_new-mdbend_bnd->num),NMEM_MIN);
    mdbend_bnd->j1 =(int*)crealloc(&(mdbend_bnd->j1)[1], 
        build_intra->nbend_bnd_max*sizeof(int),"residue_bond")-1;
    mdbend_bnd->j2 =(int*)crealloc(&(mdbend_bnd->j2)[1],
        build_intra->nbend_bnd_max*sizeof(int),"residue_bond")-1;
    mdbend_bnd->j3 =(int*)crealloc(&(mdbend_bnd->j3)[1],
        build_intra->nbend_bnd_max*sizeof(int),"residue_bond")-1;
    mdbend_bnd->jtyp=(int*)crealloc(&(mdbend_bnd->jtyp)[1],
        build_intra->nbend_bnd_max*sizeof(int),"residue_bond")-1;
  }/*endif*/

  /*-----------------------------------------------------------------------*/
  /* iii) Fill the bend_bnd list */

  ibend_bnd = mdbend_bnd->num;
  for(i=1;i<=resbond_parse->nres_bond;i++){
    /*a) 1 back from residue 1 */
    for(j=1;j<=resbond_parse->resbond_prm[i].natm_res1_1back;j++){
      ibend_bnd++;
      mdbend_bnd->j1[ibend_bnd]= 
        resbond_parse->resbond_prm[i].iatm_res2_site;
      mdbend_bnd->j2[ibend_bnd]=
        resbond_parse->resbond_prm[i].iatm_res1_site;
      mdbend_bnd->j3[ibend_bnd]=
        resbond_parse->resbond_prm[i].iatm_res1_1back[j];
    }/*endfor*/
    /*b) 1 back from residue 2 */
    for(j=1;j<=resbond_parse->resbond_prm[i].natm_res2_1back;j++){
      ibend_bnd++;
      mdbend_bnd->j1[ibend_bnd]=
        resbond_parse->resbond_prm[i].iatm_res1_site;
      mdbend_bnd->j2[ibend_bnd]=
        resbond_parse->resbond_prm[i].iatm_res2_site;
      mdbend_bnd->j3[ibend_bnd]=
        resbond_parse->resbond_prm[i].iatm_res2_1back[j];
    }/*endfor*/
  }/*endfor*/

  /*-----------------------------------------------------------------------*/
  /* iv) Check the type */

  for(i=mdbend_bnd->num+1;i<=nbend_bnd_new;i++){
    ibend_bnd = i;
    itype1 = atommaps->iatm_atm_typ[mdbend_bnd->j1[ibend_bnd]];
    itype2 = atommaps->iatm_atm_typ[mdbend_bnd->j2[ibend_bnd]];
    itype3 = atommaps->iatm_atm_typ[mdbend_bnd->j3[ibend_bnd]];
    strcpy(build_intra->cbend_bnd_typ_now->atm1,atommaps->atm_typ[itype1]);
    strcpy(build_intra->cbend_bnd_typ_now->atm2,atommaps->atm_typ[itype2]);
    strcpy(build_intra->cbend_bnd_typ_now->atm3,atommaps->atm_typ[itype3]);
    strcpy(build_intra->cbend_bnd_typ_now->label,"");
    itype = (mdbend_bnd->ntyp)+1;
    for(j=1;j<=mdbend_bnd->ntyp;j++){
      if((strcasecmp(build_intra->cbend_bnd_typ[j].atm1,
              build_intra->cbend_bnd_typ_now->atm1)==0)
          &&(strcasecmp(build_intra->cbend_bnd_typ[j].atm2,
              build_intra->cbend_bnd_typ_now->atm2)==0)
          &&(strcasecmp(build_intra->cbend_bnd_typ[j].atm3,
              build_intra->cbend_bnd_typ_now->atm3)==0)
          &&(strcasecmp(build_intra->cbend_bnd_typ[j].label,
              build_intra->cbend_bnd_typ_now->label)==0)
        ){itype=j;}
      if((strcasecmp(build_intra->cbend_bnd_typ[j].atm1,
              build_intra->cbend_bnd_typ_now->atm3)==0)
          &&(strcasecmp(build_intra->cbend_bnd_typ[j].atm2,
              build_intra->cbend_bnd_typ_now->atm2)==0)
          &&(strcasecmp(build_intra->cbend_bnd_typ[j].atm3,
              build_intra->cbend_bnd_typ_now->atm1)==0)
          &&(strcasecmp(build_intra->cbend_bnd_typ[j].label,
              build_intra->cbend_bnd_typ_now->label)==0)
        ){itype=j;}
    }/*endfor*/
    if(itype>build_intra->nbend_bnd_typ_max){
      build_intra->nbend_bnd_typ_max += NMEM_MIN;
      build_intra->cbend_bnd_typ      = (CBEND *) 
        crealloc(&(build_intra->cbend_bnd_typ)[1],
            build_intra->nbend_bnd_typ_max*sizeof(CBEND),"residue_bond")-1;
    }/*endif*/
    if(itype==(mdbend_bnd->ntyp)+1){
      mdbend_bnd->ntyp+=1;
      strcpy(build_intra->cbend_bnd_typ[itype].atm1,
          build_intra->cbend_bnd_typ_now->atm1);
      strcpy(build_intra->cbend_bnd_typ[itype].atm2,
          build_intra->cbend_bnd_typ_now->atm2);
      strcpy(build_intra->cbend_bnd_typ[itype].atm3,
          build_intra->cbend_bnd_typ_now->atm3);
      strcpy(build_intra->cbend_bnd_typ[itype].label,
          build_intra->cbend_bnd_typ_now->label);
    }/*endif*/
    mdbend_bnd->jtyp[ibend_bnd] = itype;
  }/*endfor*/
  mdbend_bnd->num = nbend_bnd_new;

  /*=======================================================================*/
  /* III) Form the res-torss */ 
  /*-----------------------------------------------------------------------*/
  /* i) count new torsion interactions */

  ntors_pow_new = mdtors->npow;
  for(i=1;i<=resbond_parse->nres_bond;i++){
    /*------------------------------------------------*/
    /* a) at bond site on res 1 and two away on res 2 */
    for(j=1;j<=resbond_parse->resbond_prm[i].natm_res1_1back;j++){
      for(k=1;k<=resbond_parse->resbond_prm[i].natm_res1_2back[j];k++){
        ntors_pow_new++;
      }/*endfor*/
    }/*endfor*/
    /*------------------------------------------------*/
    /* b) at bond site on res 2 and two away on res 1 */
    for(j=1;j<=resbond_parse->resbond_prm[i].natm_res2_1back;j++){
      for(k=1;k<=resbond_parse->resbond_prm[i].natm_res2_2back[j];k++){
        ntors_pow_new++;
      }/*endfor*/
    }/*endfor*/
    /*------------------------------------------------*/
    /* c) one away on bond site on res 1 and res 2 */
    for(j=1;j<=resbond_parse->resbond_prm[i].natm_res2_1back;j++){
      for(k=1;k<=resbond_parse->resbond_prm[i].natm_res1_1back;k++){
        ntors_pow_new++;
      }/*endfor*/
    }/*endfor*/
    /*------------------------------------------------*/
    /* d) improper if bond site on res 1 has three bonds */
    if(resbond_parse->resbond_prm[i].natm_res1_1back==2) ntors_pow_new++;
    /*------------------------------------------------*/
    /* e) improper if bond site on res 2 has three bonds */
    if(resbond_parse->resbond_prm[i].natm_res2_1back==2) ntors_pow_new++;
  }/*endfor*/

  /*-----------------------------------------------------------------------*/
  /* ii) Malloc them */

  if(ntors_pow_new > build_intra->ntors_pow_max){
    build_intra->ntors_pow_max += MAX((ntors_pow_new-mdtors->npow),
        NMEM_MIN);
    mdtors->j1_pow =(int*)crealloc(&(mdtors->j1_pow)[1],
        build_intra->ntors_pow_max*sizeof(int),"residue_bond")-1;
    mdtors->j2_pow =(int*)crealloc(&(mdtors->j2_pow)[1],
        build_intra->ntors_pow_max*sizeof(int),"residue_bond")-1;
    mdtors->j3_pow =(int*)crealloc(&(mdtors->j3_pow)[1],
        build_intra->ntors_pow_max*sizeof(int),"residue_bond")-1;
    mdtors->j4_pow =(int*)crealloc(&(mdtors->j4_pow)[1],
        build_intra->ntors_pow_max*sizeof(int),"residue_bond")-1;
    mdtors->jtyp_pow=(int*)crealloc(&(mdtors->jtyp_pow)[1],
        build_intra->ntors_pow_max*sizeof(int),"residue_bond")-1;
  }/*endif*/

  /*-----------------------------------------------------------------------*/
  /* iii) Fill the vectors */

  itors = mdtors->npow;
  for(i=1;i<=resbond_parse->nres_bond;i++){
    /*------------------------------------------------*/
    /* a)at bond site on res 1 and two away on res 2 */
    for(j=1;j<=resbond_parse->resbond_prm[i].natm_res1_1back;j++){
      for(k=1;k<=resbond_parse->resbond_prm[i].natm_res1_2back[j];k++){
        itors++;
        mdtors->j1_pow[itors] = 
          resbond_parse->resbond_prm[i].iatm_res2_site;
        mdtors->j2_pow[itors] = 
          resbond_parse->resbond_prm[i].iatm_res1_site;
        mdtors->j3_pow[itors] = 
          resbond_parse->resbond_prm[i].iatm_res1_1back[j];
        mdtors->j4_pow[itors] = 
          resbond_parse->resbond_prm[i].iatm_res1_2back[j][k];
        mdtors->jtyp_pow[itors] = 0;
      }/*endfor*/
    }/*endfor*/
    /*------------------------------------------------*/
    /* b)at bond site on res 2 and two away on res 1 */
    for(j=1;j<=resbond_parse->resbond_prm[i].natm_res2_1back;j++){
      for(k=1;k<=resbond_parse->resbond_prm[i].natm_res2_2back[j];k++){
        itors++;
        mdtors->j1_pow[itors] = 
          resbond_parse->resbond_prm[i].iatm_res1_site;
        mdtors->j2_pow[itors] = 
          resbond_parse->resbond_prm[i].iatm_res2_site;
        mdtors->j3_pow[itors] = 
          resbond_parse->resbond_prm[i].iatm_res2_1back[j];
        mdtors->j4_pow[itors] = 
          resbond_parse->resbond_prm[i].iatm_res2_2back[j][k];
        mdtors->jtyp_pow[itors] = 0;
      }/*endfor*/
    }/*endfor*/
    /*------------------------------------------------*/
    /* c) one away on bond site on res 1 and res 2 */
    for(j=1;j<=resbond_parse->resbond_prm[i].natm_res2_1back;j++){
      for(k=1;k<=resbond_parse->resbond_prm[i].natm_res1_1back;k++){
        itors++;
        mdtors->j1_pow[itors] = 
          resbond_parse->resbond_prm[i].iatm_res2_1back[j];
        mdtors->j2_pow[itors] = 
          resbond_parse->resbond_prm[i].iatm_res2_site;
        mdtors->j3_pow[itors] = 
          resbond_parse->resbond_prm[i].iatm_res1_site;
        mdtors->j4_pow[itors] = 
          resbond_parse->resbond_prm[i].iatm_res1_1back[k];
        mdtors->jtyp_pow[itors] = 0;
      }/*endfor*/
    }/*endfor*/
    /*--------------------------------------------------*/
    /* d) improper if bond site on res 1 has three bonds */
    map[1] = resbond_parse->resbond_prm[i].improp_res1[1];
    map[2] = resbond_parse->resbond_prm[i].improp_res1[2];
    map[3] = resbond_parse->resbond_prm[i].improp_res1[3];
    map[4] = resbond_parse->resbond_prm[i].improp_res1[4];
    if((resbond_parse->resbond_prm[i].natm_res1_1back==2)&&
        (map[1]+map[2]+map[3]+map[4])>0){
      itors++;
      mdtors->nimpr++;
      ind[1] = resbond_parse->resbond_prm[i].iatm_res1_site;
      ind[2] = resbond_parse->resbond_prm[i].iatm_res1_1back[1];
      ind[3] = resbond_parse->resbond_prm[i].iatm_res1_1back[2];
      ind[4] = resbond_parse->resbond_prm[i].iatm_res2_site;
      mdtors->j1_pow[itors] = ind[map[1]];
      mdtors->j2_pow[itors] = ind[map[2]];
      mdtors->j3_pow[itors] = ind[map[3]];
      mdtors->j4_pow[itors] = ind[map[4]];
      mdtors->jtyp_pow[itors] = -1;
#ifdef DEBUG
      PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
      PRINTF("WARNING: adding improper torsion for atoms %d %d %d %d\n",
          mdtors->j1_pow[itors],mdtors->j2_pow[itors],
          mdtors->j3_pow[itors],mdtors->j4_pow[itors]);
      PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
      FFLUSH(stdout);
#endif
    }/*endif*/
    /*--------------------------------------------------*/
    /* e) improper if bond site on res 2 has three bonds */
    map[1] = resbond_parse->resbond_prm[i].improp_res2[1];
    map[2] = resbond_parse->resbond_prm[i].improp_res2[2];
    map[3] = resbond_parse->resbond_prm[i].improp_res2[3];
    map[4] = resbond_parse->resbond_prm[i].improp_res2[4];
    if((resbond_parse->resbond_prm[i].natm_res2_1back==2)&&
        (map[1]+map[2]+map[3]+map[4])>0){
      itors++;
      mdtors->nimpr++;
      ind[1] = resbond_parse->resbond_prm[i].iatm_res2_site;
      ind[2] = resbond_parse->resbond_prm[i].iatm_res2_1back[1];
      ind[3] = resbond_parse->resbond_prm[i].iatm_res2_1back[2];
      ind[4] = resbond_parse->resbond_prm[i].iatm_res1_site;
      mdtors->j1_pow[itors] = ind[map[1]];
      mdtors->j2_pow[itors] = ind[map[2]];
      mdtors->j3_pow[itors] = ind[map[3]];
      mdtors->j4_pow[itors] = ind[map[4]];
      mdtors->jtyp_pow[itors] = -1;
#ifdef DEBUG
      PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
      PRINTF("WARNING: adding improper torsion for atoms %d %d %d %d\n",
          mdtors->j1_pow[itors],mdtors->j2_pow[itors],
          mdtors->j3_pow[itors],mdtors->j4_pow[itors]);
      PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
      FFLUSH(stdout);
#endif
    }/*endif*/
  }/*endfor: over the resbonds*/

  /*-----------------------------------------------------------------------*/
  /* iv) Check the types */

  for(itors=mdtors->npow+1;itors<=ntors_pow_new;itors++){
    itype1 = atommaps->iatm_atm_typ[mdtors->j1_pow[itors]];
    itype2 = atommaps->iatm_atm_typ[mdtors->j2_pow[itors]];
    itype3 = atommaps->iatm_atm_typ[mdtors->j3_pow[itors]];
    itype4 = atommaps->iatm_atm_typ[mdtors->j4_pow[itors]];
    strcpy(build_intra->ctors_typ_now->atm1,atommaps->atm_typ[itype1]);
    strcpy(build_intra->ctors_typ_now->atm2,atommaps->atm_typ[itype2]);
    strcpy(build_intra->ctors_typ_now->atm3,atommaps->atm_typ[itype3]);
    strcpy(build_intra->ctors_typ_now->atm4,atommaps->atm_typ[itype4]);
    if(mdtors->jtyp_pow[itors]==0){
      strcpy(build_intra->ctors_typ_now->label,"");
    }
    if(mdtors->jtyp_pow[itors]==-1){
      strcpy(build_intra->ctors_typ_now->label,"improper");
    }    

    itype = (mdtors->ntyp_pow)+1;
    for(j=1;j<=mdtors->ntyp_pow;j++){
      if((strcasecmp(build_intra->ctors_typ_pow[j].atm1,
              build_intra->ctors_typ_now->atm1)==0)
          &&(strcasecmp(build_intra->ctors_typ_pow[j].atm2,
              build_intra->ctors_typ_now->atm2)==0)
          &&(strcasecmp(build_intra->ctors_typ_pow[j].atm3,
              build_intra->ctors_typ_now->atm3)==0)
          &&(strcasecmp(build_intra->ctors_typ_pow[j].atm4,
              build_intra->ctors_typ_now->atm4)==0)
          &&(strcasecmp(build_intra->ctors_typ_pow[j].label,
              build_intra->ctors_typ_now->label)==0)
        ){itype=j;}
      if((strcasecmp(build_intra->ctors_typ_pow[j].atm1,
              build_intra->ctors_typ_now->atm4)==0)
          &&(strcasecmp(build_intra->ctors_typ_pow[j].atm2,
              build_intra->ctors_typ_now->atm3)==0)
          &&(strcasecmp(build_intra->ctors_typ_pow[j].atm3,
              build_intra->ctors_typ_now->atm2)==0)
          &&(strcasecmp(build_intra->ctors_typ_pow[j].atm4,
              build_intra->ctors_typ_now->atm1)==0)
          &&(strcasecmp(build_intra->ctors_typ_pow[j].label,
              build_intra->ctors_typ_now->label)==0)
        ){itype=j;}
    }/*endfor*/
    if(itype>build_intra->ntors_typ_pow_max){
      build_intra->ntors_typ_pow_max += NMEM_MIN;
      build_intra->ctors_typ_pow      = (CTORS *) 
        crealloc(&(build_intra->ctors_typ_pow)[1],
            build_intra->ntors_typ_pow_max*sizeof(CTORS),"residue_bond")-1;
    }/*endif*/
    if(itype==(mdtors->ntyp_pow)+1){
      mdtors->ntyp_pow+=1;
      strcpy(build_intra->ctors_typ_pow[itype].atm1,
          build_intra->ctors_typ_now->atm1);
      strcpy(build_intra->ctors_typ_pow[itype].atm2,
          build_intra->ctors_typ_now->atm2);
      strcpy(build_intra->ctors_typ_pow[itype].atm3,
          build_intra->ctors_typ_now->atm3);
      strcpy(build_intra->ctors_typ_pow[itype].atm4,
          build_intra->ctors_typ_now->atm4);
      strcpy(build_intra->ctors_typ_pow[itype].label,
          build_intra->ctors_typ_now->label);
    }/*endif*/
    mdtors->jtyp_pow[itors] = itype;
  }/*endfor*/
  mdtors->npow = ntors_pow_new;

  /*=======================================================================*/
  /* III) Form the res-onfo if set */ 
  if(resbond_parse->ionfo==1){    
    /*-----------------------------------------------------------------------*/
    /* i) Count the onfos */

    nonfo_new = mdonfo->num;
    for(i=1;i<=resbond_parse->nres_bond;i++){
      /*------------------------------------------------*/
      /* a) at bond site on res 1 and two away on res 2 */
      for(j=1;j<=resbond_parse->resbond_prm[i].natm_res1_1back;j++){
        for(k=1;k<=resbond_parse->resbond_prm[i].natm_res1_2back[j];k++){
          nonfo_new++;
        }/*endfor*/
      }/*endfor*/
      /*------------------------------------------------*/
      /* b) at bond site on res 2 and two away on res 1 */
      for(j=1;j<=resbond_parse->resbond_prm[i].natm_res2_1back;j++){
        for(k=1;k<=resbond_parse->resbond_prm[i].natm_res2_2back[j];k++){
          nonfo_new++;
        }/*endfor*/
      }/*endfor*/
      /*------------------------------------------------*/
      /* c) on away on bond site on res 1 and res 2 */
      for(j=1;j<=resbond_parse->resbond_prm[i].natm_res2_1back;j++){
        for(k=1;k<=resbond_parse->resbond_prm[i].natm_res1_1back;k++){
          nonfo_new++;
        }/*endfor*/
      }/*endfor*/
    }/*endfor*/

    /*-----------------------------------------------------------------------*/
    /* ii) Malloc them */

    if(nonfo_new > build_intra->nonfo_max){
      build_intra->nonfo_max += MAX((nonfo_new-mdonfo->num),
          NMEM_MIN);
      mdonfo->j1 =(int *) crealloc(&(mdonfo->j1)[1],
          build_intra->nonfo_max*sizeof(int),"residue_bond")-1;
      mdonfo->j2 =(int *) crealloc(&(mdonfo->j2)[1],
          build_intra->nonfo_max*sizeof(int),"residue_bond")-1;
      mdonfo->jtyp =(int *) crealloc(&(mdonfo->jtyp)[1],
          build_intra->nonfo_max*sizeof(int),"residue_bond")-1;
    }/*endif*/

    /*-----------------------------------------------------------------------*/
    /* iii) Fill the onfo list */

    ionfo = mdonfo->num;
    for(i=1;i<=resbond_parse->nres_bond;i++){
      /*------------------------------------------------*/
      /* a)at bond site on res 1 and two away on res 2 */
      for(j=1;j<=resbond_parse->resbond_prm[i].natm_res1_1back;j++){
        for(k=1;k<=resbond_parse->resbond_prm[i].natm_res1_2back[j];k++){
          ionfo++;
          mdonfo->j1[ionfo] = 
            resbond_parse->resbond_prm[i].iatm_res2_site;
          mdonfo->j2[ionfo] = 
            resbond_parse->resbond_prm[i].iatm_res1_2back[j][k];
        }/*endfor*/
      }/*endfor*/
      /*------------------------------------------------*/
      /* b)at bond site on res 2 and two away on res 1 */
      for(j=1;j<=resbond_parse->resbond_prm[i].natm_res2_1back;j++){
        for(k=1;k<=resbond_parse->resbond_prm[i].natm_res2_2back[j];k++){
          ionfo++;
          mdonfo->j1[ionfo] = 
            resbond_parse->resbond_prm[i].iatm_res1_site;
          mdonfo->j2[ionfo] = 
            resbond_parse->resbond_prm[i].iatm_res2_2back[j][k];
        }/*endfor*/
      }/*endfor*/
      /*------------------------------------------------*/
      /* c) on away on bond site on res 1 and res 2 */
      for(j=1;j<=resbond_parse->resbond_prm[i].natm_res2_1back;j++){
        for(k=1;k<=resbond_parse->resbond_prm[i].natm_res1_1back;k++){
          ionfo++;
          mdonfo->j1[ionfo] = 
            resbond_parse->resbond_prm[i].iatm_res2_1back[j];
          mdonfo->j2[ionfo] = 
            resbond_parse->resbond_prm[i].iatm_res1_1back[k];
        }/*endfor*/
      }/*endfor*/
    }/*endfor*/
    /*-----------------------------------------------------------------------*/
    /* iv) Compare the types */

    for(ionfo=mdonfo->num+1;ionfo<=nonfo_new;ionfo++){
      itype1 = atommaps->iatm_atm_typ[mdonfo->j1[ionfo]];
      itype2 = atommaps->iatm_atm_typ[mdonfo->j2[ionfo]];
      strcpy(build_intra->confo_typ_now->atm1,atommaps->atm_typ[itype1]);
      strcpy(build_intra->confo_typ_now->atm2,atommaps->atm_typ[itype2]);
      strcpy(build_intra->confo_typ_now->label,"");
      itype = (mdonfo->ntyp)+1;
      for(j=1;j<=mdonfo->ntyp;j++){
        if((strcasecmp(build_intra->confo_typ[j].atm1,
                build_intra->confo_typ_now->atm1)==0)
            &&(strcasecmp(build_intra->confo_typ[j].atm2,
                build_intra->confo_typ_now->atm2)==0)
            &&(strcasecmp(build_intra->confo_typ[j].label,
                build_intra->confo_typ_now->label)==0)
          ) {itype=j;}
        if((strcasecmp(build_intra->confo_typ[j].atm1,
                build_intra->confo_typ_now->atm2)==0)
            &&(strcasecmp(build_intra->confo_typ[j].atm2,
                build_intra->confo_typ_now->atm1)==0)
            &&(strcasecmp(build_intra->confo_typ[j].label,
                build_intra->confo_typ_now->label)==0)
          ){itype=j;}
      }/*endfor*/
      if(itype>build_intra->nonfo_typ_max){
        build_intra->nonfo_typ_max += NMEM_MIN;
        build_intra->confo_typ = (CBOND *) 
          crealloc(&(build_intra->confo_typ)[1],
              build_intra->nonfo_typ_max*sizeof(CBOND),"residue_bond")-1;
      }/*endif*/
      if(itype==(mdonfo->ntyp)+1){
        mdonfo->ntyp++;
        strcpy(build_intra->confo_typ[itype].atm1,
            build_intra->confo_typ_now->atm1);
        strcpy(build_intra->confo_typ[itype].atm2,
            build_intra->confo_typ_now->atm2);
        strcpy(build_intra->confo_typ[itype].label,
            build_intra->confo_typ_now->label);
      }/*endif*/
      mdonfo->jtyp[ionfo] = itype;
    }/*endfor*/
    mdonfo->num = nonfo_new;
    /*-----------------------------------------------------------------------*/
  }/*endif:onfos on*/

  /*=======================================================================*/
}/*end routine*/
/*==========================================================================*/







