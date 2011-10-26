#include "standard_include.h"
#include "../../../include/Atoms.h"
#include "../../../include/debug_flags.h"
#include "../class_defs/allclass_cp.h"
#include "../class_defs/allclass_gen.h"
#include "../class_defs/allclass_mdatoms.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../class_defs/CP_OPERATIONS/class_cprspaceion.h"

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CPRSPACEION::CP_getionforce(const int natm,FastAtoms *atoms,int myid, int nproc,
                                 double *pot_ewd_ret, double *vself_ret, double *vbgr_ret,
				 double *recip_corr, PSSCRATCH *pscratch)
//============================================================================
  { // Begin Function
//----------------------------------------------------------------------------
// Local Variables

  MDATOMS      *mdatoms      = MDATOMS::get();
  GENERAL_DATA *general_data = GENERAL_DATA::get();

#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_mdatoms.h"

  double p  = 0.3614;
  double e1 = 0.2041422096422003; double  e2 = 0.1997535956961481;
  double e3 = 0.2213176596405576; double  e4 = 0.03360430734640255;
  double e5 = 0.4732592578721755; double  e6 =-0.509078520069735;
  double e7 = 0.6772631491947646; double  e8 =-0.369912979092217;
  double e9 = 0.06965131976970335;

  double de1 = 1.0*e1; double de2 = 2.0*e2; double de3 = 3.0*e3;
  double de4 = 4.0*e4; double de5 = 5.0*e5; double de6 = 6.0*e6;
  double de7 = 7.0*e7; double de8 = 8.0*e8; double de9 = 9.0*e9;


  int natm_piny  = mdclatoms_info->natm_tot; 
  double alp_ewd = genewald->alp_ewd;
  double ecut    = genewald->ecut;

  double *hmati  = gencell->hmati;
  double *hmat   = gencell->hmat;
  double  vol    = gencell->vol;
  int    iperd   = gencell->iperd;

  double *q      = atoms->q;
  double *x      = atoms->x;
  double *y      = atoms->y;
  double *z      = atoms->z;
  double *fx     = atoms->fx;
  double *fy     = atoms->fy;
  double *fz     = atoms->fz;

//============================================================================

#ifdef _CP_DEBUG_ATM_FORC_
  double *fxt = (double *)cmalloc(natm*sizeof(double),"debug ions");
  double *fyt = (double *)cmalloc(natm*sizeof(double),"debug ions");
  double *fzt = (double *)cmalloc(natm*sizeof(double),"debug ions");
  for(int i=0;i<natm;i++){
     fxt[i]=0.0;
     fyt[i]=0.0;
     fzt[i]=0.0;
  }//endfor
#endif

//============================================================================
// Primitively Parallelize the computation

   double pot_ewd  = 0.0;
   if(natm>1 && myid<natm-1){   
     int natm1 = natm-1;
     int ndiv  = MIN(natm1,nproc);
     int nsiz  = natm1/ndiv;
     int nrem  = (natm1 % ndiv);
     int ist   = nsiz*myid + MIN(nrem,myid);
     nsiz      = (myid<nrem ? nsiz+1 : nsiz);
     int iend  = ist+nsiz;

     double talp2    = 2.0*alp_ewd*alp_ewd;
     double palp     = p*alp_ewd;
     double alen     = sqrt(hmat[1]*hmat[1]+hmat[2]*hmat[2]+hmat[3]*hmat[3]);
     double blen     = sqrt(hmat[4]*hmat[4]+hmat[5]*hmat[5]+hmat[6]*hmat[6]);
     double clen     = sqrt(hmat[7]*hmat[7]+hmat[8]*hmat[8]+hmat[9]*hmat[9]);
     double hmat_min = 0.0;
     if(iperd==3){hmat_min = 0.5*MIN3(alen,blen,clen);}
     if(iperd==2){hmat_min = 0.5*MIN(alen,blen);}
     if(iperd==1){hmat_min = 0.5*alen;}

     for(int iatm = ist; iatm < iend; iatm++){
       for(int jatm = iatm+1; jatm < natm; jatm++){

         double dx = x[iatm]-x[jatm];
         double dy = y[iatm]-y[jatm];
         double dz = z[iatm]-z[jatm];
         if(iperd==3){
          double da = dx*hmati[1]+dy*hmati[4]+dz*hmati[7];
          double db = dx*hmati[2]+dy*hmati[5]+dz*hmati[8];
          double dc = dx*hmati[3]+dy*hmati[6]+dz*hmati[9];
          da -= NINT(da);
          db -= NINT(db);
          dc -= NINT(dc);
  	  dx  = da*hmat[1]+db*hmat[4]+dc*hmat[7];
          dy  = da*hmat[2]+db*hmat[5]+dc*hmat[8];
          dz  = da*hmat[3]+db*hmat[6]+dc*hmat[9];
	 }//endif
         if(iperd==2){
          double da = dx*hmati[1]+dy*hmati[4];
	  double db = dx*hmati[2]+dy*hmati[5];
          da -= NINT(da);
          db -= NINT(db);
  	  dx  = da*hmat[1]+db*hmat[4];
          dy  = da*hmat[2]+db*hmat[5];
          dz  = da*hmat[3]+db*hmat[6];
	 }//endif
         if(iperd==1){
	  double da = dx*hmati[1];
          da -= NINT(da);
  	  dx  = da*hmat[1];
	 }//endif
         double r2 = dx*dx+dy*dy+dz*dz;
         double r  = sqrt(r2);
         int igo = 0;
         double vnow=0.0,dvnow=0.0;
         if(r<=hmat_min && iperd!=0){
    	    igo = 1;
            double qij    = q[iatm]*q[jatm];
            double ralp   = r * alp_ewd;
            double eee    = exp(-ralp*ralp);
            double tt     = 1.0/(1.0+p*ralp);
            double gerfc  = ((((((((e9*tt+e8)*tt+e7)*tt+e6)*tt+e5)*tt
                            +e4)*tt+e3)*tt+e2)*tt+e1)*tt*eee;
            double dgerfc = ((((((((de9*tt+de8)*tt+de7)*tt+de6)*tt+de5)*tt
                            +de4)*tt+de3)*tt+de2)*tt+de1)*tt*tt*eee*palp
                            +talp2*gerfc*r;
            vnow          = qij * gerfc/r;
            dvnow         = (gerfc/r2 + dgerfc/r)*qij/r;
	 }//endif
         if(iperd==0){
 	    igo = 1;
            double qij = q[iatm]*q[jatm];
            vnow       = qij/r;
            dvnow      = qij/(r*r2);
	 }//endif
         if(igo==1){
            pot_ewd += vnow;
#ifdef _CP_DEBUG_ATM_FORC_
            fxt[iatm] += dx*dvnow;
            fyt[iatm] += dy*dvnow;
            fzt[iatm] += dz*dvnow;
            fxt[jatm] -= dx*dvnow;
            fyt[jatm] -= dy*dvnow;
            fzt[jatm] -= dz*dvnow;
#endif
#ifndef _CP_DEBUG_ATM_FORC_
            fx[iatm] += dx*dvnow;
            fy[iatm] += dy*dvnow;
            fz[iatm] += dz*dvnow;
            fx[jatm] -= dx*dvnow;
            fy[jatm] -= dy*dvnow;
            fz[jatm] -= dz*dvnow;
#endif
         }//endif
       }//endfor : jatm
       CmiNetworkProgress();
     }//endfor : iatm
   }//endif : there is work to do 
  pot_ewd_ret[0] = pot_ewd;


//============================================================================

  recip_corr[0] = 0.0;
  if(iperd<3 && iperd>0){
    atm_recip_corr(natm,x,y,z,fx,fy,fz,q,vol,hmat,hmati,iperd,recip_corr,
                   myid,nproc, pscratch);
  }//endif

//============================================================================
// Self and background ewald terms

  if(myid==0 && iperd !=0){
    double q2sum = 0.0;
    double qsum  = 0.0;
    for(int i=0;i<natm_piny;i++){
      q2sum += q[i]*q[i];
      qsum  += q[i];
    }//endfor
    double vself = -(q2sum*alp_ewd)/sqrt(M_PI);
    double vbgr  = -(0.5*qsum*qsum*M_PI)/(alp_ewd*alp_ewd*vol);

    double kcut = sqrt(2.0*ecut);
    double arg  = kcut/(2.0*alp_ewd);
    double eee  = exp(-arg*arg);
    double tt   = 1.0/(1.0+p*arg);
    double self_erfc = ((((((((e9*tt+e8)*tt+e7)*tt+e6)*tt+e5)*tt
                             +e4)*tt+e3)*tt+e2)*tt+e1)*tt*eee;
    double self_erf  = 1.0-self_erfc;

    vself_ret[0] = vself*self_erf;
    vbgr_ret[0]  = vbgr;
  }//endif

//============================================================================

#ifdef _CP_DEBUG_ATM_FORC_
   if(myid==0){
     FILE *fp = fopen("atom_rspace_only_forc.out","w");
       for(int i=0;i<natm;i++){
         fprintf(fp,"%d %g %g %g\n",i+1,fxt[i],fyt[i],fzt[i]);
       }//endfor
     fclose(fp);
     CkPrintf("Atom forces written to atom_forc.out don't contain\n");
     CkPrintf("the real space contribution. Those are in rspace_only.out\n");
     CkPrintf("Currently only proc 0 rspace_only forces are given.\n");
   }//endif
   cfree(fxt,"debug ions");
   cfree(fyt,"debug ions");
   cfree(fzt,"debug ions");
#endif

//============================================================================
  } // End function 
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CPRSPACEION::atm_recip_corr(int natm,double *x, double *y, double *z,
                                 double *fx, double *fy, double *fz,double *q,
                                 double vol,double *hmat,double *hmati,int iperd,
                                 double *vnow_ret,int myid, int nproc, PSSCRATCH *pscratch){
//==========================================================================

  GENERAL_DATA *general_data = GENERAL_DATA::get();
  CP *cp                     = CP::get();
#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_cp.h"
PSNONLOCAL *nonlocal = &(cppseudo->nonlocal);

   int ncorr_b         = genewald->ncorr_b;
   int ncorr_c         = genewald->ncorr_c;
   int *ka_corr_b      = genewald->ka_corr_b;
   int *kb_corr_b      = genewald->kb_corr_b;
   int *kc_corr_b      = genewald->kc_corr_b;
   int *ka_corr_c      = genewald->ka_corr_c;
   int *kb_corr_c      = genewald->kb_corr_c;
   int *kc_corr_c      = genewald->kc_corr_c;
   double *clus_corr_c = genewald->kernel_corr_c;
   double *clus_corr_b = genewald->kernel_corr_b;

   complex *ei_inc = pscratch->ei_inc+1;
   complex *h      = pscratch->ti_inc+1;
   double *sa      = pscratch->x;
   double *sb      = pscratch->y;
   double *sc      = pscratch->z;

   double *ei_incr = reinterpret_cast<double*> (ei_inc);
   double *hr      = reinterpret_cast<double*> (h);
   double *helr    = &hr[0];
   double *heli    = &hr[natm];
   double *cossc   = &ei_incr[0];
   double *sinsc   = &ei_incr[natm];

   double tpi      = 2.0*M_PI;
   double rvol     = 1.0/vol;

   double vnow = 0.0;

  if(iperd<1 || iperd>2){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Incorrect perdiodicty in atm_recip_corr\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }//endif

//==========================================================================
// Add more kc when ka and kb are small  : always

if(myid<ncorr_c){   

   // helpful vectors for increments along c
   for(int i=0;i<natm;++i){
     sa[i]      = x[i]*hmati[1] + y[i]*hmati[4] + z[i]*hmati[7];
     sb[i]      = x[i]*hmati[2] + y[i]*hmati[5] + z[i]*hmati[8];
     sc[i]      = x[i]*hmati[3] + y[i]*hmati[6] + z[i]*hmati[9];
     double arg = sc[i]*tpi;
     cossc[i]   = cos(arg);
     sinsc[i]   = sin(arg);
   }//endfor

   // cheesy parallel decomposition
   int ndiv   = MIN(ncorr_c,nproc);
   int nsiz   = ncorr_c/ndiv;
   int nrem   = (ncorr_c % ndiv);
   int ist    = nsiz*myid + MIN(nrem,myid);
   nsiz       = (myid<nrem ? nsiz+1 : nsiz);
   int iend   = ist+nsiz;

   int ka_old = ka_corr_c[ist]-1;
   int kb_old = kb_corr_c[ist]-1;
   int kc_old = kc_corr_c[ist];

   // The comp loop
   for(int ic=ist;ic<iend;ic++){
     int ka = ka_corr_c[ic];
     int kb = kb_corr_c[ic];
     int kc = kc_corr_c[ic];
     double aka = (double)ka;
     double akb = (double)kb;
     double akc = (double)kc;
     if(ka_old!=ka || kb_old !=kb ||kc_old!=kc-1){
       for(int i=0;i<natm;++i){
         double arg = tpi*(aka*sa[i]+akb*sb[i]+akc*sc[i]);
         helr[i]    = cos(arg);
         heli[i]    = sin(arg);
       }//endfor
     }//endif
     double sumr = 0.0;
     double sumi = 0.0;
     for(int i=0;i < natm;i++){
        double qi = q[i];
        sumr += helr[i]*qi;
        sumi += heli[i]*qi;
     }//endfor
     double smag   = sumr*sumr + sumi*sumi;
     double preg   = clus_corr_c[ic]*rvol; // assumes exp(-g^2/falp2)/g2 = 0
     vnow += smag*preg;
     sumr *= (preg*2.0);
     sumi *= (preg*2.0);
     double gx = (aka*hmati[1] + akb*hmati[2] + akc*hmati[3]);
     double gy = (aka*hmati[4] + akb*hmati[5] + akc*hmati[6]);
     double gz = (aka*hmati[7] + akb*hmati[8] + akc*hmati[9]);
     for(int i = 0;i < natm;++i){
       double fdiff  = (heli[i]*sumr - helr[i]*sumi)*q[i];
       fx[i] += gx*fdiff;
       fy[i] += gy*fdiff;
       fz[i] += gz*fdiff;
     }//endfor
     for(int i = 0;i<natm;++i){
        double helrt   = helr[i];
        double helit   = heli[i];
        helr[i] = helrt*cossc[i] - helit*sinsc[i];
        heli[i] = helit*cossc[i] + helrt*sinsc[i];
     }//endfor
     ka_old = ka;
     kb_old = kb;
     kc_old = kc;
   }//endfor : ic

}//endif : this proc has work to do

//==========================================================================
// Add more kb when ka and kc are small  : for perd==1

if(myid<ncorr_b && iperd==1){   

   // helpful vectors for increments along b
   for(int i=0;i<natm;++i){
     double arg = sb[i]*tpi;
     cossc[i]   = cos(arg);
     sinsc[i]   = sin(arg);
   }//endfor

   // cheesy parallel decomposition
   int ndiv   = MIN(ncorr_b,nproc);
   int nsiz   = ncorr_b/ndiv;
   int nrem   = (ncorr_b % ndiv);
   int ist    = nsiz*myid + MIN(nrem,myid);
   nsiz       = (myid<nrem ? nsiz+1 : nsiz);
   int iend   = ist+nsiz;

   int ka_old = ka_corr_b[ist]-1;
   int kb_old = kb_corr_b[ist]-1;
   int kc_old = kc_corr_b[ist];

   // The comp loop
   for(int ic=ist;ic<iend;ic++){
     int ka = ka_corr_b[ic];
     int kb = kb_corr_b[ic];
     int kc = kc_corr_b[ic];
     double aka = (double)ka;
     double akb = (double)kb;
     double akc = (double)kc;
     if(ka_old!=ka || kb_old !=kb-1 ||kc_old!=kc){
       for(int i=0;i<natm;++i){
         double arg = tpi*(aka*sa[i]+akb*sb[i]+akc*sc[i]);
         helr[i]    = cos(arg);
         heli[i]    = sin(arg);
       }//endfor
     }//endif
     double sumr = 0.0;
     double sumi = 0.0;
     for(int i=0;i < natm;i++){
        double qi = q[i];
        sumr += helr[i]*qi;
        sumi += heli[i]*qi;
     }//endfor
     double smag   = sumr*sumr + sumi*sumi;
     double preg   = clus_corr_b[ic]*rvol; // assumes exp(-g^2/falp2)/g2 = 0
     vnow += smag*preg;
     sumr *= (preg*2.0);
     sumi *= (preg*2.0);
     double gx = (aka*hmati[1] + akb*hmati[2] + akc*hmati[3]);
     double gy = (aka*hmati[4] + akb*hmati[5] + akc*hmati[6]);
     double gz = (aka*hmati[7] + akb*hmati[8] + akc*hmati[9]);
     for(int i = 0;i < natm;++i){
       double fdiff  = (heli[i]*sumr - helr[i]*sumi)*q[i];
       fx[i] += gx*fdiff;
       fy[i] += gy*fdiff;
       fz[i] += gz*fdiff;
     }//endfor
     for(int i = 0;i < natm;++i){
        double helrt   = helr[i];
        double helit   = heli[i];
        helr[i] = helrt*cossc[i] - helit*sinsc[i];
        heli[i] = helit*cossc[i] + helrt*sinsc[i];
     }//endfor
     ka_old = ka;
     kb_old = kb;
     kc_old = kc;
   }//endfor : ic

}//endif : this proc has work to do

//============================================================================
// Return the correction energy

   vnow_ret[0] = vnow;

//============================================================================
  } // End function 
//============================================================================






