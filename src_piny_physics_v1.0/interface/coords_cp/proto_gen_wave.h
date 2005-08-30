#include "../../include/class_defs/CP/class_cpylm_cons.h"

 void splin_btrans(int nsplin,double **gpsi0,double **gpsi1,
                   double **gpsi2,double **gpsi3,
                   double gmin,double gmax,double *dg,
                   double *gpsi00,int *pn_ang,char *fname_ps,int atm);

 void get_ylm(double, double, double, double,
              double * ,double *,CPYLM_CONS *ylm_cons);

 void bess_trans(double *v_rphi,int nr,double dr,double *r,
                 double *fv_rphi,int nsplin_g,double *g,
                 int iang,double *gzero);

 void get_gpsi(double g,int nsplin,
               double **gpsi0,double **gpsi1,double **gpsi2,double **gpsi3,
               double *gpsi_now,double gmin,double dg,int n_ang);

 void  fit_spline(double *c0i,double *c1i,double *c2i,
                  double *c3i,double *xi,int nsplin);

