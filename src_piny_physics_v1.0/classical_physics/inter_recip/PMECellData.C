#include "PMECellData.h"

PMECellData::PMECellData () {
}

PMECellData::~PMECellData () {
  clear ();
}

void PMECellData::reinitialize (int nAtoms, int nLen,int interpOrder,
    int gridSizeX, int gridSizeY, int gridSizeZ) {
  clear();
  initialize (nAtoms, nLen, interpOrder,
      gridSizeX, gridSizeY, gridSizeZ);
}

void PMECellData::initialize (int nAtoms, int nLen, int interpOrder,
    int gridSizeX, int gridSizeY, int gridSizeZ) {
  natm = nAtoms;
  nlen = nLen;
  n_interp = interpOrder;

  ngrid_a = gridSizeX;
  ngrid_b = gridSizeY;
  ngrid_c = gridSizeZ;

  int nreal = ngrid_a*ngrid_b*ngrid_c;

  iatemp = (int *)malloc(natm*sizeof(int))-1; // size [natm]
  ibtemp = (int *)malloc(natm*sizeof(int))-1; // size [natm]
  ictemp = (int *)malloc(natm*sizeof(int))-1; // size [natm]


  igrid_a   = cmall_int_mat(1,n_interp,1,nlen,"PMECellData::initialize"); // size [n_interp*nlen]
  igrid_b   = cmall_int_mat(1,n_interp,1,nlen,"PMECellData::initialize"); // size [n_interp*nlen]
  igrid_c   = cmall_int_mat(1,n_interp,1,nlen,"PMECellData::initialize"); // size [n_interp*nlen]
  igrid_now = cmall_int_mat(1,n_interp,1,nlen,"PMECellData::initialize"); // size [n_interp*nlen]

  aj  = (double *)malloc(n_interp*sizeof(double))-1; // size [n_interp]
  rn  = (double *)malloc(n_interp*sizeof(double))-1; // size [n_interp]
  rn1 = (double *)malloc(n_interp*sizeof(double))-1; // size [n_interp]

  frac_a = (double *)malloc(natm*sizeof(double))-1; // size [natm]
  frac_b = (double *)malloc(natm*sizeof(double))-1; // size [natm]
  frac_c = (double *)malloc(natm*sizeof(double))-1; // size [natm]

  qgrid_rspace = (double *)malloc(nreal*sizeof(double))-1; // size [nreal]

  qgrid_now   = cmall_mat(1,n_interp,1,nlen,"PMECellData::initialize"); // size [n_interp*nlen]
  ua          = cmall_mat(1,n_interp,1,nlen,"PMECellData::initialize"); // size [n_interp*nlen]
  ub          = cmall_mat(1,n_interp,1,nlen,"PMECellData::initialize"); // size [n_interp*nlen]
  uc          = cmall_mat(1,n_interp,1,nlen,"PMECellData::initialize"); // size [n_interp*nlen]
  mn_a        = cmall_mat(1,n_interp,1,nlen,"PMECellData::initialize"); // size [n_interp*nlen]
  mn_b        = cmall_mat(1,n_interp,1,nlen,"PMECellData::initialize"); // size [n_interp*nlen]
  mn_c        = cmall_mat(1,n_interp,1,nlen,"PMECellData::initialize"); // size [n_interp*nlen]
  dmn_a       = cmall_mat(1,n_interp,1,nlen,"PMECellData::initialize"); // size [n_interp*nlen]
  dmn_b       = cmall_mat(1,n_interp,1,nlen,"PMECellData::initialize"); // size [n_interp*nlen]
  dmn_c       = cmall_mat(1,n_interp,1,nlen,"PMECellData::initialize"); // size [n_interp*nlen]

  precompute();
}

void PMECellData::clear () {

  if (NULL != iatemp) { free(&iatemp[1]); iatemp = NULL;}
  if (NULL != ibtemp) { free(&ibtemp[1]); ibtemp = NULL;}
  if (NULL != ictemp) { free(&ictemp[1]); ictemp = NULL;}

  cfree_int_mat(igrid_a,1,n_interp,1,nlen); // size [n_interp*nlen]
  cfree_int_mat(igrid_b,1,n_interp,1,nlen); // size [n_interp*nlen]
  cfree_int_mat(igrid_c,1,n_interp,1,nlen); // size [n_interp*nlen]
  cfree_int_mat(igrid_now,1,n_interp,1,nlen); // size [n_interp*nlen]

  igrid_a   = NULL;
  igrid_b   = NULL;
  igrid_c   = NULL;
  igrid_now = NULL;

  if (NULL != aj) { free(&aj[1]); aj = NULL;}
  if (NULL != rn) { free(&rn[1]); rn = NULL;}
  if (NULL != rn1) { free(&rn1[1]); rn1 = NULL;}

  if (NULL != frac_a) { free(&frac_a[1]); frac_a = NULL;}
  if (NULL != frac_b) { free(&frac_b[1]); frac_b = NULL;}
  if (NULL != frac_c) { free(&frac_c[1]); frac_c = NULL;}

  if (NULL != qgrid_rspace) { free(&qgrid_rspace[1]); qgrid_rspace = NULL;}

  cfree_mat(qgrid_now,1,n_interp,1,nlen); // size [n_interp*nlen]
  cfree_mat(ua,1,n_interp,1,nlen); // size [n_interp*nlen]
  cfree_mat(ub,1,n_interp,1,nlen); // size [n_interp*nlen]
  cfree_mat(uc,1,n_interp,1,nlen); // size [n_interp*nlen]
  cfree_mat(mn_a,1,n_interp,1,nlen); // size [n_interp*nlen]
  cfree_mat(mn_b,1,n_interp,1,nlen); // size [n_interp*nlen]
  cfree_mat(mn_c,1,n_interp,1,nlen); // size [n_interp*nlen]
  cfree_mat(dmn_a,1,n_interp,1,nlen); // size [n_interp*nlen]
  cfree_mat(dmn_b,1,n_interp,1,nlen); // size [n_interp*nlen]
  cfree_mat(dmn_c,1,n_interp,1,nlen); // size [n_interp*nlen]

  qgrid_now = NULL;

  ua = NULL;
  ub = NULL;
  uc = NULL;

  mn_a = NULL;
  mn_b = NULL;
  mn_c = NULL;

  dmn_a = NULL;
  dmn_b = NULL;
  dmn_c = NULL;

  natm = 0;
  nlen = 0;
  n_interp = 0;

  ngrid_a = 0;
  ngrid_b = 0;
  ngrid_c = 0;
}

void PMECellData::pup (PUP::er& p) {
  p|natm;
  p|nlen;
  p|n_interp;

  p|ngrid_a;
  p|ngrid_b;
  p|ngrid_c;

  int nreal = ngrid_a*ngrid_b*ngrid_c;

  if (natm > 0) {
    pup1d_int(p, &iatemp, natm); // size [natm]
    pup1d_int(p, &ibtemp, natm); // size [natm]
    pup1d_int(p, &ictemp, natm); // size [natm]
  }

  if (0 < n_interp) {  
    pup1d_dbl(p, &aj, n_interp);  // size [n_interp]
    pup1d_dbl(p, &rn, n_interp);  // size [n_interp]
    pup1d_dbl(p, &rn1, n_interp);  // size [n_interp]
  }

  if (natm > 0) {
    pup1d_dbl(p, &frac_a, natm); // size [natm]
    pup1d_dbl(p, &frac_b, natm); // size [natm]
    pup1d_dbl(p, &frac_c, natm); // size [natm]
  }

  if (0 < nreal) {
    pup1d_dbl(p, &qgrid_rspace, nreal); // size [ngrid_a*ngrid_b*ngrid_c]
  }

  if (0 < n_interp && 0 < nlen) {
    pup2d_int (p, &igrid_a, n_interp, nlen); // size [n_interp*nlen]
    pup2d_int (p, &igrid_b, n_interp, nlen); // size [n_interp*nlen]
    pup2d_int (p, &igrid_c, n_interp, nlen); // size [n_interp*nlen]

    pup2d_int (p, &igrid_now, n_interp, nlen); // size [n_interp*nlen]

    pup2d_dbl (p, &qgrid_now, n_interp, nlen); // size [n_interp*nlen]
    pup2d_dbl (p, &ua, n_interp, nlen); // size [n_interp*nlen]
    pup2d_dbl (p, &ub, n_interp, nlen); // size [n_interp*nlen]
    pup2d_dbl (p, &uc, n_interp, nlen); // size [n_interp*nlen]

    pup2d_dbl (p, &mn_a, n_interp, nlen); // size [n_interp*nlen]
    pup2d_dbl (p, &mn_b, n_interp, nlen); // size [n_interp*nlen]
    pup2d_dbl (p, &mn_c, n_interp, nlen); // size [n_interp*nlen]

    pup2d_dbl (p, &dmn_a, n_interp, nlen); // size [n_interp*nlen]
    pup2d_dbl (p, &dmn_b, n_interp, nlen); // size [n_interp*nlen]
    pup2d_dbl (p, &dmn_c, n_interp, nlen); // size [n_interp*nlen]
  }
}

void PMECellData::precompute () {
  // could be precomputed but costs nothing
  for(int j=1;j<=n_interp;j++){
    aj[j] = (double) (j-1);
    rn[j] = (double) (j);
    if(j > 1){rn1[j] = 1.0/((double)(j-1));}
  }//endfor
  rn1[1] = 0.0;
#ifdef PME_DEBUG
  for(int j=1;j<=n_interp;j++){
    CmiPrintf ("aj[%d]=%e, rn[%d]=%e, rn1[%d]=%e\n",j,aj[j],j,rn[j],j,rn1[j]);
  }
#endif  
}

void PMECellData::computeMN (Vector* _pos, double* _hmati) {

  Vector* pos   = _pos-1; // arrays should start from 1
  double* hmati = _hmati-1; // arrays should start from 1
  //==============================================================================
  //==============================================================================
  // I) Construct some useful constants                                       
  //==============================================================================
  const double    grid_a   = (double) ngrid_a;
  const double    grid_b   = (double) ngrid_b;
  const double    grid_c   = (double) ngrid_c;
  const int       ngrid_bc = ngrid_b*ngrid_c;  

  //==============================================================================
  //==============================================================================
  // II)  Get the scaled imaged particle coordinates
  //      for a general box
  //==============================================================================
  for(int iatm = 1;iatm <= natm;++iatm){
    double atemp = pos[iatm].x*hmati[1]
      + pos[iatm].y*hmati[4]
      + pos[iatm].z*hmati[7];
    double btemp = pos[iatm].x*hmati[2]
      + pos[iatm].y*hmati[5]
      + pos[iatm].z*hmati[8];
    double ctemp = pos[iatm].x*hmati[3]
      + pos[iatm].y*hmati[6]
      + pos[iatm].z*hmati[9];
    // atemp,btemp and ctemp are between zero and 1
    atemp = atemp - NINT((atemp-0.5));
    btemp = btemp - NINT((btemp-0.5));
    ctemp = ctemp - NINT((ctemp-0.5));
    //atemp,btemp and ctemp are between zero and grida,gridb,gridc
    atemp *= grid_a;
    btemp *= grid_b;
    ctemp *= grid_c;
    // compute the fractional remainder (u) from closest grid point
    // so that frac_a, frac_b and frac_c are between zero and 1
    iatemp[iatm] = (int) (atemp);
    ibtemp[iatm] = (int) (btemp);
    ictemp[iatm] = (int) (ctemp);
    frac_a[iatm] = atemp - (double) (iatemp[iatm]);
    frac_b[iatm] = btemp - (double) (ibtemp[iatm]);
    frac_c[iatm] = ctemp - (double) (ictemp[iatm]);
  }//endfor

  //==============================================================================
  //==============================================================================
  // III) Calculate the Cardinal B spline interpolation functions of the
  //      charge weighted density on the real space grid
  //==============================================================================

  // put the atoms on the grid in ``chunks'' of length nlen
  // nlen is cache optimization parameter
  for(int iatm=1;iatm<=natm;iatm+=nlen){
    // clean up for natm % nlen not equal to zero
    int iatm1 = iatm-1;
    int iend = MIN(natm,iatm1+nlen);
    int nnow = iend-iatm1;

    //------------------------------------------------------------------------------
    // A) Using current fraction, find the grid points on which M_n is non-zero 
    //    n_interp-1 is the order of spline. and M_n(a) Mn(b) and Mn(c) are
    //    each non-zero on n_interp points of the grid.
    //------------------------------------------------------------------------------

    for(int j=1;j<=n_interp;j++){
      for(int i=1;i<=nnow;i++){
        const int itmp = i+iatm1;
        ua[j][i]  = frac_a[itmp] + aj[j];
        ub[j][i]  = frac_b[itmp] + aj[j];
        uc[j][i]  = frac_c[itmp] + aj[j];
        int j2    = j-2;
        int ia    = iatemp[itmp] - j2;
        int ib    = ibtemp[itmp] - j2;
        int ic    = ictemp[itmp] - j2;
        ia        = (ia>0 ? ia:ngrid_a+ia);  // range 1-ngrid_a
        ib        = (ib>0 ? ib:ngrid_b+ib);  // range 1-ngrid_b
        ic        = (ic>0 ? ic:ngrid_c+ic);  // range 1-ngrid_c
        ia        = (ia<=ngrid_a ? ia:ia-ngrid_a);  // range 1-ngrid_a
        ib        = (ib<=ngrid_b ? ib:ib-ngrid_b);  // range 1-ngrid_b
        ic        = (ic<=ngrid_c ? ic:ic-ngrid_c);  // range 1-ngrid_c
        // The local mapping of the non-zero points in qgrid.
        // assumed that I have some part of the array starting igrid_str
        // and all the particle you sent fit. In scalar igrid_str=0.
        // The igrid variable are scratch storage used to optimized.
        igrid_a[j][i] = (ia-1)*ngrid_bc;
        igrid_b[j][i] = (ib-1)*ngrid_c;
        igrid_c[j][i] = ic;
      }//endfor
    }//endfor
#ifdef PME_DEBUG    
    for(int j=1;j<=n_interp;j++){
      for(int i=1;i<=nnow;i++){
        CmiPrintf ("igrid[%d][%d]= %d, %d, %d\n", j, i, igrid_a[j][i],igrid_b[j][i],igrid_c[j][i]);
      }
    }
#endif    
    //------------------------------------------------------------------------------
    // B) Initialize M2 and get the Mn's using the recursion relation            
    //    Note: M_n is defined on 0 to n. Since frac between 0 and 1, the       
    //         calculation is performed in an order that takes advantage of it 
    //------------------------------------------------------------------------------
    for(int i=1;i<=nnow;i++){
      mn_a[1][i] = 1.0 - fabs(ua[1][i]-1.0);
      mn_b[1][i] = 1.0 - fabs(ub[1][i]-1.0);
      mn_c[1][i] = 1.0 - fabs(uc[1][i]-1.0);
      mn_a[2][i] = 1.0 - fabs(ua[2][i]-1.0);
      mn_b[2][i] = 1.0 - fabs(ub[2][i]-1.0);
      mn_c[2][i] = 1.0 - fabs(uc[2][i]-1.0);
    }//endfor

    for(int j=3;j<=n_interp;j++){
      for(int i=1;i<=nnow;i++){
        mn_a[j][i]   = 0.0;
        mn_b[j][i]   = 0.0;
        mn_c[j][i]   = 0.0;
      }//endfor
    }//endfor

    // recursive loop
    for(int n=3;n<=n_interp;n++){
      for(int j=n;j>=2;j--){
        const int j1 = j-1;
        // nnow steps of each recursion
        for(int i=1;i<=nnow;i++){
          const double mn_a_tmp = (ua[j][i]*mn_a[j][i]
              + (rn[n]-ua[j][i])*mn_a[j1][i])*rn1[n];
          const double mn_b_tmp = (ub[j][i]*mn_b[j][i]
              + (rn[n]-ub[j][i])*mn_b[j1][i])*rn1[n];
          const double mn_c_tmp = (uc[j][i]*mn_c[j][i]
              + (rn[n]-uc[j][i])*mn_c[j1][i])*rn1[n];
          mn_a[j][i] = mn_a_tmp;
          mn_b[j][i] = mn_b_tmp;
          mn_c[j][i] = mn_c_tmp;
        }//end for: i
      }//end for: j
      for(int i=1;i<=nnow;i++){
        mn_a[1][i] = ua[1][i]*mn_a[1][i]*rn1[n];
        mn_b[1][i] = ub[1][i]*mn_b[1][i]*rn1[n];
        mn_c[1][i] = uc[1][i]*mn_c[1][i]*rn1[n];
      }//endfor

      // derivative of the mn to get the forces
      if(n==(n_interp-1)){
        for(int i=1;i<=nnow;i++){
          dmn_a[1][i] = mn_a[1][i];
          dmn_b[1][i] = mn_b[1][i];
          dmn_c[1][i] = mn_c[1][i];
        }//endfor
        for(int j=2;j<=n_interp;j++){
          int j1 = j-1;
          for(int i=1;i<=nnow;i++){
            dmn_a[j][i] = mn_a[j][i] - mn_a[j1][i];
            dmn_b[j][i] = mn_b[j][i] - mn_b[j1][i];
            dmn_c[j][i] = mn_c[j][i] - mn_c[j1][i];
          }//endfor
        }//endfor
      }//endif : n=n_interp-1
    }//end for: n
  } // end for: iatm
}

void PMECellData::computeRSpaceGrid (Vector* _pos, double* _q, double* _hmati) {

  computeMN(_pos, _hmati);

  // array's should start from 1
  //Vector* pos = _pos-1;
  double* q = _q-1;
  //double* hmati = _hmati-1;

  //-------------------------------------------------------------------------- 
  // Stuff Mn's on the grid:Sum is ordered to remove recursive dependencies 
  // in qgrid. igrid_now is unique for fixed i       
  const int ngrid   = ngrid_a*ngrid_b*ngrid_c;
  for(int i=1;i<=ngrid;i++){
    qgrid_rspace[i]=0.0;
  }//endfor

  for(int iatm=1;iatm<=natm;iatm+=nlen){
    // clean up for natm % nlen not equal to zero
    int iatm1 = iatm-1;
    int iend  = MIN(natm,iatm1+nlen);
    int nnow  = iend-iatm1;
    for(int jc=1;jc<=n_interp;jc++){
      for(int jb=1;jb<=n_interp;jb++){

        // compute index into 1d array [][][] -> []
        // compute the overall weight Mn(x)*Mn(y)*Mn(z)
        for(int ja=1;ja<=n_interp;ja++){
          for(int i=1;i<=nnow;i++){
            int itmp     = i+iatm1;
            igrid_now[ja][i] = igrid_a[ja][i]+igrid_b[jb][i]+igrid_c[jc][i];
            qgrid_now[ja][i] = mn_a[ja][i]*mn_b[jb][i]*mn_c[jc][i]*q[itmp];
          }
        }//endfor

        // loop which optimizations badly
        // due to cache coherency if nnow is on the inside
        // but for each for particle, except for periodic Bc.
        // igrid_now moves linearly through memory with ja.
        for(int i=1;i<=nnow;i++){
          for(int ja=1;ja<=n_interp;ja++){
            qgrid_rspace[igrid_now[ja][i]] += qgrid_now[ja][i];
          }
        }// i,ja
      }
    }//end for: jb,jc
  }//end for: iatm

#ifdef PINY_PME_DEBUG
  for (int i=1; i<ngrid; i++) {
    CmiPrintf ("qgrid_rspace[%d]=%e\n",i,qgrid_rspace[i]);
  }
#endif
}

double* PMECellData::getRSpaceGrid () {
  return &(qgrid_rspace[1]);
}

void PMECellData::calcForces (Vector* _force, double* _q, double* _hmati) {
  const double    grid_a   = (double) ngrid_a;
  const double    grid_b   = (double) ngrid_b;
  const double    grid_c   = (double) ngrid_c;

  // arrays should start from 1
  Vector* force = _force-1;
  double* q     = _q-1;
  double* hmati = _hmati-1;

  for(int iatm=1;iatm<=natm;iatm+=nlen){
    // clean up for natm % nlen not equal to zero
    int iatm1 = iatm-1;
    int iend  = MIN(natm,iatm1+nlen);
    int nnow  = iend-iatm1;

    // loop over for cache opt
    for(int jc=1;jc<=n_interp;jc++){
      for(int jb=1;jb<=n_interp;jb++){

        // pull out a portion of the qgrid that you need
        for(int ja=1;ja<=n_interp;ja++){
          for(int i=1;i<=nnow;i++){
            igrid_now[ja][i] = igrid_a[ja][i]+igrid_b[jb][i]+igrid_c[jc][i];
          }
        }//endfor
        for(int ja=1;ja<=n_interp;ja++){
          for(int i=1;i<=nnow;i++){
            qgrid_now[ja][i] = qgrid_rspace[igrid_now[ja][i]];
          }
        }//endfor

        // use qgrid, the mn's and the dmn's to get the forces on the particles.
        for(int ja=1;ja<=n_interp;ja++){
          for(int i=1;i<=nnow;i++){
            const int itmp        = i+iatm1;
            const double atemp    = dmn_a[ja][i]* mn_b[jb][i]* mn_c[jc][i]*grid_a;
            const double btemp    = mn_a[ja][i]*dmn_b[jb][i]* mn_c[jc][i]*grid_b;
            const double ctemp    = mn_a[ja][i]* mn_b[jb][i]*dmn_c[jc][i]*grid_c;
            const double qgrid_dx = atemp*hmati[1]+btemp*hmati[2]+ctemp*hmati[3];
            const double qgrid_dy = atemp*hmati[4]+btemp*hmati[5]+ctemp*hmati[6];
            const double qgrid_dz = atemp*hmati[7]+btemp*hmati[8]+ctemp*hmati[9];
            const double scale    = q[itmp]*qgrid_now[ja][i];
            force[itmp].x        -= (qgrid_dx*scale);
            force[itmp].y        -= (qgrid_dy*scale);
            force[itmp].z        -= (qgrid_dz*scale);
          }
        }//endfor : atoms, ja
      }
    }//end for: jb,jc
  }//end for: iatm
}
