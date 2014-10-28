
//==========================================================================
// 1 week: 
//   Glenn will make this work in scalar on a small system which
//   will be read in. particle positions, charges, grid size,
//   and the order of the interpolation. Write some neater notes

// parallelization scheme:
//        identify type of chares
//           do you want to parallelize the FFT by planes.
//           gather either the grid OR compute the grid of each
//           patch and communicate to the plane chares.
//           Interpolation spreads the particle to N^3 points
//           One dimensional example
//                        |    o          | goes to 
//                        |    ooooooo    | 
//           Spherically cutoff FFT helps ALOT in scalar and parallel
//
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

double ewald3d_recip_pme(
    int natm,             // number of atoms
    const VECTOR& pos,    // atom positions
    VECTOR& force,        // atom forces
    double *q,            // atom charges
    double *pvten,        // pressure tensor
    int iperd,            // periodicity
    double *hmati,        // inverse cell matrix
    double vol,           // cell volume
    double alp_ewd,       // ewald convergence parameter
    int n_interp,         // pme interpolation level
    int ngrid_a,          // grid pts along a in real space
    int ngrid_b,          // grid pts along b in real space
    int ngrid_c,          // grid pts along c in real space
    int igrid_strt,       // parallel offset index
    int nktot,            // pts in spherically cutoff g-space
    double ecut_pme,      // spherical cutoff in g-space
    int *ka,              // ka of g-space vector
    int *kb,              // kb of g-space vector
    int *kc,              // kc of g-space vector
    PME_CONS_SCR *pme_cons_scr
)

//==========================================================================*/
{//begin routine 
  //==========================================================================
  //  Strip out the scratch memory and constants from pme_cons_scr

  double pvten_tmp[10];

  int *iatemp            = pme_cons_scr->iatemp;     // size natm
  int *ibtemp            = pme_cons_scr->ibtemp;     // size natm
  int *ictemp            = pme_cons_scr->ictemp;     // size natm
  double *aj             = pme_cons_scr->aj;         // size n_interp
  double *rn             = pme_cons_scr->rn;         // size n_interp
  double *rn1            = pme_cons_scr->rn1;        // size n_interp
  double *frac_a         = pme_cons_scr->frac_a;     // size natm
  double *frac_b         = pme_cons_scr->frac_b;     // size natm
  double *frac_c         = pme_cons_scr->frac_c;     // size natm
  double *qgrid          = pme_cons_scr->qgrid;      // size fft r-space
  double *qgrid_tmp_re   = pme_cons_scr->qgrid_tmp_re;// size fft g-space
  double *qgrid_tmp_im   = pme_cons_scr->qgrid_tmp_im;// size fft g-space

  double *bweight_tot    = pme_cons_scr->bweight_tot;// size fft g-space
  // precomputed const

  int nlen               = pme_cons_scr->nlen;       // scratch size
  int **igrid_a          = pme_cons_scr->igrid_a;    // n_interp x nlen
  int **igrid_b          = pme_cons_scr->igrid_b;
  int **igrid_c          = pme_cons_scr->igrid_c;
  int **igrid_now        = pme_cons_scr->igrid_now;
  double **qgrid_now     = pme_cons_scr->qgrid_now;
  double **ua            = pme_cons_scr->ua;
  double **ub            = pme_cons_scr->ub;
  double **uc            = pme_cons_scr->uc;
  double **mn_a          = pme_cons_scr->mn_a;
  double **mn_b          = pme_cons_scr->mn_b;
  double **mn_c          = pme_cons_scr->mn_c;
  double **dmn_a         = pme_cons_scr->dmn_a;
  double **dmn_b         = pme_cons_scr->dmn_b;
  double **dmn_c         = pme_cons_scr->dmn_c;

  //==========================================================================
  // I) Construct some useful constants                                       

  const int nrecip_grid   = 2*ngrid_a*ngrid_b*ngrid_c;
  const double    grid_a  = (double) ngrid_a;
  const double    grid_b  = (double) ngrid_b;
  const double    grid_c  = (double) ngrid_c;
  const int    ngrid_bc   = ngrid_b*ngrid_c;
  const int    ngrid_ab   = ngrid_a*ngrid_b;

  // could be precomputed
  for(int j=1;j<=n_interp;j++){
    aj[j] = (double) (j-1);
    rn[j] = (double) (j);
    if(j > 1){rn1[j] = 1.0/((double)(j-1));}
  }//endfor
  rn1[1] = 0.0;

  //==========================================================================
  // II) Get the scaled imaged particle coordinates
  //     for a general box

  for(int iatm = 0;iatm < natm;++iatm){
    double atemp = pos[iatm].x*hmati[1]
      + pos[iatm].y*hmati[4]
      + pos[iatm].z*hmati[7];
    double btemp = pos[iatm].x*hmati[2]
      + pos[iatm].y*hmati[5]
      + pos[iatm].z*hmati[8];
    double ctemp = pos[iatm].x*hmati[3]
      + pos[iatm].y*hmati[6]
      + pos[iatm].z*hmati[9];
#ifdef ORTHO_BOX_EXAMPLE
    double atemp = pos[iatm].x/Lx
      double btemp = pos[iatm].y/Ly
      double ctemp = pos[iatm].z/Lz;
#endif
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
    iatemp[i] = (int) (pos[i].x);
    ibtemp[i] = (int) (pos[i].y);
    ictemp[i] = (int) (pos[i].z);
    frac_a[i] = pos[i].x - (double) (iatemp[i]);
    frac_b[i] = pos[i].y - (double) (ibtemp[i]);
    frac_c[i] = pos[i].z - (double) (ictemp[i]);
  }//endfor

  //==========================================================================
  // III) Calculate the Cardinal B spline interpolation functions of the      
  //     charge weighted density on the real space grid                       

  // In scalar this would be of size ngrid_a*ngrid_b*ngrid_c
  // grid[i][j][j] expressed as a 1d array.

  for(int i=0;i<nrecip_grid;i++){
    qgrid[i]=0.0;
  }//endfor

  // put the atoms on the grid in ``chunks'' of length nlen
  // nlen is cache optimization parameter

  for(int iatm=0;iatm<natm;iatm+=nlen){
    // clean up for natm % nlen not equal to zero
    int iend = MIN(natm,iatm+nlen);
    int nnow = iend-iatm;
    //-------------------------------------------------------------------------- 
    // A) Using current fraction, find the grid points on which M_n is non-zero 
    //    n_interp-1 is the order of spline. and M_n(a) Mn(b) and Mn(c) are
    //    each non-zero on n_interp points of the grid.

    for(int j=1;j<=n_interp;j++){
      for(int i=0;i<nnow;i++){
        const int itmp = i+iatm;
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
        igrid_a[j][i] = 2*ia - 2 - igrid_str;
        igrid_b[j][i] = 2*(ib - 1)*ngrid_a;
        igrid_c[j][i] = 2*(ic - 1)*ngrid_ab;
      }//endfor
    }//endfor

    //-------------------------------------------------------------------------- 
    // B) Initialize M2 and get the Mn's using the recursion relation            
    //    Note: M_n is defined on 0 to n. Since frac between 0 and 1, the       
    //          calculation is performed in an order that takes advantage of it 

    for(int i=0;i<nnow;i++){
      mn_a[1][i] = 1.0 - fabs(ua[1][i]-1.0);
      mn_b[1][i] = 1.0 - fabs(ub[1][i]-1.0);
      mn_c[1][i] = 1.0 - fabs(uc[1][i]-1.0);
      mn_a[2][i] = 1.0 - fabs(ua[2][i]-1.0);
      mn_b[2][i] = 1.0 - fabs(ub[2][i]-1.0);
      mn_c[2][i] = 1.0 - fabs(uc[2][i]-1.0);
    }//endfor

    for(int j=3;j<=n_interp;j++){
      for(int i=0;i<nnow;i++){
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
        for(int i=0;i<nnow;i++){
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
      for(int i=0;i<nnow;i++){
        mn_a[1][i] = ua[1][i]*mn_a[1][i]*rn1[n];
        mn_b[1][i] = ub[1][i]*mn_b[1][i]*rn1[n];
        mn_c[1][i] = uc[1][i]*mn_c[1][i]*rn1[n];
      }//endfor 

    }//end for: n

    //-------------------------------------------------------------------------- 
    // C) Stuff Mn's on the grid:Sum is ordered to remove recursive dependencies 
    //                           in qgrid. igrid_now is unique for fixed i       

    for(int jc=1;jc<=n_interp;jc++){
      for(int jb=1;jb<=n_interp;jb++){

        // compute index into 1d array [][][] -> []
        // compute the overall weight Mn(x)*Mn(y)*Mn(z)
        for(int ja=1;ja<=n_interp;ja++){
          for(int i=0;i<nnow;i++){
            int itmp     = i+iatm;
            igrid_now[ja][i] = igrid_a[ja][i]+igrid_b[jb][i]+igrid_c[jc][i];
            qgrid_now[ja][i] = mn_a[ja][i]*mn_b[jb][i]*mn_c[jc][i]*q[itmp];
          }}//endfor

        // loop which optimizations badly
        // due to cache coherency if nnow is on the inside
        // but for each for particle, except for periodic Bc.
        // igrid_now moves linearly through memory with ja. 
        for(int i=0;i<nnow;i++){
          for(int ja=1;ja<=n_interp;ja++){
            qgrid[igrid_now[ja][i]] += qgrid_now[ja][i];
          }}// i,ja

      }}//end for: jb,jc

  }//end for: iatm

  //==========================================================================
  // IV) Fourier Transform qgrid  (3d parallel FFT)

  fft3d_to_gspace(qgrid,qgrid_tmp_re,qgrid_tmp_im); 

  //==========================================================================
  // VI) Compute the potential energy on the spherically cutoff grid.         

  const double falp2 = 4.0*alp_ewd*alp_ewd;
  const double tpi   = 2.0*M_PI;
  const double pivol = vol/(4.0*M_PI);
  const double rvol  = 1.0/vol;  

  for(int i=1;i<=9;i++){pvten_tmp[i]=0.0;}
  double pot_pme = 0.0;

  // in parallel my part of reciprocal space
  // in scalar its all of it.
  for(int i=0;i < nktot; ++i) {
    const double aka  = tpi*( (double) kastr[i] );
    const double akb  = tpi*( (double) kbstr[i] );
    const double akc  = tpi*( (double) kcstr[i] );
    const double xk   = (aka*hmati[1] + akb*hmati[2] + akc*hmati[3]);
    const double yk   = (aka*hmati[4] + akb*hmati[5] + akc*hmati[6]);
    const double zk   = (aka*hmati[7] + akb*hmati[8] + akc*hmati[9]);
    const double g2   = xk*xk + yk*yk + zk*zk;
    if(0.5*g2 < ecut_pme && g2 != 0.0){
      const double preg = exp((-g2/falp2))/(g2*pivol);
      const double smag = (qgrid_tmp_re[i]*qgrid_tmp_re[i]
          +qgrid_tmp_im[i]*qgrid_tmp_im[i])*bweight_tot[i];
      const double prep = -2.0*preg*smag*((g2/falp2)+1.0)/g2;
      pot_pme      += (smag*preg);
      pvten_tmp[1] += prep*xk*xk;
      pvten_tmp[5] += prep*yk*yk;
      pvten_tmp[9] += prep*zk*zk;
      pvten_tmp[2] += prep*xk*yk;
      pvten_tmp[3] += prep*xk*zk;
      pvten_tmp[6] += prep*yk*zk;
      qgrid_tmp_re[i] *= (preg*bweight_tot[i]);
      qgrid_tmp_im[i] *= (preg*bweight_tot[i]);
    }else{
      qgrid_tmp_re[i] = 0.0;
      qgrid_tmp_im[i] = 0.0;
    }//endif
  }//endfor

  pvten_tmp[4]  = pvten_tmp[2];
  pvten_tmp[7]  = pvten_tmp[3];
  pvten_tmp[8]  = pvten_tmp[6];
  pvten_tmp[1] += pot_pme;
  pvten_tmp[5] += pot_pme;
  pvten_tmp[9] += pot_pme;

  for(int i=1;i<=9;i++){
    pvten[i] += pvten_tmp[i];
  }//endfor

  //==========================================================================
  // VII) Fourier Transform qgrid back to get the forces on the particles

  fft3d_to_gspace(qgrid,qgrid_tmp_re,qgrid_tmp_im); 

  //==========================================================================
  // VIII) Calculate the force                                                

  // because things weren't done for all the particles
  // there is some recompute. On many procs with natm<nlen
  // vaoid the recompute. If you are on 1 proc with tons
  // of memory you can avoid the recompute.

  for(int iatm=0;iatm<natm;iatm+=nlen){
    int iend = MIN(natm,iatm+nlen);
    int nnow = iend-iatm;
    //-------------------------------------------------------------------------- 
    // A) Using current fraction, find the grid points on which M_n is non-zero 

    // exactly above
    for(int j=1;j<=n_interp;j++){
      for(int i=0;i<nnow;i++){
        const int itmp = i+iatm;
        ua[j][i]    = frac_a[itmp] + aj[j];
        ub[j][i]    = frac_b[itmp] + aj[j];
        uc[j][i]    = frac_c[itmp] + aj[j];
        int j2      = j-2;
        int ia      = iatemp[itmp] - j2;
        int ib      = ibtemp[itmp] - j2;
        int ic      = ictemp[itmp] - j2;
        ia        = (ia>0 ? ia:ngrid_a+ia);  // range 1-ngrid_a
        ib        = (ib>0 ? ib:ngrid_b+ib);  // range 1-ngrid_b
        ic        = (ic>0 ? ic:ngrid_c+ic);  // range 1-ngrid_c
        ia        = (ia<=ngrid_a ? ia:ia-ngrid_a);  // range 1-ngrid_a
        ib        = (ib<=ngrid_b ? ib:ib-ngrid_b);  // range 1-ngrid_b
        ic        = (ic<=ngrid_c ? ic:ic-ngrid_c);  // range 1-ngrid_c
        igrid_a[j][i] = 2*ia - 2 - igrid_str;
        igrid_b[j][i] = 2*(ib - 1)*ngrid_a;
        igrid_c[j][i] = 2*(ic - 1)*ngrid_ab;
      }//endfor
    }//endfor

    //-------------------------------------------------------------------------- 
    // B) Initialize M2 and get the Mn's using the recursion relation            
    //    Note: M_n is defined on 0 to n. Since frac between 0 and 1, the       
    //          calculation is performed in an order that takes advantage of it 

    for(int i=0;i<nnow;i++){
      mn_a[1][i] = 1.0 - fabs(ua[1][i]-1.0);
      mn_b[1][i] = 1.0 - fabs(ub[1][i]-1.0);
      mn_c[1][i] = 1.0 - fabs(uc[1][i]-1.0);
      mn_a[2][i] = 1.0 - fabs(ua[2][i]-1.0);
      mn_b[2][i] = 1.0 - fabs(ub[2][i]-1.0);
      mn_c[2][i] = 1.0 - fabs(uc[2][i]-1.0);
    }//endfor
    for(int j=3;j<=n_interp;j++){
      for(int i=0;i<nnow;i++){
        mn_a[j][i]   = 0.0;
        mn_b[j][i]   = 0.0;
        mn_c[j][i]   = 0.0;
      }//endfor
    }//endfor

    for(int n=3;n<=n_interp;n++){

      for(int j=n;j>=2;j--){
        int j1 = j-1;
        for(int i=0;i<nnow;i++){
          const double mn_a_tmp = (ua[j][i]*mn_a[j][i]
              +(rn[n]-ua[j][i])*mn_a[j1][i])*rn1[n];
          const double mn_b_tmp = (ub[j][i]*mn_b[j][i]
              +(rn[n]-ub[j][i])*mn_b[j1][i])*rn1[n];
          const double mn_c_tmp = (uc[j][i]*mn_c[j][i]
              +(rn[n]-uc[j][i])*mn_c[j1][i])*rn1[n];
          mn_a[j][i] = mn_a_tmp;
          mn_b[j][i] = mn_b_tmp;
          mn_c[j][i] = mn_c_tmp;
        }//end for: i
      }//end for: j
      for(int i=0;i<nnow;i++){
        mn_a[1][i] = ua[1][i]*mn_a[1][i]*rn1[n];
        mn_b[1][i] = ub[1][i]*mn_b[1][i]*rn1[n];
        mn_c[1][i] = uc[1][i]*mn_c[1][i]*rn1[n];
      }//endfor 

      // derivative of the mn to get the forces
      if(n==(n_interp-1)){
        for(int i=0;i<nnow;i++){
          dmn_a[1][i] = mn_a[1][i];
          dmn_b[1][i] = mn_b[1][i];
          dmn_c[1][i] = mn_c[1][i];
        }//endfor
        for(int j=2;j<=n_interp;j++){    
          int j1 = j-1;
          for(int i=0;i<nnow;i++){
            dmn_a[j][i] = mn_a[j][i] - mn_a[j1][i];
            dmn_b[j][i] = mn_b[j][i] - mn_b[j1][i];
            dmn_c[j][i] = mn_c[j][i] - mn_c[j1][i];
          }//endfor
        }//endfor
      }//endif : n=n_interp-1

    }//end for: n

    //-------------------------------------------------------------------------- 
    // C) Calculate the force                                                    

    // loop over for cache opt
    for(int jc=1;jc<=n_interp;jc++){
      for(int jb=1;jb<=n_interp;jb++){

        // pull out a portion of the qgrid that you need
        for(int ja=1;ja<=n_interp;ja++){
          for(int i=0;i<nnow;i++){
            igrid_now[ja][i] = igrid_a[ja][i]+igrid_b[jb][i]+igrid_c[jc][i];
          }}//endfor
        for(int ja=1;ja<=n_interp;ja++){
          for(int i=0;i<nnow;i++){
            qgrid_now[ja][i] = qgrid[igrid_now[ja][i]];
          }}//endfor

        // use qgrid, the mn's and the dmn's to get the forces on the particles.
        for(int ja=1;ja<=n_interp;ja++){
          for(int i=0;i<nnow;i++){
            const int itmp        = i+iatm;
            const double atemp    = dmn_a[ja][i]* mn_b[jb][i]* mn_c[jc][i]*grid_a;
            const double btemp    =  mn_a[ja][i]*dmn_b[jb][i]* mn_c[jc][i]*grid_b;
            const double ctemp    =  mn_a[ja][i]* mn_b[jb][i]*dmn_c[jc][i]*grid_c;
            const double qgrid_dx = atemp*hmati[1]+btemp*hmati[2]+ctemp*hmati[3];
            const double qgrid_dy = atemp*hmati[4]+btemp*hmati[5]+ctemp*hmati[6];
            const double qgrid_dz = atemp*hmati[7]+btemp*hmati[8]+ctemp*hmati[9];
            const double scale    = q[itmp]*qgrid_now[ja][i];
            force[itmp].x        -= (qgrid_dx*scale);
            force[itmp].y        -= (qgrid_dy*scale);
            force[itmp].z        -= (qgrid_dz*scale); 
          }}//endfor : atoms, ja
      }}//end for: jb,jc

  }//end for: iatm

  return pot_pme;
  //----------------------------------------------------------------------- 
}//end routine
//=======================================================================
