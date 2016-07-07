#ifndef SYSINFO_H
#define SYSINFO_H

#include <cstdio>
#include <cstdlib>

class SYSINFO{
 public:

  int nspin;       // number of spin
  int nkpt;        // number of k points in each spin channel
  int nstate;      // total number of states
  int nocc;        // total number of valence (filled) state
  int nunocc;      // total number of conduction (empty) state
  
  double a1[3], a2[3], a3[3];  // lattice vectors (atomic unit)
  double b1[3], b2[3], b3[3];  // reciprocal lattice vectors
  double shift[3];             // shift vector for Vcoulb(q+G) when q=0 and G=0
  
  double vol;      // volume of the simulation box (atomic unit)
  double alat;     // lattice constant (atomic unit)
  
  double *kwt;
  double **kcart;  // cartesian coordinates
  double **qcart;  // cartesian coordinates
  double **kvec;   // crystal coordinates (reciprocal basis)
  double **qvec;   // crystal coordinates (reciprocal basis)
  
  int *npwk;       // number of planewaves at each k point
  int nfftDen[3];  // number of dense fft grid (for density)

/*
  // constructor
  SYSINFO(){
    read_sysinfo();
    cartesian_to_crystal();
    calc_qvec();
    calc_vol();
#ifdef DEBUG
    system_class_out();
#endif
    }
  }
*/
  //===========================
  // member functions
  //===========================

  void construct(){
    read_sysinfo();
    cartesian_to_crystal();
    calc_qvec();
    calc_vol();
#ifdef DEBUG
    system_class_out();
#endif
  }
  // read it from file
  void read_sysinfo();
  void cartesian_to_crystal();
  void calc_qvec();
  void calc_vol();
  
  // say what you have
  void system_class_out(){
    FILE *fp; fp = fopen("SYSINFO_class_out","w");
    fprintf(fp,"number of spin: %d\n", nspin);
    fprintf(fp,"number of k point: %d\n", nkpt);
    fprintf(fp,"number of total states: %d\n",nstate);
    fprintf(fp,"lattice vectors:\n");
    fprintf(fp,"a1:  %lg  %lg  %lg \n",a1[0],a1[1],a1[2]);
    fprintf(fp,"a2:  %lg  %lg  %lg \n",a2[0],a2[1],a2[2]);
    fprintf(fp,"a3:  %lg  %lg  %lg \n",a3[0],a3[1],a3[2]);
    fprintf(fp,"reciprocal lattice vectors:\n");
    fprintf(fp,"b1:  %lg  %lg  %lg \n",b1[0],b1[1],b1[2]);
    fprintf(fp,"b2:  %lg  %lg  %lg \n",b2[0],b2[1],b2[2]);
    fprintf(fp,"b3:  %lg  %lg  %lg \n",b3[0],b3[1],b3[2]);
    fprintf(fp,"lattice constant:  %lg\n",alat);
    fprintf(fp,"volume of the cell: %lg\n",vol);
    fprintf(fp,"k points in inverse bohr unit (cartesian coordinates) and weight\n");
    for (int i=0; i<nkpt; i++){
      fprintf(fp,"%lg  %lg  %lg  %lg\n",kcart[i][0],kcart[i][1],kcart[i][2],kwt[i]);
    }
    fprintf(fp,"k points in crystal coordinates\n");
    for (int i=0; i<nkpt; i++){
      fprintf(fp,"%lg  %lg  %lg\n",kvec[i][0],kvec[i][1],kvec[i][2]);
    }
    fprintf(fp,"q points in crystal coordinates\n");
    for (int i=0; i<nkpt; i++){
      fprintf(fp,"%lg  %lg  %lg\n",qvec[i][0],qvec[i][1],qvec[i][2]);
    }

    fprintf(fp,"number of planewaves at each k points:\n");
    for (int i=0; i<nkpt; i++){
      fprintf(fp,"kpoint %d:  %d \n",i, npwk[i]);
    }
    fclose(fp);
  }//end class_out function
  
};


#endif
