#include "debug_gsize.h"

void calc_gsize(STATES* psi, SYSINFO sys, int this_kIndex){

  if (this_kIndex==0 || this_kIndex==1 || this_kIndex == 3){
    char fileName[200];
    sprintf(fileName,"KineticEnergy_kIndex_%d.dat",this_kIndex);
    FILE *fp = fopen(fileName,"w");
    printf("Calculate Kinetic energy at %d\n", this_kIndex);
    printf("k vector is:  %lg  %lg  %lg\n",sys.kvec[this_kIndex][0], sys.kvec[this_kIndex][1], sys.kvec[this_kIndex][2]);
  
    // number of data
    int numData = psi->ndata;

    for (int i=0; i<numData; i++){
      double ga, gb, gc; // k+G vector
      ga = psi->gvec.ig[i] + sys.kvec[this_kIndex][0];
      gb = psi->gvec.jg[i] + sys.kvec[this_kIndex][1];
      gc = psi->gvec.kg[i] + sys.kvec[this_kIndex][2];

      double gx, gy, gz; // k+G in cartesian coordinates
      gx = ga*sys.b1[0] + gb*sys.b2[0] + gc*sys.b3[0];
      gy = ga*sys.b1[1] + gb*sys.b2[1] + gc*sys.b3[1];
      gz = ga*sys.b1[2] + gb*sys.b2[2] + gc*sys.b3[2];

      gx *= 2*M_PI/sys.alat;
      gy *= 2*M_PI/sys.alat;
      gz *= 2*M_PI/sys.alat;

      double Ekin = 0;
      Ekin = 0.5 * ( gx*gx + gy*gy + gz*gz );
      fprintf(fp, "%d %d %d   %lg \n", psi->gvec.ig[i], psi->gvec.jg[i], psi->gvec.kg[i], Ekin);
    }//end i loop
    fclose(fp);
  }//end if

}
