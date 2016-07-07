// hardcoded nfft value!!!
#include <cstdlib>
#include "print_util.h"
// 1D index to 3D G index 
void PdataIndexG(int thisIndex, int (&gidx)[3]){

  int nfft[3] = {12,12,12};
  int ijk = 0;
  bool match = false;
  for(int i=0; i<nfft[0]; i++){
    for(int j=0; j<nfft[1]; j++){
      for(int k=0; k<nfft[2]; k++){
        if(ijk == thisIndex){
          gidx[0] = i;
          gidx[1] = j;
          gidx[2] = k;
          match = true;
          break;
        }
        ijk += 1; //increment 1D index
      }//end k loop
      if(match){break;}
    } //end j loop
    if(match){break;}
  }//end i loop

  for (int i=0; i<3; i++){
    if (gidx[i] >= nfft[i]/2){
      gidx[i] -= nfft[i];
    }
  }

}



/*
// writing file to compare our values to Berkeley GW
void compare_BGW(CMATRIX* epsInv, GSPACE* g, int iq, char* fromFile){
  FILE *fp = fopen(fromFile,"r");
  int numdata;
  // read the first line, which is the number of data written in the file
  fscanf(fp,"%d",&numdata);
  // read data
  complex bgw[numdata];
  int ga1[numdata];
  int gb1[numdata];
  int gc1[numdata];
  int ga2[numdata];
  int gb2[numdata];
  int gc2[numdata];
  for (int i=0; i<numdata; i++){
    fscanf(fp,"%d %d %d %d %d %d %lg %lg", &ga1[i], &gb1[i], &gc1[i], 
                &ga2[i], &gb2[i], &gc2[i], &bgw[i].re, &bgw[i].im);
  }
  fclose(fp);

  char writeFile[200];
  sprintf(writeFile, "EpsInv_q%d_comp.dat", iq);
  fp = fopen(writeFile,"w");
  // let's compare data
  for (int r=0; r < epsInv->nrow; r++){
    for (int c=0; c < epsInv->ncol; c++){
      int igr = g->ig[r];
      int jgr = g->jg[r];
      int kgr = g->kg[r];
      int igc = g->ig[c];
      int jgc = g->jg[c];
      int kgc = g->kg[c];
      // let's find
      for (int n=0; n < numdata; n++){
        if(igr == ga1[n] && jgr == gb1[n] && kgr == gc1[n] &&
           igc == ga2[n] && jgc == gb2[n] && kgc == gc2[n]){
           fprintf(fp," (%d, %d, %d)   (%d, %d, %d)   Yale:   %lg   %lg    BGW:   %lg   %lg\n",
                   igr, jgr, kgr, igc, jgc, kgc, epsInv->get(r,c).re, epsInv->get(r,c).im, bgw[n].re, bgw[n].im);
	   break;
        }//endif
      }//end searching for the matching data
    }//end for the column
  }//end for the row
}



// writing file to compare our values to Berkeley GW
void compareP_BGW(CMATRIX* P, int nfft[3], int iq, char* fromFile){
  FILE *fp = fopen(fromFile,"r");
  int numdata;
  // read the first line, which is the number of data written in the file
  fscanf(fp,"%d",&numdata);
  // read data
  complex bgw[numdata];
  int ga1[numdata];
  int gb1[numdata];
  int gc1[numdata];
  int ga2[numdata];
  int gb2[numdata];
  int gc2[numdata];
  for (int i=0; i<numdata; i++){
    double dummy1, dummy2;
    fscanf(fp,"%d %d %d %lg %d %d %d %lg %lg %lg", &ga1[i], &gb1[i], &gc1[i], &dummy1,
	     &ga2[i], &gb2[i], &gc2[i], &dummy2, &bgw[i].re, &bgw[i].im);
  }
  fclose(fp);

printf("reading done. Comparison\n");

  char writeFile[200];
  sprintf(writeFile, "P_q%d_comp.dat", iq);
  fp = fopen(writeFile,"w");
  // let's compare data
  int r=0, c=0; // row and column conter
  for (int ir=0; ir<nfft[0]; ir++){
    for (int jr=0; jr<nfft[1]; jr++){
      for (int kr=0; kr<nfft[2]; kr++){
	  c = 0;
        for (int ic=0; ic<nfft[0]; ic++){
	  for (int jc=0; jc<nfft[1]; jc++){
	    for (int kc=0; kc<nfft[2]; kc++){
	      // let's find
              int g1=ir, g2=jr, g3=kr, g4=ic, g5=jc, g6=kc;
              if (g1 >= nfft[0]/2 ){ g1 -= nfft[0];}
              if (g2 >= nfft[1]/2 ){ g2 -= nfft[1];}
              if (g3 >= nfft[2]/2 ){ g3 -= nfft[2];}
              if (g4 >= nfft[0]/2 ){ g4 -= nfft[0];}
              if (g5 >= nfft[1]/2 ){ g5 -= nfft[1];}
              if (g6 >= nfft[2]/2 ){ g6 -= nfft[2];}
 
	      for (int n=0; n < numdata; n++){
	        if(g1 == ga1[n] && g2 == gb1[n] && g3 == gc1[n] &&
		   g4 == ga2[n] && g5 == gb2[n] && g6 == gc2[n]){
		   fprintf(fp," (%d, %d, %d)   (%d, %d, %d)   Yale:   %lg   %lg    BGW:   %lg   %lg\n",
                   g1, g2, g3, g4, g5, g6, P->get(r,c).re, P->get(r,c).im, bgw[n].re, bgw[n].im);
		   //break;
		 }//endif
	       }//end for
	       c += 1; // increase column counter
	     }//end kc
	   }//end jc
	 }//end ic
	 r += 1; // increase row counter
       }//end kr
     }//end jr
   }//end ir

}//end function
*/
