/* Set size of the FFT */

#include "get_fftsize.h"
#include <cstdlib>
#include <cmath>
#include <iostream>

void wfn_fftsize(USRINPUT usrin, SYSINFO sys, int (&nfft)[3]){

  // input EcutFFT value is in Rydberg unit, so divide by 2 
  double Ecut = usrin.EcutFFT*0.5;

  // set temporary value for nfft
  int maxga, maxgb, maxgc;
  maxga = sys.nfftDen[0];
  maxgb = sys.nfftDen[1];
  maxgc = sys.nfftDen[2];

  double xk, yk, zk, tryme;
  const double factor = 2*M_PI/sys.alat;

  for (int ga=1; ga<maxga; ga++){
    xk = double(ga)*sys.b1[0]*factor;
    yk = double(ga)*sys.b1[1]*factor;
    zk = double(ga)*sys.b1[2]*factor;
    tryme = 0.5 * ( xk*xk + yk*yk + zk*zk );
    if ( Ecut <= tryme ){
      nfft[0] = ga;
      break;
    }
  }

  for (int gb=1; gb<maxgb; gb++){
    xk = double(gb)*sys.b2[0]*factor;
    yk = double(gb)*sys.b2[1]*factor;
    zk = double(gb)*sys.b2[2]*factor;
    tryme = 0.5 * ( xk*xk + yk*yk + zk*zk );
    if ( Ecut <= tryme ){
      nfft[1] = gb;
      break;
    }
  }

  for (int gc=1; gc<maxgc; gc++){
    xk = double(gc)*sys.b3[0]*factor;
    yk = double(gc)*sys.b3[1]*factor;
    zk = double(gc)*sys.b3[2]*factor;
    tryme = 0.5 * ( xk*xk + yk*yk + zk*zk );
    if ( Ecut <= tryme ){
      nfft[2] = gc;
      break;
    }
  }

  for (int i=0; i<3; i++){
   nfft[i] = nfft[i]*2 + 1;
  }


  int nrad_in = 200;
  int nrad;
  int k;
  int krad[201];
  set_radix(nrad_in, &nrad, krad);
  for (int j = 0; j < 3; j++) {
    for (k = 1; k <= nrad; k++) {
      if (krad[k] > nfft[j]) {
        break;
      }
    }
    nfft[j] = krad[k];
  }


#ifdef DEBUG
  printf("FFT grid for wavefunction: %d %d %d\n",nfft[0], nfft[1], nfft[2]);
#endif

}//end function




// set the optimal number of fft grids
void set_radix(int nrad_in,int *nrad_ret, int *krad)

  /*==========================================================================*/
{/*begin routine */
  /*-------------------------------------------------------------------------*/

  int nrad=179;
  (*nrad_ret)  = nrad;

  if(nrad_in<nrad){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    printf("Internal Error in hardcoded radix size array.\n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    exit(1);
  }/*endif*/

  krad[1]   = 4;
  krad[2]   = 6;
  krad[3]   = 8;
  krad[4]   = 10;
  krad[5]   = 12;
  krad[6]   = 14;
  krad[7]   = 16;
  krad[8]   = 18;
  krad[9]   = 20;
  krad[10]  = 22;
  krad[11]  = 24;
  krad[12]  = 28;
  krad[13]  = 30;
  krad[14]  = 32;
  krad[15]  = 36;
  krad[16]  = 40;
  krad[17]  = 42;
  krad[18]  = 44;
  krad[19]  = 48;
  krad[20]  = 56;
  krad[21]  = 60;
  krad[22]  = 64;
  krad[23]  = 66;
  krad[24]  = 70;
  krad[25]  = 72;
  krad[26]  = 80;
  krad[27]  = 84;
  krad[28]  = 88;
  krad[29]  = 90;
  krad[30]  = 96;
  krad[31]  = 112;
  krad[32]  = 112;
  krad[33]  = 120;
  krad[34]  = 128;
  krad[35]  = 128;
  krad[36]  = 132;
  krad[37]  = 140;
  krad[38]  = 144;
  krad[39]  = 154;
  krad[40]  = 160;
  krad[41]  = 168;
  krad[42]  = 176;
  krad[43]  = 180;
  krad[44]  = 192;
  krad[45]  = 198;
  krad[46]  = 210;
  krad[47]  = 220;
  krad[48]  = 224;
  krad[49]  = 240;
  krad[50]  = 252;
  krad[51]  = 256;
  krad[52]  = 264;
  krad[53]  = 280;
  krad[54]  = 288;
  krad[55]  = 308;
  krad[56]  = 320;
  krad[57]  = 330;
  krad[58]  = 336;
  krad[59]  = 352;
  krad[60]  = 360;
  krad[61]  = 384;
  krad[62]  = 396;
  krad[63]  = 420;
  krad[64]  = 440;
  krad[65]  = 448;
  krad[66]  = 462;
  krad[67]  = 480;
  krad[68]  = 504;
  krad[69]  = 512;
  krad[70]  = 528;
  krad[71]  = 560;
  krad[72]  = 576;
  krad[73]  = 616;
  krad[74]  = 630;
  krad[75]  = 640;
  krad[76]  = 660;
  krad[77]  = 672;
  krad[78]  = 704;
  krad[79]  = 720;
  krad[80]  = 768;
  krad[81]  = 770;
  krad[82]  = 792;
  krad[83]  = 840;
  krad[84]  = 880;
  krad[85]  = 896;
  krad[86]  = 924;
  krad[87]  = 960;
  krad[88]  = 990;
  krad[89]  = 1008;
  krad[90]  = 1024;
  krad[91]  = 1056;
  krad[92]  = 1120;
  krad[93]  = 1152;
  krad[94]  = 1232;
  krad[95]  = 1260;
  krad[96]  = 1280;
  krad[97]  = 1320;
  krad[98]  = 1344;
  krad[99]  = 1386;
  krad[100] = 1408;
  krad[101] = 1440;
  krad[102] = 1536;
  krad[103] = 1540;
  krad[104] = 1584;
  krad[105] = 1680;
  krad[106] = 1760;
  krad[107] = 1792;
  krad[108] = 1848;
  krad[109] = 1920;
  krad[110] = 1980;
  krad[111] = 2016;
  krad[112] = 2048;
  krad[113] = 2112;
  krad[114] = 2240;
  krad[115] = 2304;
  krad[116] = 2310;
  krad[117] = 2464;
  krad[118] = 2520;
  krad[119] = 2560;
  krad[120] = 2640;
  krad[121] = 2688;
  krad[122] = 2772;
  krad[123] = 2816;
  krad[124] = 2880;
  krad[125] = 3072;
  krad[126] = 3080;
  krad[127] = 3168;
  krad[128] = 3360;
  krad[129] = 3520;
  krad[130] = 3584;
  krad[131] = 3696;
  krad[132] = 3840;
  krad[133] = 3960;
  krad[134] = 4032;
  krad[135] = 4096;
  krad[136] = 4224;
  krad[137] = 4480;
  krad[138] = 4608;
  krad[139] = 4620;
  krad[140] = 4928;
  krad[141] = 5040;
  krad[142] = 5120;
  krad[143] = 5280;
  krad[144] = 5376;
  krad[145] = 5544;
  krad[146] = 5632;
  krad[147] = 5760;
  krad[148] = 6144;
  krad[149] = 6160;
  krad[150] = 6336;
  krad[151] = 6720;
  krad[152] = 6930;
  krad[153] = 7040;
  krad[154] = 7168;
  krad[155] = 7392;
  krad[156] = 7680;
  krad[157] = 7920;
  krad[158] = 8064;
  krad[159] = 8192;
  krad[160] = 8448;
  krad[161] = 8960;
  krad[162] = 9216;
  krad[163] = 9240;
  krad[164] = 9856;
  krad[165] = 10080;
  krad[166] = 10240;
  krad[167] = 10560;
  krad[168] = 10752;
  krad[169] = 11088;
  krad[170] = 11264;
  krad[171] = 11520;
  krad[172] = 12288;
  krad[173] = 12320;
  krad[174] = 12672;
  krad[175] = 13440;
  krad[176] = 13860;
  krad[177] = 14080;
  krad[178] = 14336;
  krad[179] = 14784;

  /*--------------------------------------------------------------------------*/
}/*end routine */
/*==========================================================================*/
