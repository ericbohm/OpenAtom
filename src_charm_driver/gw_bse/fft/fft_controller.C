#include "fft_controller.h"
#include "fft_routines.h"
#include "controller.h"

static CmiNodeLock fft_plan_lock;
void init_plan_lock() {
  fft_plan_lock = CmiCreateLock();
}

FFTController::FFTController() {
  first_time = true;
  in_pointer = out_pointer = NULL;

  geps = new GSPACE();
  

  // TODO: A group dependency could probably solve this better
  contribute(CkCallback(CkReductionTarget(Controller, fftControllerReady), controller_proxy));
}

void FFTController::do_fftw() {
  fftw_execute(plan);
}

void FFTController:: calc_vcoulb(double* qvec, double* b1, double* b2, double * b3, double shift[3], double alat, double vol, int nkpt, int iq){
//CkPrintf("\nStarted calculate vcoulb\n");fflush(stdout);
    double* vcoulb;
    vcoulb = new double [geps->ng];

    double gx, gy, gz;
    double gq[3];
    const double fact = 4*PI/vol/nkpt;

    for (int i=0; i<geps->ng; i++) {


        if (iq==0){
            gx = geps->ig[i] + shift[0];
            gy = geps->jg[i] + shift[1];
            gz = geps->kg[i] + shift[2];
        }
        else{
            gx = geps->ig[i] + qvec[0];
            gy = geps->jg[i] + qvec[1];
            gz = geps->kg[i] + qvec[2];
        }

        vcoulb[i] = 0;
        for (int j=0; j<3; j++) {
            gq[j] =  gx*b1[j] + gy*b2[j] + gz*b3[j];
            gq[j] *= 2*PI/alat;

            vcoulb[i] += gq[j]*gq[j];
        }
        vcoulb[i] = 1/vcoulb[i];
        vcoulb[i] *= fact;
    }



    std::vector<double> vcoulb_v;
    vcoulb_v.resize(geps->ng);
    for(int i=0;i<geps->ng;i++)
      vcoulb_v[i] = vcoulb[i];
//CkPrintf("\nCalculated vcoulb\n");fflush(stdout);
    controller_proxy.got_vcoulb(vcoulb_v);
}

void FFTController::get_geps(double epsCut, double* qvec, double* b1, double* b2, double * b3, 
                              double alat, int nfft[3]){

  int ndata = nfft[0]*nfft[1]*nfft[2];
  bool accept[ndata];
  int *gx, *gy, *gz;
    
  gx = new int [ndata];
  gy = new int [ndata];
  gz = new int [ndata];

  fftidx_to_gidx(gx, gy, gz, nfft);

//Values would need to be sent to Pmatrix geps

   double gxtmp, gytmp, gztmp;
    double vtmp[3];
    double Ekin;
    int count = 0;
      
    for (int i=0; i<ndata; i++) { //can't be 0?
        gxtmp = gx[i] + qvec[0];
        gytmp = gy[i] + qvec[1];
        gztmp = gz[i] + qvec[2];//iq was removed assuming 0 index - might be wrong, since we have only one node
        /* transfer to cartesian unit to calculate energy */
        Ekin = 0;
        for (int j=0; j<3; j++) {
            vtmp[j] = gxtmp*b1[j] + gytmp*b2[j] + gztmp*b3[j];
            vtmp[j] *= 2*PI/alat;
            Ekin += 0.5 * vtmp[j] * vtmp[j];
        }
        
        if (Ekin <= epsCut) {
            accept[i] = true;
            count += 1;
        }
        else{
            accept[i] = false;
        }

    }

    CkPrintf("[FFT CONTROLLER] Dimension of epsilon matrix = %d\n", count);
    // set values
    geps->ng = count;
    geps->ig = new int [count];
    geps->jg = new int [count];
    geps->kg = new int [count];
   
    int j=0;
   
    for (int i=0; i<ndata; i++) {
        if (accept[i]) {
            geps->ig[j] = gx[i];
            geps->jg[j] = gy[i];
            geps->kg[j] = gz[i];
            j += 1;
        }
    }
   
    if ( j!= count ) {
        CkPrintf(" Oops. Error when reducing gspace!!!");
    }
   
    delete[] gx;
    delete[] gy;
    delete[] gz;

    std::vector<int> accept_v;
    accept_v.resize(ndata);
    for(int i=0;i<ndata;i++)
      if(accept[i])
        accept_v[i] = 1;
      else
        accept_v[i] = 0;

    int epsilon_size = count;
    controller_proxy.got_geps(accept_v, epsilon_size);
}

void FFTController::destroy_fftw_stuff() {
  fftw_destroy_plan(plan);
  fftw_free(in_pointer);
  fftw_free(out_pointer);
  in_pointer = out_pointer = NULL;
}

void FFTController::setup_fftw_3d(int nfft[3], int direction) {
  // check for some dumb input values
  if (nfft[0]<=0 || nfft[1] <=0 || nfft[2] <=0) {
    CkPrintf("setup_fftw_3d routine received illegal value for nfft. \
              nfft should be positive number.");
    CkExit();
  }
  if (!(direction==-1 || direction==1)) {
    CkPrintf("setup_fftw_3d routine received illegal value for direction. \
              FFTW direction must either 1 or -1");
    CkExit();
  }

  // if not the first time and any parameters (size, direction) mismatch
  // we need to destroy old plans to set up new ones
  if (!first_time && (direction != old_direction ||
      nfft[0] != old_nfft[0] ||
      nfft[1] != old_nfft[1] ||
      nfft[2] != old_nfft[2])) {
    destroy_fftw_stuff();
    first_time = true;
  }

  // if first time, we need to set up
  if(first_time) {
    const int ndata = nfft[0]*nfft[1]*nfft[2];
    in_pointer = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*ndata);
    out_pointer = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*ndata);

    CmiLock(fft_plan_lock);
    plan = fftw_plan_dft_3d(nfft[0], nfft[1], nfft[2],
        in_pointer, out_pointer, direction, FFTW_ESTIMATE);
    CmiUnlock(fft_plan_lock);

  }

  // now the old value changes to new value
  first_time = false;
  old_nfft[0] = nfft[0];
  old_nfft[1] = nfft[1];
  old_nfft[2] = nfft[2];
  old_direction = direction;
}
