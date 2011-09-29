//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime>

//--------------------------------------------------------------------------

typedef struct complex{
    double  re;
    double  im;   
};
#define PRINTF printf
#define EXIT(N) {exit(N);}

int main();
void flip_data_set(int , int *, int *, int *,complex *);

//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
int main(){
//==========================================================================

  int nktot,nktot2,n1,n2,n3;
  int *kx,*ky,*kz;
  complex *data;
  char fname[1024];
  char directory[1024];
  FILE *fp;

//==========================================================================

  int nstate = 128;
  int ibinary=0;
  sprintf(directory,"Spin.0_Kpt.0_Bead.0_Temper.0");

  PRINTF("I am converting %d doublePacked states to singlePack form\n",nstate);
  PRINTF("I am assuming doublePacked states are stored in the directory %s\n",directory);
  PRINTF("I will write the singlePacked states in the directory %s_flipped\n",directory);
  PRINTF("I am using the binary option %d : (0/1) = (ascii/binary)\n",ibinary);
  PRINTF("The user can change the defaults on lines 41-43 of the source\n");
  PRINTF("The motivated user can write a nice input file and check in\n");

//==========================================================================
// Read in the header information and malloc some memory

  if(ibinary==0){
   sprintf(fname,"%s/state1.out",directory);
   fp = fopen(fname,"r");
   fscanf(fp,"%d %d %d %d",&nktot,&n1,&n2,&n3);
  }else{
    fp = fopen(fname,"rb");
    int n=1;
    fread(&(nktot),sizeof(int),n,fp);
  }//endif
  printf("nktot = %d\n",nktot);
  fclose(fp);

  nktot2 = 2*nktot-1;
  data = (complex *)malloc(nktot2*sizeof(complex));    
  kx   = (int *)malloc(nktot2*sizeof(int));
  ky   = (int *)malloc(nktot2*sizeof(int));
  kz   = (int *)malloc(nktot2*sizeof(int));

//==========================================================================
// Read in the each state, flip it, and write the flipped guy back out

  for(int is=0;is<nstate;is++){
//---------------------------------------------------------------------------
// Read in the state
    sprintf(fname,"%s/state%d.out",directory,is + 1);
    if(ibinary==0){
      fp = fopen(fname,"r");
      printf("Reading file: %s\n",fname);
      fscanf(fp,"%d %d %d %d",&nktot,&n1,&n2,&n3);
      for(int i =0;i<nktot;i++){
        fscanf(fp,"%lf %lf %d %d %d",&data[i].re,&data[i].im,&kx[i],&ky[i],&kz[i]);    
      }//endfor
    }else{
      fp = fopen(fname,"rb");
      int n=1;
      fread(&(nktot),sizeof(int),n,fp);
      fread(&(n1),sizeof(int),n,fp);
      fread(&(n2),sizeof(int),n,fp);
      fread(&(n3),sizeof(int),n,fp);
      for(int i=0;i<nktot;i++) {
        int x,y,z;
        double re,im;
        fread(&(re),sizeof(double),n,fp);
        fread(&(im),sizeof(double),n,fp);
        fread(&(x),sizeof(int),n,fp);
        fread(&(y),sizeof(int),n,fp);
        fread(&(z),sizeof(int),n,fp);
        data[i].re = re; data[i].im=im;
        kx[i]    = x;
        ky[i]    = y;
        kz[i]    = z; 
      }//endfor 
    }//endif
    fclose(fp);
//-----------------------------------------------------------------------------------
// Flip the state
    flip_data_set(nktot,kx,ky,kz,data);

    sprintf(fname,"%s_flipped/state%d.out",directory,is + 1);
//---------------------------------------------------------------------------
// Write the flipped state out
    if(ibinary==0){
      fp = fopen(fname,"w");
      printf("Writing file: %s\n",fname);
      fprintf(fp,"%d %d %d %d\n",nktot2,n1,n2,n3);
      for(int i =0;i<nktot2;i++){
        fprintf(fp,"%lf %lf %d %d %d\n",data[i].re,data[i].im,kx[i],ky[i],kz[i]);    
      }//endfor
    }else{
      fp = fopen(fname,"wb");
      int n=1;
      fwrite(&(nktot2),sizeof(int),n,fp);
      fwrite(&(n1),sizeof(int),n,fp);
      fwrite(&(n2),sizeof(int),n,fp);
      fwrite(&(n3),sizeof(int),n,fp);
      double sum=0.0;
      for(int i=0;i<nktot2;i++) {
        int x,y,z;
        double re,im;
        re = data[i].re; im = data[i].im; 
        sum += re*re + im*im;
        x = kx[i];
        y = ky[i];
        z = kz[i];
        fwrite(&(re),sizeof(double),n,fp);
        fwrite(&(im),sizeof(double),n,fp);
        fwrite(&(x),sizeof(int),n,fp);
        fwrite(&(y),sizeof(int),n,fp);
        fwrite(&(z),sizeof(int),n,fp);
      }//endfor 
      //printf("%.10g \n",sum);
    }//endif
    fclose(fp);
  }//endfor; is=states

//==========================================================================

  return 1;

//--------------------------------------------------------------------------
  }//end main
//==========================================================================



//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Take the funny piny output and reorder it to full complex beast
//  kx,ky,kz had better be malloced nktot*2 - 1 but contain nktot data
//  at input
//==========================================================================

void flip_data_set(int nktot, int *kx, int *ky, int *kz,complex *data)

//==========================================================================
    {//begin routine 
//==========================================================================
// Count half plane kx=0 of piny data : check piny data

  int nplane0 = 0;
  for(int i=0;i<nktot;i++){
    if(kx[i]==0){nplane0++;}
  }//endfor
  int *kxt = (int *)malloc(nplane0*sizeof(int));
  int *kyt = (int *)malloc(nplane0*sizeof(int));
  int *kzt = (int *)malloc(nplane0*sizeof(int));

  complex *datat;
  datat = (complex *)malloc(nplane0*sizeof(complex));

  for(int i=0;i<nplane0-1;i++){
    kxt[i]= kx[i];
    kyt[i]= ky[i];
    kzt[i]= kz[i];
    datat[i]= data[i];
    if(kx[i]!=0){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("Error.1 while flipping piny dblpack data set\n");
      PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT(1);
    }//endif
  }//endfor
  kxt[(nplane0-1)]   = kx[(nktot-1)];
  kyt[(nplane0-1)]   = ky[(nktot-1)];
  kzt[(nplane0-1)]   = kz[(nktot-1)];
  datat[(nplane0-1)] = data[(nktot-1)];

  if(kx[(nktot-1)]!=0 || ky[(nktot-1)] || kz[(nktot-1)]){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Error.2 while flipping piny dblpack data set\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
    EXIT(1);
  }//endif

//==========================================================================
// Expand the data set

  for(int i=nktot-2;i>=0;i--){
    kx[(i+nplane0)]   = kx[i];
    ky[(i+nplane0)]   = ky[i];
    kz[(i+nplane0)]   = kz[i];
    data[(i+nplane0)] = data[i];
  }//endfor

//==========================================================================
// Create the bottom half of plane zero by symmetry : 

  int i1 = 0;
  for(int i=0;i<nplane0-1;i++){
    int ind = nplane0-i-2;
    kx[i]      =  kxt[ind];
    ky[i]      = -kyt[ind];
    kz[i]      = -kzt[ind];
    data[i].re =  datat[ind].re;
    data[i].im = -datat[ind].im;
    if(kx[i]!=0 || ky[i]<ky[i1]){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("Error.3 while flipping piny dblpack data set %d %d %d %d %d\n",kx[i],ky[i],kz[i],i,ind);
      PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT(1);
    }//endif
    i1 = i;
  }//endfor
  kx[(nplane0-1)]   = kxt[(nplane0-1)];
  ky[(nplane0-1)]   = kyt[(nplane0-1)];
  kz[(nplane0-1)]   = kzt[(nplane0-1)];
  data[(nplane0-1)] = datat[(nplane0-1)];

  if(nktot>=nplane0+1){
    if(kx[nplane0]!= 0 || ky[nplane0]!=0 || kz[nplane0]!=1){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("Error.4 while flipping piny dblpack data set\n");
      PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT(1);
    }//endif
  }//endfor

#ifdef _GLENN_DBG_FLIP
  int nnn = MIN(nplane0+3,nktot);
  for(int i=0;i<nnn;i++){
    PRINTF(" %d : %d %d %d \n",i,kx[i],ky[i],kz[i]);
  }
#endif

//==========================================================================
// Exit

  free(kxt);
  free(kyt);
  free(kzt);
  free(datat);

//==========================================================================
// Now we have full planes in the upper half space. We can flip them over

  int nnow = nktot + nplane0 - 1;
  int nktot2 = 2*nktot - 1;

  datat = (complex *)malloc(nnow*sizeof(complex));
  kxt = (int *)malloc(nnow*sizeof(int));
  kyt = (int *)malloc(nnow*sizeof(int));
  kzt = (int *)malloc(nnow*sizeof(int));
  for(int i=0;i<nnow;i++){
    kxt[i]   = kx[i];
    kyt[i]   = ky[i];
    kzt[i]   = kz[i];
    datat[i] = data[i];
  }//endfor

  // The top half + plane0 is now in place!!!
  int ioff = nktot2 - nnow;
  for(int i=0,j=ioff;i<nnow;i++,j++){
    kx[j]    = kxt[i];
    ky[j]    = kyt[i];
    kz[j]    = kzt[i];
    data[j]  = datat[i];
  }

  // Get the bottom half using hermitian conjg sym!
  for(int i=0,j=nktot2-1;i<nnow-2*nplane0+1;i++,j--){
    kx[i]       =  -kx[j];
    ky[i]       =  -ky[j];
    kz[i]       =  -kz[j];
    data[i].re =  data[j].re;
    data[i].im = -data[j].im;
  }

//==========================================================================
    }//end routine
//==========================================================================
