/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*   This program converts piny_md configs to angstroms                     */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#define MAXWORD 50
#define IBM_ESSL_OFF
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define NINT(X) ( (int) ((X)>=0.0 ? ((X)+0.5):((X)-0.5)) )
#define MAX(A,B) (((A)>(B))?(A):(B))
#define MIN(A,B) (((A)<(B))?(A):(B))

void readtoendofline(FILE *);
FILE *cfopen(char [],char *);
void init_kvec(double *,double *,double *,double,int *,int *,int *,int *);
void make_map(double ,double ,double *,int ,int ,int *,int *,int *);
void *cmalloc(size_t );
void countkvec3d(int *,double ,int *,double *);
void setkvec3d(int ,double ,int *,double *,int *, int *, int *);
void calc_cutoff(int , double *,double *,int ,int *, int *, double *, double );
void radixme(int *, int *, int *);
void gethinv(double *, double *, double *, int );

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
 int main (int argc, char *argv[])
/*==========================================================================*/
   {/* begin routine */ 
/*==========================================================================*/
/*   Local Variables */

  int i,is,ic,iii,n;

  int ncoef;                /* # of coef    */
  int ncoef_old;            /* # of coef    */
  int ncoef_new;            /* # of coef    */
  int nstate_up, nstate_dn; /* # of states  */
  int izero,itime;
  int *map;                 /* map          */
  int  kmaxv_new[4];        /* kvec limit   */
  int  kmaxv_old[4];        /* kvec limit   */
  int ibinary;
  int csize=MAXWORD;

  double ecut_old,ecut_new; /* cutoffs      */
  double occ_up,occ_dn;     /* occupation # */
  double hmat[10];          /* h-matrix     */
  double hmati[10];         /* h^-1 matrix  */
  double deth;              /* Volume       */

  float *cre,*cim;          /* coef vectors */
  float zero=0.0;

  char input_file[MAXWORD];
  char output_file[MAXWORD];
  char binary_opt[MAXWORD];
  char c_array1[MAXWORD];
  char c_array2[MAXWORD];
  char c_array3[MAXWORD];
  char c_array4[MAXWORD];
  char c_array5[MAXWORD];
  char c_array6[MAXWORD];
  char c_array7[MAXWORD];
  char c_array8[MAXWORD];
  char c_array9[MAXWORD];
  char c_array10[MAXWORD];

  FILE *fp_in,*fp_out;      /* Files        */
  FILE *fp_info;            /* Files        */

/*==========================================================================*/
/* 0) Get the info                                                          */

   fp_info = cfopen(argv[1],"r");
   fscanf(fp_info,"%s",&binary_opt);                readtoendofline(fp_info);
   fscanf(fp_info,"%s %lf",&input_file,&ecut_old);  readtoendofline(fp_info);
   fscanf(fp_info,"%s %lf",&output_file,&ecut_new); readtoendofline(fp_info);
   for(i=0;i<=2;i++){
     fscanf(fp_info,"%lf %lf %lf",&hmat[1+i],&hmat[4+i],&hmat[7+i]);
     readtoendofline(fp_info);
   }/*endfor*/
   fclose(fp_info);

   for(i=1;i<=9;i++){
     hmat[i] /= 0.529177;
   }
   gethinv(hmat,hmati,&deth,3);    

   ibinary = -1;
   if(strcasecmp(binary_opt,"binary_off")==0){ibinary=0;}
   if(strcasecmp(binary_opt,"binary_on") ==0){ibinary=1;}

   if(ibinary==-1){
     printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     printf("Chaos reigns within.\n");
     printf("Reflect, repent, and retype.\n");
     printf("Binary_opt is on or off.\n");
     printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);
     exit(1);
   }/*endif*/

/*==========================================================================*/
/* I) Initialize                                                            */

   init_kvec(&ecut_old,&ecut_new,hmati,deth,&ncoef_old,&ncoef_new,
             kmaxv_old,kmaxv_new);

   if(ncoef_old>ncoef_new){
     printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     printf("Chaos reigns within.\n");
     printf("Reflect, repent, and retype.\n");
     printf("More old coefs found than new.\n");
     printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);
     exit(1);
   }/*endif*/

   map = (int *)cmalloc(ncoef_new*sizeof(int))-1;
   make_map(ecut_old,ecut_new,hmati,ncoef_old,ncoef_new,kmaxv_new,kmaxv_old,
            map);

/*==========================================================================*/
/* II) Carefully open/output the input file:                                */

  if(ibinary==1){

    fp_in   = cfopen(input_file,"rb");
    fp_out  = cfopen(output_file,"wb");

  }else{

    fp_in   = cfopen(input_file,"r");
    fp_out  = cfopen(output_file,"w");

  }/*endif:binary*/

/*==========================================================================*/
/* III) Read/write the number of coefs/states                               */

  if(ibinary==1){

    fread(c_array1,sizeof(char),csize,fp_in);
    fread(c_array2,sizeof(char),csize,fp_in);
    fread(c_array3,sizeof(char),csize,fp_in);
    fread(c_array4,sizeof(char),csize,fp_in);
    fread(c_array5,sizeof(char),csize,fp_in);
    fread(c_array6,sizeof(char),csize,fp_in);
    fread(c_array7,sizeof(char),csize,fp_in);  
    n = 1;
    fread(&ncoef,sizeof(int),n,fp_in);
    fread(&izero,sizeof(int),n,fp_in);
    fread(&nstate_up,sizeof(int),n,fp_in);
    fread(&nstate_dn,sizeof(int),n,fp_in);
    fread(c_array8,sizeof(char),csize,fp_in);
    fread(c_array9,sizeof(char),csize,fp_in);
    fread(&itime,sizeof(int),n,fp_in);
    fread(c_array10,sizeof(char),csize,fp_in); 

    fwrite(c_array1,sizeof(char),csize,fp_out);
    fwrite(c_array2,sizeof(char),csize,fp_out);
    fwrite(c_array3,sizeof(char),csize,fp_out);
    fwrite(c_array4,sizeof(char),csize,fp_out);
    fwrite(c_array5,sizeof(char),csize,fp_out);
    fwrite(c_array6,sizeof(char),csize,fp_out);
    fwrite(c_array7,sizeof(char),csize,fp_out);  
    n = 1;
    fwrite(&ncoef_new,sizeof(int),n,fp_out);
    fwrite(&izero,sizeof(int),n,fp_out);
    fwrite(&nstate_up,sizeof(int),n,fp_out);
    fwrite(&nstate_dn,sizeof(int),n,fp_out);
    fwrite(c_array8,sizeof(char),csize,fp_out);
    fwrite(c_array9,sizeof(char),csize,fp_out);
    fwrite(&itime,sizeof(int),n,fp_out);
    fwrite(c_array10,sizeof(char),csize,fp_out); 

  }else{

    readtoendofline(fp_in);
    fprintf(fp_out," \n");

    fscanf(fp_in,"%d %d %d %d",&ncoef,&iii,&nstate_up,&nstate_dn);
           readtoendofline(fp_in);
    fprintf(fp_out,"%d 0 %d %d lda restart_pos 1 norb_off\n",
                       ncoef_new,nstate_up,nstate_dn);
  }/*endif:binary*/

  if(ncoef!=ncoef_old){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The Tao that is read,\n");
    printf("Is not the true Tao, until,\n");
    printf("You specify the correct cutoff:\n");
    printf("   ncoef_info %d ncoef_dump %d\n",ncoef_old,ncoef);
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }/*endif*/

/*==========================================================================*/
/* IV) Read/write the occupations numbers                                  */


  if(ibinary==1){

    n = 1;
    fread (c_array1,sizeof(char),csize,fp_in);
    fread (c_array2,sizeof(char),csize,fp_in);
    fwrite(c_array1,sizeof(char),csize,fp_out);
    fwrite(c_array2,sizeof(char),csize,fp_out);

    for(i=1;i<=nstate_up;i++){
      fread(&occ_up,sizeof(double),n,fp_in);
      fread(&occ_dn,sizeof(double),n,fp_in);
      fwrite(&occ_up,sizeof(double),n,fp_out);
      fwrite(&occ_dn,sizeof(double),n,fp_out);
    }/*endfor*/

  }else{

    readtoendofline(fp_in);
    fprintf(fp_out," \n");

    for(i=1;i<=nstate_up;i++){
      fscanf(fp_in,"%lf %lf",&occ_up,&occ_dn);readtoendofline(fp_in);
      fprintf(fp_out,"%g %g\n",occ_up,occ_dn);
    }/*endfor*/

  }/*endif:binary*/

/*==========================================================================*/
/* V) Read each state and write it out zero filled                          */

  if(ibinary==1){

    fread (c_array1,sizeof(char),csize,fp_in);
    fread (c_array2,sizeof(char),csize,fp_in);
    fwrite(c_array1,sizeof(char),csize,fp_out);
    fwrite(c_array2,sizeof(char),csize,fp_out);
    
  }else{

    readtoendofline(fp_in);
    fprintf(fp_out," \n");

  }/*endif:binary*/

  cre     = (float *)cmalloc(ncoef_old*sizeof(float))-1;
  cim     = (float *)cmalloc(ncoef_old*sizeof(float))-1;
			    
  for(is=1;is<=nstate_up;is++){

    if(ibinary==1){

      n = 1;
      for(i=1;i<=ncoef;i++){
        fread(&(cre[i]),sizeof(float),n,fp_in);
        fread(&(cim[i]),sizeof(float),n,fp_in);
      }/*endfor*/
      ic = 0;
      for(i=1;i<=ncoef_new;i++){
        if(map[i]==1 && ic <= ncoef){
          ic++;
          fwrite(&(cre[ic]),sizeof(float),n,fp_out);
          fwrite(&(cim[ic]),sizeof(float),n,fp_out);
        }else{
          fwrite(&zero,sizeof(float),n,fp_out);
          fwrite(&zero,sizeof(float),n,fp_out);
        }/*endif*/
      }/*endfor*/


    }else{

      for(i=1;i<=ncoef;i++){
        fscanf(fp_in,"%f %f",&cre[i],&cim[i]);readtoendofline(fp_in);
      }/*endfor*/
      ic = 0;
      for(i=1;i<=ncoef_new;i++){
        if(map[i]==1 && ic <= ncoef){
          ic++;
          fprintf(fp_out,"%g %g\n",cre[ic],cim[ic]);
        }else{
          fprintf(fp_out,"0.0 0.0\n");
        }/*endif*/
      }/*endfor*/

    }/*endif:binary*/

  }/*endfor:states*/

/*==========================================================================*/
/* V) Close the files                                                       */

  fclose(fp_in);
  fclose(fp_out);

/*--------------------------------------------------------------------------*/
   }/*end routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void init_kvec(double *ecut_old,double *ecut_new,double *hmati,double deth,
               int *ncoef_old,int *ncoef_new,int *kmaxv_old,int *kmaxv_new)

/*==========================================================================*/
/*             Begin routine                                                */
{/* begin routine */
/*==========================================================================*/

  int kmax_cp[4];
  int nktot,iii;
  double ecut_now;

  int cp_on    = 1;
  int kmax_ewd = 0;

/*==========================================================================*/
/* I) Old k-vectors */

   calc_cutoff(kmax_ewd,&ecut_now,ecut_old,cp_on,kmax_cp,kmaxv_old,hmati,
               deth);
   countkvec3d(&nktot,ecut_now,kmax_cp,hmati);
   *ncoef_old = nktot+1;

/*==========================================================================*/
/* II) New k-vectors */

   calc_cutoff(kmax_ewd,&ecut_now,ecut_new,cp_on,kmax_cp,kmaxv_new,hmati,
               deth);
   countkvec3d(&nktot,ecut_now,kmax_cp,hmati);
   *ncoef_new = nktot+1;

/*--------------------------------------------------------------------------*/
   }/*end routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void make_map(double ecut_old,double ecut_new,double *hmati,
              int ncoef_old,int ncoef_new,int *kmaxv_new,
              int *kmaxv_old,int *map)

/*==========================================================================*/
/*             Begin routine                                                */
{/* begin routine */
/*==========================================================================*/
/*  Local Variables */
  
  int i,ic;
  int nktot;
  int *kastr,*kbstr,*kcstr;
  int *kastr_old,*kbstr_old,*kcstr_old;
  int kmax_cp[4];
  double aka,akb,akc,tpi;
  double xk,yk,zk,try;

/*==========================================================================*/
/* I) Make the new kvectors */

   kastr     = (int *)cmalloc(ncoef_new*sizeof(int))-1;
   kbstr     = (int *)cmalloc(ncoef_new*sizeof(int))-1;
   kcstr     = (int *)cmalloc(ncoef_new*sizeof(int))-1;

   kastr_old = (int *)cmalloc(ncoef_old*sizeof(int))-1;
   kbstr_old = (int *)cmalloc(ncoef_old*sizeof(int))-1;
   kcstr_old = (int *)cmalloc(ncoef_old*sizeof(int))-1;

   nktot      = ncoef_new-1;
   kmax_cp[1] = kmaxv_new[1]/2;
   kmax_cp[2] = kmaxv_new[2]/2;
   kmax_cp[3] = kmaxv_new[3]/2;
   setkvec3d(nktot,ecut_new,kmax_cp,hmati,kastr,kbstr,kcstr);

   nktot      = ncoef_old-1;
   kmax_cp[1] = kmaxv_old[1]/2;
   kmax_cp[2] = kmaxv_old[2]/2;
   kmax_cp[3] = kmaxv_old[3]/2;
   setkvec3d(nktot,ecut_old,kmax_cp,hmati,kastr_old,kbstr_old,kcstr_old);

/*==========================================================================*/
/* II) Make the map */

   tpi = 2.0*M_PI;
   ic  = 1;
   for(i=1;i<=ncoef_new-1;i++){
 
     map[i] = 0;
     if(kastr_old[ic]==kastr[i]&&
        kbstr_old[ic]==kbstr[i]&&
        kcstr_old[ic]==kcstr[i]){
       map[i] = 1;
     }/*endif*/ 
     ic    += map[i];

   }/*endfor*/

   map[ncoef_new] = 1;

   if(ic!=ncoef_old){
     printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     printf("Chaos reigns within.\n");
     printf("Reflect, repent, and retype.\n");
     printf("All old coefs not found. %d %d\n",ic,ncoef_old);
     printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);
     exit(1);
   }/*endif*/

/*==========================================================================*/
/* Free the memory */

   free(&kastr[1]);
   free(&kbstr[1]);
   free(&kcstr[1]);
   free(&kastr_old[1]);
   free(&kbstr_old[1]);
   free(&kcstr_old[1]);

/*--------------------------------------------------------------------------*/
   }/*end routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* cfopen:  Function to open files and check existence (careful open)       */
/*==========================================================================*/

FILE *cfopen(char file_name[],char *mode)

/*==========================================================================*/
{  /* begin routine */
  FILE *fp;
  
  if(file_name==NULL){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    printf("With searching comes loss\n");
    printf("And the presence of absence:\n");
    printf("Requesting a file with a null pointer  !\n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }

  if(strlen(file_name)==0){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    printf("With searching comes loss\n");
    printf("And the presence of absence:\n");
    printf("File pointer with no filename!\n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }

  if(mode[0] == 'w'){
    if((fp=fopen(file_name,"w")) == NULL){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      printf("With searching comes loss\n");
      printf("And the presence of absence:\n");
      printf("File \"%s\" not found \n",file_name);
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }
    return fp;
  } else {
    if((fp=fopen(file_name,"r")) == NULL){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      printf("ERROR: can't open \"%s\" for reading (exiting)\n",
                      file_name);
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }
    return fp;
  }
/*--------------------------------------------------------------------------*/
  }/* end routine */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* readtoendofline: Function to read to end of line in read_coord files     */
/*==========================================================================*/
void readtoendofline(FILE *fp){
  int eol,ch;
  eol = (int )'\n';
  ch = eol+1;
  while(ch!=eol){ch=fgetc(fp);}
  if(ch==EOF){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      printf("ABORTED effort:\n");
      printf("Close all that you have.\n");
      printf("You ask way too much.\n");
      printf("Unexpected end of file reached\n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
  }/*endif*/  
}/* end routine */    
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* cmalloc: Careful malloc                                                  */
/*==========================================================================*/
void *cmalloc(size_t len)
{ /* begin routine */
  void *mem_ptr;
  double request;
  
  if(len == 0) return NULL;
  mem_ptr = malloc(len);
  if(mem_ptr == NULL) {
   request = ((double) len)*1.0e-6;
   printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
   printf("Dude, like you've just requested %g MBytes\n",request);
   printf("of memory -- get real\n\n");
   printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
   fflush(stdout);
   exit(1);
  }

  memset(mem_ptr,0,len);

  return mem_ptr;

} /* end routine */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void calc_cutoff(int kmax_ewd, double *ecut,double *ecut_cp,int cp_on,
                 int *kmax_cp, int *kmaxv, double *hmatik, double deth)

/*==========================================================================*/
/*               Begin subprogram:                                          */
      {/*begin routine*/
/*==========================================================================*/
/*               Local variable declarations:                               */

  int iii,nk1,nk2,nk3;
  double rtwoth,tpi,rvol23;
  double d1,d2,d3;
  double try1,try2,try3;
  double temp1,temp2,temp3;

/*==========================================================================*/
/* II) Calculate Ewald Cutoff */

   rtwoth = -(2./3.);
   rvol23 = pow(deth,rtwoth);
   tpi = M_PI * 2.0;
   *ecut = M_PI * .5 * M_PI * (double) (kmax_ewd * kmax_ewd)*rvol23;

/*==========================================================================*/
/* III) Compare Ewald to CP cutoff  */

   if (cp_on == 1) {
     if (*ecut_cp * .5 < *ecut) {
	printf("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");
	printf("warning: ewald cutoff greater than cp cutoff\n");
	printf("%g vs %g\n",*ecut*2.,*ecut_cp);
        printf("using the ewald value\n");
        printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
     }/*endif*/
     d1 = *ecut, d2 = *ecut_cp * .5;
     *ecut = MAX(d1,d2);
   }/*endif*/
   *ecut_cp = *ecut;

/*==========================================================================*/
/* III) Adjust shape of reciprocal space                                    */

   d1 = hmatik[1];  d2 = hmatik[4];   d3 = hmatik[7];
   try1 = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
   d1 = hmatik[2];  d2 = hmatik[5];   d3 = hmatik[8];
   try2 = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
   d1 = hmatik[3];  d2 = hmatik[6];   d3 = hmatik[9];
   try3 = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
   temp1 = sqrt(*ecut * .5) / (M_PI * try1);
   temp2 = sqrt(*ecut * .5) / (M_PI * try2);
   temp3 = sqrt(*ecut * .5) / (M_PI * try3);
   if(cp_on==1){
      kmax_cp[1] = NINT(temp1);
      kmax_cp[2] = NINT(temp2);
      kmax_cp[3] = NINT(temp3);
      radixme(&(kmax_cp[1]),&(kmax_cp[2]),&(kmax_cp[3]));
      kmaxv[1] = kmax_cp[1] << 1;
      kmaxv[2] = kmax_cp[2] << 1;
      kmaxv[3] = kmax_cp[3] << 1;
      nk1 = 4*(kmax_cp[1]+1);
      nk2 = 4*(kmax_cp[2]+1);
      nk3 = 4*(kmax_cp[3]+1);
      printf("Using FFT size %d %d %d\n",nk1,nk2,nk3);
   }else{
      kmaxv[1] = (int) ( 2.0*temp1);
      kmaxv[2] = (int) ( 2.0*temp2);
      kmaxv[3] = (int) ( 2.0*temp3);
    }/*endif*/

/*-------------------------------------------------------------------------*/
     } /* end routine */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Count the number of k vectors on large grid */
/*==========================================================================*/

void countkvec3d(int *nktot,double ecut,int *kmaxv,double *hmatik)

/*==========================================================================*/
/*       Begin routine */
{/*begin routine */
/*==========================================================================*/
/* Local variables */

  int iii,icount;
  int i, kbmin, kcmin, kbmax, kcmax, kamax, ka, kb, kc;
  double xk, yk, zk;
  double aka, akb, akc;
  double tpi,try;

/*==========================================================================*/
/* Count the kvectors */

  tpi = 2.0*M_PI;
  icount = 0;
  kamax = kmaxv[1];

/*********************************/

   for (i = 1; i <= kamax; ++i) {
    aka = (double) i;
    xk = aka * hmatik[1] * tpi;
    yk = aka * hmatik[4] * tpi;
    zk = aka * hmatik[7] * tpi;
    try = (xk * xk + yk * yk + zk * zk) * .5;
    if (try > ecut ) {
      break;
    }
  }

  kamax = i - 1;

/***********************************/

  for (ka = 0; ka <= kamax; ++ka) {
    aka = (double) ka;
    kbmin = -kmaxv[2];
    if (ka == 0) {
      kbmin = 0;
    }
    for (i = kbmin; i <= 0; ++i) {
      akb = (double) i;
      xk = (aka * hmatik[1] + akb * hmatik[2]) * tpi;
      yk = (aka * hmatik[4] + akb * hmatik[5]) * tpi;
      zk = (aka * hmatik[7] + akb * hmatik[8]) * tpi;
      try = (xk * xk + yk * yk + zk * zk) * .5;
      if (try <= ecut ) {
	break;
      }
    }

/*********************************/

    kbmin = i;
    for (i = 1; i <= kmaxv[2]; ++i) {
      akb = (double) i;
      xk = (aka * hmatik[1] + akb * hmatik[2]) * tpi;
      yk = (aka * hmatik[4] + akb * hmatik[5]) * tpi;
      zk = (aka * hmatik[7] + akb * hmatik[8]) * tpi;
      try = (xk * xk + yk * yk + zk * zk) * .5;
      if (try > ecut ) {
	break;
      }
    }
    
    kbmax = i - 1;
    for (kb = kbmin; kb <= kbmax; ++kb) {
      
      akb = (double) kb;
      kcmin = -kmaxv[3];
      if (ka == 0 && kb == 0) {
	kcmin = 1;
      }
      for (i = kcmin; i <= 0; ++i) {
	akc = (double) i;
	xk = (aka * hmatik[1] + akb * hmatik[2] + akc * hmatik[3]) * tpi;
	yk = (aka * hmatik[4] + akb * hmatik[5] + akc * hmatik[6]) * tpi;
	zk = (aka * hmatik[7] + akb * hmatik[8] + akc * hmatik[9]) * tpi;
	try = (xk * xk + yk * yk + zk * zk) * .5;
	if (try <= ecut ) {
	  break;
	}
      }
/*********************************/

      kcmin = i;
      for (i = 1; i <= kmaxv[3]; ++i) {
	akc = (double) i;
	xk = (aka * hmatik[1] + akb * hmatik[2] + akc * hmatik[3]) * tpi;
	yk = (aka * hmatik[4] + akb * hmatik[5] + akc * hmatik[6]) * tpi;
	zk = (aka * hmatik[7] + akb * hmatik[8] + akc * hmatik[9]) * tpi;
	try = (xk * xk + yk * yk + zk * zk) * .5;
	if (try > ecut ) {
	  break;
	}
      }

      kcmax = i - 1;
      akc = (double) kcmin;
      for (kc = kcmin; kc <= kcmax; ++kc) {
	++icount;
      }
    }
  }
  *nktot = icount;
/*--------------------------------------------------------------------------*/
} /* countkvec3d */
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* This subroutine determines the allowed spherically truncated             */
/* half space k (.i.e. g) vectors given a cutoff. it is used by             */
/* both the ewald and cp modules.                                           */
/* Sets up the k vectors for a given cutoff and shape                       */
/*==========================================================================*/

void setkvec3d(int nktot,double ecut,int *kmaxv,double *hmatik,
		int *kastore, int *kbstore, int *kcstore)

/*==========================================================================*/
/*       Begin routine */
{/*begin routine */
/*==========================================================================*/
/* Local variables */

    int iii;
    int i1, i2, i3;

    int i, kbmin, kcmin, kbmax, kcmax, kamax, ka, kb, kc;
    double xk, yk, zk;
    int icount;
    double aka, akb, akc;
    double tpi;
    double try;

/*==========================================================================*/
/* Count the K-vectors */

    tpi = 2.0*M_PI;
    icount = 0;

/*=============================*/

    kamax = kmaxv[1];
    i1 = kamax;
    for (ka = 0; ka <= i1; ++ka) {
      aka = (double) ka;
      kbmin = -kmaxv[2];
      if (ka == 0) {
	kbmin = 0;
      }
      for (i = kbmin; i <= 0; ++i) {
	akb = (double) i;
	xk = (aka * hmatik[1] + akb * hmatik[2]) * tpi;
	yk = (aka * hmatik[4] + akb * hmatik[5]) * tpi;
	zk = (aka * hmatik[7] + akb * hmatik[8]) * tpi;
	try = (xk * xk + yk * yk + zk * zk) * .5;
	if (try <= ecut ) {
	  break;
	}
      }
      kbmin = i;
      i2 = kmaxv[2];
      for (i = 1; i <= i2; ++i) {
	akb = (double) i;
	xk = (aka * hmatik[1] + akb * hmatik[2]) * tpi;
	yk = (aka * hmatik[4] + akb * hmatik[5]) * tpi;
	zk = (aka * hmatik[7] + akb * hmatik[8]) * tpi;
	try = (xk * xk + yk * yk + zk * zk) * .5;
	if (try > ecut ) {
	  break;
	}
      }

/*=============================*/

      kbmax = i - 1;
      i2 = kbmax;
      for (kb = kbmin; kb <= i2; ++kb) {
	akb = (double) kb;
	kcmin = -kmaxv[3];
	if (ka == 0 && kb == 0) {
	  kcmin = 1;
	}
	for (i = kcmin; i <= 0; ++i) {
	  akc = (double) i;
	  xk = (aka * hmatik[1] + akb * hmatik[2] + akc * hmatik[3]) * tpi;
	  yk = (aka * hmatik[4] + akb * hmatik[5] + akc * hmatik[6]) * tpi;
	  zk = (aka * hmatik[7] + akb * hmatik[8] + akc * hmatik[9]) * tpi;
	  try = (xk * xk + yk * yk + zk * zk) * .5;
	  if (try <= ecut ) {
	    break;
	  }
	}

/*=============================*/

	kcmin = i;
	i3 = kmaxv[3];
	for (i = 1; i <= i3; ++i) {
	  akc = (double) i;
	  xk = (aka * hmatik[1] + akb * hmatik[2] + akc * hmatik[3]) * tpi;
	  yk = (aka * hmatik[4] + akb * hmatik[5] + akc * hmatik[6]) * tpi;
	  zk = (aka * hmatik[7] + akb * hmatik[8] + akc * hmatik[9]) * tpi;
	  try = (xk * xk + yk * yk + zk * zk) * .5;
	  if (try > ecut ) {
	    break;
	  }
	}
	kcmax = i - 1;
	i3 = kcmax;
	for (kc = kcmin; kc <= i3; ++kc) {
	  ++icount;
	  aka = (double) ka;
	  akb = (double) kb;
	  akc = (double) kc;
	  xk = (aka * hmatik[1] + akb * hmatik[2] + akc * hmatik[3]) * tpi;
	  yk = (aka * hmatik[4] + akb * hmatik[5] + akc * hmatik[6]) * tpi;
	  zk = (aka * hmatik[7] + akb * hmatik[8] + akc * hmatik[9]) * tpi;
	  kastore[icount] = ka;
	  kbstore[icount] = kb;
	  kcstore[icount] = kc;
	}
      }
    }
    if(nktot!=icount){
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Mismatch number of kvectors\n");
        printf("%d vs %d\n",icount,nktot);
        printf("Contact technical support\n");
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
    }
    ++icount;
    kastore[icount] = 0;
    kbstore[icount] = 0;
    kcstore[icount] = 0;

/*--------------------------------------------------------------------------*/
 } /* setkvec3d */
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* This routine makes sure that the k1,k2,k3 for the fft                    */
/*  satisfy the radix condition                                             */
/*==========================================================================*/

void radixme(int *kmax1, int *kmax2, int *kmax3)

/*==========================================================================*/
/* Calculate the quantity to be radicized */
/* it written with the plus one to get rid */
/* of the stupid annoying edge vectors. */
/* the factor of 4 appears because the kmax's are the */
/* maximum k vector along a direction. the normal fft */
/* grid is therefore (2*(kmax1+1))(2*(kmax2+1))(2*(kmax3+1)). */
/* the additional factor of two comes from the fact the density */
/* needs to be defined on twice as fine an fft grid */
/* (4*(kmax1+1))(4*(kmax2+1))(4*(kmax3+1)). */
/*==========================================================================*/
/*       Begin routine */
{/*begin routine */
/*==========================================================================*/
/* Local variables */

  int i1, i2;
  int krad[100], nrad, i, k1, k2, k3, iii, jjj, kkk;

/*==========================================================================*/
  k1 = (*kmax1 + 1) << 2;
  k2 = (*kmax2 + 1) << 2;
  k3 = (*kmax3 + 1) << 2;

  nrad = 50;

/* radix condition IBM_ESSL (2^h)(3^i)(5^j)(7^k)(11^m)                */
/* radix condition HP_VECLIB SGI_COMPLIB   (2^h)(3^i)(5^j)             */
#ifdef IBM_ESSL
  krad[0] = 4;
  krad[1] = 8;
  krad[2] = 12;
  krad[3] = 16;
  krad[4] = 20;
  krad[5] = 24;
  krad[6] = 32;
  krad[7] = 36;
  krad[8] = 40;
  krad[9] = 48;
  krad[10] = 60;
  krad[11] = 64;
  krad[12] = 72;
  krad[13] = 80;
  krad[14] = 96;
  krad[15] = 112;
  krad[16] = 112;
  krad[17] = 120;
  krad[18] = 128;
  krad[19] = 144;
  krad[20] = 160;
  krad[21] = 180;
  krad[22] = 192;
  krad[23] = 220;
  krad[24] = 224;
  krad[25] = 240;
  krad[26] = 252;
  krad[27] = 288;
  krad[28] = 308;
  krad[29] = 320;
  krad[30] = 336;
  krad[31] = 360;
  krad[32] = 384;
  krad[33] = 420;
  krad[34] = 440;
  krad[35] = 480;
  krad[36] = 504;
  krad[37] = 512;
  krad[38] = 560;
  krad[39] = 576;
  krad[40] = 616;
  krad[41] = 640;
  krad[42] = 660;
  krad[43] = 720;
  krad[44] = 768;
  krad[45] = 792;
  krad[46] = 880;
  krad[47] = 896;
  krad[48] = 960;
  krad[49] = 1008;
#else
  krad[0]=   2;
  krad[1]=   4;
  krad[2]=   6;
  krad[3]=   8;
  krad[4]=  10;
  krad[5]=  12;
  krad[6]=  14;
  krad[7]=  16;
  krad[8]=  18;
  krad[9]=  20;
  krad[10]=  22;
  krad[11]=  24;
  krad[12]=  28;
  krad[13]=  30;
  krad[14]=  32;
  krad[15]=  36;
  krad[16]=  40;
  krad[17]=  42;
  krad[18]=  44;
  krad[19]=  48;
  krad[20]=  56;
  krad[21]=  60;
  krad[22]=  64;
  krad[23]=  66;
  krad[24]=  70;
  krad[25]=  72;
  krad[26]=  80;
  krad[27]=  84;
  krad[28]=  88;
  krad[29]=  90;
  krad[30]=  96;
  krad[31]= 110;
  krad[32]= 112;
  krad[33]= 120;
  krad[34]= 126;
  krad[35]= 128;
  krad[36]= 132;
  krad[37]= 140;
  krad[38]= 144;
  krad[39]= 154;
  krad[40]= 160;
  krad[41]= 168;
  krad[42]= 176;
  krad[43]= 180;
  krad[44]= 192;
  krad[45]= 198;
  krad[46]= 210;
  krad[47]= 220;
  krad[48]= 224;
  krad[49]= 240;
  krad[50]= 252;
  krad[51]= 256;
  krad[52]= 264;
  krad[53]= 280;
  krad[54]= 288;
  krad[55]= 308;
  krad[56]= 320;
  krad[57]= 330;
  krad[58]= 336;
  krad[59]= 352;
  krad[60]= 360;
  krad[61]= 384;
  krad[62]= 396;
  krad[63]= 420;
  krad[64]= 440;
  krad[65]= 448;
  krad[66]= 462;
  krad[67]= 480;
  krad[68]= 504;
  krad[69]= 512;
  krad[70]= 528;
  krad[71]= 560;
  krad[72]= 576;
  krad[73]= 616;
  krad[74]= 630;
  krad[75]= 640;
  krad[76]= 660;
  krad[77]= 672;
  krad[78]= 704;
  krad[79]= 720;
  krad[80]= 768;
  krad[81]= 770;
  krad[82]= 792;
  krad[83]= 840;
  krad[84]= 880;
  krad[85]= 896;
  krad[86]= 924;
  krad[87]= 960;
  krad[88]= 990;
  krad[89]=1008;
  krad[90]=1024;
  krad[91]=1056;
  krad[92]=1120;
  krad[93]=1152;
  krad[94]=1232;
  krad[95]=1260;
  krad[96]=1280;
  krad[97]=1320;
  krad[98]=1344;
  krad[99]=1386;
#endif

  i1 = nrad;
  for (i = 1; i <= i1; ++i) {
    if (krad[i - 1] > k1) {
      i2 = i - 1;
      iii = MAX(i2,1);
      jjj = krad[(i - 1)];
      kkk = krad[(iii - 1)];
      if (k1 > kkk && k1 < jjj) {
	jjj = (i2 = k1 - krad[(i - 1)], abs(i2));
	kkk = (i2 = k1 - krad[(iii - 1)], abs(i2));
	if (jjj < kkk) {
	  k1 = krad[(i - 1)];
	} else {
	  k1 = krad[(iii - 1)];
	}
      }
    }
    if (krad[(i - 1)] > k2) {
      i2 = i - 1;
      iii = MAX(i2,1);
      jjj = krad[(i - 1)];
      kkk = krad[(iii - 1)];
      if (k2 > kkk && k2 < jjj) {
	jjj = (i2 = k2 - krad[(i - 1)], abs(i2));
	kkk = (i2 = k2 - krad[(iii - 1)], abs(i2));
	if (jjj < kkk) {
	  k2 = krad[(i - 1)];
	} else {
	  k2 = krad[(iii - 1)];
	}
      }
    }
    if (krad[(i - 1)] > k3) {
      i2 = i - 1;
      iii = MAX(i2,1);
      jjj = krad[(i - 1)];
      kkk = krad[(iii - 1)];
      if (k3 > kkk && k3 < jjj) {
	jjj = (i2 = k3 - krad[(i - 1)], abs(i2));
	kkk = (i2 = k3 - krad[(iii - 1)], abs(i2));
	if (jjj < kkk) {
	  k3 = krad[(i - 1)];
	} else {
	  k3 = krad[(iii - 1)];
	}
      }
    }
  }
  *kmax1 = k1 / 4 - 1;
  *kmax2 = k2 / 4 - 1;
  *kmax3 = k3 / 4 - 1;

/*--------------------------------------------------------------------------*/
} /* radixme */
/*==========================================================================*/



/*===============================================================*/
/*  Inverse of a 3x3 */
/*===============================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===============================================================*/
void gethinv(double *hmat, double *hmati, double *deth, int iperd)
{
  double vol;
  int i;
  /* gets inverse, hmati, of the iperd dimensional matrix hmat */
  /* (stored as a 3 x 3) */
  for(i=1;i<=9;i++){hmati[i]=0.0;}
  if (iperd == 3) {
    vol = (hmat[1] * (hmat[5] * hmat[9] - hmat[8] * hmat[6]) + 
	   hmat[4] * (hmat[8] * hmat[3] - hmat[2] * hmat[9]) + 
	   hmat[7] * (hmat[2] * hmat[6] - hmat[5] * hmat[3]));
    *deth = vol;
    hmati[1] = (hmat[5] * hmat[9] - hmat[8] * hmat[6]) / vol;
    hmati[5] = (hmat[1] * hmat[9] - hmat[7] * hmat[3]) / vol;
    hmati[9] = (hmat[1] * hmat[5] - hmat[4] * hmat[2]) / vol;
    hmati[4] = (hmat[7] * hmat[6] - hmat[4] * hmat[9]) / vol;
    hmati[2] = (hmat[3] * hmat[8] - hmat[2] * hmat[9]) / vol;
    hmati[7] = (hmat[4] * hmat[8] - hmat[7] * hmat[5]) / vol;
    hmati[3] = (hmat[2] * hmat[6] - hmat[3] * hmat[5]) / vol;
    hmati[8] = (hmat[7] * hmat[2] - hmat[8] * hmat[1]) / vol;
    hmati[6] = (hmat[3] * hmat[4] - hmat[6] * hmat[1]) / vol;
  }
  if (iperd == 2) {
    vol = hmat[1] * hmat[5] - hmat[4] * hmat[2];
    hmati[1] = hmat[5] / vol;
    hmati[5] = hmat[1] / vol;
    hmati[4] = -hmat[4] / vol;
    hmati[2] = -hmat[2] / vol;
        
    /* never really in two dimensions so take hmati(3,3) as follows */
    
    hmati[9] = 1. / hmat[9];
    *deth = vol * hmat[9];
  }

  if(iperd==1 || iperd ==0) {
   *deth = hmat[1]*hmat[5]*hmat[9];
   hmati[1] = 1.0/hmat[1];
   hmati[5] = 1.0/hmat[5];
   hmati[9] = 1.0/hmat[9];
  }
  if(*deth==0.0){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");    
    printf("The present volume is zero.                 \n");
    printf("If this is not an error in your input data, \n");
    printf("contact technical support                  \n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  } /*endif*/

} /* gethinv */
/*===============================================================*/

