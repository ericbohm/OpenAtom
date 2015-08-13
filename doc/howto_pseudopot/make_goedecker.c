/*==================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==================================================================*/
/*                                                                  */
/*                Create Goedecker Pseudo files                     */
/*                    Yeah, Dawn!                                   */
/*                                                                  */
/*==================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*==================================================================*/
/*               Structures and Typedefs                            */
/*==================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==================================================================*/

#define MAX_CHAR  50

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct{int nc,lmax;
               double z,zv,rloc,*c,erf_alp;
               int *nrad_l;
               double **h,*rh;
              } GOEDECKER;

typedef struct{int nr;
               double rmax,dr;
               double *vloc;
               double *vlong;
               double **p0,**p1,**p2,**p3;
              } PSEUDO;

void readtoendofline(FILE *);
void read_goedecker(GOEDECKER *,char *);
void print_goedecker(GOEDECKER *);
void make_off_diagonal_h(GOEDECKER *);
void make_goedecker_local(GOEDECKER *,PSEUDO *);
void check_goedecker_local(GOEDECKER *);
void make_goedecker_projectors(GOEDECKER *,PSEUDO *);
void check_goedecker_projectors(GOEDECKER *,PSEUDO *);
void check_goedecker_gspace(GOEDECKER *);
void pseudo_initial(PSEUDO *,int *);
void make_gamma_function(double *);
void write_goedecker_file(GOEDECKER *,PSEUDO *);
void check_nonlocal_channels(GOEDECKER *,PSEUDO *);

/*==================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==================================================================*/



/*==================================================================*/
/*               Main program                                       */
/*==================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==================================================================*/
   int main(void)
/*==================================================================*/
  {/*begin routine*/
/*==================================================================*/

   GOEDECKER *goedecker;
   PSEUDO *pseudo;
   char *gfile_name;

/*------------------------------------------------------------------*/

   pseudo    = (PSEUDO *) calloc(1,sizeof(PSEUDO));

   goedecker = (GOEDECKER *) calloc(1,sizeof(GOEDECKER));
   gfile_name = (char *) calloc(MAX_CHAR,sizeof(char ));

/*------------------------------------------------------------------*/
/* Read in parameter file for Goedecker Pseudo potential            */

   printf("Enter the name of the Goedecker pseudo potential file \n");
   fflush(stdin);
   scanf("%s",gfile_name);

   printf("Reading file %s \n",gfile_name);

   read_goedecker(goedecker,gfile_name);
   make_off_diagonal_h(goedecker);
   print_goedecker(goedecker);

/*------------------------------------------------------------------*/
/* Create Pseudo Potential    */

   pseudo_initial(pseudo,goedecker->nrad_l);

   make_goedecker_local(goedecker,pseudo);
   check_goedecker_local(goedecker);
   make_goedecker_projectors(goedecker,pseudo);

   check_goedecker_projectors(goedecker,pseudo);

   check_goedecker_gspace(goedecker);

   write_goedecker_file(goedecker,pseudo);

   check_nonlocal_channels(goedecker,pseudo);

/*==================================================================*/
    return 0;
    }/*end routine*/
/*==================================================================*/



/*==================================================================*/
/*               Read in the parameters                             */
/*==================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==================================================================*/
void read_goedecker(GOEDECKER *goedecker,char *gfile_name)
/*==================================================================*/
  {/*begin routine*/
/*==================================================================*/
  int i,j;
  int nrad_mall,nrad_mall2;
  int lmax,lmax_mall;
  FILE *fp_in;
  char *cdum,*cdum1;

   cdum = (char *) calloc(10,sizeof(char));
   cdum1 = (char *) calloc(10,sizeof(char));

  if( (fp_in = fopen(gfile_name,"r")) == NULL){
     fprintf(stderr,"Error opening %s \n",gfile_name);
     exit(1);
  } 

  readtoendofline(fp_in); 
  fscanf(fp_in,"%lf",&(goedecker->z));readtoendofline(fp_in);
  fscanf(fp_in,"%lf",&(goedecker->zv));readtoendofline(fp_in);

  fscanf(fp_in,"%d",&(goedecker->lmax)); readtoendofline(fp_in);
  fscanf(fp_in,"%lf",&(goedecker->rloc)); readtoendofline(fp_in);
  fscanf(fp_in,"%d",&(goedecker->nc));

  /* Zero c1-c4 */
   goedecker->c = (double *) malloc(4*sizeof(double )) -1;
   goedecker->c[1] = 0.0;
   goedecker->c[2] = 0.0;
   goedecker->c[3] = 0.0;
   goedecker->c[4] = 0.0;

  switch(goedecker->nc) {
   case 1: fscanf(fp_in,"%lf",&(goedecker->c[1])); 
           readtoendofline(fp_in);
           break;
   case 2: fscanf(fp_in,"%lf %lf",&(goedecker->c[1]),&(goedecker->c[2]));
           readtoendofline(fp_in);
           break;
   case 3: fscanf(fp_in,"%lf %lf %lf",&(goedecker->c[1]),&(goedecker->c[2]),
                                  &(goedecker->c[3]));
           readtoendofline(fp_in);
           break;
   case 4: fscanf(fp_in,"%lf %lf %lf %lf",&(goedecker->c[1]),
                                          &(goedecker->c[2]),
                                          &(goedecker->c[3]),
                                          &(goedecker->c[4]));
           readtoendofline(fp_in);
           break;
  }/*end switch*/

/*---------------------------------------------------------------*/
/*---------------------------------------------------------------*/
/* Read in projection/wavefunction parameters                    */
/* Number of radial channels for given angular momentum channel  */

    lmax = goedecker->lmax;
    lmax_mall = ( lmax > 4 ? lmax : 4);
    

    goedecker->nrad_l = (int *) calloc(lmax_mall,sizeof(int));
    goedecker->rh     = (double *) calloc(lmax_mall,sizeof(double));
    goedecker->h      = (double **) calloc(lmax_mall,sizeof(double *));

    for(i=0; i<=lmax; i++){
      fscanf(fp_in,"%d",&(goedecker->nrad_l[i])); readtoendofline(fp_in);
      nrad_mall  = goedecker->nrad_l[i];
      nrad_mall2 = nrad_mall*nrad_mall; 

      if(nrad_mall > 0){
        fscanf(fp_in,"%lf",&(goedecker->rh[i])); readtoendofline(fp_in);

        goedecker->h[i] = (double *) calloc(nrad_mall2,sizeof(double));

        for(j=0; j< nrad_mall2; j+=(nrad_mall+1) ){
          fscanf(fp_in,"%lf",&(goedecker->h[i][j])); readtoendofline(fp_in);
        }/*endfor*/

      }/*endif*/

    }/*endfor*/

  
   fclose(fp_in);
   free(cdum);
   free(cdum1);

/*==================================================================*/
   }/*end routine*/
/*==================================================================*/



/*==================================================================*/
/*               Write out the parameters to the screen             */
/*==================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==================================================================*/
void print_goedecker(GOEDECKER *goedecker)
/*==================================================================*/
   {/*begin routine*/
/*==================================================================*/
   int nrad_l,nrad_l2;
   int i,j;

   printf("z  %lg \n",goedecker->z);
   printf("zv %lg \n",goedecker->zv);
   printf("lmax %d \n",goedecker->lmax);
   printf("rloc %lg \n",goedecker->rloc);
   printf("nc %d \n",goedecker->nc);
   printf("c1 %lg c2 %lg c3 %lg c4 %lg \n",goedecker->c[1],goedecker->c[2],
                                           goedecker->c[3],goedecker->c[4]);

   for(i=0; i< (goedecker->lmax+1); i++){
    printf("number of radial channels for channel %d  %d \n",i,goedecker->nrad_l[i]); 
     if(goedecker->nrad_l[i] > 0){
       printf("rvalue for channel %d  %lg \n",i,goedecker->rh[i]); 
       nrad_l  = goedecker->nrad_l[i];
       nrad_l2 = nrad_l*nrad_l;
       printf("nrad square %d \n",nrad_l2);
       for(j=0; j< nrad_l2; j++){
          printf("%lg ",goedecker->h[i][j]); 
       }
        printf("\n");
     }/*endif*/ 
   }/*endfor*/

/*==================================================================*/
   }/*end routine*/
/*==================================================================*/


/*==================================================================*/
/*               Skip information at the end of a line              */
/*==================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==================================================================*/
 void readtoendofline(FILE *fp)
/*==================================================================*/
  {/*begin routine*/
/*==================================================================*/
  int eol = (int) '\n';
  int ch;

   ch = eol + 1;
   while( ch != eol && ch != EOF ){ ch = fgetc(fp);}
   if( ch == EOF){
      fprintf(stderr,"@@@@@@@@@@@@@@@@@@@@ ERROR @@@@@@@@@@@@@@@@@@@@@ \n");
      fprintf(stderr,"   END OF FILE REACHED \n");
      fprintf(stderr,"@@@@@@@@@@@@@@@@@@@@ ERROR @@@@@@@@@@@@@@@@@@@@@ \n");
      exit(1);
   }/*endif*/

/*==================================================================*/
   }/*end routine*/
/*==================================================================*/


/*==================================================================*/
/*               Initialize arrays                                  */
/*==================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==================================================================*/
void pseudo_initial(PSEUDO *pseudo,int *nrad_l)
/*==================================================================*/
  {/*begin routine*/
/*==================================================================*/
   int i;
   int npts;

   pseudo->nr = 1000;
   pseudo->rmax = 15.0;
   pseudo->dr   = pseudo->rmax/((double)pseudo->nr);

   npts = pseudo->nr;

   pseudo->vloc  = (double *) calloc(npts,sizeof(double));
   pseudo->vlong = (double *) calloc(npts,sizeof(double));

 /* Allocate array for S projectors */
   if(nrad_l[0] > 0){
     pseudo->p0 = (double **) calloc(nrad_l[0],sizeof(double*));
      for(i=0; i< nrad_l[0]; i++){
       pseudo->p0[i] = (double *) calloc(npts,sizeof(double));
      }/*endfor*/
   }/*endif*/

 /* Allocate array for P projectors */
   if(nrad_l[1] > 0){
     pseudo->p1 = (double **) calloc(nrad_l[1],sizeof(double*));
      for(i=0; i< nrad_l[1]; i++){
       pseudo->p1[i] = (double *) calloc(npts,sizeof(double));
      }/*endfor*/
   }/*endif*/

 /* Allocate array for D projectors */
   if(nrad_l[2] > 0){
     pseudo->p2 = (double **) calloc(nrad_l[2],sizeof(double*));
      for(i=0; i< nrad_l[2]; i++){
       pseudo->p2[i] = (double *) calloc(npts,sizeof(double));
      }/*endfor*/
   }/*endif*/

 /* Allocate array for F projectors */
   if(nrad_l[3] > 0){
     pseudo->p3 = (double **) calloc(nrad_l[3],sizeof(double*));
      for(i=0; i< nrad_l[3]; i++){
       pseudo->p3[i] = (double *) calloc(npts,sizeof(double));
      }/*endfor*/
   }/*endif*/


/*==================================================================*/
   }/*end routine*/
/*==================================================================*/


/*==================================================================*/
/*               Local Pseudopotential                              */
/*==================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==================================================================*/
void make_goedecker_local(GOEDECKER *goedecker,PSEUDO *pseudo)
/*==================================================================*/
  {/*begin routine*/
/*==================================================================*/
  int i;

  double p=0.3614;
  double e1 = 0.2041422096422003, e2 = 0.1997535956961481;
  double e3 = 0.2213176596405576, e4 = 0.03360430734640255;
  double e5 = 0.4732592578721755, e6 =-0.509078520069735;
  double e7 = 0.6772631491947646, e8 =-0.369912979092217;
  double e9 = 0.06965131976970335;
  double r,r2,eee,tt;
  double alp,alpc,ralp;
  double gerfc,gerf;
  double gauss,ralpc,ralpc2,ralpc4,ralpc6;
  double zv   = goedecker->zv;
  double *c   = goedecker->c;
  int    nr    = pseudo->nr; 
  double rmax  = pseudo->rmax;
  double *vloc  = pseudo->vloc;
  double *vlong = pseudo->vlong;
  double dr    = pseudo->dr;

/*--------------------------------------------------------*/

    goedecker->erf_alp = 1.0/(goedecker->rloc*sqrt(2.0));

   alp = 1.0/(goedecker->rloc*sqrt(2.0));
   printf("zv %lg alp %lg\n",zv,alp);
    vlong[0] = 0.0;
    for(i=1; i< nr;i++){
      r      = ((double) i)*dr;
      ralp   = r * alp;
      eee    = exp(-ralp*ralp);
      tt     = 1.0/(1.0+p*ralp);
      gerfc  = ((((((((e9*tt+e8)*tt+e7)*tt+e6)*tt+e5)*tt
                      +e4)*tt+e3)*tt+e2)*tt+e1)*tt*eee;
      gerf     = 1.0 - gerfc;
      vlong[i] = -zv * gerf/r;

    }/*endfor*/

/* SHORT RANGE PIECE */
    alpc = 1.0/goedecker->rloc;

    for(i=0; i< nr;i++){
      r      = ((double) i)*dr;
      ralpc  = (r*alpc);
      ralpc2 = ralpc*ralpc;
      ralpc4 = ralpc2*ralpc2;
      ralpc6 = ralpc4*ralpc2;
      gauss = exp(-ralpc2/2.0);
      vloc[i] = (gauss*(c[1] + c[2]*ralpc2 + c[3]*ralpc4 + c[4]*ralpc6));
    }/*endfor*/

/*==================================================================*/
   }/*end routine*/
/*==================================================================*/


/*==================================================================*/
/*               Radial Channel Matrix                              */
/*==================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==================================================================*/
void make_off_diagonal_h(GOEDECKER *goedecker)
/*==================================================================*/
   {/*begin routine*/
/*==================================================================*/
  int *nrad_l = goedecker->nrad_l;
  double **h  = goedecker->h;

  double hs11,hs22,hs33,hs12,hs13,hs23;
  double hp11,hp22,hp33,hp12,hp13,hp23;
  double hd11,hd22,hd33,hd12,hd13,hd23;

 /* Create off diagonal elements for S Radial channel */ 
  switch(nrad_l[0]){
    case 0: break;
    case 1: break;
    case 2:
          hs11 = h[0][0];
          hs22 = h[0][3];
          hs12 = -(1.0/2.0)*sqrt(3.0/5.0)*hs22;

          h[0][1] = hs12;
          h[0][2] = hs12;
        break;
    case 3:
          hs11 = h[0][0];
          hs22 = h[0][4];
          hs33 = h[0][8];

          hs12 = -(1.0/2.0)*sqrt(3.0/5.0)*hs22;
          hs13 =  (1.0/2.0)*sqrt(5.0/21.0)*hs33;
          hs23 = -(1.0/2.0)*sqrt(100.0/63.0)*hs33;
         
          h[0][1] = hs12;
          h[0][3] = hs12;

          h[0][2] = hs13;
          h[0][6] = hs13;

          h[0][5] = hs23;
          h[0][7] = hs23;

        break;
  }/*end switch*/
  
 /* Create off diagonal elements for P Radial channel */ 
  switch(nrad_l[1]){
    case 0: break;
    case 1: break;
    case 2:
          hp11 = h[1][0];
          hp22 = h[1][3];
 
          hp12 = -(1.0/2.0)*sqrt(5.0/7.0)*hp22;

          h[1][1] = hp12;
          h[1][2] = hp12;
         break;
    case 3:
          hp11 = h[1][0];
          hp22 = h[1][4];
          hp33 = h[1][8];

          hp12 = -(1.0/2.0)*sqrt(5.0/7.0)*hp22;
          hp13 =  (1.0/6.0)*sqrt(35.0/11.0)*hp33; 
          hp23 = -(1.0/6.0)*(14.0/sqrt(11.0))*hp33;

          h[1][1] = hp12;
          h[1][3] = hp12;

          h[1][2] = hp13;
          h[1][6] = hp13;

          h[1][5] = hp23;
          h[1][7] = hp23;
         break;

  }/*end switch*/

 /* Create off diagonal elements for D Radial channel */ 
  switch(nrad_l[2]){
    case 0: break;
    case 1: break;
    case 2: 
       hd11 = h[2][0];       
       hd22 = h[2][3];       

       hd12 = -(1.0/2.0)*sqrt(7.0/9.0)*hd22;
    
       h[2][1] = hd12;
       h[2][2] = hd12;

       break;
    case 3: 
        hd11 = h[2][0];
        hd22 = h[2][4];
        hd33 = h[2][8];

        hd12 = -(1.0/2.0)*sqrt(7.0/9.0)*hd22;     
        hd13 =  (1.0/2.0)*sqrt(63.0/143.0)*hd33;
        hd23 = -(1.0/2.0)*(18.0/sqrt(143.0))*hd33;

        h[2][1] = hd12;
        h[2][3] = hd12;

        h[2][2] = hd13;
        h[2][6] = hd13;

        h[2][5] = hd23;
        h[2][7] = hd23;
       
       break;
  }/*end switch*/

/*==================================================================*/
  }/*end routine*/
/*==================================================================*/



/*==================================================================*/
/*               Radial Projectors                                  */
/*==================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==================================================================*/
void make_goedecker_projectors(GOEDECKER *goedecker,PSEUDO *pseudo)
/*==================================================================*/
   {/*begin routine*/
/*==================================================================*/
  int i,j;
  double eee,gauss;
  double rtw,value,rtgamma;
  double r,rs,rp,rd,ralp,r_pow,rs_pow,rp_pow,rd_pow;
  double gamma[20];

  double **p0  = pseudo->p0; /*S projectors*/
  double **p1  = pseudo->p1; /*P projectors*/
  double **p2  = pseudo->p2; /*D projectors*/
  double **p3  = pseudo->p3; /*F projectors*/
  int    nr    = pseudo->nr; /*number of points*/
  double dr    = pseudo->dr;  /* r spacing*/

  int  *nrad_l = goedecker->nrad_l;  /* number of radial channels for l */
  double *rh   = goedecker->rh;      /* radius rs (0) rp (1) rd (2) */

/*------------------------------------------------------*/
  rtw  = sqrt(2.0);
  make_gamma_function(&(gamma[0]));

/*--------------------------------*/
/* Make the S channel Projectors  */
/*--------------------------------*/
  for(i=0; i< nrad_l[0] ; i++){
      rs = rh[0];
    for(j=0; j< nr; j++){
      r = ((double) j)*dr; 

      ralp = r/rs;
      eee  = ralp*ralp/2.0;
     gauss = exp(-eee);

     value = ((double)(2*i));
     r_pow = pow(r,value); 
   
     value  = ((double)(4*i + 3))/2.0;
     rs_pow = pow(rs,value);

     rtgamma = sqrt(gamma[(2*i+1)]);

     p0[i][j] =  r*(rtw*r_pow*gauss)/(rs_pow*rtgamma); /*PINY wants r*Pr*/
      if(j==0){printf("p0 %lg \n",p0[i][j]);}
    }/*endfor*/
  }/*endfor*/

/*--------------------------------*/
/* Make the P channel Projectors  */
/*--------------------------------*/
  for(i=0; i< nrad_l[1] ; i++){
      rp = rh[1];
    for(j=0; j< nr; j++){
      r = ((double) j)*dr; 

      ralp = r/rp;
      eee  = ralp*ralp/2.0;
     gauss = exp(-eee);

     value = ((double)(2*i + 1));
     r_pow = pow(r,value); 
   
     value  = ((double)(4*i + 5))/2.0;
     rp_pow = pow(rp,value);

     rtgamma = sqrt(gamma[(2*i+2)]);

     p1[i][j] =  r*(rtw*r_pow*gauss)/(rp_pow*rtgamma); /* PINY wants r*Pr*/

    }/*endfor*/
  }/*endfor*/


/*--------------------------------*/
/* Make the D channel Projectors  */
/*--------------------------------*/
  for(i=0; i< nrad_l[2] ; i++){
      rd = rh[2];
    for(j=0; j< nr; j++){
      r = ((double) j)*dr; 

      ralp = r/rd;
      eee  = ralp*ralp/2.0;
     gauss = exp(-eee);

     value = ((double)(2*i + 2));
     r_pow = pow(r,value); 
   
     value  = ((double)(4*i + 7))/2.0;
     rd_pow = pow(rd,value);

     rtgamma = sqrt(gamma[(2*i+3)]);

     p2[i][j] =  r*(rtw*r_pow*gauss)/(rd_pow*rtgamma); /*PINY wants r*Pr */
    }/*endfor*/
  }/*endfor*/

/*--------------------------------*/
/* Make the F channel Projectors  */
/*--------------------------------*/
  for(i=0; i< nrad_l[3] ; i++){
      rd = rh[3];
    for(j=0; j< nr; j++){
      r = ((double) j)*dr; 

      ralp = r/rd;
      eee  = ralp*ralp/2.0;
     gauss = exp(-eee);

     value = ((double)(2*i + 3));
     r_pow = pow(r,value); 
   
     value  = ((double)(4*i + 9))/2.0;
     rd_pow = pow(rd,value);

     rtgamma = sqrt(gamma[(2*i+4)]);

     p3[i][j] =  r*(rtw*r_pow*gauss)/(rd_pow*rtgamma); /*PINY wants r*Pr */
    }/*endfor*/
  }/*endfor*/

/*==================================================================*/
   }/*end routine*/
/*==================================================================*/



/*==================================================================*/
/*               Gamma Function                                     */
/*==================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==================================================================*/
 void  make_gamma_function(double *gamma)
/*==================================================================*/
   {/*begin routine*/
/*==================================================================*/
   int i;
   double pi,rpi;

   pi  = M_PI;
   rpi = sqrt(pi);

/* gamma[i] = gamma( i + 1/2)  eg  gamma[1] = gamma(3/2) */

   gamma[0] = rpi;

   for(i=1; i<19; i++){
     gamma[i] = (((double) 2*i - 1)/2.0)*gamma[(i-1)];
   }/*endfor*/

/*==================================================================*/
   }/*end routine*/
/*==================================================================*/



/*==================================================================*/
/*               Write out the PSEUDO                              */
/*==================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==================================================================*/
void write_goedecker_file(GOEDECKER *goedecker,PSEUDO *pseudo)
/*==================================================================*/
  {/*begin routine*/
/*==================================================================*/
  FILE *fp_out; 
  int i,j,lval;
  int lmax    = goedecker->lmax;
  int *nrad_l = goedecker->nrad_l;
  double zv   = goedecker->zv;
  double alp  = goedecker->erf_alp;
  double **h  = goedecker->h;
  double **p0 = pseudo->p0;
  double **p1 = pseudo->p1;
  double **p2 = pseudo->p2;
  double **p3 = pseudo->p3;
  int nrad_l_sq;

  int nr       = pseudo->nr;
  double rmax  = pseudo->rmax;
  double *vloc = pseudo->vloc;

/*------------------------------------------*/
  if( (fp_out = fopen("PSEUDO.OUT","w")) == NULL){
     printf("error opening PSEUDO.OUT \n");
     exit(1);
  }/*endif*/

   fprintf(fp_out,"%d %lf %d \n",nr,rmax,lmax);
   fprintf(fp_out,"%d %d %d %d \n",nrad_l[0],nrad_l[1],nrad_l[2],nrad_l[3]);
   fprintf(fp_out,"%lg %lg  0.0  1.0 \n",zv,alp);
   fprintf(fp_out,"0.0  1.0 \n");  /* zpol  gamma */

/* Write out Radial Channel Matrices */
/*   fprintf(fp_out,"RADIAL CHANNEL MATRICES \n"); */
   for(i=0; i< 4; i++){
     nrad_l_sq = nrad_l[i]*nrad_l[i];
     for(j=0; j< nrad_l_sq; j++){ 
       fprintf(fp_out,"%12.8lg  ",h[i][j]);
       if((j+1)%nrad_l[i] == 0){fprintf(fp_out,"\n");}
     }/*endfor*/
   }/*endfor*/

/* Write out Projectors for Angular Momentum Channels */
/* S CHANNELS */
/*   fprintf(fp_out,"S PROJECTORS \n"); */
   for(i=0; i< nrad_l[0]; i++){
     for(j=0; j< nr; j++){
        fprintf(fp_out,"%12.8lg %12.8lg\n",p0[i][j],p0[i][j]);
     } /*endfor*/
   }/*endfor*/

/* P CHANNELS */
/*   fprintf(fp_out,"P PROJECTORS \n"); */
   for(i=0; i< nrad_l[1]; i++){
     for(j=0; j< nr; j++){
        fprintf(fp_out,"%12.8lg %12.8lg \n",p1[i][j],p1[i][j]);
     }/*endfor*/
   }/*endfor*/
/* D CHANNELS */
/*   fprintf(fp_out,"D PROJECTORS \n"); */
   for(i=0; i< nrad_l[2]; i++){
     for(j=0; j< nr; j++){
        fprintf(fp_out,"%12.8lg %12.8lg \n",p2[i][j],p2[i][j]);
     }/*endfor*/ 
   }/*endfor*/
/* F CHANNELS */
/*   fprintf(fp_out,"F PROJECTORS \n"); */
   for(i=0; i< nrad_l[3]; i++){
     for(j=0; j< nr; j++){
        fprintf(fp_out,"%12.8lg %12.8lg \n",p3[i][j],p3[i][j]);
     }/*endfor*/ 
   }/*endfor*/
  

/* Write out Short range local potential */

/*   fprintf(fp_out,"VLOC \n"); */
   for(j=0; j< nr; j++){
     fprintf(fp_out,"%12.8lg %12.8lg \n",vloc[j],vloc[j]);
   }/*endfor*/

  fclose(fp_out);

/*==================================================================*/
   }/*end routine*/
/*==================================================================*/


/*==================================================================*/
/*               Check PINY versus analytical results               */
/*==================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==================================================================*/
void check_goedecker_projectors(GOEDECKER *goedecker,PSEUDO *pseudo)
/*==================================================================*/
   {/*begin routine*/
/*==================================================================*/
  double **p0 = pseudo->p0;
  double **p1 = pseudo->p1;
  double **p2 = pseudo->p2;
  double **p3 = pseudo->p3;
  double   dr = pseudo->dr;
  double   nr = pseudo->nr;

  int *nrad_l = goedecker->nrad_l;
  int i,j;
  double total;

/* CHECK S PROJECTORS NORMALIZATION */
  for(i=0; i< nrad_l[0]; i++){
     total = 0.0;
    for(j=0; j< nr; j++){
      total += (p0[i][j]*p0[i][j]*dr);
    }/*endfor*/
     printf("S radidal channel %d  normalized %lg \n",i,total);
  }/*endfor*/

/* CHECK P PROJECTORS NORMALIZATION */
  for(i=0; i< nrad_l[1]; i++){
     total = 0.0;
    for(j=0; j< nr; j++){
      total += (p1[i][j]*p1[i][j]*dr);
    }/*endfor*/
     printf("P radidal channel %d  normalized %lg \n",i,total);
  }/*endfor*/

/* CHECK D PROJECTORS NORMALIZATION */
  for(i=0; i< nrad_l[2]; i++){
    total = 0.0;
    for(j=0; j< nr; j++){
      total += (p2[i][j]*p2[i][j]*dr);
    }/*endfor*/
    printf("D radidal channel %d  normalized %lg \n",i,total);
  }/*endfor*/

/* CHECK D PROJECTORS NORMALIZATION */
  for(i=0; i< nrad_l[3]; i++){
     total = 0.0;
    for(j=0; j< nr; j++){
      total += (p3[i][j]*p3[i][j]*dr);
    }/*endfor*/
     printf("F radidal channel %d  normalized %lg \n",i,total);
  }/*endfor*/

/*==================================================================*/
    }/*end routine*/
/*==================================================================*/



/*==================================================================*/
/*               Check PINY versus analytical results               */
/*==================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==================================================================*/
void check_goedecker_gspace(GOEDECKER *goedecker)
/*==================================================================*/
   {/*begin routine*/
/*==================================================================*/
  FILE *fp_out;
  int *nrad_l = goedecker->nrad_l;
  double  *rh = goedecker->rh;
  double rs,rp,rd,rf;
  double pi_fact;
  double pi;
  double eee,gauss,pre_fact,poly;
  double g,g2,g3,g4;
  double grs,grs2,grs4,rs3;
  double grp,grp2,grp4,rp5;
  double grd,grd2,rd7;
  double grf,grf2,rf9;
  int i,ir,ig;

  double **p0_g,**p1_g,**p2_g,**p3_g;

  double dg_spl = 0.00160278;
  double gmin_spl = 0.392699 ;

/*------------------------------------------------------*/
   if( nrad_l[0] > 0){
    p0_g = (double **) calloc(nrad_l[0],sizeof(double *));
     for(i=0; i< nrad_l[0]; i++){
       p0_g[i] = (double *) calloc(10,sizeof(double ));
     }/*endfor*/
   }/*endif*/

   if( nrad_l[1] > 0){
    p1_g = (double **) calloc(nrad_l[1],sizeof(double *));
     for(i=0; i< nrad_l[1]; i++){
       p1_g[i] = (double *) calloc(10,sizeof(double ));
     }/*endfor*/
   }/*endif*/

   if( nrad_l[2] > 0){
    p2_g = (double **) calloc(nrad_l[2],sizeof(double *));
     for(i=0; i< nrad_l[2]; i++){
       p2_g[i] = (double *) calloc(10,sizeof(double ));
     }/*endfor*/
   }/*endif*/

   if( nrad_l[3] > 0){
    p3_g = (double **) calloc(nrad_l[3],sizeof(double *));
     for(i=0; i< nrad_l[3]; i++){
       p3_g[i] = (double *) calloc(10,sizeof(double ));
     }/*endfor*/
   }/*endif*/

   pi      = M_PI;
   pi_fact = pow(pi,(5.0/4.0));

   if( (fp_out = fopen("proj_gspace","w")) == NULL){
     printf("Error opening proj_gspace \n");
     exit(1);
   }/*endif*/

/* S Channel projectors in G space */
   fprintf(fp_out,"S projectors \n");
  for(ir=1; ir <= nrad_l[0]; ir++){
     rs    = rh[0];   
     rs3   = rs*rs*rs;
    for(ig=0; ig < 10; ig++){
     g     = dg_spl*((double) (ig)) + gmin_spl;
     grs   = g*rs;
     grs2  = grs*grs;
     grs4  = grs2*grs2;
     eee   = -grs2/2.0; 
     gauss = exp(eee); 
    if(ir == 1){
      pre_fact = 4.0*sqrt((2.0*rs3))*pi_fact;
      p0_g[(ir-1)][ig] = pre_fact*gauss;
      fprintf(fp_out,"proj 1 %lg  \n",p0_g[(ir-1)][ig]);
    }/*endif 1 channel*/
    if(ir == 2){
      pre_fact = 8.0*sqrt((2.0/15.0)*rs3)*pi_fact;
      poly  = (3.0 - grs2);
      p0_g[(ir-1)][ig] = pre_fact*poly*gauss;
      fprintf(fp_out,"proj 2 %lg  \n",p0_g[(ir-1)][ig]);
    }/*endif 2 channel */
    if(ir == 3){
      pre_fact = (16.0/3.0)*sqrt((2.0/105.0)*rs3)*pi_fact; 
      poly  = (15.0 - 10.0*grs2 + grs4); 
      p0_g[(ir-1)][ig] = pre_fact*poly*gauss;
      fprintf(fp_out,"proj 3 %lg  \n",p0_g[(ir-1)][ig]);
    }/*endif 3 channel */
   }/*endfor g vector*/
  }/*endfor radial channels*/

/* P Channel projectors in G space */
   fprintf(fp_out," P PROJECTORS \n");
  for(ir=1; ir <= nrad_l[1]; ir++){
     rp    = rh[1];   
     rp5   = pow(rp,5.0);
   for(ig=0; ig < 10; ig++){
     g     = dg_spl*((double) (ig)) + gmin_spl;
     grp   = g*rp;
     grp2  = grp*grp;
     grp4  = grp2*grp2; 
     eee   =  -grp2/2.0;
     gauss = exp(eee);
    if(ir == 1){
     pre_fact = 8.0*sqrt(rp5/3.0)*pi_fact;
     p1_g[(ir-1)][ig] = pre_fact*g*gauss;
      fprintf(fp_out,"proj 1 %lg  \n",p1_g[(ir-1)][ig]);
    }/*endif 1 channel */
    if(ir == 2){
     pre_fact = 16.0*sqrt(rp5/105.0)*pi_fact; 
     poly  = (5.0 - grp2);
     p1_g[(ir-1)][ig] = pre_fact*poly*g*gauss;
      fprintf(fp_out,"proj 2 %lg  \n",p1_g[(ir-1)][ig]);
    }/*endif 1 channel */
    if(ir == 3){
     pre_fact = (32.0/3.0)*sqrt(rp5/1155.0)*pi_fact;
     poly  = (35.0 - 14.0*grp2 + grp4);
     p1_g[(ir-1)][ig] = pre_fact*poly*g*gauss;
      fprintf(fp_out,"proj 3 %lg  \n",p1_g[(ir-1)][ig]);
    }/*endif 3 channel */

   }/*endfor g vector */
  }/*endfor radial channels */


/* D Channel projectors in G space */
   fprintf(fp_out," D PROJECTORS \n");
  for(ir=1; ir <= nrad_l[2]; ir++){
     rd    = rh[2];   
     rd7   = pow(rd,7.0);
    for(ig=0; ig < 10; ig++){
     g     = dg_spl*((double) (ig)) + gmin_spl;
     g2    = g*g;
     grd   = g*rd;
     grd2  = grd*grd;
     eee   =  -grd2/2.0;
     gauss = exp(eee);
    if(ir == 1){
      pre_fact = 8.0*sqrt(2.0*rd7/15.0)*pi_fact;
      p2_g[(ir-1)][ig] = pre_fact*g2*gauss; 
      fprintf(fp_out,"proj 1 %lg  \n",p2_g[(ir-1)][ig]);
    }/*endif*/
    if(ir == 2){
      pre_fact = (16.0/3.0)*sqrt(2.0*rd7/105.0)*pi_fact; 
      poly  = (7.0 - grd2); 
      p2_g[(ir-1)][ig] = pre_fact*poly*g2*gauss; 
      fprintf(fp_out,"proj 2 %lg  \n",p2_g[(ir-1)][ig]);
    }/*endif*/

    }/*endfor g vector */
  }/*endfor radial channels */

/* F Channel projectors in G space */
   fprintf(fp_out," F PROJECTORS \n");
  for(ir=1; ir <= nrad_l[3]; ir++){
     rf    = rh[3];   
     rf9   = pow(rf,9.0);
    for(ig=0; ig < 10; ig++){
      g    = dg_spl*((double) (ig)) + gmin_spl;
      g3   = g*g*g; 
      grf  = g*rf;
      grf2 = grf*grf;
      eee  = -grf2/2.0;
     gauss = exp(eee);
     if(ir == 1){
       pre_fact = 16.0*sqrt(rf9/105.0)*pi_fact;
       p3_g[(ir-1)][ig] = pre_fact*g3*gauss;
      fprintf(fp_out,"proj 1 %lg  \n",p3_g[(ir-1)][ig]);
     }/*endif*/

    }/*endfor g vector */
  }/*endfor radial channel*/

/* free memory */
   if( nrad_l[0] > 0){ free(p0_g); }
   if( nrad_l[1] > 0){ free(p1_g); }
   if( nrad_l[2] > 0){ free(p2_g); }
   if( nrad_l[3] > 0){ free(p3_g); }

/*==================================================================*/
   }/*end routine*/
/*==================================================================*/



/*==================================================================*/
/*               Check PINY versus analytical results               */
/*==================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==================================================================*/
void check_goedecker_local(GOEDECKER *goedecker)
/*==================================================================*/
  {/*begin routine*/
/*==================================================================*/
  FILE *fp_out;
  int ig,ic;
  int      nc  = goedecker->nc;
  double   *c  = goedecker->c;
  double rloc  = goedecker->rloc; 
  double   zv  = goedecker->zv;
  double erf_alp = goedecker->erf_alp;
  double dg_spl = 0.001602784;
  double gmin_spl = 0.39269908 ;
  double g,g2,gr,gr2,gr4,gr6,rloc3;
  double vloc;
  double pi,pi_fact,eee,gauss;
  double poly; 

   fp_out = fopen("local_gspace","w");

     pi      = M_PI;
     pi_fact = sqrt(8.0*pi*pi*pi); 
     rloc3   = rloc*rloc*rloc;

  for(ig=0; ig < 10; ig++){
     vloc = 0.0;
       g  = dg_spl*((double) (ig)) + gmin_spl;
       g2 = g*g;
      gr  = g*rloc;
     gr2  = gr*gr;
     gr4  = gr2*gr2;
     gr6  = gr2*gr2*gr2;
     eee  = -gr2/2.0;
    gauss = exp(eee);
     vloc += (-4.0*pi*zv*gauss/g2);  
    for(ic = 1; ic <= nc; ic++){
      if(ic == 1){
        poly = c[ic];
        vloc += pi_fact*rloc3*gauss*poly; 
      }/*endif*/
      if(ic == 2){
        poly = c[ic]*(3.0 - gr2);
        vloc += pi_fact*rloc3*gauss*poly; 
      }/*endif*/
      if(ic == 3){
        poly = c[ic]*(15.0 - 10.0*gr2 + gr4);
        vloc += pi_fact*rloc3*gauss*poly; 
      }/*endif*/
      if(ic == 4){
        poly = c[ic]*(105.0 - 105.0*gr2 + 21.0*gr4 - gr6);
        vloc += pi_fact*rloc3*gauss*poly; 
      }/*endif*/
    }/*endfor constants*/
     fprintf(fp_out,"vloc %12.8lg \n",vloc);
  }/*endfor g vector */

  fclose(fp_out);

/*==================================================================*/
   }/*end routine*/
/*==================================================================*/



/*==================================================================*/
/*               Check PINY versus analytical results               */
/*==================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==================================================================*/
void check_nonlocal_channels(GOEDECKER *goedecker,PSEUDO *pseudo)
/*==================================================================*/
   {/*begin routine*/
/*==================================================================*/
  int iang,irad1,irad2,ih,ir;
  double value1,value2,value;
  double value_chan,value_tot,val_loc,val_loc_tot;
  int lmax       = goedecker->lmax;
  int *nrad_l    = goedecker->nrad_l;
  double **h     = goedecker->h;   /* (l,n) */
  int nr         = pseudo->nr;
  double dr      = pseudo->dr;
  double **p0    = pseudo->p0;
  double **p1    = pseudo->p1;
  double **p2    = pseudo->p2;
  double **p3    = pseudo->p3;
  double *vloc   = pseudo->vloc;
  double *vlong  = pseudo->vlong;
  double *psi_s,*psi_p,*psi_d,*psi_f; 

   value_tot = 0.0;
   val_loc_tot = 0.0;

/*------------*/
/* S channel  */
/*------------*/

 if(lmax >=0){
  printf("S CHANNEL \n");
      ih = 0;
      value_chan = 0.0;
   for(irad1 = 0; irad1 < nrad_l[0]; irad1++){
     psi_s = p0[0];
    for(irad2 = 0; irad2 < nrad_l[0]; irad2++){
     value1 = 0.0;
     value2 = 0.0;

      for(ir=0; ir < nr; ir++){
        value1 += ((psi_s[ir]*p0[irad1][ir])*dr);
        value2 += ((psi_s[ir]*p0[irad2][ir])*dr); 
      }/*endfor ir*/

      value = value1*value2;
      printf("%d %d value %lg value*h %lg \n",
             irad1+1,irad2+1,value,value*h[0][ih]);
      value = value1*value2*h[0][ih];
      value_chan += value;
     ih++;
    }/*endfor*/
   }/*endfor*/
   val_loc = 0.0;
   for(ir=0; ir < nr; ir++){
     val_loc += ((psi_s[ir]*psi_s[ir]*(vloc[ir]+vlong[ir]))*dr);
   }/*endfor ir*/
   printf("S channel total %lg local %lg\n",value_chan,val_loc); 
   val_loc_tot += val_loc;
   value_tot   += value_chan;
  }
/*------------*/
/* P channel  */
/*------------*/
if(lmax>=1){
  printf("P CHANNEL \n");
      ih = 0;
      value_chan = 0.0;
   for(irad1 = 0; irad1 < nrad_l[1]; irad1++){
     psi_p = p1[0];
    for(irad2 = 0; irad2 < nrad_l[1]; irad2++){
     value1 = 0.0;
     value2 = 0.0;
      for(ir=0; ir < nr; ir++){
        value1 += ((psi_p[ir]*p1[irad1][ir])*dr);
        value2 += ((psi_p[ir]*p1[irad2][ir])*dr); 
      }/*endfor ir*/
      value = value1*value2;
    printf("%d %d value %lg value*h %lg \n",
              irad1+1,irad2+1,value,value*h[1][ih]);
      value = value1*value2*h[1][ih];
      value_chan += value;
      ih++;
    }/*endfor*/
   }/*endfor*/
   val_loc = 0.0;
   for(ir=0; ir < nr; ir++){
     val_loc += ((psi_p[ir]*psi_p[ir]*(vloc[ir]+vlong[ir]))*dr);
   }/*endfor ir*/

   value_chan *= 3.0;
   val_loc    *= 3.0;
   printf("P channel total %lg local %lg\n",value_chan,val_loc); 
   val_loc_tot += val_loc;
   value_tot   += value_chan;
}
/*------------*/
/* D channel  */
/*------------*/
if(lmax>=2){
  printf("D CHANNEL \n");
   ih = 0;
   value_chan = 0.0;
   for(irad1 = 0; irad1 < nrad_l[2]; irad1++){
     psi_d = p2[0];
    for(irad2 = 0; irad2 < nrad_l[2]; irad2++){
     value1 = 0.0;
     value2 = 0.0;
      for(ir=0; ir < nr; ir++){
        value1 += ((psi_d[ir]*p2[irad1][ir])*dr); 
        value2 += ((psi_d[ir]*p2[irad2][ir])*dr); 
      }/*endfor ir*/
      value = value1*value2;
      printf("%d %d value %lg value*h %lg \n",
           irad1+1,irad2+1,value,value*h[2][ih]);
      value = value1*value2*h[2][ih];
      value_chan += value;
      ih++;
    }/*endfor*/
   }/*endfor*/
   val_loc = 0.0;
   for(ir=0; ir < nr; ir++){
     val_loc += ((psi_d[ir]*psi_d[ir]*(vloc[ir]+vlong[ir]))*dr);
   }/*endfor ir*/

   val_loc    *= 5.0;
   value_chan *= 5.0;
   printf("D channel total %lg local %lg \n",value_chan,val_loc); 

   value_tot   += value_chan;
   val_loc_tot += val_loc;
}
/*------------*/
/* F channel  */
/*------------*/
if(lmax>=3){
  printf("F CHANNEL \n");
      value_chan = 0.0;
      ih = 0;
   for(irad1 = 0; irad1 < nrad_l[3]; irad1++){
     psi_f = p3[0];
    for(irad2 = 0; irad2 < nrad_l[3]; irad2++){
     value1 = 0.0;
     value2 = 0.0;
      for(ir=0; ir < nr; ir++){
        value1 += ((psi_f[ir]*p3[irad1][ir])*dr); 
        value2 += ((psi_f[ir]*p3[irad2][ir])*dr); 
      }/*endfor ir*/
      value = value1*value2;
      printf("%d %d value %lg value*h %lg \n",
              irad1+1,irad2+1,value,value*h[3][ih]);
      value = value1*value2*h[3][ih];
      value_chan += value;
      ih++;
    }/*endfor*/
   }/*endfor*/
   val_loc = 0.0;
   for(ir=0; ir < nr; ir++){
     val_loc += ((psi_f[ir]*psi_f[ir]*(vloc[ir]+vlong[ir]))*dr);
   }/*endfor ir*/

   val_loc    *= 7.0;
   value_chan *= 7.0;
   printf("F channel total %lg local %lg\n",value_chan,val_loc); 
   value_tot   += value_chan;
   val_loc_tot += val_loc;
}
   printf("Total non-local %lg Total local %lg\n",value_tot,val_loc_tot);

/*==================================================================*/
   }/*end routine*/
/*==================================================================*/
