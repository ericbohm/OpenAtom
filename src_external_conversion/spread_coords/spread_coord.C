//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//
//  This program takes classical input and make the PIMD ring polymers
//  To run this program : executable inputfile
//
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Standard include files

#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cmath>

//==========================================================================
// PINY constants, inline functions and typedefs

#define BOLTZ 315777.0
#define BOHR  0.529177

#define MAXWORD   80
#define MAXLINE  100
typedef char NAME[MAXWORD];  

#define PRINTF printf
#define FFLUSH fflush
#define EXIT(N) {exit(N);}
#define PRINT_LINE_STAR {PRINTF("==============================================================================\n");}
#define PRINT_LINE_DASH {PRINTF("------------------------------------------------------------------------------\n");}

#define MAX(A,B) (((A)>(B))?(A):(B))
#define MIN(A,B) (((A)<(B))?(A):(B))

//==========================================================================
// Function prototypes

int main (int , char *[]);
void read_natm_tot(int , int *, char *);
void read_coord(int , int , double *, double *, double *,NAME *, double *, char *);
void spread_coord(int ,double *,double *,double *,double , double ,
    double , double *, double *, long *);
void stage_1_part(double *,double *,double ,double,double ,int ,long *,int *);
void gaussran(int , long *, double *);
double altRandom(long *);
void readtoendofline(FILE *);
double dsum1(int ,double *);

//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//  Main Program : Controller
//==========================================================================
int main (int argc, char *argv[]){
  //==========================================================================
  // Local variables

  int i,j;
  int natm_tot;
  int pi_beads;
  int istart;
  double rcut, T_ext;
  double *x,*y,*z,*mass;
  double *xp, *yp, *zp, *pos, *gauss;
  double hmat[10];

  long seed;
  double dseed;

  NAME *atm_typ;
  char start_typ[MAXWORD];
  char mass_file[MAXWORD];
  char fnameIn[MAXLINE];
  char dnameOut[MAXLINE];
  char fnameOut[MAXLINE];
  char fnameNow[MAXLINE];
  FILE *fp;

  //==========================================================================
  // Tell everyone what you are doing

  PRINTF("\n");
  PRINT_LINE_STAR
    PRINTF("Creating PIMD coordinate input files from a single classical starting pt\n");
  PRINT_LINE_DASH

    //=========================================================================
    //             Check for input file                                 

    if(argc < 2) {
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("No input file specified\n");
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }/*endif*/

  //==========================================================================
  // Get some input from the user

  PRINTF("\nReading input parameters from %s\n",argv[1]);
  fp = fopen(argv[1],"r");
  fscanf(fp,"%d %lf %s",&pi_beads,&dseed,start_typ); readtoendofline(fp);
  fscanf(fp,"%lf %lf",&rcut,&T_ext); readtoendofline(fp);
  fscanf(fp,"%s %s %s",fnameIn,dnameOut,fnameOut); readtoendofline(fp);
  fscanf(fp,"%s",mass_file); 
  fclose(fp);
  PRINTF("Finished reading input parameters from %s\n\n",argv[1]);

  PRINTF("Spreading out to %d beads at T=%gK ",pi_beads,T_ext);
  PRINTF("keeping the ring polymer diameters < %gA\n\n",rcut);

  seed = (long) dseed;
  rcut /= BOHR;

  istart = 0;
  if(strcasecmp(start_typ,"initial")    == 0){istart = 1;}
  if(strcasecmp(start_typ,"restart_pos")== 0){istart = 2;}

  if(istart==0){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Error in start_typ, %s, only initial and restart_pos allowed.\n",start_typ);
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }//endif

  if(pi_beads<=1){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Number of beads, %d, must be greater than 1\n",pi_beads);
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }//endif

  //==========================================================================
  // Get the number of atoms and malloc

  read_natm_tot(istart,&natm_tot,fnameIn);

  x       = (double *)malloc(natm_tot*sizeof(double))-1;
  y       = (double *)malloc(natm_tot*sizeof(double))-1;
  z       = (double *)malloc(natm_tot*sizeof(double))-1;
  mass    = (double *)malloc(natm_tot*sizeof(double))-1;
  atm_typ = (NAME *)  malloc(natm_tot*sizeof(NAME))  -1;

  xp      = (double *)malloc(pi_beads*sizeof(double))-1;
  yp      = (double *)malloc(pi_beads*sizeof(double))-1;
  zp      = (double *)malloc(pi_beads*sizeof(double))-1;
  pos     = (double *)malloc(pi_beads*sizeof(double))-1;
  gauss   = (double *)malloc(pi_beads*sizeof(double))-1;

  //==========================================================================
  // Read in the coordinates in specified format

  read_coord(istart,natm_tot,x,y,z,atm_typ,hmat,fnameIn);

  //==========================================================================
  // Set up the masses by reading from file

  // PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
  // PRINTF("The atom masses are hard coded for water on lines 159-163\n");
  // PRINTF("The program needs an upgrade for general use\n");
  // PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n\n");    
  // for(i=1;i<=natm_tot;i+=3){
  //   mass[(i  )] = 1822.0*16.0; // oxygen
  //   mass[(i+1)] = 1836.0;      // hydrogen
  //   mass[(i+2)] = 1836.0;      // hydrogen
  // }//endfor

  fp = fopen(mass_file,"r");
  fscanf(fp,"%d",&i);  readtoendofline(fp);
  if (i!=natm_tot) {
    PRINTF("\nNaughty!  Your number of atoms in the mass file does not match the natm_tot\n\n");
    PRINTF("From masses file %d  while natm_tot = %d\n",i,natm_tot);
    EXIT(1);
  }
  for(i=1;i<=natm_tot;i++){
    fscanf(fp,"%lg",&mass[i]); readtoendofline(fp);
    mass[i] *= 1822.0;  // covert from amu to atomic units
  }//endfor
  fclose(fp);
  
  //==========================================================================
  // Write the coordinate files in initial restart mode

  PRINTF("Writing the PIMD coord input files to %s/Bead.I_Temper.0/%s : I=0...%d\n",
      dnameOut,fnameOut,pi_beads-1);

  for(j=1;j<=pi_beads;j++){
    sprintf (fnameNow, "%s/Bead.%d_Temper.0/%s",dnameOut,(j-1),fnameOut);
    fp = fopen(fnameNow,"w");
    fprintf(fp,"%d 1 1\n",natm_tot);
    fclose(fp);
  }//endfor

  for(i=1;i<=natm_tot;i++){
    xp[1] = x[i]; yp[1] = y[i]; zp[1] = z[i];
    spread_coord(pi_beads,xp,yp,zp,mass[i],T_ext,rcut,pos,gauss,&seed);
    for(j=1;j<=pi_beads;j++){
      sprintf (fnameNow, "%s/Bead.%d_Temper.0/%s",dnameOut,(j-1),fnameOut);
      fp = fopen(fnameNow,"a");
      fprintf(fp,"%.10g %.10g %.10g\n",xp[j]*BOHR,yp[j]*BOHR,zp[j]*BOHR);
      fclose(fp);
    }//endfor
  }//endfor

  for(j=1;j<=pi_beads;j++){
    sprintf (fnameNow, "%s/Bead.%d_Temper.0/%s",dnameOut,(j-1),fnameOut);
    fp = fopen(fnameNow,"a");
    for(i=0;i<3;i++){
      fprintf(fp,"%.10g %.10g %.10g\n",
          hmat[(i+1)]*BOHR, hmat[(i+4)]*BOHR,hmat[(i+7)]*BOHR);
    }//endfor
    fclose(fp);
  }//endfor

  PRINTF("Finished writing the PIMD coord input files to %s/Bead.I_Temper.0/%s : I=0...%d\n\n",
      dnameOut,fnameOut,pi_beads-1);

  //==========================================================================
  // Tell everyone you are done

  PRINT_LINE_DASH
    PRINTF("Finished creating PIMD coordinate input files from a single classical starting pt\n");
  PRINT_LINE_STAR
    PRINTF("\n");

  return 0;

  //--------------------------------------------------------------------------
}//end routine
//==========================================================================



//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Manage the sampling of the free particle propagator in 3D
//==========================================================================
void spread_coord(int pi_beads,double *x,double *y,double *z,double mass, double T_ext,
    double rcut, double *pos, double *gauss, long *seed)
  //======================================================================
{//begin routine 
  //======================================================================
  //               Local variable declarations                            

  int ip,ierr,ierr_now;

  //========================================================================
  //  II)Spread out classical coords to path integral coords:
  //     Place centroid of chain at classical position
  //     pos has centroid at 0

  ierr = 0;

  stage_1_part(pos,gauss,mass,T_ext,rcut,pi_beads,seed,&ierr_now);
  ierr += ierr_now;
  for(ip=2;ip<=pi_beads;ip++){x[ip] = pos[ip] + x[1];}
  x[1] = pos[1] + x[1];

  stage_1_part(pos,gauss,mass,T_ext,rcut,pi_beads,seed,&ierr_now);
  ierr += ierr_now;
  for(ip=2;ip<=pi_beads;ip++){y[ip] = pos[ip] + y[1];}
  y[1] = pos[1] + y[1];

  stage_1_part(pos,gauss,mass,T_ext,rcut,pi_beads,seed,&ierr_now);
  ierr += ierr_now;
  for(ip=2;ip<=pi_beads;ip++){z[ip] = pos[ip] + z[1];}
  z[1] = pos[1] + z[1];

  //========================================================================
  // IV) Error mesage 

  if(ierr>0){
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    PRINTF("%d coordinates outside spread cutoff\n",ierr);
    PRINTF("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }//endif

  //========================================================================
}// end routine 
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//   Directly sample the free particle propagator (1D)
//==========================================================================
void stage_1_part(double *x,double *gaussx,double mass,double T_ext,
    double rcut,int pi_beads,long *seed,int *ierr_ret)
  //======================================================================
{  //begin routine 
  //======================================================================
  //               Local variable declarations                            

  int i,nmove,not_done,ierr,done,iii;

  double tau;
  double endlx,endhx,prer,prerp,sig2,sig,ameanx,dist,xcm;
  double dpi_beads = (double)pi_beads;

  //======================================================================

  if(pi_beads==1){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Number of beads, %d, must be greater than 1\n",pi_beads);
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }//endif

  //======================================================================

  nmove = pi_beads-1;
  tau   = BOLTZ/(T_ext*dpi_beads);

  //======================================================================
  //    Spread out the particle into a cyclic chain of diameter < rcut

  //----------------------------------------------------------
  // Put 1st bead is at the origin : Get ready to move bead 2

  x[1]  = 0.0;  endlx = 0.0; endhx = 0.0;
  prer  = (double)(nmove)/(double)(nmove+1);
  prerp = 1./(double)(nmove+1);
  sig2  = prer*tau/mass;
  sig   = sqrt(sig2);
  ameanx = prer*endlx + prerp*endhx;

  //--------------------------
  // Move beads 2 to pi_beads

  ierr = 0;  iii = 1;
  gaussran(nmove,seed,gaussx);
  for(i=1;i<=nmove;i++){
    not_done=0;done=0;
    while(not_done<=5&&done==0){ // try to keep diameter small
      x[(i+1)] =  gaussx[iii]*sig + ameanx;
      iii++;
      dist = fabs(x[i+1]);
      if(iii>nmove){gaussran(nmove,seed,gaussx);iii=1;}
      if(dist<=rcut){done++;}else{not_done++;}
    }//endwhile
    if(not_done>5){ierr++;}
    if(i<nmove){// get ready to move the next bead
      prer  = (double)(nmove-i)/(double)(nmove+1-i);
      prerp = 1./(double)(nmove+1-i);
      sig2  = prer*tau/mass;
      sig   = sqrt(sig2);
      ameanx = prer*x[i+1] + prerp*endhx;
    }//endif
  }//endfor

  //======================================================================
  // Place the centroid and not the 1st bead at the origin

  xcm = dsum1(pi_beads,x)/(double)(pi_beads);
  for(i=1;i<=pi_beads;i++){x[i] = x[i] - xcm;}

  //======================================================================
  // Return the error flag

  *ierr_ret = ierr;

  //========================================================================
}// end routine 
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//   Little utility routine
//==========================================================================
double dsum1(int n,double *x){
  double add = 0.0;
  for(int i=1;i<=n;i++){add += x[i];}
  return add;
}//end routine
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//  Read in the total number of atoms 
//==========================================================================
void read_natm_tot(int istart, int *natm_tot, char *dnamei)
  //==========================================================================
{// begin routine
  //==========================================================================

  int istart_now;
  int natm_tot_now;
  int pi_beads_now;
  int itime_dump;

  char restart_type_now[MAXWORD];
  char restart_type_spec[MAXWORD];

  FILE *fp_dnamei;

  //==========================================================================

  PRINTF("Reading the number of atoms from %s\n",dnamei);
  fp_dnamei = fopen((const char *)dnamei,"r");

  //========================================================================
  //  I)Read in header                                                    

  //-------------------
  // A) Type 1 start : 
  if(istart==1){
    if(fscanf(fp_dnamei,"%d %d %d",&natm_tot_now,&istart_now,
          &pi_beads_now)!=3){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("Error reading start type and number of atoms \n");
      PRINTF("in file %s\n",dnamei);
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }//endif
    readtoendofline(fp_dnamei);
  }//endif
  //---------------------
  // B) Type 2,3,4 start 
  if(istart>1){
    if(fgetc(fp_dnamei)==EOF){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("Error reading header information \n");
      PRINTF("in file %s\n",dnamei);
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }//endif
    readtoendofline(fp_dnamei);

    if(fscanf(fp_dnamei,"%d %s %d %d",&natm_tot_now,restart_type_now,
          &itime_dump,&pi_beads_now)!=4){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("Error reading start type and number of atoms \n");
      PRINTF("in file %s\n",dnamei);
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    }//endif
    readtoendofline(fp_dnamei);

    istart_now = 0;
    if(strcasecmp(restart_type_now,"initial") == 0){        istart_now = 1;}
    if(strcasecmp(restart_type_now,"restart_pos") == 0){    istart_now = 2;}
    if(strcasecmp(restart_type_now,"restart_posvel") == 0){ istart_now = 3;}
    if(strcasecmp(restart_type_now,"restart_all") == 0){    istart_now = 4;}

    if(istart==1){strcpy(restart_type_spec,"initial");}
    if(istart==2){strcpy(restart_type_spec,"restart_pos");}
    if(istart==3){strcpy(restart_type_spec,"restart_posvel");}
    if(istart==4){strcpy(restart_type_spec,"restart_all");}

    if(istart_now==0){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("Start up %s option in ",restart_type_now);
      PRINTF("user specified coordinate file %s\n",dnamei);
      PRINTF("not supported. Supported genre: \n");
      PRINTF("initial, restart_pos, restart_posvel, restart_all-> \n");     
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    } // endif 
  } // endif: istart 

  //--------------------- 
  // C) General checks   

  if(istart_now < istart ) {
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Start up option, %s, in ",restart_type_now);
    PRINTF("user specified coordinate file %s\n",dnamei);
    PRINTF("Incompatible with class setup,%s\n",restart_type_spec);
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  } // endif 

  //==========================================================================

  natm_tot[0] = natm_tot_now;
  fclose(fp_dnamei);

  PRINTF("Finished reading the number of atoms from %s\n\n",dnamei);

  //==========================================================================
}//end routine
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//  Read in atom positions and atom types 
//==========================================================================
void read_coord(int istart, int natm_tot, double *x, double *y, double *z,
    NAME *atm_typ, double *hmat, char *dnamei)
  //==========================================================================
{// begin routine
  //==========================================================================

  int i;

  int istart_now;
  int natm_tot_now;
  int pi_beads_now;
  int itime_dump;
  int imol_num_now;

  double h1,h2,h3;

  char restart_type_spec[MAXWORD];
  char restart_type_now[MAXWORD];
  char atm_typ_now[MAXWORD];
  char res_typ_now[MAXWORD];
  char mol_typ_now[MAXWORD];

  char line[MAXLINE];

  FILE *fp_dnamei;

  //==========================================================================

  PRINTF("Reading the atom positions and system cell from %s\n",dnamei);
  fp_dnamei = fopen((const char *)dnamei,"r");

  //========================================================================
  //  I)Read in header                                                    

  //-------------------
  // A) Type 1 start : 
  if(istart==1){
    if(fscanf(fp_dnamei,"%d %d %d",&natm_tot_now,&istart_now,
          &pi_beads_now)!=3){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("Error reading start type and number of atoms \n");
      PRINTF("in file %s\n",dnamei);
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }//endif
    readtoendofline(fp_dnamei);
  }//endif
  //---------------------
  // B) Type 2,3,4 start 
  if(istart>1){
    if(fgetc(fp_dnamei)==EOF){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("Error reading header information \n");
      PRINTF("in file %s\n",dnamei);
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }//endif
    readtoendofline(fp_dnamei);

    if(fscanf(fp_dnamei,"%d %s %d %d",&natm_tot_now,restart_type_now,
          &itime_dump,&pi_beads_now)!=4){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("Error reading start type and number of atoms \n");
      PRINTF("in file %s\n",dnamei);
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    }//endif
    readtoendofline(fp_dnamei);

    istart_now = 0;
    if(strcasecmp(restart_type_now,"initial") == 0){        istart_now = 1;}
    if(strcasecmp(restart_type_now,"restart_pos") == 0){    istart_now = 2;}
    if(strcasecmp(restart_type_now,"restart_posvel") == 0){ istart_now = 3;}
    if(strcasecmp(restart_type_now,"restart_all") == 0){    istart_now = 4;}

    if(istart==1){strcpy(restart_type_spec,"initial");}
    if(istart==2){strcpy(restart_type_spec,"restart_pos");}
    if(istart==3){strcpy(restart_type_spec,"restart_posvel");}
    if(istart==4){strcpy(restart_type_spec,"restart_all");}

    if(istart_now==0){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("Start up %s option in ",restart_type_now);
      PRINTF("user specified coordinate file %s\n",dnamei);
      PRINTF("not supported. Supported genre: \n");
      PRINTF("initial, restart_pos, restart_posvel, restart_all-> \n");     
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    } // endif 
  } // endif: istart 

  //--------------------- 
  // C) General checks   
  if(istart_now < istart ) {
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Start up option, %s, in ",restart_type_now);
    PRINTF("user specified coordinate file %s\n",dnamei);
    PRINTF("Incompatible with class setup,%s\n",restart_type_spec);
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  } // endif 

  if(natm_tot_now != natm_tot ) {
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Incorrect number of atoms %d versus %d\n",natm_tot, natm_tot_now);
    PRINTF("user specified coordinate file %s\n",dnamei);
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  } // endif 

  //========================================================================
  // IV)istart = 1 (initial)                                              

  if(istart == 1) {
    //-----------------------------------------------------------------------
    // A)Atm positions                                                  
    for(i=1;i<=natm_tot;i++){
      if(fscanf(fp_dnamei,"%lf %lf %lf",
            &(x[i]),
            &(y[i]),
            &(z[i])) != 3) {
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("Error while reading in the %d atom coordinate\n",i);
        PRINTF("in file \"%s\"\n",dnamei);
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        FFLUSH(stdout);
        EXIT(1);
      }//endif
      readtoendofline(fp_dnamei);
      x[i] /= BOHR;      y[i] /= BOHR;      z[i] /= BOHR;
    }//endfor:atoms
    //------------------------------------------------------------------
    //B) Read the box
    for(i=0;i<3;i++){
      if(fscanf(fp_dnamei,"%lf %lf %lf",&h1,&h2,&h3) != 3){
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("Error while reading in the %d cell vector \n",i+1);
        PRINTF("in file \"%s\"\n",dnamei);
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        FFLUSH(stdout);
        EXIT(1);
      }//endif
      hmat[(i+1)]=h1/BOHR;  hmat[(i+4)]=h2/BOHR;  hmat[(i+7)]=h3/BOHR;
      readtoendofline(fp_dnamei);
    }//endfor
  }//endif: start=1 

  //========================================================================
  // V)istart = 2 (restart_pos)                                             

  if(istart >= 2) {
    if(fgets(line,MAXLINE,fp_dnamei)==NULL){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("EOF before particle coordinates \n");
      PRINTF("in file \"%s\"\n",dnamei);
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }//endif
    for(i=1;i<=natm_tot;i++){
      if(fscanf(fp_dnamei,"%lf %lf %lf %s %s %s %d",
            &(x[i]),
            &(y[i]),
            &(z[i]),
            atm_typ_now,res_typ_now,mol_typ_now,
            &imol_num_now) != 7) {
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("Error while reading in the %d atom coordinate\n",i);
        PRINTF("in file \"%s\"\n",dnamei);
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        FFLUSH(stdout);
        EXIT(1);
      }//endif
      readtoendofline(fp_dnamei);
      strcpy(atm_typ[i],atm_typ_now);
    }//endfor
    PRINTF("Reading in cell shape information\n");
    if(fgets(line,MAXLINE,fp_dnamei)==NULL){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("EOF before cell vectors \n");
      PRINTF("in file \"%s\"\n",dnamei);
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }//endif
    for(i=0;i<3;i++){
      if(fscanf(fp_dnamei,"%lf %lf %lf",&h1,&h2,&h3) != 3){
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("Error reading in cell vector %d\n",i);
        PRINTF("in file \"%s\"\n",dnamei);
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        FFLUSH(stdout);
        EXIT(1);
      }// endif
      hmat[(i+1)]=h1;  hmat[(i+4)]=h2;  hmat[(i+7)]=h3;
      readtoendofline(fp_dnamei);
    }// endfor
  }//endif: start>=2

  //==========================================================================
  // Done

  fclose(fp_dnamei);
  PRINTF("Finished reading the atom positions and system cell from %s\n\n",dnamei);

  //--------------------------------------------------------------------------
}//end routine
//==========================================================================



//=================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=================================================================
//  Generate Uniform Random Numbers
//=================================================================
#define MODULUS_R    2147483647 // DON'T CHANGE THIS VALUE       
#define MULTIPLIER_R 48271      // DON'T CHANGE THIS VALUE       
//=================================================================
// Random returns a pseudo-random real number uniformly distributed 
// between 0.0 and 1.0. 
//=================================================================
double altRandom(long *seed){
  long t;
  const long Q = MODULUS_R / MULTIPLIER_R;
  const long R = MODULUS_R % MULTIPLIER_R;

  t = MULTIPLIER_R * (seed[0] % Q) - R * (seed[0] / Q);
  if(t > 0){
    seed[0] = t;
  }else {
    seed[0] = t + MODULUS_R;
  }//endif
  return ((double) seed[0] / MODULUS_R);
}
//=================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Generate Gaussian Random numbers 
//===============================================================
void gaussran(int nran, long *seed, double *gauss)
  //========================================================================
{//begin routine
  //========================================================================
  //             Local variable declarations                                
  int i,loop;
  double twopi,rad2,al,r,phi,arg;
  //========================================================================
  // I) Constants 

  twopi = 2.0*M_PI;
  rad2  = sqrt(2.0);
  loop  = nran/2;

  //========================================================================
  // II) Make nran (or nran-1 if nran odd) Gaussian random numbers          
  for(i=1;i<=loop;i++){
    //------------------------------------------------------------------------
    // A) uniform random numbers in r and phi 
    r   = altRandom(seed); 
    phi = altRandom(seed); 
    r   = MAX(r,1e-30);
    r   = MIN(r,1.0);
    //------------------------------------------------------------------------
    // B) Gaussify in x and y
    al  = sqrt(-log(r))*rad2;
    arg = twopi*phi;
    gauss[2*i-1] = al*cos(arg);
    gauss[2*i]   = al*sin(arg);
  }//endfor
  //========================================================================
  // III) Make one more if nran is odd 

  if((nran % 2)!=0){
    r   = altRandom(seed); 
    phi = altRandom(seed); 
    r   = MAX(r,1e-30);
    r   = MIN(r,1.0);
    arg = twopi*phi;
    al  = sqrt(-log(r))*rad2;
    gauss[nran] = al*cos(arg);
  }//endif
  //------------------------------------------------------------------------
}//end routine
//========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// readtoendofline: Function to read to end of line in read_coord files     
//==========================================================================
void readtoendofline(FILE *fp){
  int eol,ch;
  eol = (int )'\n';
  ch = eol+1;
  while(ch!=eol){ch=fgetc(fp);}
  if(ch==EOF){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("ERROR: Unexpected end of file reached          \n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }//endif
}// end routine 
//==========================================================================
