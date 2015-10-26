//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//
//                  GWBSE simulation options                                   
//
//            class definition for GWBSE input options
//                                                                          
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
#ifndef _GWBSEOPTS_
#define _GWBSEOPTS_


class GWBSEOPTS{

  //----------------
  public:
    bool gwbse_on;               // Opt: Is GW-BSE of any type ``on''
    bool write_epsmatinv;        // Opt: write epsilon inverse matrix
    bool read_epsmatinv;         // Opt: read epsilon inverse matrix (if it skipps epsilon)
    bool doublePack;             // Opt: gamma point real wavefunctions if on
    
    int nstate;                  // Num: totan number of states in wavefunction files
    int nocc;                    // Num: number of occupied state
    int nunocc;                  // Num: number of unoccupied state
    int nkpt;                    // Num: number of k points
    int nspin;                   // Num: number of spin
    int coulb_trunc_opt;         // Num: coulomb truncation option
                                 // 0-no truncation 1-wire 2-sheet 3-molecule
    int ibinary_opt;             // Num: binary option 1, 2, 3
    double latt;                 // lattice parameter
    double a1[3], a2[3], a3[3];  // Num: lattice vectors (read it from lattice.dat)
    double b1[3], b2[3], b3[3];  // Num: reciprocal lattice vectors (read it from lattice.dat)
    double **kvec;               // Num: klists (read it from klist.dat)
    double *kwt;                 // Num: k weights
    double **qvec;               // Num: qlists (after reading klist, qvec will be calculated)

    char fileEpsMatInv[255];         // Chr: file name for inverse epsilon matrix
    char fileRho[255];               // Chr: file name for density


    
    //----------------
    //con-destruct:
    GWBSEOPTS(){
      gwbse_on        = false;
      write_epsmatinv = false;
      read_epsmatinv  = false;
      doublePack      = false;
      nstate          = 0;
      nocc            = 0;               
      nunocc          = 0;              
      nkpt            = 0;
      nspin           = 0;
      coulb_trunc_opt = 0;
      ibinary_opt     = 0;
    };
    ~GWBSEOPTS(){};

#ifdef PUP_ON
    //----------------
    //pupping :
    void pup(PUP::er &p){
      p | gwbse_on;
      p | write_epsmatinv;
      p | read_epsmatinv;
      p | doublePack;
      //pupping ints
      p | nstate;
      p | nocc;
      p | nunocc;
      p | nkpt;
      p | nspin;
      p | coulb_trunc_opt;
      p | ibinary_opt;
      //pupping dbles;
      p | latt;

      if(p.isUnpacking()) {
         kvec   = new double* [nkpt];
	 kwt = new double [nkpt];
         qvec = new double* [nkpt];
         for(int i=0; i<nkpt ; i++){
           kvec[i] = new double [3];
           qvec[i] = new double [3];
          }
      }
      //pupping 1d dble arrays;
      PUParray(p,a1,3);
      PUParray(p,a2,3);
      PUParray(p,a3,3);
      PUParray(p,b1,3);
      PUParray(p,b2,3);
      PUParray(p,b3,3);
      PUParray(p,kwt,nkpt);
      //pupping 2d dble arrars;
      //pup2d_dbl(p,&kvec,nkpt,3,"gwbseopts");
      for(int i=0;i<nkpt;i++){
        PUParray(p,kvec[i],3);
        PUParray(p,qvec[i],3);
      }
#ifdef _PARALLEL_DEBUG_        
      if (p.isUnpacking())
        state_class_out ();
#endif           
    } // end pup
#endif
    
    void state_class_out(){
      char fileName [255];
      sprintf (fileName, "%d_gwbseopts.out", CKMYPE());
      FILE *fp;  fp = fopen(fileName,"w");
      // ints
      fprintf(fp,"bwgse_on %d\n",gwbse_on);
      fprintf(fp,"write_epsmatinv %d\n",write_epsmatinv);
      fprintf(fp,"read_epsmatinv %d\n",read_epsmatinv);
      fprintf(fp,"doublePack %d\n",doublePack);
      fprintf(fp,"nstate %d\n",nstate);
      fprintf(fp,"nocc %d\n",nocc);
      fprintf(fp,"nunocc %d\n",nunocc);
      fprintf(fp,"nkpt %d\n",nkpt);
      fprintf(fp,"nspin %d\n", nspin);
      fprintf(fp,"coulb_trunc_opt %d\n",coulb_trunc_opt);
      fprintf(fp,"ibinary_opt %d\n",ibinary_opt);

      fprintf(fp,"lattice vectors \n");
      fprintf(fp,"%lg   %lg   %lg\n",a1[0],a1[1],a1[2]);
      fprintf(fp,"%lg   %lg   %lg\n",a2[0],a2[1],a2[2]);
      fprintf(fp,"%lg   %lg   %lg\n",a3[0],a3[1],a3[2]);
      fprintf(fp,"reciprocal lattice vectors \n");
      fprintf(fp,"%lg   %lg   %lg\n",b1[0],b1[1],b1[2]);
      fprintf(fp,"%lg   %lg   %lg\n",b2[0],b2[1],b2[2]);
      fprintf(fp,"%lg   %lg   %lg\n",b3[0],b3[1],b3[2]);
      fprintf(fp,"kpoint list and their weight \n");
      for (int i=0; i<nkpt; i++){
	  fprintf(fp,"%lg   %lg   %lg   %lg\n",kvec[i][0],kvec[i][1],kvec[i][2],kwt[i]);
      }
      fprintf(fp,"qpoint list \n");
      for (int i=0; i<nkpt; i++){
	  fprintf(fp,"%lg   %lg   %lg\n",qvec[i][0],qvec[i][1],qvec[i][2]);
      }


      // dbles
      
      fclose(fp);
    }// end routine 

}; // GWBSEOPTS

#ifdef PUP_ON
PUPmarshall(GWBSEOPTS);
#endif


#endif

//==========================================================================
