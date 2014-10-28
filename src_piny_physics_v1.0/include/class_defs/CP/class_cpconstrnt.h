//==========================================================================
//  CP constrnt stuff 

#ifndef _CPCONSTRNT_
#define _CPCONSTRNT_

class CPCONSTRNT{

  //----------------
  public:
    double c_tolshake;          // Num: PW coef shake tolerence        
    double c_tolratl;           // Num: PW coef rattle tolerence       
    double c_tolnorb;           // Num: PW coef norb tolerence         

    //----------------
    //con-destruct:
    CPCONSTRNT(){
      c_tolshake = 0.0;
      c_tolratl  = 0.0;
      c_tolnorb  = 0.0;
    };
    ~CPCONSTRNT(){};

#ifdef PUP_ON
    //----------------
    //pupping
    void pup(PUP::er &p){
      //pupping dbles
      p | c_tolshake;
      p | c_tolratl;
      p | c_tolnorb;
#ifdef _PARALLEL_DEBUG_        
      if (p.isUnpacking())
        state_class_out ();
#endif           
    } // end pup
#endif

    void state_class_out(){

      char fileName [255];
      sprintf (fileName, "%d_cpconstrnt.state", CkMyPe());
      FILE *fp;  fp = fopen(fileName,"w");

      fprintf(fp,"c_tolshake %g\n",c_tolshake);
      fprintf(fp,"c_tolratl %g\n",c_tolratl);
      fprintf(fp,"c_tolnorb %g\n",c_tolnorb);
      fclose(fp);

    }// end routine

}; //CPCONSTRNT

#ifdef PUP_ON
PUPmarshall(CPCONSTRNT);
#endif

#endif


//==========================================================================
