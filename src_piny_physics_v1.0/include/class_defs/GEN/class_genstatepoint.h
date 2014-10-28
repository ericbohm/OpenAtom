//==========================================================================
//                  State point data                                        
//             {Variables needed for mem allocation:}                       

#ifndef _GENSTATEPOINT_
#define _GENSTATEPOINT_

class GENSTATEPOINT{

  //----------------
  public:
    double pext,t_ext, stens_ext;  // Num: External press-temp, surface tens  

    //----------------
    //con-destruct:
    GENSTATEPOINT(){};
    ~GENSTATEPOINT(){};

#ifdef PUP_ON
    //----------------
    //pupping
    void pup(PUP::er &p){
      //pupping dbles
      p | pext;
      p | t_ext;
      p | stens_ext;
#ifdef _PARALLEL_DEBUG_        
      if (p.isUnpacking())
        state_class_out ();
#endif      
    } // end pup
#endif

    void state_class_out(){

      char fileName [255];
      sprintf (fileName, "%d_genstatepoint.state", CkMyPe());
      FILE *fp;  fp = fopen(fileName,"w");

      fprintf(fp,"pext %g\n",pext);
      fprintf(fp,"t_ext %g\n",t_ext);
      fprintf(fp,"stens_ext %g\n",stens_ext);
      fclose(fp);
    }// end routine



}; //STATEPOINT;

#ifdef PUP_ON
PUPmarshall(GENSTATEPOINT);
#endif

#endif

//==========================================================================
