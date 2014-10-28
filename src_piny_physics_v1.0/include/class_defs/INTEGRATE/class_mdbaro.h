//==========================================================================
//                  Andersen-Hoover NPT (NPT_I) info                        
//             {Variables needed for mem allocation:                        
//                                                  }                       
//                                                                          
//==========================================================================

#ifndef _MDBARO_
#define _MDBARO_

class MDBARO {

  //-----------------------------------------------------------------------
  public:
    int iopt;                   // Opt: npt_iso on/off
    int len_nhc;                // Num: length of NHC                  
    double x_lnv,v_lnv;         // Num: log(vol),dlog(vol)/dt          
    double v_lnv_glob;          // Num: dlog(vol)/dt                   
    double v_lnv_g;             // Num: dlog(vol)/dt                   
    double f_lnv_p,f_lnv_v;     // Num: d^2log(vol)/dt^2               
    double vol;                 // Num: volume                         
    double mass_lnv;            // Num: Mass of log(vol)               
    double c2_lnv;              // Num: Useful constant                
    double x_lnv_o;             // Num: old log(vol)                   
    double v_lnv_g_wght;        // Num: weighted v_lnv_g               
    double area;                // Num: area                           

    double *x_vol_nhc,*v_vol_nhc;// Lst: Volume NHCs Lth:len_nhc        
    double *f_vol_nhc,*mass_vol_nhc,*gkt_vol;
    double *hmato;              // Num: old cell matrix  :Lth 9        

    //-----------------------------------------------------------------------
    //con-destruct:

    MDBARO(){
      iopt    = 0;
      len_nhc = 0;
      hmato   = (double *) cmalloc(9*sizeof(int),"constr:mdbaro")-1;
    }
    ~MDBARO(){};

    //-----------------------------------------------------------------------
#ifdef PUP_ON
    //pupping
    void pup(PUP::er &p){
      //pupping ints
      p | iopt;
      p | len_nhc;
      //pupping dbles
      p | x_lnv;
      p | v_lnv;
      p | v_lnv_glob;
      p | v_lnv_g;
      p | f_lnv_p;
      p | f_lnv_v;
      p | vol;
      p | mass_lnv;
      p | c2_lnv;
      p | x_lnv_o;
      p | v_lnv_g_wght;
      p | area;
      //pupping dbl  arrays
      if(iopt>0){
        pup1d_dbl(p,&hmato,9);
        if(len_nhc>0){
          pup1d_dbl(p,&x_vol_nhc,len_nhc);
          pup1d_dbl(p,&v_vol_nhc,len_nhc);
          pup1d_dbl(p,&f_vol_nhc,len_nhc);
          pup1d_dbl(p,&mass_vol_nhc,len_nhc);
          pup1d_dbl(p,&gkt_vol,len_nhc);
        }//endif      
      }//endif      
#ifdef _PARALLEL_DEBUG_        
      if (p.isUnpacking())
        state_class_out ();
#endif     
    }// end pup
#endif

    //-----------------------------------------------------------------------
    void state_class_out(){

      int i;
      char fileName [255];
      sprintf (fileName, "%d_mdbaro.state", CkMyPe());
      FILE *fp; fp = fopen(fileName,"w"); 

      // int 
      fprintf(fp,"len_nhc %d\n",len_nhc);
      // dbles
      fprintf(fp,"x_lnv %g\n",x_lnv);
      fprintf(fp,"v_lnv %g\n",v_lnv);
      fprintf(fp,"v_lnv_glob %g\n",v_lnv_glob);
      fprintf(fp,"v_lnv_g %g\n",v_lnv_g);
      fprintf(fp,"f_lnv_p %g\n",f_lnv_p);
      fprintf(fp,"f_lnv_v %g\n",f_lnv_v);
      fprintf(fp,"vol %g\n",vol);
      fprintf(fp,"mass_lnv %g\n",mass_lnv);
      fprintf(fp,"c2_lnv %g\n",c2_lnv);
      fprintf(fp,"x_lnv_o %g\n",x_lnv_o);
      fprintf(fp,"v_lnv_g_wght %g\n",v_lnv_g_wght);
      fprintf(fp,"area %g\n",area);
      // dbl  arrays
      if(iopt>0){
        for(i=1;i<=len_nhc;i++){fprintf(fp,"mass_vol_nhc[%d] %g\n",i,
            mass_vol_nhc[i]);}
        for(i=1;i<=len_nhc;i++){fprintf(fp,"gkt_vol[%d] %g\n",i,gkt_vol[i]);}
      }// endif
      fclose(fp);
    } // end routine

    //-----------------------------------------------------------------------
}; // MDBARO;
//==========================================================================

#ifdef PUP_ON
PUPmarshall(MDBARO);
#endif

#endif
//==========================================================================
