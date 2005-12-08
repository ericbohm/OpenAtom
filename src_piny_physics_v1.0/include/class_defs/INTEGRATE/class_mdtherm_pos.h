//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                                                                          
//                      class_mdtherm_pos.h                                 
//                                                                          
//    Class definition for atom thermostat positions, velocities, forces    
//                                                                          
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================


#ifndef _MDTHERM_POS_
#define _MDTHERM_POS_

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================


class MDTHERM_POS {
 public:
  int num_nhc;                    // Num: Number of thermostats 
  int len_nhc;                    // Num: Length of thermostat chains 
  double **x_nhc,**v_nhc;         // Lst: Atm NHC pos,vel,for;         
                                  // Lth: num_nhc x len_nhc            

//--------------------------------------------------------------------
// Default constructor/destructor

  MDTHERM_POS(){
    num_nhc = 0;
    len_nhc = 0;
  }
  ~MDTHERM_POS(){}

//--------------------------------------------------------------------
#ifdef PUP_ON
  void pup(PUP::er &p){

    // PUP ints

    p | num_nhc;
    p | len_nhc;

  // PUP arrays

    pup2d_dbl(p,&x_nhc,len_nhc,num_nhc,"mdtherm_pos");
    pup2d_dbl(p,&v_nhc,len_nhc,num_nhc,"mdtherm_pos");


  } // end pack/unpack 
#endif

}; // end class definition 
//--------------------------------------------------------------------

#ifdef PUP_ON
PUPmarshall(MDTHERM_POS);
#endif

//--------------------------------------------------------------------
#endif
//==========================================================================

