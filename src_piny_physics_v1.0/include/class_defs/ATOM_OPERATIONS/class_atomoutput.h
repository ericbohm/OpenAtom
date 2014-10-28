//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                          Atom Integration                                    
//==========================================================================

#ifndef _ATOMOUTPUT_
#define _ATOMOUTPUT_

class ATOMOUTPUT{

  //---------------------------------------------------------------------------
  public:

    //---------------------------------------------------------------------------
    //con-destruct:

    ATOMOUTPUT(){};
    ~ATOMOUTPUT(){};

    //---------------------------------------------------------------------------
    // Functions

    static void initialize_piny_output();
    static void write_gen_header(int ,int ,NAME ,char *);
    static void write_gen_header_cp(int ,int ,NAME ,char *);

    static void ctrl_piny_output(int ,int ,int ,int, int , Atom *, AtomNHC *,int*,int,int,int);
    static void write_atom_output_dump(int , int ,int , int,Atom *,AtomNHC *, char *);
    static void write_atom_output_conf(int , int , int , Atom *,char *);

    //---------------------------------------------------------------------------
}; //ATOMINTEGRATE
//==========================================================================

#endif
