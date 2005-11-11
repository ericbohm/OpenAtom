#ifndef _RunDescriptor_
#define _RunDescriptor_

//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
// Index logic for lines of constant x,y in gspace.
// Each line is broken up into two parts (0...-X) and (-X...-1)
//===================================================================================

struct RunDescriptor {
  int x, y, z;   // g-space starting coordinate of the line in g-space
  int position;  // starting position of run in the serialized array
  int length;    // length of the run (# of pts in this line
  int inc;
  int mychar;
  RunDescriptor(int px, int py, int pz, int pposition, int plength, int pinc, 
		int inchar) 
    :   x(px), y (py), z (pz), position(pposition), length(plength),inc( pinc),
       mychar(inchar) {}
  RunDescriptor(int px, int py, int pz, int pposition, int plength, int pinc) 
    :   x(px), y (py), z (pz), position(pposition), length(plength),inc( pinc)
  {mychar=0;}
  RunDescriptor(){}
  RunDescriptor(int i){}
  RunDescriptor(const RunDescriptor &r) {x = r.x; y = r.y; z = r.z;
                position = r.position; length = r.length; inc = r.inc; 
		mychar=r.mychar;}
  void operator=(const RunDescriptor &r) {x = r.x; y = r.y; z = r.z; 
                 position = r.position; length = r.length; inc = r.inc;
		 mychar=r.mychar;
  }
  void pup(PUP::er &p) {
    p|x; 
    p|y; 
    p|z; 
    p|position; 
    p|length; 
    p|inc;
    p|mychar;
  }
};

//===================================================================================

#endif
