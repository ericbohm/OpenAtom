#ifndef _atoms_h_
#define _atoms_h_
class Atom {

 public:
  Atom(double x_, double y_, double z_, double q_, double m_) {
     x = x_;   y = y_;   z = z_;  q = q_; m = m_; 
    vx = 0.0; vy = 0.0; vz = 0.0;
    fx = 0.0; fy = 0.0; fz = 0.0;
  }
  Atom(Atom &a) { x = a.x;   y = a.y;   z = a.z;  q = a.q; m = a.m; 
                 vx = a.vx; vy = a.vy; vz = a.vz;
                 fx = a.fx; fy = a.fy; fz = a.fz;
  }
  void Init(double x_, double y_, double z_, double q_, double m_) {
    x = x_; y = y_; z = z_; q = q_; m = m_; 
    vx = 0.0; vy = 0.0; vz = 0.0;
    fx = 0.0; fy = 0.0; fz = 0.0;
  }
  Atom() {}

  double x, y, z, q;
  double fx,fy,fz;
  double vx,vy,vz;
  double m;

  void pup(PUP::er &p) {
     p|x;  p|y;  p|z;  p|q; p|m; 
     p|fx; p|fy; p|fz; 
     p|vx; p|vy; p|vz;
  }//routine

};

PUPmarshall(Atom);
#endif
