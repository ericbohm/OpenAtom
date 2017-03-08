#ifndef __MAIN_H__
#define __MAIN_H__

class Main : public CBase_Main {

 public:

  /// Constructors ///
  Main(CkArgMsg* msg);
  Main(CkMigrateMessage* msg);
  
  /// Entry Methods ///
  void done();
};

class EpsMap : public CkArrayMap {
    int _x, _y, pe, rem, quo;
  public:
    EpsMap(){
      _x = _y = 7;
      int total = _x*_y;
      rem = total%CkNumNodes();
      quo = total/CkNumNodes();
    }
    ~EpsMap(){}
    int procNum(int, const CkArrayIndex &idx) {
      int count = idx.data()[0]*_y +idx.data()[1];
      int pe = 0;
      if(count < quo+rem) return pe;
      int sub = count - (quo+rem);
      pe = sub/quo +1;
      return pe;
    }
};


#endif //__MAIN_H__
