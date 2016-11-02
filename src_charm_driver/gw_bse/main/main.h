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


#endif //__MAIN_H__
