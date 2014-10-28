#include "CPcharmParaInfo.decl.h"
#include "CPcharmParaInfo.h"

#ifndef PARAINFO_GROUP
#define PARAINFO_GROUP

class CPcharmParaInfoGrp: public Group 
{
  public:
    CPcharmParaInfoGrp(CkMigrateMessage *m)  {}
    CPcharmParaInfoGrp(CPcharmParaInfo &s)   { 
      for(int i=0;i<11;i++){
        CkPrintf("before group Rcom copy constructore %d : %d %d \n",i,s.RCommPkg[i].num_recv_tot,s.RCommPkg[i].num_send_tot);
      }
      cpcharmParaInfo = new CPcharmParaInfo(s); 
      for(int i=0;i<11;i++){
        CkPrintf("cpcharmParaInfo->RCommPkg %d : %d %d \n",i,cpcharmParaInfo->RCommPkg[i].num_recv_tot,cpcharmParaInfo->RCommPkg[i].num_send_tot);
      }

    }
    ~CPcharmParaInfoGrp()                    { delete cpcharmParaInfo; }
    CPcharmParaInfo *cpcharmParaInfo;
};

#endif // PARAINFO_GROUP

