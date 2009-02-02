#include "CPcharmParaInfo.decl.h"
#include "CPcharmParaInfo.h"

#ifndef PARAINFO_GROUP
#define PARAINFO_GROUP

class CPcharmParaInfoGrp: public Group 
{
    public:
        CPcharmParaInfoGrp(CkMigrateMessage *m)  {}
        CPcharmParaInfoGrp(CPcharmParaInfo &s)   { cpcharmParaInfo = new CPcharmParaInfo(s); }
        ~CPcharmParaInfoGrp()                    { delete cpcharmParaInfo; }
        CPcharmParaInfo *cpcharmParaInfo;
};

#endif // PARAINFO_GROUP

