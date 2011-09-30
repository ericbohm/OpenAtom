/** \file atomMessages.h
 */
#ifndef ATOMMESSAGES_H
#define ATOMMESSAGES_H


class AtomMsg: public CMessage_AtomMsg 
{
    public:
        int nsize;
        int natmStr,natmEnd;
        double *data;  
};


//class AtomXYZMsg: public CkMcastBaseMsg, public CMessage_AtomXYZMsg
class AtomXYZMsg: public CMessage_AtomXYZMsg
{
    public:
        double *x;  
        double *y;  
        double *z;  
	int index;
};
#endif //ATOMMESSAGES_H
