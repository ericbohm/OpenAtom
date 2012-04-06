#ifndef STREAMED_CHUNK_H
#define STREAMED_CHUNK_H

#include "charm++.h"

struct streamedChunk
{
    const static int sz = 20;
    // The indices of the sending GSpace chare
    short i, j, k;
    // The total number of data elements that this GSpace chare will send to the receiver
    short numDatums;
    // The sequence number of this chunk in the series of chunks being sent between this send-recv pair
    short chunkSeqNum;
    // Store either data or a pointer
    complex data[sz];

    streamedChunk()
        : i(0), j(0), k(0)
        , numDatums(0)
        , chunkSeqNum(-1)
    {}

    streamedChunk(const short _i, const short _j, const short _k, int nDatum, int seqNum, complex *buf = NULL)
        : i(_i), j(_j), k(_k)
        , numDatums(nDatum)
        , chunkSeqNum(chunkSeqNum)
    {
        if (buf)
            memcpy(data, buf, sz*sizeof(complex));
    }
};
PUPbytes(streamedChunk)

#endif // STREAMED_CHUNK_H

