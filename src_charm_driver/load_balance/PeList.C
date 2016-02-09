/** \file PeList.C
 *
 */


#include "charm++.h"
#include "PeList.h"
#include <algorithm>
#include "TopoManager.h"
extern TopoManager *topoMgr;
extern Config config;

#ifndef USE_ROUND_ROBIN_INSTANCE_MAP
static inline int IN_THIS_INSTANCE(int globalCount, int instanceIndex) {
  if(globalCount >= config.numPesPerInstance * instanceIndex &&
    globalCount < config.numPesPerInstance * (instanceIndex + 1)) {
    return 1;
  } else return 0;
}
#else
static inline int IN_THIS_INSTANCE(int globalCount, int instanceIndex) {
  if(globalCount % config.numInstances == instanceIndex) {
    return 1;
  } else return 0;
}
#endif

PeList::PeList(int ndims, const int bdims[10], int instanceIndex) {
  isSet = 1; //needed because there are PE tht may not be part of a box
  isList = 0; //has to be a vector since we want to jump to a loc
  isSetOnly = 0;
  LisV = 0;
    
#if USE_BITVECTOR
  if(isSet) {
    allPE.clear();
    allPE.resize(config.numPes, false);
    true_size = 0;
  }
#endif
  int pe_added = 0;
  int size = config.numPes;
  vectorList.reserve(size);

  if(config.simpleTopo == 1) {
    CkAssert(topoMgr != NULL);

    if(ndims == 3) {
      int globalCount = 0;
      int boxes[4], dims[4], order[4];
      TopoManager_getDims(dims);
      dims[3] *= CkMyNodeSize();
      printf("Dims are %d %d %d %d\n", dims[0], dims[1], dims[2], dims[3]);
      if(dims[0] >= dims[1] && dims[0] >= dims[2]) {
        order[0] = 0;
        if(dims[1] >= dims[2]) {
          order[1] = 1; order[2] = 2;
        } else {
          order[1] = 2; order[2] = 1;
        }
      } else if(dims[1] >= dims[0] && dims[1] >= dims[2]) {
        order[0] = 1;
        if(dims[0] >= dims[2]) {
          order[1] = 0; order[2] = 2;
        } else {
          order[1] = 2; order[2] = 0;
        }
      } else if(dims[2] >= dims[0] && dims[2] >= dims[1]) {
        order[0] = 2;
        if(dims[0] >= dims[1]) {
          order[1] = 0; order[2] = 1;
        } else {
          order[1] = 1; order[2] = 0;
        }
      } else {
        printf("Weird math...none of the dimension seems to be largest\n");
        order[0] = 0; order[1] = 1; order[2] = 2;
      }
      printf("Order is %d %d %d\n", order[0], order[1], order[2]);

      for(boxes[ order[0] ] = 0; boxes[ order[0] ] < dims[ order[0] ]; boxes[ order[0] ] += bdims[ order[0] ]) {
        for(boxes[ order[1] ] = 0; boxes[ order[1] ] < dims[ order[1] ]; boxes[ order[1] ]  += bdims[ order[1] ]) {
          for(boxes[ order[2] ] = 0; boxes[ order[2] ] < dims[ order[2] ]; boxes[ order[2] ] += bdims[ order[2] ]) {

            int begin[3], end[3], curLoc[4];
            for(int i = 0; i < 3; i++) {
              begin[i] = boxes[i];
              end[i] = std::min((begin[i] + bdims[i]), dims[i]);
            }

            for(curLoc[ order[0] ] = begin[ order[0] ]; curLoc[ order[0] ] < end[ order[0] ]; curLoc[ order[0] ] ++) {
              for(curLoc[ order[1] ] = begin[ order[1] ]; curLoc[ order[1] ] < end[ order[1] ]; curLoc[ order[1] ] ++) {
                for(curLoc[ order[2] ] = begin[ order[2] ]; curLoc[ order[2] ] < end[ order[2] ]; curLoc[ order[2] ] ++) {
                  for(curLoc[3] = 0; curLoc[3] < dims[3]; curLoc[3]++) {
                    int pe = -1;
                    TopoManager_getPeRank(&pe, curLoc);
                    if(pe > 0 || (!config.excludePE0 && pe == 0)) {
                      if(IN_THIS_INSTANCE(globalCount, instanceIndex)) {
#if USE_BITVECTOR
                        allPE[pe] = true;
                        true_size++;
#else
                        setList.insert(pe);
#endif
                        vectorList.push_back(pe);
                        pe_added++;
                      }
                    }
                    if(pe >= 0) {
                      globalCount++;
                    }
                  }
                }
              }
            }
          }
        }
      }
    } else if(ndims == 5) {
      int globalCount = 0;
      int boxes[6], dims[6], order[6];
      TopoManager_getDims(dims);
      dims[5] *= CkMyNodeSize();
      printf("Dims are %d %d %d %d %d %d\n", dims[0], dims[1], dims[2], dims[3],
          dims[4], dims[5]);

      //find largest two index
      int largest = 0, val = dims[0];
      for(int i = 1; i < 4; i++) {
        if(val < dims[i]) {
          val = dims[i];
          largest = i;
        }
      }
      order[0] = largest;

      if(order[0] != 0) { 
        largest = 0; val = dims[0];
      } else {
        largest = 1; val = dims[1];
      }

      for(int i = 1; i < 4; i++) {
        if(i != order[0]) {
          if(val < dims[i]) {
            val = dims[i];
            largest = i;
          }
        }
      }
      order[1] = largest;

      int loc = 2;
      for(int i = 0; i < 4; i++) {
        if(order[0] != i && order[1] != i) {
          order[loc++] = i;
        }
      }
      order[4] = 4;
      printf("Order is %d %d %d %d %d\n", order[0], order[1], order[2], order[3],
          order[4]);

      for(boxes[ order[0] ] =  0 ; boxes[ order[0] ] < dims[ order[0] ]; boxes[ order[0] ] += bdims[ order[0] ]) {
        for(boxes[ order[1] ] = 0; boxes[ order[1] ] < dims[ order[1] ]; boxes[ order[1] ]  += bdims[ order[1] ]) {
          for(boxes[ order[2] ] = 0; boxes[ order[2] ] < dims[ order[2] ]; boxes[ order[2] ] += bdims[ order[2] ]) {
            for(boxes[ order[3] ] = 0; boxes[ order[3] ] < dims[ order[3] ]; boxes[ order[3] ] += bdims[ order[3] ]) {
              for(boxes[4] = 0; boxes[4] < dims[4]; boxes[4] += bdims[4]) {

                int begin[5], end[5], curLoc[6];
                for(int i = 0; i < 5; i++) {
                  begin[i] = boxes[i];
                  end[i] = std::min((begin[i] + bdims[i]), dims[i]);
                }

                for(curLoc[ order[0] ] = begin[ order[0] ]; curLoc[ order[0] ] < end[ order[0] ]; curLoc[ order[0] ]++) {
                  for(curLoc[ order[1] ] = begin[ order[1] ]; curLoc[ order[1] ] < end[ order[1] ]; curLoc[ order[1] ]++) {
                    for(curLoc[ order[2] ] = begin[ order[2] ]; curLoc[ order[2] ] < end[ order[2] ]; curLoc[ order[2] ]++) {
                      for(curLoc[ order[3] ] = begin[ order[3] ]; curLoc[ order[3] ] < end[ order[3] ]; curLoc[ order[3] ]++) {
                        for(curLoc[4] = begin[4]; curLoc[4] < end[4]; curLoc[4]++) {
                          for(curLoc[5] = 0; curLoc[5] < dims[5]; curLoc[5]++) {
                            int pe = -1;
                            TopoManager_getPeRank(&pe, curLoc);
                            if(pe > 0 || (!config.excludePE0 && pe == 0)) {
                              if(IN_THIS_INSTANCE(globalCount, instanceIndex)) {
#if USE_BITVECTOR
                                allPE[pe] = true;
                                true_size++;
#else
                                setList.insert(pe);
#endif
                                vectorList.push_back(pe);
                                pe_added++;
                              }
                            }
                            if(pe >= 0) {
                              globalCount++;
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    } else {
      CkAbort("Unsupported ndims\n");
    }
  
  } else {
    CkAbort("You shouldn't be calling this function\n");
  }
  reset();
}

PeList::PeList() {
  CkAbort("Default constructor not implemented\n");
}

void PeList::sortSource(int srcPe, int pushToList) {
  CkAssert(srcPe >= 0);
  CkAssert(srcPe < config.numPes);
  int pe_dims[10];
  TopoManager_getPeCoordinates(srcPe, pe_dims);
  sortSource(pe_dims, pushToList);
}

void PeList::sortSource(int *srcPe, int pushToList) {
  if(isSetOnly) {
    CkAbort("Cannot sort a list which is setOnly\n");
  }
  //first ensure, we have an updated vector
  std::vector< int > sortIdx, tempVec;
  if(isList && !LisV) {
    if(config.simpleTopo && !pushToList) { 
      //store to a tempArray
      tempVec.reserve(list.size());
      std::copy(list.begin(), list.end(), std::back_inserter(tempVec));
      vectorList.clear();
      vectorList.resize(list.size());
    } else {
      vectorList.clear();
      vectorList.reserve(list.size());
      std::copy(list.begin(), list.end(), std::back_inserter(vectorList));
    }
  } else {
    if(config.simpleTopo && !pushToList) { 
      //store to a tempArray
      tempVec.resize(vectorList.size());
      memcpy(&tempVec[0], &vectorList[0], vectorList.size() * sizeof(int));
    }
  }
  if(config.simpleTopo) {
    sortIdx.resize(vectorList.size());
    for(int i = 0; i < sortIdx.size(); i++) {
      sortIdx[i] = i;
    }
    if(pushToList) {
      topoMgr->sortRanksByHops(srcPe, &vectorList[0], &sortIdx[0], vectorList.size());
    } else {
      topoMgr->sortRanksByHops(srcPe, &tempVec[0], &sortIdx[0], tempVec.size());
    }
    if(pushToList) {
      list.clear();
      for(int i = 0; i < sortIdx.size(); i++) {
        list.push_back(vectorList[sortIdx[i]]);
      }
      isList = 1;
      LisV = 1;
    } else {
      for(int i = 0; i < sortIdx.size(); i++) {
        vectorList[i] = tempVec[sortIdx[i]];
      }
      if(isList) {
        isList = 0;
        LisV = 1;
      }
    }
  } else {
    std::sort(vectorList.begin(), vectorList.end());
    if(pushToList) {
      list.clear();
      for(int i = 0; i < vectorList.size(); i++) {
        list.push_back(vectorList[i]);
      }
      isList = 1;
      LisV = 1;
    } else {
      if(isList) {
        isList = 0;
        LisV = 1;
      }
    }
  }
  reset();
}

int PeList::minDist(int srcPe) {
  CkAssert(srcPe >= 0);
  CkAssert(srcPe < config.numPes);
  int pe_dims[10];
  TopoManager_getPeCoordinates(srcPe, pe_dims);
  return minDist(pe_dims);
}

int PeList::minDist(int *srcPe)
{
  if(isSetOnly) {
    CkAbort("Cannot find minDist in a list which is setOnly\n");
    return -1;
  }
  if(config.simpleTopo) {
    if(isList) {
      std::list< int >::iterator it = list.begin();
      offsetList = 0;
      int count = 0;
      currList = it;
      int minHops = topoMgr->getHopsBetweenRanks(srcPe, *it);
      int minPE = *it;
      if(minHops == 0) return minPE;
      it++;
      while(it != list.end()) {
        count++;
        int nowHops = topoMgr->getHopsBetweenRanks(srcPe, *it);
        if(nowHops < minHops) {
          minHops = nowHops;
          minPE = *it;
          currList = it;
          offsetList = count;
          if(minHops == 0) return minPE;
        }
        it++;
      }
      return minPE;
    } else {
      std::vector< int >::iterator it = vectorList.begin();
      int minHops = topoMgr->getHopsBetweenRanks(srcPe, *it);
      int minPE = *it;
      if(minHops == 0) return minPE;
      it++;
      while(it != vectorList.end()) {
        int nowHops = topoMgr->getHopsBetweenRanks(srcPe, *it);
        if(nowHops < minHops) {
          minHops = nowHops;
          minPE = *it;;
          if(minHops == 0) return minPE;
        }
        it++;
      }
      return minPE;
    }
  } else {
    // if not a topology distance is meaningless we just pick the first element
    if(isList) {
      std::list< int >::iterator it = list.begin();
      currList = it;
      offsetList = 0;
      return *it;
    } else {
      return vectorList[0];
    }
  }
}

int PeList::findNext() {
  int pe;
  CkAssert(isSetOnly == 0);
  if(size() == 0) {
    reset();
    return -1;
  }
  if(isList) {
    if(!dontIncr) {
      currList++;
      offsetList++;
      if(currList == list.end()) {
        currList = list.begin();
        offsetList = 0;
      }
    } else {
      dontIncr = false;
    }
    pe = *currList;
  } else {
    if(!dontIncr) {
      offsetVector++;
      if(offsetVector >= vectorList.size()) {
        offsetVector = 0;
      }
    } else {
      dontIncr = false;
    }
    pe = vectorList[offsetVector];
  }
  return pe;
}

PeList *PeList::distributeAcrossPelist(int numElements, int listoffset)
{
  if(isSetOnly) {
    CkAbort("Cannot distributeAcrossPelist in a list which is setOnly\n");
  }
  int newListSize = (size() > numElements) ? numElements : size();
  PeList *outlist = new PeList(1, 1, 0);
  int stride = size() / numElements + 1;
  int mod = size() % numElements;
  if(isList) {
    std::list< int >::iterator it = list.begin();
    for(int i = 0; i < listoffset; i++) it++;
    for(int i = 0; i < newListSize; i++) {
      if(it == list.end()) it = list.begin();
      outlist->checkAndAdd(*it);
      if(i == mod) stride--;
      for(int j = 0; j < stride; j++) {
        if(it == list.end()) it = list.begin();
        it++;
      }
    }
  } else {
    for(int i = 0; i < newListSize; i++, listoffset += stride) {
      outlist->checkAndAdd(vectorList[listoffset % size()]);
      if(i == mod) stride--;
    }
  }
  outlist->reset();
  return outlist;
}


