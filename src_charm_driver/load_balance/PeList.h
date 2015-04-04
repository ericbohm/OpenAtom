/** \ingroup mapping
 *
 */
//@{
#ifndef _PELIST_H
#define _PELIST_H

#include <cmath>
#include "configure.h"
#include "TopoManager.h"
#include "cklists.h"
#include <set>
#include <vector>
#include <list>

#define USE_BITVECTOR 1
extern Config                            config;

#if USE_BITVECTOR
static inline void copy_if(std::vector< bool > & from, std::vector< int > & to) {
  for(int i = 0; i < from.size(); i++) {
    if(from[i]) {
      to.push_back(i);
    }
  }
}
static inline void copy_if(std::vector< bool > & from, std::list< int > & to) {
  for(int i = 0; i < from.size(); i++) {
    if(from[i]) {
      to.push_back(i);
    }
  }
}
#endif

/* rewriting PeList - dropping sorting functionality */
class PeList 
{
  public:
#if USE_BITVECTOR
    std::vector< bool > allPE;
    int true_size;
#else
    std::set< int > setList;
#endif
    std::vector< int > vectorList;
    std::list< int > list;
    std::list< int >::iterator currList;
    int offsetList, offsetVector;;
    int isSet, isList;
    int isSetOnly;
    int LisV;
    bool dontIncr;

    PeList();  //default constructor

    PeList(int ndims, const int bdims[10], int instanceIndex);
    
    PeList(int _isSet, int _isList, int _size) {
      isSet = _isSet;
      isList = _isList;
      if(!isList) vectorList.reserve(_size);
#if USE_BITVECTOR
      if(isSet) {
        allPE.clear();
        allPE.resize(config.numPes, false);
        true_size = 0;
      }
#endif
      for(int i = 0; i < _size; i++) {
        if(isSet) {
#if USE_BITVECTOR
          allPE[i] = true;
          true_size++;
#else
          setList.insert(i);
#endif
        }
        if(isList) {
          list.push_back(i);
        } else { 
          vectorList.push_back(i);
        }
      }
      LisV = 0;
      isSetOnly = 0;
      reset();
    }

    PeList(int _isSet, int _isList, CkVec <int> inList) {
      isSet = _isSet;
      isList = _isList;
      vectorList.reserve(inList.size());
#if USE_BITVECTOR
      if(isSet) {
        allPE.clear();
        allPE.resize(config.numPes, false);
        true_size = 0;
      }
#endif
      for(int i = 0; i < inList.size(); i++) {
        if(isSet) {
#if USE_BITVECTOR
          allPE[inList[i]] = true;
          true_size++;
#else
          setList.insert(inList[i]);
#endif
        }
        if(isList) {
          list.push_back(inList[i]);
        } else {
          vectorList.push_back(inList[i]);
        }
      }
      LisV = 0;
      isSetOnly = 0;
      reset();
    }

    PeList(int _isSet, int _isList, PeList &inList) {
      isSet = _isSet;
      isList = _isList;
      isSetOnly = 0;
      
      if(inList.isSetOnly) {
        CkAbort("Cannot create from a list with has only set\n");
      }
      
#if USE_BITVECTOR
      if(isSet) {
        allPE.clear();
        allPE.resize(config.numPes, false);
        true_size = 0;
      }
#endif

      if(!isList) vectorList.reserve(inList.size());
      append(inList);

      LisV = 0;
      reset();
    }

    PeList(int _isSet, int _isList, int _size, int *a) {
      isSet = _isSet;
      isList = _isList;
      if(!isList) vectorList.reserve(_size);
#if USE_BITVECTOR
      if(isSet) {
        allPE.clear();
        allPE.resize(config.numPes, false);
        true_size = 0;
      }
#endif
      for(int i = 0; i < _size; i++) {
        if(isSet) {
#if USE_BITVECTOR
          allPE[a[i]] = true;
          true_size++;
#else
          setList.insert(a[i]);
#endif
        }
        if(isList) {
          list.push_back(a[i]);
        } else {
          vectorList.push_back(a[i]);
        }
      }
      LisV = 0;
      isSetOnly = 0;
      reset();
    }	

    PeList(int _isSet, int _isList, PeList &inList, int start, int _size) {
      CkAssert(start < inList.size());
      CkAssert(start + _size <= inList.size());
     
      if(inList.isSetOnly) {
        CkAbort("Cannot create from a list with has only set\n");
      }
      if(inList.isList && !LisV) {
        CkAbort("Cannot insert a PeList from an offset if not in vector form\n");
      }

      isSet = _isSet;
      isList = _isList;
      if(!isList) vectorList.reserve(_size);
#if USE_BITVECTOR
      if(isSet) {
        allPE.clear();
        allPE.resize(config.numPes, false);
        true_size = 0;
      }
#endif
      while(_size > 0) {
        if(isSet) {
#if USE_BITVECTOR
          allPE[inList.vectorList[start]] = true;
          true_size++;
#else
          setList.insert(inList.vectorList[start]);
#endif
        }
        if(isList) {
          list.push_back(inList.vectorList[start]);
        } else {
          vectorList.push_back(inList.vectorList[start]);
        }
        start++; 
        _size--;
      }
      LisV = 0;
      isSetOnly = 0;
      reset();
    }

    inline int size() {
      if(isSetOnly) {
#if USE_BITVECTOR
        return true_size;
#else 
        return setList.size();
#endif
      } else if(isList) {
        return list.size();
      } else {
        return vectorList.size();
      }
    }

    inline int count() { 
      if(isSetOnly) {
        CkAbort("Cannot get a count from a list with has only set\n");
      }
      if(isList) {
        return (list.size() - offsetList);
      } else {
        return (vectorList.size() - offsetVector);
      }
    }

    inline bool exists(int target) {
      if(isSet) {
#if USE_BITVECTOR
        return allPE[target];
#else
        return (setList.count(target) == 1);
#endif
      } else {
        CkAbort("Cannot query when set is not maintained\n");
        return false;
      }
    }

    void reset() {
      currList = list.begin();
      offsetList = 0;
      offsetVector = 0;
      dontIncr = true;
    } 

    void trimUsed(int ignoreSet = 0) {
      if((!ignoreSet && isSet) || !isList) {
        CkAbort("Trim used not supported for set and vectors\n");
      }

      if(!ignoreSet) {
        std::list< int >::iterator it = list.begin();
        while(it != currList) {
#if USE_BITVECTOR
          allPE[*it] = false;
          true_size--;
#else
          setList.erase(*it);
#endif
          it++;
        }
      }

      list.erase(list.begin(), currList);
      LisV = 0;
      reset();
    }

    void append(PeList &inList) {
      int currSize = size();
#if USE_BITVECTOR
      bool exists;
#else
      std::pair<std::set< int >::iterator,bool> exists;
#endif
      if(inList.isSetOnly) {
        CkAbort("Cannot append using a list with has only set\n");
      }
      
      if(inList.isList == 1) {
        std::list< int >::iterator it = inList.list.begin();
        while(it != inList.list.end()) {
#if USE_BITVECTOR
          exists = false;
          if(isSet) {
            exists = !(allPE[*it]);
            allPE[*it] = true;
            if(exists) true_size++;
          }
          if((!isSet || exists) && !isSetOnly) {
#else
          exists.second = false;
          if(isSet) {
            exists = setList.insert(*it);
          }
          if((!isSet || exists.second) && !isSetOnly) {
#endif
            if(isList) {
              list.push_back(*it);
            } else {
              vectorList.push_back(*it);
            }
          }
          it++;
        }
      } else {
        std::vector< int >::iterator it = inList.vectorList.begin();
        while(it != inList.vectorList.end()) {
#if USE_BITVECTOR
          exists = false;
          if(isSet) {
            exists = !(allPE[*it]);
            allPE[*it] = true;
            if(exists) true_size++;
          }
          if((!isSet || exists) && !isSetOnly) {
#else
          exists.second = false;
          if(isSet) {
            exists = setList.insert(*it);
          }
          if((!isSet || exists.second) && !isSetOnly) {
#endif
            if(isList) {
              list.push_back(*it);
            } else {
              vectorList.push_back(*it);
            }
          }
          it++;
        }
      }
      LisV = 0;
      if(currSize == 0) reset();
    }

    inline void deleteCurrent() {
      if(!isList) {
        CkAbort("Cannot delete when list is not maintained\n");
      }
      if(isSet) {
#if USE_BITVECTOR
        allPE[*currList] = false;
        true_size--;
#else
        setList.erase(*currList);
#endif
      }
      //printf("Erase %d %d %d %d\n", size(), offsetList, *currList, currList == list.begin());
      currList = list.erase(currList);
      if(currList == list.end()) {
        currList = list.begin();
        offsetList = 0;
      }
      LisV = 0;
      dontIncr = true;
    }

    inline void addOne(int pe) {
      int currSize = size();
      if(isSet) {
#if USE_BITVECTOR
        allPE[pe] = true;
        true_size++;
#else
        setList.insert(pe);
#endif
      }
      if(!isSetOnly) {
        if(isList) {
          list.push_back(pe);
        } else {
          vectorList.push_back(pe);
        }
      }
      LisV = 0;
      if(currSize == 0) reset();
    }

    inline void checkAndAdd(int pe) {
      int currSize = size();
      if(isSet) {
#if USE_BITVECTOR
        bool exists = !allPE[pe];
        allPE[pe] = true;
        if(exists) true_size++;
        if(exists && !isSetOnly) {
#else
        std::pair<std::set< int >::iterator,bool> exists;
        exists = setList.insert(pe);
        if(exists.second && !isSetOnly) {
#endif
          if(isList) {
            list.push_back(pe);
          } else {
            vectorList.push_back(pe);
          }
        }
        LisV = 0;
      } else {
        CkAbort("Cannot check and add when set is not maintained\n");
      }
      if(currSize == 0) reset();
    }

    int findNext();       // return next available, 
    void sortSource(int srcPe, int pushToList);
    void sortSource(int *srcPe, int pushToList);
    int minDist(int srcPe);
    int minDist(int *srcPe);

    PeList &operator=(PeList &inList) {
      isList = inList.isList;
      isSet = inList.isSet;
      LisV = inList.LisV;
      isSetOnly = inList.isSetOnly;
    
#if USE_BITVECTOR
      allPE = inList.allPE;
      true_size = inList.true_size;
#else
      setList = inList.setList;
#endif
      vectorList = inList.vectorList;
      list = inList.list;
  
      reset();
      return *this;
    }

    inline void addListorVector(int _isList) {
      if(!isSetOnly && (_isList == isList)) return;
      if(isSetOnly) {
        isList = _isList;
        isSetOnly = 0;
        list.clear();
        vectorList.clear();
        if(isList) {
#if USE_BITVECTOR
          copy_if(allPE, list);
#else
          std::copy(setList.begin(), setList.end(), std::back_inserter(list));
#endif
        } else {
#if USE_BITVECTOR
          vectorList.reserve(true_size);
          copy_if(allPE, vectorList);
#else
          vectorList.reserve(setList.size());
          std::copy(setList.begin(), setList.end(), std::back_inserter(vectorList));
#endif
        }
        LisV = 0;
        reset();
        return;
      }
      if(_isList) {
        addList();
      } else {
        addVector();
      }
    }

    inline void addList() {
      if(isSetOnly) {
        CkAbort("Cannot add a list for a list which only has set\n");
      }
      if(isList || LisV) {
        isList = 1;
        reset();
      }
      isList = 1;
      list.clear();
      std::copy(vectorList.begin(), vectorList.end(), std::back_inserter(list));
      reset();
    }

    inline void addVector() {
      if(isSetOnly) {
        CkAbort("Cannot add a vector for a list which only has set\n");
      }
      if(!isList || LisV) {
        isList = 0;
        reset();
      }
      isList = 0;
      vectorList.clear();
      vectorList.reserve(list.size());
      std::copy(list.begin(), list.end(), std::back_inserter(vectorList));
      reset();
    }

    // need to rebuild your sortIdx after set union
    inline PeList &operator+(PeList &inList) { 
      // make array large enough for both, paste together
      if(!isSet) {
        CkAbort("Cannot add to a list with no set\n");
      }
      append(inList);
      return *this; 
    }

    inline void deleteOne(int key) {
      if(!isSetOnly) {
        CkAbort("Cannot delete an element from a list which is not setOnly\n");
      }
#if USE_BITVECTOR
      if(allPE[key] == true) {
        allPE[key] = false;
        true_size--;
      }
#else
      setList.erase(key);
#endif
    }

    // need to rebuild your sortIdx after unary set difference
    inline void deleteList(PeList &inList, int _setOnly, int _isList) {
      if(!isSet) {
        CkAbort("Cannot mass delete from a list with no set\n");
      }
      if(inList.isSetOnly) {
        CkAbort("Cannot mass delete using a list with is setOnly\n");
      }

      if(inList.size() == 0) return;

      if(inList.isList) {
        std::list< int >::iterator it = inList.list.begin();
        while(it != inList.list.end()) {
#if USE_BITVECTOR
          if(allPE[*it]) {
            allPE[*it] = false;
            true_size--;
          }
#else
          setList.erase(*it);
#endif
          it++;
        }
      } else {
        std::vector< int >::iterator it = inList.vectorList.begin();
        while(it != inList.vectorList.end()) {
#if USE_BITVECTOR
          if(allPE[*it]) {
            allPE[*it] = false;
            true_size--;
          }
#else
          setList.erase(*it);
#endif
          it++;
        }
      }

      list.clear();
      vectorList.clear();
      isSetOnly = _setOnly;
      if(!isSetOnly) {
        isList = _isList;
        if(isList) {
#if USE_BITVECTOR
          copy_if(allPE, list);
#else
          std::copy(setList.begin(), setList.end(), std::back_inserter(list));
#endif
        } else {
#if USE_BITVECTOR
          vectorList.reserve(true_size);
          copy_if(allPE, vectorList);
#else
          vectorList.reserve(setList.size());
          std::copy(setList.begin(), setList.end(), std::back_inserter(vectorList));
#endif
        }
        LisV = 0;
      }
      reset();
    }
};


/**
 * Hacky solution to passing a PeList to GSpace(0,0) for use in paircalc mapping
 * without actually having to pup the arrays in a PeList
 *
 * Written without any global knowledge of the mapping code
 */
class PeListFactory
{
  public:
    /// Default constructor
    PeListFactory()
    {
      CkAbort("I hate that we end up in this function");
    }

    PeListFactory(int _size): 
      useDefault(0), size(_size) {}

    PeListFactory(int _ndims, int _bdims[10], int _instanceIndex) {
      ndims = _ndims;
      useDefault = 1;
      instanceIndex = _instanceIndex;
      for(int i = 0; i < ndims; i++) {
        bdims[i] = _bdims[i];
      }
    }

    /// Return an appropriately constructed PeList
    PeList* operator() () const
    {
      if (useDefault == 0)
        return new PeList(1, 0, size);
      else if(useDefault == 1) 
        return new PeList(ndims, bdims, instanceIndex);

      return NULL;
    }

  private:
    // Should I return a default PeList or a topoAware PeList
    int useDefault;
    // Parameters used to create a boxy PeList. Refer PeList code for more details
    int size, ndims;
    int instanceIndex, bdims[10];
};
PUPbytes(PeListFactory)

  //@}
#endif
