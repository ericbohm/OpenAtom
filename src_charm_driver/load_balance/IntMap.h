/** \file IntMap.h
 *  Author: Eric J Bohm
 *  Date Created: June 4th, 2006
 *
 * Simple abstraction to replace CkHashtable with a straight array of
 * integers when we know we can perfectly hash the keys into an int
 * and only want an int back.  Meaning storage is just int
 * map[numKeys] and lookup is a constant time array offset.  Array
 * offset is an int multiply and int add.  Should be significantly
 * speedier than CkHashtable for the conditions under which it is
 * appropriate.
 *
 * 
 * IntMap4 is for 4D paircalc
 *
 * IntMap2 is the more standard 2D
 *
 * This is a much less flexible object than CkHashtable and probably
 * only suitable for use within procmaps.
 */

#ifndef _INTMAP_H_
#define _INTMAP_H_


class IntMap4 {
 private:
  int ****Map;
  int keyWmax;
  int keyXmax;
  int keyYmax;
  int keyZmax;
  int keyStep;
 public:

    IntMap4(int keyW, int keyX, int keyY, int keyZ, int step) 
           : keyWmax(keyW), keyXmax(keyX) keyYmax(keyY), keyZmax(keyZ), keyStep(step)
      {
	Map=new int***[keyWmax];
	for(int w=0;w<keyWmax;w++)
	  {
	    Map[w]= new int**[keyXmax];
	    for(int x=0;x<keyXmax;x++)
	      {
		Map[w][x]= new int*[keyYmax];
		for(int y=0;y<keyYmax;y++)
		  Map[w][x][y]= new int[keyZmax];
	      }
	
	  }
      }
    void buildMap(int keyW=1, int keyX=1, int keyY=1, int keyZ=1, int step=1)
      {
	CkAssert(keyW>0);
	CkAssert(keyX>0);
	CkAssert(keyY>0);
	CkAssert(keyZ>0);
	keyWmax=keyW;
	keyXmax=keyX;
	keyYmax=keyY;
	keyZmax=keyZ;
	keyStep=step;
	Map=new int***[keyWmax];
	for(int w=0;w<keyWmax;w++)
	  {
	    Map[w]= new int**[keyXmax];
	    for(int x=0;x<keyXmax;x++)
	      {
		Map[w][x]= new int*[keyYmax];
		for(int y=0;y<keyYmax;y++)
		  Map[w][x][y]= new int[keyZmax];
	      }
	  }
      }
    void pup(PUP::er &p)
      {
	  p|keyWmax;
	  p|keyXmax;
	  p|keyYmax;
	  p|keyZmax;
	  p|keyStep;
	  if(p.isUnpacking())
	    Map=new int***[keyWmax];
	  for(int w=0;w<keyWmax;w++)
	    {
	      if(p.isUnpacking())
		Map[w]= new int**[keyXmax];
	      for(int x=0;x<keyXmax;x++)
		{
		  if(p.isUnpacking())
		    Map[w][x]= new int*[keyYmax];
		  for(int y=0;y<keyYmax;y++)
		    {
		      if(p.isUnpacking())
			Map[w][x][y]= new int[keyZmax];
		      PUParray(p,Map[w][x][y],keyZmax);
		    }
		}
	    }
      }
    inline int getWmax(){return(keyWmax);}
    inline int getXmax(){return(keyXmax);}
    inline int getYmax(){return(keyYmax);}
    inline int getZmax(){return(keyZmax);}
    inline int get(int W, int X, int Y, int Z)  {
      /*
      CkAssert(W<keyWmax);
      CkAssert(X/keyStep<keyXmax);
      CkAssert(Y/keyStep<keyYmax);
      CkAssert(Z<keyZmax);
      */
      return(Map[W][X/keyStep][Y/keyStep][Z]);
    }
    //    inline int &put(int W, int X, int Y, int Z){return(&(Map[W][X/keyStep][Y/keyStep][Z]));}
    inline void set(int W, int X, int Y, int Z, int value){
      CkAssert(W<keyWmax);
      CkAssert(X/keyStep<keyXmax);
      CkAssert(Y/keyStep<keyYmax);
      CkAssert(Z<keyZmax);
      CkAssert(CkNumPes()>value);
      Map[W][X/keyStep][Y/keyStep][Z]=value;
    }
    void dump()
      {
	for(int w=0;w<keyWmax;w++)
	  for(int x=0;x<keyXmax;x++)
	    for(int y=0;y<keyYmax;y++)
	      for(int z=0;z<keyZmax;z++)
		CkPrintf("%d %d %d %d %d \n",w,x,y,z, get(w,x,y,z));
      }
    IntMap4(){keyWmax=0;keyXmax=0; keyYmax=0, keyZmax=0; keyStep=1; Map=NULL;}
};

class IntMap2 {
 private:
  int **Map;
  int keyXmax;
  int keyYmax;
  int keyStep;
 public:
    IntMap2(){keyXmax=0; keyYmax=0; Map=NULL;}
    IntMap2(int keyX, int keyY): keyXmax(keyX), keyYmax(keyY) 
      {
	Map= new int*[keyXmax];
	for(int x=0;x<keyXmax;x++)
	  Map[x]= new int[keyYmax];
	
      }
    void buildMap(int keyX=1, int keyY=1)
      {
	CkAssert(keyX>0);
	CkAssert(keyY>0);
	keyXmax=keyX;
	keyYmax=keyY;
	Map= new int*[keyXmax];
	for(int x=0;x<keyXmax;x++)
	  Map[x]= new int[keyYmax];
      }
    void pup(PUP::er &p)
      {
	  p|keyXmax;
	  p|keyYmax;
	  if(p.isUnpacking())
	    Map=new int*[keyXmax];
	  for(int x=0;x<keyXmax;x++)
	    {
	      if(p.isUnpacking())
		Map[x]= new int[keyYmax];
	      PUParray(p,Map[x], keyYmax);
	    }
      }
    inline int getXmax(){return(keyXmax);}
    inline int getYmax(){return(keyYmax);}
    inline int get(int X, int Y)  {
      /*
      CkAssert(X<keyXmax);
      CkAssert(Y<keyYmax);
      */
      return(Map[X][Y]);
    }
    //    inline &int put(int X, int Y){&(Map[X][Y]);}
    inline void set(int X, int Y, int value){
      CkAssert(CkNumPes()>value);
      CkAssert(X<keyXmax);
      CkAssert(Y<keyYmax);
      Map[X][Y]=value;
    }
    void dump()
      {
	  for(int x=0;x<keyXmax;x++)
	    for(int y=0;y<keyYmax;y++)
		CkPrintf("%d %d %d \n",x,y, get(x,y));
      }
};


#endif
