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
  int *stepTable;
 public:

    IntMap4(int keyW, int keyX, int keyY, int keyZ, int step) 
           : keyWmax(keyW), keyXmax(keyX), keyYmax(keyY), keyZmax(keyZ), keyStep(step)
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
	stepTable= new int [keyXmax*keyStep];
	for(int s=0;s<keyXmax*keyStep;s++)
	  stepTable[s]=s/keyStep;
      }
    ~IntMap4(){
      if(Map!=NULL)
	{
	  for(int w=0;w<keyWmax;w++)
	    {
	      for(int x=0;x<keyXmax;x++)
		{
		  for(int y=0;y<keyYmax;y++)
		    delete [] Map[w][x][y];
		  delete [] Map[w][x];
		}
	      delete [] Map[w];
	    }
	  delete [] Map;
	  Map=NULL;
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
	stepTable= new int [keyXmax*keyStep];
	for(int s=0;s<keyXmax*keyStep;s++)
	  stepTable[s]=s/keyStep;
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
	  if(p.isUnpacking())
	    stepTable= new int[keyXmax*keyStep];
	  PUParray(p,stepTable,keyXmax*keyStep);

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
#define USE_INT_MAP_MATH
#ifdef USE_INT_MAP_MATH
       return(Map[W][X/keyStep][Y/keyStep][Z]);
#else
       return(Map[W][stepTable[X]][stepTable[Y]][Z]);
#endif
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
		CkPrintf("%d %d %d %d %d \n",w,x,y,z, get(w,x*keyStep,y*keyStep,z));
      }
    IntMap4(){keyWmax=0;keyXmax=0; keyYmax=0, keyZmax=0; keyStep=1; Map=NULL;}
};

class IntMap3 {
 private:
  int ***Map;
  int keyXmax;
  int keyYmax;
  int keyZmax;

 public:

    IntMap3(int keyX, int keyY, int keyZ, int step) 
           :  keyXmax(keyX), keyYmax(keyY), keyZmax(keyZ)
      {
	Map=new int**[keyXmax];
	for(int x=0;x<keyXmax;x++)
	  {
	    Map[x]= new int*[keyYmax];
	    for(int y=0;y<keyYmax;y++)
	      {
		Map[x][y]= new int[keyZmax];
		memset(Map[x][y],-1,keyZmax*sizeof(int));
	      }
	  }
      }
    ~IntMap3(){
      if(Map!=NULL)
	{
	  for(int x=0;x<keyXmax;x++)
	    {
	      for(int y=0;y<keyYmax;y++)
		delete [] Map[x][y];
	      delete [] Map[x];
	    }
	  Map=NULL;
	}
      
    }

    void buildMap(int keyX=1, int keyY=1, int keyZ=1, int step=1)
      {
	CkAssert(keyX>0);
	CkAssert(keyY>0);
	CkAssert(keyZ>0);
	keyXmax=keyX;
	keyYmax=keyY;
	keyZmax=keyZ;
	Map=new int**[keyXmax];
	for(int x=0;x<keyXmax;x++)
	  {
	    Map[x]= new int*[keyYmax];
	    for(int y=0;y<keyYmax;y++)
	      {
		Map[x][y]= new int[keyZmax];
		memset(Map[x][y],-1,keyZmax*sizeof(int));
	      }
	  }
      }
    void pup(PUP::er &p)
      {
	  p|keyXmax;
	  p|keyYmax;
	  p|keyZmax;
	  if(p.isUnpacking())
	    Map=new int**[keyXmax];
	  for(int x=0;x<keyXmax;x++)
	    {
	      if(p.isUnpacking())
		Map[x]= new int*[keyYmax];
	      for(int y=0;y<keyYmax;y++)
		{
		  if(p.isUnpacking())
		    Map[x][y]= new int[keyZmax];
		  PUParray(p,Map[x][y],keyZmax);
		}
	    }
      }
    inline int getXmax(){return(keyXmax);}
    inline int getYmax(){return(keyYmax);}
    inline int getZmax(){return(keyZmax);}
    inline int get(int X, int Y, int Z)  {
      return(Map[X][Y][Z]);
    }
    //    inline int &put(int W, int X, int Y, int Z){return(&(Map[W][X/keyStep][Y/keyStep][Z]));}
    inline void set(int X, int Y, int Z, int value){
      CkAssert(X<keyXmax);
      CkAssert(Y<keyYmax);
      CkAssert(Z<keyZmax);
      CkAssert(CkNumPes()>value);
      Map[X][Y][Z]=value;
    }
    int getCentroid();
    void dump()
      {
	  for(int x=0;x<keyXmax;x++)
	    for(int y=0;y<keyYmax;y++)
	      for(int z=0;z<keyZmax;z++)
		CkPrintf("%d %d %d %d \n",x,y,z, get(x,y,z));
      }
    IntMap3(){keyXmax=0; keyYmax=0, keyZmax=0;  Map=NULL;}
};

class IntMap2on2 {
 private:
  int **Map;
  int keyXmax;
  int keyYmax;
 public:
    IntMap2on2(){keyXmax=0; keyYmax=0; Map=NULL;}
    ~IntMap2on2(){
      if(Map!=NULL)
	{
	  for(int x=0;x<keyXmax;x++)
	    delete [] Map[x];
	  delete [] Map;
	  Map=NULL;
	}
    }
    IntMap2on2(int keyX, int keyY): keyXmax(keyX), keyYmax(keyY) 
      {
	Map= new int*[keyXmax];
	for(int x=0;x<keyXmax;x++)
	  {
	    Map[x]= new int[keyYmax];
	    memset(Map[x],-1,keyYmax*sizeof(int));
	  }
	
      }
    void buildMap(int keyX=1, int keyY=1)
      {
	CkAssert(keyX>0);
	CkAssert(keyY>0);
	keyXmax=keyX;
	keyYmax=keyY;
	Map= new int*[keyXmax];
	for(int x=0;x<keyXmax;x++)
	  {
	    Map[x]= new int[keyYmax];
	    memset(Map[x],-1,keyYmax*sizeof(int));
	  }
	
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
    int getCentroid();
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


class IntMap2on1 {
 private:
  int *Map;
  int keyXmax;
  int keyYmax;
  int keymax;
 public:
    ~IntMap2on1(){
      if(keymax>0 && Map!=NULL)
	{
	  delete [] Map;
	  Map=NULL;
	}
    }
    IntMap2on1(){keyXmax=0; keyYmax=0; keymax=0;Map=NULL;}
    IntMap2on1(int keyX, int keyY): keyXmax(keyX), keyYmax(keyY) 
      {
	CkAssert(keyXmax>0);
	CkAssert(keyYmax>0);
	CkAssert(keyXmax>keyYmax);
	keymax=keyXmax*keyYmax;
	Map= new int[keymax];
      }
    void buildMap(int keyX=1, int keyY=1)
      {
	CkAssert(keyX>0);
	CkAssert(keyY>0);
	keyXmax=keyX;
	keyYmax=keyY;
	CkAssert(keyXmax>keyYmax);
	keymax=keyXmax*keyYmax;
	Map= new int[keyXmax];
      }
    void pup(PUP::er &p)
      {
	  p|keyXmax;
	  p|keyYmax;
	  p|keymax;
	  if(p.isUnpacking())
	    Map=new int[keymax];
	  PUParray(p,Map,keymax);

      }
    inline int getXmax(){return(keyXmax);}
    inline int getYmax(){return(keyYmax);}
    inline int getmax(){return(keymax);}
    inline int get(int X, int Y)  {
      /*
      CkAssert(X<keyXmax);
      CkAssert(Y<keyYmax);
      */
      return(Map[Y*keyXmax +X]);
    }
    //    inline &int put(int X, int Y){&(Map[X][Y]);}
    inline void set(int X, int Y, int value){
      CkAssert(CkNumPes()>value);
      CkAssert(X<keyXmax);
      CkAssert(Y<keyYmax);
      Map[Y*keyXmax +X]=value;
    }
    void dump()
      {
	  for(int x=0;x<keyXmax;x++)
	    for(int y=0;y<keyYmax;y++)
		CkPrintf("%d %d %d \n",x,y, get(x,y));
      }
};


#endif
