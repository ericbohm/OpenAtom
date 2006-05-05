/** \file FindProcessor.C
 *  Author: Abhinav S Bhatele
 *  Date Created: November 30th, 2005
 *  Topologies supported: 3D Mesh  (on BlueGene/L)
 * 			  3D Torus (on Bluegene/L): Co-processor mode
 * 			  3D Torus (on Bluegene/L): Virtual NOde mode
 */

#include "FindProcessor.h"
#include <charm++.h>

int distance=0;

FindProcessor::FindProcessor()
{
	for(int i=0; i<3;i++)
	{	
		start[i]=0;
		next[i]=0;
	}
}

FindProcessor::FindProcessor(int a[])
{
	for(int i=0; i<3;i++)
		start[i]=a[i];
}


void FindProcessor::findNext(int a[])
{
	if(a[0]==0 && a[1]==0 && a[2]==0)
	{
		next[2]=1;
		cout<<"-------------------------\n";
		cout<<"Distance "<<next[2]<<" starting\n";
		cout<<"-------------------------\n";
		printing(0, 0, 1);
	}	
	else
	{
		if(a[0]<=0 && a[1]<=0 && a[2]<=0)
		{
			for(int i=0; i<3;i++)
				a[i]=-1*a[i];
			
			if(a[2]>0)
			{
				next[0]=a[0]; 
				next[1]=a[1]+1;
				next[2]=a[2]-1;
				printing(next[0], next[1], next[2]);
			}
			else
			{
				if(a[1]>0)
				{
					next[0]=a[0]+1; 
					next[1]=0;
					next[2]=a[1]-1;
					printing(next[0], next[1], next[2]);
				}
				else
				{
					next[0]=0;
					next[1]=0;
					next[2]=a[0]+1;
					cout<<"-------------------------\n";
					cout<<"Distance "<<next[2]<<" starting\n";
					cout<<"-------------------------\n";
					printing(next[0], next[1], next[2]);
				}	
			}	
		}
		else
		{
			if(a[2]>0)
			{
				next[0]=a[0]; 
				next[1]=a[1];
				next[2]=a[2]*-1;
				printing(next[0], next[1], next[2]);
			}	
			else 
			{
				if(a[1]>0)
				{
					next[0]=a[0]; 
					next[1]=a[1]*-1;
					next[2]=a[2]*-1;
					printing(next[0], next[1], next[2]);
				}	
				else
				{
					next[0]=a[0]*-1; 
					next[1]=a[1]*-1;
					next[2]=a[2]*-1;
					printing(next[0], next[1], next[2]);
				}
			}
		}
	}
}

int FindProcessor::findNextInMIter(int a[])
{
        if(a[0]==0 && a[1]==0 && a[2]==0)
        {
                next[2]=1;
                //cout<<"-------------------------\n";
                //cout<<"Distance "<<next[2]<<" starting\n";
                //cout<<"-------------------------\n";
                //printing(0, 0, 1);
        }
        else
        {
                if(a[2]>0)
                {
                        next[0]=a[0];
                        next[1]=a[1]+1;
                        next[2]=a[2]-1;
                        if(next[0]>=nopZ || next[1]>=nopY || next[2]>=nopX)
                        {
                          for(int i=0;i<3;i++)
                            start[i]=next[i];
                          return 2;
                        }
                        //printing(next[0], next[1], next[2]);
                }
		else
                {
                        if(a[1]>0)
                        {
                                next[0]=a[0]+1;
                                next[1]=0;
                                next[2]=a[1]-1;
                                if(next[0]>=nopZ || next[1]>=nopY || next[2]>=nopX)
                                {
                                  for(int i=0;i<3;i++)
                                    start[i]=next[i];
                                  return 2;
                                }
                                //printing(next[0], next[1], next[2]);
                        }
                        else
                        {
                                next[0]=0;
                                next[1]=0;
                                next[2]=a[0]+1;
                                if(next[0]>=nopZ || next[1]>=nopY || next[2]>=nopX)
                                {
                                  for(int i=0;i<3;i++)
                                    start[i]=next[i];
                                  return 2;
                                }
                                //cout<<"-------------------------\n";
                                //cout<<"Distance "<<next[2]<<" starting\n";
                                //cout<<"-------------------------\n";
                                //printing(next[0], next[1], next[2]);
                        }
                }

        }
        count++;
        return 1;
}

int FindProcessor::findNextInMesh(int a[])
{
	int ret = findNextInMIter(a);
	if(count == nopX*nopY*nopZ)
	{
        	cout<<"-------------------------\n";
        	cout<<"No more processors left\n";
        	cout<<"-------------------------\n";
        	CkAbort("inconsistent no. of chares and processors\n");
        	return 0;
	}
	if(ret==1)
		return ret;
	else
	{
		ret=findNextInMIter(start);
		while(ret==2)
		{
			ret=findNextInMIter(start);
		}
		return ret;
	}
}

int FindProcessor::findNextIter(int a[])
{
	int newa[3];
	int negXL=0, negYL=0, negZL=0;
	int posXL=0, posYL=0, posZL=0;
	
	if(nopX%2==0)
		negXL = nopX/2;
	else
		negXL = nopX/2 + 1;
	if(nopY%2==0)
		negYL = nopY/2;
	else
		negYL = nopY/2 + 1;
	if(nopZ%2==0)
		negZL = nopZ/2;
	else
		negZL = nopZ/2 + 1;
	//cout<<"actual -> "<<negXL<<" "<<negYL<<" "<<negZL<<"\n";
	/*int negXL=(int)(ceil(nopX/2.0));
	int negYL=(int)(ceil(nopY/2.0));
	int negZL=(int)(ceil(nopZ/2.0));*/
	
	posXL=nopX/2 + 1;
	posYL=nopY/2 + 1;
	posZL=nopZ/2 + 1;
	
	if(count==nopX*nopY*nopZ)
	{
		cout<<"-------------------------\n";
		cout<<"No more processors left\n";
		cout<<"-------------------------\n";
		CkAbort("inconsistent no. of chares and processors\n");
		return 0;
	}

	if(a[0]==0 && a[1]==0 && a[2]==0)
	{
		next[2]=1;
		distance=1;
		if(next[2]>=posXL)
                {
                        for(int i=0;i<3;i++)
			{
				newa[i]=next[i];
                                start[i]=next[i];
			}
			//printing_sp(start[0], start[1], start[2]);
                        return 2;
                }
                else
                {
			//cout<<"-------------------------\n";
			//cout<<"Distance "<<distance<<" starting\n";
			//cout<<"-------------------------\n";
			//printing(0, 0, 1);
		}
	}	
	else
	{
		if(a[0]<=0 && a[1]<=0 && a[2]<=0)
		{
			for(int i=0; i<3;i++)
				a[i]=-1*a[i];
			
			if(a[2]>0)
			{
				next[0]=a[0]; 
				next[1]=a[1]+1;
				next[2]=a[2]-1;
				if(abs(next[0])>=posZL || abs(next[1])>=posYL || abs(next[2])>=posXL)
				{
					for(int i=0;i<3;i++)
					{
						newa[i]=next[i];
						start[i]=next[i];
					}
					//printing_sp(start[0], start[1], start[2]);
					return 2;
				}
				else
				{
					if(next[2]<0)
						next[2]=next[2]+nopX;
					if(next[1]<0)
						next[1]=next[1]+nopY;
					if(next[0]<0)
						next[0]=next[0]+nopZ;
				}
			}
			else
			{
				if(a[1]>0)
				{
					next[0]=a[0]+1; 
					next[1]=0;
					next[2]=a[1]-1;
					if(abs(next[0])>=posZL || abs(next[1])>=posYL || abs(next[2])>=posXL)
					{
						for(int i=0;i<3;i++)
						{
							start[i]=next[i];
							newa[i]=next[i];
						}
						//printing_sp(start[0], start[1], start[2]);
						return 2;
					}
					else
					{
						if(next[2]<0)
							next[2]=next[2]+nopX;
						if(next[1]<0)
							next[1]=next[1]+nopY;
						if(next[0]<0)
							next[0]=next[0]+nopZ;
					}
				}
				else
				{
					next[0]=0;
					next[1]=0;
					next[2]=a[0]+1;
					if(abs(next[0])>=posZL || abs(next[1])>=posYL || abs(next[2])>=posXL)
					{
						distance=distance+1;
						//cout<<"-------------------------\n";
						//cout<<"Distance "<<distance<<" starting\n";
						//cout<<"-------------------------\n";
						for(int i=0;i<3;i++)
						{
							start[i]=next[i];
							newa[i]=next[i];
						}
						//printing_sp(start[0], start[1], start[2]);
						return 2;
					}
					else
					{
						if(next[2]<0)
							next[2]=next[2]+nopX;
						if(next[1]<0)
							next[1]=next[1]+nopY;
						if(next[0]<0)
							next[0]=next[0]+nopZ;
					}
				}	
			}	
		}
		else
		{
			if(a[2]>0)
			{
				next[0]=a[0]; 
				next[1]=a[1];
				next[2]=a[2]*-1;
				if(compare(next[0], posZL, negZL) || compare(next[1], posYL, negYL) || compare(next[2], posXL, negXL))
				{
					for(int i=0;i<3;i++)
					{
						start[i]=next[i];
						newa[i]=next[i];
					}
					//printing_sp(start[0], start[1], start[2]);
					return 2;
				}
				else
				{
					if(next[2]<0)
						next[2]=next[2]+nopX;
					if(next[1]<0)
						next[1]=next[1]+nopY;
					if(next[0]<0)
						next[0]=next[0]+nopZ;
				}
			}	
			else 
			{
				if(a[1]>0)
				{
					next[0]=a[0]; 
					next[1]=a[1]*-1;
					next[2]=a[2]*-1;
					if(compare(next[0], posZL, negZL) || compare(next[1], posYL, negYL) || compare(next[2], posXL, negXL))
					{
						for(int i=0;i<3;i++)
						{
							start[i]=next[i];
							newa[i]=next[i];
						}
						//printing_sp(start[0], start[1], start[2]);
						return 2;
					}
					else
					{
						if(next[2]<0)
							next[2]=next[2]+nopX;
						if(next[1]<0)
							next[1]=next[1]+nopY;
						if(next[0]<0)
							next[0]=next[0]+nopZ;
					}
				}	
				else
				{
					next[0]=a[0]*-1; 
					next[1]=a[1]*-1;
					next[2]=a[2]*-1;
					if(compare(next[0], posZL, negZL) || compare(next[1], posYL, negYL) || compare(next[2], posXL, negXL))
					{
						for(int i=0;i<3;i++)
						{
							start[i]=next[i];
							newa[i]=next[i];
						}
						//printing_sp(start[0], start[1], start[2]);
						return 2;
					}
					else
					{
						if(next[2]<0)
							next[2]=next[2]+nopX;
						if(next[1]<0)
							next[1]=next[1]+nopY;
						if(next[0]<0)
							next[0]=next[0]+nopZ;
					}
				}
			}
		}
	//printing(next[0], next[1], next[2]);
	}
	count++; //add near all other prints
	return 1;
}

int FindProcessor::findNextInTorus(int a[])
{
    int ret = findNextIter(a);
    while(ret==2)
    {
        ret=findNextIter(start);
    }
    return ret;
}

int FindProcessor::findNextInTorusV(int t, int a[])
{
	int ret;
	if(t==0)
	{
		w=1;
		//printing(w, next[0], next[1], next[2]);
		return 1;
	}
	else
	{
		w=0;
		if(count==nopX*nopY*nopZ) {
                  cout<<"-------------------------\n";
		  cout<<"No more processors left\n";
		  cout<<"-------------------------\n";
		  CkAbort("inconsistent no. of chares and processors\n");
		  return 0;
                }
                
                ret = findNextIter(a);
                while(ret==2)
                {
                    ret=findNextIter(start);
                }
                return ret;
	}	
}

int FindProcessor::compare(int n, int a, int b)
{
	if(n>=0)
		if(n>=a)
			return 1;
		else
			return 0;
	else
		if(abs(n)>=b)
			return 1;
		else
			return 0;
}

/* Various functions which can be used for printing the processor coordinates
 * printSome(int) - print the no. of processors the user wishes
 * printing(int, int, int) - print for bluegene CO mode 
 * printing_sp(int, int, int) - for debugging purposes
 * printing(int, int, int, int) - print for bluegene VN mode 
 */

void FindProcessor::printSome(int n)
{
	int val;
	int array[3];
	if(option==1)
		findNext(start);
	if(option==2)
		findNextInMesh(start);
	if(option==3)
	{
		//count=1;
		if(start[2]>nopX/2)
			array[2]=start[2]-nopX;
		else
			array[2]=start[2];
		if(start[1]>nopY/2)
			array[1]=start[1]-nopY;
		else
			array[1]=start[1];
		if(start[0]>nopZ/2)
			array[0]=start[0]-nopZ;
		else
			array[0]=start[0];
		//cout<<"check: ";
		//printing(array[0],array[1],array[2]);
		val=findNextInTorus(array);
	}
	if(option==4)
	{
		if(start[2]>nopX/2)
			array[2]=start[2]-nopX;
		else
			array[2]=start[2];
		if(start[1]>nopY/2)
			array[1]=start[1]-nopY;
		else
			array[1]=start[1];
		if(start[0]>nopZ/2)
			array[0]=start[0]-nopZ;
		else
			array[0]=start[0];
		//cout<<"check: ";
		//printing(array[0],array[1],array[2]);
		val=findNextInTorusV(w, array);
	}
	while(n>0)
	{
		for(int i=0;i<3;i++)
			start[i]=next[i];
		if(option==1)
			findNext(start);
		if(option==2)
			findNextInMesh(start);
		if(option==3)
		{
			if(val!=0)
			{
				if(start[2]>nopX/2)
					array[2]=start[2]-nopX;
				else
					array[2]=start[2];
				if(start[1]>nopY/2)
					array[1]=start[1]-nopY;
				else
					array[1]=start[1];
				if(start[0]>nopZ/2)
					array[0]=start[0]-nopZ;
				else
					array[0]=start[0];
				//cout<<"check: ";
				//printing(array[0],array[1],array[2]);
				val=findNextInTorus(array);
				/*if(count==nopX*nopY*nopZ-1)
					return;*/
			}
		}
		if(option==4)
		{
			if(val!=0)
			{
				if(start[2]>nopX/2)
					array[2]=start[2]-nopX;
				else
					array[2]=start[2];
				if(start[1]>nopY/2)
					array[1]=start[1]-nopY;
				else
					array[1]=start[1];
				if(start[0]>nopZ/2)
					array[0]=start[0]-nopZ;
				else
					array[0]=start[0];
				//cout<<"check: ";
				//printing(array[0],array[1],array[2]);
				val=findNextInTorusV(w, array);
			}
		}
		n--;
	}
}

void FindProcessor::printing(int a, int b, int c)
{
	if(count<10)
		cout<<"  "<<count<<".    ";
	else if(count<100)
		cout<<" "<<count<<".    ";
	else
		cout<<count<<".    ";
	if(a<0)
		cout<<a;
	else
		cout<<" "<<a;
	if(b<0)
		cout<<" "<<b;
	else
		cout<<"  "<<b;
	if(c<0)
		cout<<" "<<c<<"\n";
	else
		cout<<"  "<<c<<"\n";
	//count++; needs to be placed in other two functions!
}

void FindProcessor::printing_sp(int a, int b, int c)
{
	if(a<0)
		cout<<a;
	else
		cout<<" "<<a;
	if(b<0)
		cout<<" "<<b;
	else
		cout<<"  "<<b;
	if(c<0)
		cout<<" "<<c<<"\n";
	else
		cout<<"  "<<c<<"\n";
}

void FindProcessor::printing(int w, int a, int b, int c)
{
	if(count<10)
		cout<<"  "<<count<<".    ";
	else if(count<100)
		cout<<" "<<count<<".    ";
	else
		cout<<count<<".    ";
	if(w<0)
		cout<<w;
	else
		cout<<" "<<w;
	if(a<0)
		cout<<" "<<a;
	else
		cout<<"  "<<a;
	if(b<0)
		cout<<" "<<b;
	else
		cout<<"  "<<b;
	if(c<0)
		cout<<" "<<c<<"\n";
	else
		cout<<"  "<<c<<"\n";
	//count++; needs to be placed in other two functions!
}

/* int main which can be used for debugging
 * functions in this file

int main(int argc, char*argv[])
{
	FindProcessor fp=FindProcessor();
	if(strcmp(argv[1],"-3Dmesh")==0)
		fp.option=1;
	if(strcmp(argv[1],"-Bluegene")==0)
		fp.option=2;
	for(int i=2;i<5;i++)
		fp.start[i-2]=atoi(argv[i]);
	if(strcmp(argv[1],"-Torus")==0)
	{
		fp.option=3;
		fp.nopX=atoi(argv[2]);
		fp.nopY=atoi(argv[3]);
		fp.nopZ=atoi(argv[4]);
		for(int i=5;i<8;i++)
			fp.start[i-5]=atoi(argv[i]);
	}
	if(strcmp(argv[1],"-TorusV")==0)
	{
		fp.option=4;
		fp.nopX=atoi(argv[2]);
		fp.nopY=atoi(argv[3]);
		fp.nopZ=atoi(argv[4]);
		fp.w=atoi(argv[5]);
		for(int i=6;i<9;i++)
			fp.start[i-6]=atoi(argv[i]);
	}
	fp.printSome(250);
	//fp.findNext(fp.start);
	//fp.findNextInMesh(fp.start);
	//sfp.findNextInTorus(fp.start);	
} //end main

 *
 */


