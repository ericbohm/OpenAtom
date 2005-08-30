#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

void bar(int *,int );
void foo(int **,int );

int main(){

  int nlen = 10;
  int *junk;
  int i;

  foo(&junk,nlen);
  printf("\n");
  for(i=1;i<=nlen;i++){printf("in main : %d %d\n",i,junk[i]);}
  printf("\n");

}


void foo(int **junk,int nlen){
  int i;
  printf("nlen %d\n",nlen);
  scanf("%d",&i);
  *junk = (int *)malloc(nlen*sizeof(int))-1;
  for(i=1;i<=nlen;i++){(*junk)[i]=i;}
  bar(&(*junk)[1],nlen);
  scanf("%d",&i);
}


void bar(int *junk,int nlen){
  int i;
  scanf("%d",&i);
  printf("\n");
  for(i=0;i<nlen;i++){printf("in bar : %d %d\n",i,junk[i]);}
  scanf("%d",&i);
}
