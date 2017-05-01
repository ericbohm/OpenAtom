#include<stdio.h>
#include "ckcomplex.h"

bool compare_files(int ir){
  char file1[100];
  char file2[100];

  sprintf(file1,"../serial_yale/state_0_%d_%d", 0, 0);
  sprintf(file2,"../input/state_0_%d_%d", 0, 0);

	int states_size = 8;
  int size = states_size;
	FILE* fp1 = fopen(file1, "r");
	FILE* fp2 = fopen(file2, "r");
//  printf("\nOpening files %s and %s for verification\n", file1, file2);
	complex vector1[size];
	complex vector2[size];
	for (int ic=0; ic<size; ic++){
		fscanf(fp1, "%lg  %lg", &vector1[ic].re, &vector1[ic].im);
		fscanf(fp2, "%lg  %lg", &vector2[ic].re, &vector2[ic].im);
	}
  fclose(fp1);
  fclose(fp2);
  for (int ic=0; ic<size; ic++){
		printf("\n re = %lg %lg\n",vector1[ic].re,vector2[ic].re);
    if(vector1[ic].re != vector2[ic].re || vector1[ic].im != vector2[ic].im){
      printf("\nAborting due to mismatch in PMatrix row #%d\n", ir);
      return false;
    }
	}
  return true;
}

int main(){
//States: kpoint, state_index
  int size = 1728;
	int ir = 0;
//  for(int ir=0;ir<size;ir++)
//    if(!compare_files(ir)) break;
	compare_files(ir);
	return 0;
}
