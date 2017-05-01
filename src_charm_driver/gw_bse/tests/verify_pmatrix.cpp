#include<stdio.h>
#include "ckcomplex.h"

bool compare_files(int ir){
  char file1[100];
  char file2[100];

  sprintf(file1,"../serial_yale/PMatrix_row%d", ir);
  sprintf(file2,"../input/PMatrix_row%d", ir);

  int size = 1728;
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
  for (int ic=0; ic<size; ic++)
    if(vector1[ic].re != vector2[ic].re || vector1[ic].im != vector2[ic].im){
      printf("\nAborting due to mismatch in PMatrix row #%d\n", ir);
      return false;
    }
  return true;
}

int main(){
  int size = 1728;
  for(int ir=0;ir<size;ir++)
    if(!compare_files(ir)) break;

	return 0;
}
