/* 
   3 x 3 matrix and vector operation
*/

#include <iostream>
#include <cstdlib>
#include <cmath>
using namespace std;

//----------------------------
// calculate inverse matrix (double)
//----------------------------
void inverse(double m[][3], double imtrx[][3]){
    
    double a[3][3];
    
    /* compute matrix of cofactors */
    a[0][0] =  m[1][1]*m[2][2] - m[1][2]*m[2][1];
    a[1][0] = -m[1][0]*m[2][2] + m[1][2]*m[2][0];
    a[2][0] =  m[1][0]*m[2][1] - m[1][1]*m[2][0];
    a[0][1] = -m[0][1]*m[2][2] + m[0][2]*m[2][1];
    a[1][1] =  m[0][0]*m[2][2] - m[0][2]*m[2][0];
    a[2][1] = -m[0][0]*m[2][1] + m[0][1]*m[2][0];
    a[0][2] =  m[0][1]*m[1][2] - m[0][2]*m[1][1];
    a[1][2] = -m[0][0]*m[1][2] + m[0][2]*m[1][0];
    a[2][2] =  m[0][0]*m[1][1] - m[0][1]*m[1][0];
    
    /* compute determinant */
    double det;
    det = m[0][0]*a[0][0] + m[0][1]*a[1][0] + m[0][2]*a[2][0];
    
    double del = 0.00001;
    if ( abs(det) < del) {
        cout << "ERROR: Determinant of matrix is zero. Check your matrix" << endl;
    }
    
    double me;
    
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            me = a[i][j] / det;
            imtrx[i][j] = me;
        }
    }
}//-----end function


//-------------------------------
// matrix vector multiplication (double)
//-------------------------------
void matvec3(double mtrx[][3], double *v, double *a){
    
    for (int i=0; i<3; i++) {
        a[i]=0;
        for (int j=0; j<3; j++) {
            a[i] += mtrx[i][j]*v[j];
        }
    }
    
}//-----end function



//------------------------------
// dot product of two vectors (double)
//------------------------------
double dot_product(double v1[3], double v2[3]){
  double result;
  for(int i=0; i<3; i++){
    result += v1[i]*v2[i];
  }
  return result;
}//-----end function



