#include "matrixInverse.h"
#include <stdio.h>
/*
Compute the inverse of a 3x3 matrix 
*/
void matrix3Inverse(double A[3][3], double invA[3][3]){

	double BLdet = A[1][0]*A[2][1]-A[1][1]*A[2][0];
	double BMdet = A[1][2]*A[2][0]-A[1][0]*A[2][2];
	double BRdet = A[1][1]*A[2][2]-A[1][2]*A[2][1];
	double det = A[0][0]*BRdet+A[0][1]*BMdet+A[0][2]*BLdet;

	if(det==0){
//		printf("Error matrix is noninvertible! det=%E \n",det);
		int i,j;
		for(i=0; i<3; i++){
			for(j=0;j<3; j++){
				invA[i][j]=0.0;
			}
		}
	} else {

		invA[0][0]=BRdet/det;
		invA[0][1]=(A[0][2]*A[2][1]-A[2][2]*A[0][1])/det;
		invA[0][2]=(A[0][1]*A[1][2]-A[1][1]*A[0][2])/det;
		invA[1][0]=BMdet/det;
		invA[1][1]=(A[0][0]*A[2][2]-A[2][0]*A[0][2])/det;
		invA[1][2]=(A[0][2]*A[1][0]-A[1][2]*A[0][0])/det;
		invA[2][0]=BLdet/det;
		invA[2][1]=(A[0][1]*A[2][0]-A[2][1]*A[0][0])/det;
		invA[2][2]=(A[0][0]*A[1][1]-A[1][0]*A[0][1])/det;
	}
}
