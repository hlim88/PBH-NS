#include "variableDef.h"
#include "boundaryConditions.h"

//Interpolate functions even in r as r goes to zero
void interpolateEven(double *q){
	q[0] = (3311.0*q[1]-2413.0*q[2]+851.0*q[3]-122.0*q[4])/1627.0;
}

//Interpolate functions odd in r as r goes to zero
void interpolateOdd(double *q){
	q[0] = (35819.0*q[1]-16777.0*q[2]+4329.0*q[3]-488.0*q[4])/36833.0;
}

void setOriginBoundaryCondition(double *consVar){

	

	int interpLength = 5;
	double D[interpLength];
	double S[interpLength];
	double tau[interpLength];
	
	int i;
	for(i=1; i<interpLength; i++){
		D[i] = consVar[i*numVariables+CONS_DENSITY];
		S[i] =	0.5*(consVar[i*numVariables+PI] - consVar[i*numVariables+PHI]);
		tau[i] = 0.5*(consVar[i*numVariables+PI] + consVar[i*numVariables+PHI]);
	}
	
	interpolateEven(D);
	interpolateEven(tau);
	interpolateOdd(S);

	consVar[CONS_DENSITY] = D[0];
	consVar[PI] = tau[0] + S[0];
	consVar[PHI] = tau[0] - S[0];
}

/*
Use two ghost cells at the outer boundary
*/
void setOuterBoundaryCondition(double *consVar, int L){

	if(L<3) return;

	int i,j;
	for(i=L-2; i<L; i++){
		for(j=0; j<numVariables; j++){
			consVar[i*numVariables+j] = consVar[(L-3)*numVariables+j];
		}
	}

}
/* 
Sommerfeld boundary condition 
*/

void sommerfeldBoundary(double *consVar){

  const double f_falloff[3] = { 1.0, 1.0, 1.0};

  const double f_asymptotic[3] = {1.0, 0.0, 0.0}; 
  
  int i;
  for (i = 0; i<numVariables; i++){
      consVar[i] = f_falloff[i] * (consVar[i] - f_asymptotic[i]);
  }


}
