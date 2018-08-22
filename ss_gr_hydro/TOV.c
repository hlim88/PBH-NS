#define CONSTANT_PI 3.1415926535897932384626433
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "TOV.h"
#include "rungeKutta.h"

enum{
	MASS,
	PRESSURE,
	LNALPHA
};		

//The TOV differential equations
void derivative(double r, double* inputValues, double* deriv, double(*getDensity)(double)){

	double density = getDensity(inputValues[PRESSURE]);
	deriv[MASS] = 4.0*CONSTANT_PI*r*r*density;
	deriv[PRESSURE] = -1.0*(density +inputValues[PRESSURE])*(inputValues[MASS]/r+4.0*CONSTANT_PI*r*r*inputValues[PRESSURE])/(r-2.0*inputValues[MASS]);
	deriv[LNALPHA] = (inputValues[MASS]/r+4.0*CONSTANT_PI*r*r*inputValues[PRESSURE])/(r-2.0*inputValues[MASS]);
}

/*
Note by density we mean (1+specificEnergy)*restDensity.
The getDensity function takes pressure and returns density.
The value returned it the radius of the solution.
*/
double getTOVsolution(double centralPressure, double atmPressure, double (*getDensity)(double), double *pressure, double* a, double* alpha, double* r, int length){

	
	

	int dim = 3;	

	double currValue[dim];
	double derivValue[dim];
	double* currDeriv = derivValue;
	
	//Set boundary conditions
	currValue[MASS] = 0.0;
	currValue[PRESSURE] = centralPressure;
	currValue[LNALPHA] = 0.0; //Arbitrary
	derivValue[MASS] = 0.0;
	derivValue[PRESSURE] = 0.0;
	derivValue[LNALPHA] = 0.0;

	a[0]=1.0;
	pressure[0]=centralPressure;
	alpha[0] = 1.0;	


	//Integrate equations from center outwards until pressure drops to atmosphere pressure
	int i;
	for(i=1; i<length; i++){
		fourthOrderRKStep(currValue, r[i-1], currValue, currDeriv, dim, r[i]-r[i-1], &derivative, getDensity);	

		currDeriv = NULL;
		if(currValue[PRESSURE]<atmPressure){
			break;
		}
		pressure[i] = currValue[PRESSURE];
		a[i] = sqrt(r[i]/(r[i]-2.0*currValue[MASS]));
		alpha[i] = exp(currValue[LNALPHA]);
	
	}
	double totMass = currValue[MASS];
	double totRadius = r[i];

	//Go back and correct so that alpha matches Schwarschild at boundary
	double factor = sqrt(1.0-2.0*totMass/r[i])/exp(currValue[LNALPHA]);
	int j;
	for(j=i-1; j>=0; j--){
		alpha[j]*=factor;
	} 

	//Now match a Schwarschild metric for the outer radius
	for(; i<length; i++){
		pressure[i] = atmPressure;
		a[i] = sqrt(r[i]/(r[i]-2.0*totMass));
		alpha[i] = 1.0/a[i];		
	}
	
	return totRadius;
}

