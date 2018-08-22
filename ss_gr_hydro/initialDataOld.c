#include <math.h>
#include "conservedPrimative.h"
#include "initialData.h"
#include "flux.h"
#include <stdio.h>
/*
Set inital data.  All the input arrays are assumed to be uninitialized with the exception of r.
*/
void getInitialData(double *consVar, double *primVar, double *a, double *alpha, double *flux, double *fluxAlt, double *source, double *rGravity, double *rFluid, int length){

	double aFluid[length];
	double alphaFluid[length];
	double rNot = 0.1;//0.1;
	double delta =100.0;
	double rho = 0.1;
	double energy = 1.0e-6;
	double v = -0.01;
	int i;
		
	for(i=0; i<length; i++){
		a[i] = 1.0;
		alpha[i] = 1.0;
		
		primVar[i*numVariables+SPECIFIC_ENERGY]=energy;//exp(-1.0*(r[i]-r0)*(r[i]-r0)/(0.01))+energy0;
		/*if(r[i]<r0){
			primVar[i*numVariables+SPECIFIC_ENERGY]= 1.0;
		}
		*/
		primVar[i*numVariables+REST_DENSITY]= rho*exp(-1.0*(rFluid[i]-rNot)*(rFluid[i]-rNot)/(0.5*0.5));
		primVar[i*numVariables+VELOCITY]= v*(rFluid[i]/rNot);

		primativeToConserved((primVar+i*numVariables), a[i], (consVar+i*numVariables));
		floorPrimative((primVar+i*numVariables));
	}
	
	solveHamiltonianConstraint(consVar, rGravity, length, a);
	//Interpolate metric to center of cells where fluid variables are defined
	for(i=0; i<length-1; i++){
		aFluid[i] = 0.5*(a[i]+a[i+1]);
	}
	aFluid[length-1] = aFluid[length-2];
	
	getNumericalFlux(consVar, primVar, aFluid, rFluid, length, flux, fluxAlt);

	solveSlicingCondition(consVar, primVar, a, rGravity, length, alpha);
	//Interpolate metric to center of cells where fluid variables are defined
	for(i=0; i<length-1; i++){
		alphaFluid[i] = 0.5*(alpha[i]+alpha[i+1]);
	}
	alphaFluid[length-1] = alphaFluid[length-2];

	if(fluxAlt==NULL){
		getSourceArray(primVar, consVar, aFluid, alphaFluid, rFluid, length, source);
	} else {
		getSourceArrayAlt(primVar, consVar, aFluid, alphaFluid, rFluid, length, source);
	}

}

