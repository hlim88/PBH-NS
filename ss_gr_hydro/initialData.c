#include <math.h>
#include <stdio.h>
#include "conservedPrimitive.h"
#include "initialData.h"
#include "gravity.h"
#include "equationOfState.h"
#include "TOV.h"
#include "initialVelocity.h"

//Density (not rest density) as a function of pressure using EOS and polytropic relation
double getDensity(double pressure){
	if(pressure<0.0){
		pressure = ATM_PRESSURE;
	}
	double restDensity = getPolytropicRestDensity(pressure);
	double specificEnergy = getSpecificEnergy(restDensity, pressure);
	double density = restDensity*(1.0+specificEnergy);
	return density;
}


// Aug.22.2018 
// Old initial data
/*
Set inital data.  All the input arrays are assumed to be empty with the exception of r.
velocityAmp is parameter determining magnitude of velcity profile.
length is number of elements in a and alpha.  consVar and primVar have length-1 elements.
*/
double getInitialData(double centralPressure, double velocityAmp, double *consVar, double *primVar, double *a, double *alpha, double *rGravity, double *rFluid, int length){

	double pressure[length];
	double atmPressure = ATM_PRESSURE;
	double xatm = X_ATM; //Parameter determining the radius (in units of the star radius) above which the velocity is zero
	

	printf("central density=%E\n", getDensity(centralPressure));

	//Calculate TOV solution
	double starRadius = getTOVsolution(centralPressure, atmPressure, &getDensity,  pressure, a, alpha, rGravity, length);

	//Interpolate functions to center of cells where fluid variables are defined
	double aFluid[length-1];
	double alphaFluid[length-1];
	double pressureFluid[length-1];
	int i,j;
	for(i=0; i<length-1; i++){
		pressureFluid[i] = 0.5*(pressure[i]+pressure[i+1]);
	}

	double coordVelocity[length-1];
	
	for(i=0; i<length-1; i++){
		double restDensity = getPolytropicRestDensity(pressureFluid[i]);
		double specificEnergy = getSpecificEnergy(restDensity, pressureFluid[i]);
		primVar[i*numVariables+SPECIFIC_ENERGY]=specificEnergy;
		primVar[i*numVariables+REST_DENSITY]= restDensity;

		//Add a velocity profile to the solution
		double x = rFluid[i]/starRadius;
		if(x<1.0){
			coordVelocity[i] = 0.5*velocityAmp*(x*x*x-3.0*x);
		} else if(x<xatm ){	
			coordVelocity[i] = 0.5*velocityAmp*(-4.0*x*x*x+6.0*x*x*(1.0+xatm)-12.0*xatm*x+2.0*xatm*xatm*(3.0-xatm))/(xatm-1.0)/(xatm-1.0)/(xatm-1.0);
		} else {
			coordVelocity[i] = 0.0;
		}
		
	}

	//Use coordinate velocity to solve for a and alpha
	solveMetricCoordVelocity(coordVelocity, primVar, rGravity, length, a, alpha);

	
	//Interpolate metric to center of cells where fluid variables are defined
	for(i=0; i<length-1; i++){
		aFluid[i] = 0.5*(a[i]+a[i+1]);
	}



	for(i=0; i<length-1; i++){
		floorPrimitive((primVar+i*numVariables));
		primitiveToConserved((primVar+i*numVariables), aFluid[i], (consVar+i*numVariables));
	}
	

	//Because of the velocity profile we need to recalculate the metric to be consistent with the BC	
	solveHamiltonianConstraint(consVar, rGravity, length, a);
	
	//Interpolate metric to center of cells where fluid variables are defined
	for(i=0; i<length-1; i++){
		aFluid[i] = 0.5*(a[i]+a[i+1]);
	}
	
	for(i=0; i<length-1; i++){
		primitiveToConserved((primVar+i*numVariables), aFluid[i], (consVar+i*numVariables));
	}

	//Last two cells are ghost cells
	for(i=length-3; i<length-1; i++){ 
		for(j=0; j<numVariables; j++) consVar[i*numVariables +j] = consVar[(length-4)*numVariables+j];
	}

	//Impose outerboundary condition a(r_max)*alpha(r_max)=1
	alpha[length-1] = 1.0/a[length-1];

	solveSlicingCondition(consVar, primVar, a, rGravity, length, alpha);
	
	return starRadius;
}

// HL : New initial data based on HP coordinate

double getID2(double centralPressure, double velocityAmp, 
		      double *consVar, double *primVar, double *a, 
                      double *alpha, double *rGravity, double *rFluid, 
                      int length){

	double pressure[length];
	double atmPressure = ATM_PRESSURE;
	double xatm = X_ATM; //Parameter determining the radius (in units of the star radius) above which the velocity is zero
	

	printf("central density=%E\n", getDensity(centralPressure));

	//Calculate TOV solution
	double starRadius = getTOVsolution(centralPressure, atmPressure, &getDensity,  pressure, a, alpha, rGravity, length);

	//Interpolate functions to center of cells where fluid variables are defined
	double aFluid[length-1];
	double alphaFluid[length-1];
	double pressureFluid[length-1];
	int i,j;
	for(i=0; i<length-1; i++){
		pressureFluid[i] = 0.5*(pressure[i]+pressure[i+1]);
	}

	double coordVelocity[length-1];
	
	for(i=0; i<length-1; i++){
		double restDensity = getPolytropicRestDensity(pressureFluid[i]);
		double specificEnergy = getSpecificEnergy(restDensity, pressureFluid[i]);
		primVar[i*numVariables+SPECIFIC_ENERGY]=specificEnergy;
		primVar[i*numVariables+REST_DENSITY]= restDensity;

		//Add a velocity profile to the solution
		double x = rFluid[i]/starRadius;
		if(x<1.0){
			coordVelocity[i] = 0.5*velocityAmp*(x*x*x-3.0*x);
		} else if(x < xatm){	
			coordVelocity[i] = 0.5*velocityAmp*(-4.0*x*x*x+6.0*x*x*(1.0+xatm)
                                         -12.0*xatm*x+2.0*xatm*xatm*(3.0-xatm))/(xatm-1.0)
                                         /(xatm-1.0)/(xatm-1.0);
		} else {
			coordVelocity[i] = 0.0;
		}
		
	}

	//Use coordinate velocity to solve for a and alpha
	solveMetricCoordVelocity(coordVelocity, primVar, rGravity, length, a, alpha);

	
	//Interpolate metric to center of cells where fluid variables are defined
	for(i=0; i<length-1; i++){
		aFluid[i] = 0.5*(a[i]+a[i+1]);
	}



	for(i=0; i<length-1; i++){
		floorPrimitive((primVar+i*numVariables));
		primitiveToConserved((primVar+i*numVariables), aFluid[i], (consVar+i*numVariables));
	}
	

	//Because of the velocity profile we need to recalculate the metric to be consistent with the BC	
	solveHamiltonianConstraint(consVar, rGravity, length, a);
	
	//Interpolate metric to center of cells where fluid variables are defined
	for(i=0; i<length-1; i++){
		aFluid[i] = 0.5*(a[i]+a[i+1]);
	}
	
	for(i=0; i<length-1; i++){
		primitiveToConserved((primVar+i*numVariables), aFluid[i], (consVar+i*numVariables));
	}

	//Last two cells are ghost cells
	for(i=length-3; i<length-1; i++){ 
		for(j=0; j<numVariables; j++) consVar[i*numVariables +j] = consVar[(length-4)*numVariables+j];
	}

	//Impose outerboundary condition a(r_max)*alpha(r_max)=1
	alpha[length-1] = 1.0/a[length-1];

	solveSlicingCondition(consVar, primVar, a, rGravity, length, alpha);
	
	return starRadius;
}

/*
 For partial domains, this function extends the domain from r=0 to r=rMaxExt
 Then it uses getInitialData() to find ID over this extended domain and then
 copies in the required part.
 */
double getInitialDataPartialDomain(double centralPressure, double velocityAmp, double rMaxExt, 
				   double *consVar, double *primVar, double *a, double *alpha, 
   		                   double *rGravity, double *rFluid, int length){

	
	double dr = rGravity[1]-rGravity[0];
	int NrExt = (int) (rMaxExt/dr+1.5);
	double rGravityExt[NrExt];
	double aExt[NrExt];
	double alphaExt[NrExt];
	double rFluidExt[NrExt-1];
	double primVarExt[(NrExt-1)*numVariables];
	double consVarExt[(NrExt-1)*numVariables];
	
	int i,j;
	for(i=0; i<NrExt; i++) rGravityExt[i] = rMaxExt*((double) i)/(((double) NrExt)-1.0);
	for(i=0; i<NrExt-1; i++) rFluidExt[i] = 0.5*(rGravityExt[i]+rGravityExt[i+1]);

	double radius = getInitialData(centralPressure, velocityAmp, consVarExt, primVarExt, aExt, alphaExt, rGravityExt, rFluidExt, NrExt);

	int startIndex = (int) (rGravity[0]/dr+0.5);
	for(i=0; i<length; i++){
		a[i] = aExt[i+startIndex];
		alpha[i] = alphaExt[i+startIndex];
	}

	for(i=0; i<length-1; i++){
		for(j=0; j<numVariables; j++){
			consVar[i*numVariables+j]=consVarExt[(i+startIndex)*numVariables+j];
			primVar[i*numVariables+j]=primVarExt[(i+startIndex)*numVariables+j];
		}
	}

	printf("Partial ID r=%e to %e with extended r=%e and %e\n",rGravity[0],rGravity[length-1],rGravityExt[0],rGravityExt[NrExt-1]);	
	printf("aExt %e %e %e a %e %e %e\n", aExt[NrExt-1],aExt[NrExt-2],aExt[NrExt-3],a[length-1],a[length-2],a[length-3]);
	return radius;
}

