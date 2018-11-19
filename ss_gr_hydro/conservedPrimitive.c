#define VERBOSE 1
#define CONSTANT_PI 3.1415926535897932384626433
#define LAMBDA_HIGH (1.0e4)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "equationOfState.h"
#include "conservedPrimitive.h"
#include "brent.h"

/*
Determine conserved variables from primitive variables
*/
void primitiveToConserved(double *primVar, double a, double *consVar){
	
	double W = 1.0/sqrt(1.0-primVar[VELOCITY]*primVar[VELOCITY]);
	double pressure = getPressure(primVar[REST_DENSITY], primVar[SPECIFIC_ENERGY]);
	double h = 1.0+primVar[SPECIFIC_ENERGY]+pressure/primVar[REST_DENSITY];

	//DEBUG
	if(isinf(W)){
		double eps = 1.0-fabs(primVar[VELOCITY]);
		printf("W=inf and eps=%e\n",eps);
	}	
	consVar[CONS_DENSITY] = a*primVar[REST_DENSITY]*W;
	consVar[PI] = primVar[REST_DENSITY]*W*(h*W*(1.0+primVar[VELOCITY])-a)-pressure;
	consVar[PHI] = primVar[REST_DENSITY]*W*(h*W*(1.0-primVar[VELOCITY])-a)-pressure;

}


/*
 Old method to determine primitive variables from conserved variables.
 Primitive variables should already be set to some initial value which is used as a guess for the new value.
 */
/*
void conservedToPrimitiveOld(double *consVar, double a, double *primVar){
	
	 //Function that is zero when pressure is set to the correct value.
	double pressureFunction(double pressure){
		double velocity =0.5*(consVar[PI]-consVar[PHI])/(pressure+consVar[CONS_DENSITY]+0.5*(consVar[PI]+consVar[PHI]));
		double W = 1.0/sqrt(1.0-velocity*velocity);
		double specificEnergy = a*(0.5*(consVar[PI]+consVar[PHI])+consVar[CONS_DENSITY]+pressure*(1.0-W*W))/consVar[CONS_DENSITY]/W-1.0;
		double restDensity = consVar[CONS_DENSITY]/a/W;
		double func = getPressure(restDensity, specificEnergy) - pressure;
		return func;
	}
	
	//Use old pressure as intial guess	
	double pressureGuess = getPressure(primVar[REST_DENSITY], primVar[SPECIFIC_ENERGY]);
	
	if(pressureGuess<=0){
		printf("rho=%E e=%E\n",primVar[REST_DENSITY],primVar[SPECIFIC_ENERGY]);
	}
	
	//Use guess to get bracket around new value of pressure
	double pressure1, pressure2;
	int error = bracketRoot(pressureGuess, &pressure1, &pressure2, &pressureFunction);
	
	double pressure;
	if(error){
#if VERBOSE >= 2
		printf("Error converting from conservative to primitive variables.\n");
#endif
		pressure = 0.0;
	} else{ 
		//Use bracket to find root
		pressure = findRoot(pressure1, pressure2, &pressureFunction);
	}
	
	primVar[VELOCITY] =0.5*(consVar[PI]-consVar[PHI])/(pressure+consVar[CONS_DENSITY]+0.5*(consVar[PI]+consVar[PHI]));
	if(fabs(primVar[VELOCITY])>=1.0){
		if(primVar[VELOCITY]>0.0){
			primVar[VELOCITY] = (1.0-FLOOR_VALUE);
		} else {
			primVar[VELOCITY] = (-1.0+FLOOR_VALUE);
		}
#if VERBOSE >= 2
		printf("Error converting from conservative to primitive variables - superluminal.\n");
#endif
	}
	double W = 1.0/sqrt(1.0-primVar[VELOCITY]*primVar[VELOCITY]);
	primVar[SPECIFIC_ENERGY] = a*(0.5*(consVar[PI]+consVar[PHI])+consVar[CONS_DENSITY]+pressure*(1.0-W*W))/consVar[CONS_DENSITY]/W-1.0;
	primVar[REST_DENSITY] = consVar[CONS_DENSITY]/a/W;
	
	if(error){
		//	printf("Spec Energy=%E W=%E tau=%E term=%E\n",primVar[SPECIFIC_ENERGY],W,0.5*(consVar[PI]+consVar[PHI]),(0.5*(consVar[PI]+consVar[PHI])+consVar[CONS_DENSITY]+pressure*(1.0-W*W))/consVar[CONS_DENSITY]);
	}
	
}
*/


//Paramters used by enthalpyFunction()
struct enthalpyFuncParams{
	double *consVar;
	double a;
	double *primVar;
};

/*
Function that is zero when enthalpy is set to the correct value.
Used in conservedToPrimitive() for conversion of conserved variables to primitive variables.
*/
double enthalpyFunction(double H, void *vParam){
	struct enthalpyFuncParams *param = (struct enthalpyFuncParams *) vParam;
	double *consVar = param->consVar;
	double *primVar = param->primVar;
	double a = param->a;
	double LambdaHigh = LAMBDA_HIGH;
	double Lambda = 0.5*(consVar[PI]-consVar[PHI])/H;
	double func;
	
	if(fabs(Lambda)>LambdaHigh){
		double b = fabs(H/(consVar[PI]-consVar[PHI]));
		double btwo = b*b;
		double bfour = btwo*btwo;
		double beight = bfour*bfour;
		double X = 2.0*b-2.0*btwo+b*btwo-0.25*bfour*b+0.125*b*btwo*bfour; //W^-2 expanded to eighth order
		primVar[REST_DENSITY] = consVar[CONS_DENSITY]/a*sqrt(X);
		primVar[VELOCITY] = 1.0-b+0.5*btwo-0.125*bfour+bfour*btwo/16-5.0*beight/128;  //expanded to eighth order
		double onePlusVel,oneMinusVel, oneOverVelMinusOne, oneOverVelPlusOne;
		if(consVar[PI]<consVar[PHI]) {
			primVar[VELOCITY]*=-1.0;  //velocity should have same sign as S=0.5*(Pi-Phi)
			onePlusVel = b-0.5*btwo+0.125*bfour-bfour*btwo/16+5.0*beight/128;  //expanded to eighth order
			oneMinusVel = 1.0-primVar[VELOCITY];
			oneOverVelPlusOne = -1.0*b-0.5*btwo+0.125*bfour-bfour*btwo/16+5.0*beight/128;
			oneOverVelMinusOne = 1.0/primVar[VELOCITY]-1.0;
		} else {
			oneMinusVel = b-0.5*btwo+0.125*bfour-bfour*btwo/16+5.0*beight/128;  //expanded to eighth order
			onePlusVel = 1.0+primVar[VELOCITY];
			oneOverVelMinusOne = b+0.5*btwo-0.125*bfour+bfour*btwo/16-5.0*beight/128;
			oneOverVelPlusOne = 1.0/primVar[VELOCITY]+1.0;
		}
		primVar[SPECIFIC_ENERGY] = (consVar[CONS_DENSITY]+0.5*(consVar[PI]*(oneMinusVel)+consVar[PHI]*(onePlusVel)))/primVar[REST_DENSITY]-1.0;
		//func = 0.5*(consVar[PI]-consVar[PHI])/velocity-consVar[CONS_DENSITY]-0.5*(consVar[PI]+consVar[PHI])-getPressure(restDensity, specificEnergy);
		func = 0.5*(oneOverVelMinusOne)*consVar[PI]-0.5*(oneOverVelPlusOne)*consVar[PHI]-consVar[CONS_DENSITY]-getPressure(primVar[REST_DENSITY],primVar[SPECIFIC_ENERGY]);
	} else{
	
		if(fabs(Lambda)<(1.0/LambdaHigh)){
			double LambdaSq = Lambda*Lambda;
			primVar[VELOCITY]=(1.0+(-1.0+(2.0+(-5.0+(14.0-42.0*LambdaSq)*LambdaSq)*LambdaSq)*LambdaSq)*LambdaSq)*Lambda;
		} else{  
			primVar[VELOCITY]=0.5*(sqrt(1.0+4.0*Lambda*Lambda)-1.0)/Lambda;
		}
		double W = 1.0/sqrt(1.0-primVar[VELOCITY]*primVar[VELOCITY]);
		primVar[REST_DENSITY] = consVar[CONS_DENSITY]/a/W;
	
		primVar[SPECIFIC_ENERGY] = (consVar[CONS_DENSITY]+0.5*(consVar[PI]*(1.0-primVar[VELOCITY])+consVar[PHI]*(1.0+primVar[VELOCITY])))/primVar[REST_DENSITY]-1.0;
		
		if(fabs(Lambda)<(1.0/LambdaHigh)){
			func = (H*W*W-consVar[CONS_DENSITY]-0.5*(consVar[PI]+consVar[PHI]))-getPressure(primVar[REST_DENSITY], primVar[SPECIFIC_ENERGY]);
		} else{
			func = 0.5*(consVar[PI]-consVar[PHI])/primVar[VELOCITY]-consVar[CONS_DENSITY]-0.5*(consVar[PI]+consVar[PHI])-getPressure(primVar[REST_DENSITY], primVar[SPECIFIC_ENERGY]);
		}
	}
	return func;
}


/*
Determine primitive variables from conserved variables.
Primitive variables should already be set to some initial value which is used as a guess for the new value.
*/
void conservedToPrimitive(double *consVar, double a, double *primVar){

	
	double currPrimVar[numVariables];

	//Parameters used by enthalpyFunction()
	struct enthalpyFuncParams param = {consVar, a, currPrimVar};	


	//Use old enthalpy as intial guess	
	double HGuess = primVar[REST_DENSITY]*(1.0+primVar[SPECIFIC_ENERGY])+getPressure(primVar[REST_DENSITY], primVar[SPECIFIC_ENERGY]);


	//Use guess to get bracket around new value of H
	double H1, H2;
	int error = bracketRoot(HGuess, &H1, &H2, &enthalpyFunction, &param);

	double H;
	if(error){
		#if VERBOSE >= 2
		printf("Error converting from conservative to primitive variables.\n");
		#endif
		H = HGuess;
	} else{ 
		//Use bracket to find root
		H = findRoot(H1, H2, &enthalpyFunction, &param);
	}

	enthalpyFunction(H,&param);
	primVar[REST_DENSITY] = param.primVar[REST_DENSITY];
	primVar[SPECIFIC_ENERGY] = param.primVar[SPECIFIC_ENERGY];
	double velocity = param.primVar[VELOCITY];
	if(fabs(velocity)>=1.0){
		if(velocity>0.0){
			velocity = (1.0-EPS_VELOCITY);
		} else {
			velocity = (-1.0+EPS_VELOCITY);
		}
		#if VERBOSE >= 2
			printf("Error converting from conservative to primitive variables - superluminal.\n");
		#endif
	}
	primVar[VELOCITY] = velocity;
	
	
	if(isnan(primVar[REST_DENSITY]) || isnan(primVar[SPECIFIC_ENERGY]) || isnan(primVar[VELOCITY])){
		printf("rho=%E e=%E v=%E\n",primVar[REST_DENSITY],primVar[SPECIFIC_ENERGY],primVar[VELOCITY]);	
		abort();
	}
}

/*
Compute physical flux from primitive and conservative variables
*/
void getPhysicalFlux(double *primVar, double *consVar, double *flux){
	double pressure = getPressure(primVar[REST_DENSITY],primVar[SPECIFIC_ENERGY]);
	flux[CONS_DENSITY] = primVar[VELOCITY]*consVar[CONS_DENSITY];
	flux[PI] = primVar[VELOCITY]*(consVar[PI]+pressure)+pressure;
	flux[PHI] = primVar[VELOCITY]*(consVar[PHI]+pressure)-pressure;
}

/*
Compute source function
*/
void getSource(double *primVar, double *consVar, double a, double alpha, double r, double *source){
	double pressure = getPressure(primVar[REST_DENSITY], primVar[SPECIFIC_ENERGY]);
	double h = 1.0+primVar[SPECIFIC_ENERGY]+pressure/primVar[REST_DENSITY];
	double mass = 0.5*r*(1.0-1.0/a/a);

	source[CONS_DENSITY]=0.0;
	source[PI]=alpha*a*((pressure - h *primVar[REST_DENSITY])*(8.0*CONSTANT_PI*r*pressure+mass/r/r)+pressure*mass/r/r) + 2.0*pressure*alpha/a/r;
	source[PHI]=-1.0*source[PI];
}

/*
This is an alternate formulation of the above source and flux terms into two flux terms
*/

/*
Compute first physical flux for alternative form from primitive and conservative variables
*/
void getPhysicalFluxAltOne(double *primVar, double *consVar, double *flux){
	double pressure = getPressure(primVar[REST_DENSITY],primVar[SPECIFIC_ENERGY]);
	flux[CONS_DENSITY] = primVar[VELOCITY]*consVar[CONS_DENSITY];
	flux[PI] = primVar[VELOCITY]*(consVar[PI]+pressure);
	flux[PHI] = primVar[VELOCITY]*(consVar[PHI]+pressure);
}

/*
Compute second physical flux for alternative form from primitive and conservative variables
*/
void getPhysicalFluxAltTwo(double *primVar, double *consVar, double *flux){
	double pressure = getPressure(primVar[REST_DENSITY],primVar[SPECIFIC_ENERGY]);
	flux[CONS_DENSITY] = 0.0;
	flux[PI] = pressure;
	flux[PHI] = -1.0*pressure;
}

/*
Compute alternate source function
*/
void getSourceAlt(double *primVar, double *consVar, double a, double alpha, double r, double *source){
	double pressure = getPressure(primVar[REST_DENSITY], primVar[SPECIFIC_ENERGY]);
	double h = 1.0+primVar[SPECIFIC_ENERGY]+pressure/primVar[REST_DENSITY];
	double mass = 0.5*r*(1.0-1.0/a/a);

	source[CONS_DENSITY]=0.0;
	source[PI]=alpha*a*((pressure - h *primVar[REST_DENSITY])*(8.0*CONSTANT_PI*r*pressure+mass/r/r)+pressure*mass/r/r);
	//source[PI]=alpha*a*(-8.0*CONSTANT_PI*r*pressure*primVar[REST_DENSITY]*(1.0+primVar[SPECIFIC_ENERGY])+(pressure-primVar[REST_DENSITY]*primVar[SPECIFIC_ENERGY])*mass/r/r)-primVar[REST_DENSITY]*mass/r/r);
	source[PHI]=-1.0*source[PI];
}

/*
Compute right eigenvectors and eigenvalues of the Jacbian of the flux (w.r.t to conservative variables).
Return eigenvectors as matrix where the columns are the eigenvectors and eigenvalues as array.
*/
void getFluxJacobianEigen(double *primVar, double a, double A[numVariables][numVariables], double lambda[numVariables]){
	
	double W = 1.0/sqrt(1.0-primVar[VELOCITY]*primVar[VELOCITY]);
	double pressure = getPressure(primVar[REST_DENSITY], primVar[SPECIFIC_ENERGY]);
	double h = 1.0+primVar[SPECIFIC_ENERGY]+pressure/primVar[REST_DENSITY];	
	/*
	double temp = (getDerivPressureDensity(primVar[REST_DENSITY], primVar[SPECIFIC_ENERGY])+pressure/primVar[REST_DENSITY]/primVar[REST_DENSITY]*getDerivPressureEnergy(primVar[REST_DENSITY], primVar[SPECIFIC_ENERGY]))/h;
	double soundSpeed = sqrt((getDerivPressureDensity(primVar[REST_DENSITY], primVar[SPECIFIC_ENERGY])+pressure/primVar[REST_DENSITY]/primVar[REST_DENSITY]*getDerivPressureEnergy(primVar[REST_DENSITY], primVar[SPECIFIC_ENERGY]))/h);

	if(soundSpeed!=soundSpeed){
		printf("Sound speed =%E temp=%E density=%e specificEnergy=%E vel=%E \n",soundSpeed,temp, primVar[REST_DENSITY], primVar[SPECIFIC_ENERGY], primVar[VELOCITY]);
		abort();
	}
	*/
	double soundSpeed = getSoundSpeed(primVar[REST_DENSITY], primVar[SPECIFIC_ENERGY]);
	
	A[CONS_DENSITY][CONS_DENSITY] = 1.0;
	A[PI][CONS_DENSITY] = W*(1.0+primVar[VELOCITY])/a-1.0;
	A[PHI][CONS_DENSITY] = W*(1.0-primVar[VELOCITY])/a-1.0;

	A[CONS_DENSITY][PI] = 1.0;
	A[PI][PI] = W*(1.0+primVar[VELOCITY])/a*h*(1.0+soundSpeed)-1.0;
	A[PHI][PI] = W*(1.0-primVar[VELOCITY])/a*h*(1.0-soundSpeed)-1.0;
	
	A[CONS_DENSITY][PHI] = 1.0;
	A[PI][PHI] = W*(1.0+primVar[VELOCITY])/a*h*(1.0-soundSpeed)-1.0;
	A[PHI][PHI] = W*(1.0-primVar[VELOCITY])/a*h*(1.0+soundSpeed)-1.0;

	lambda[CONS_DENSITY] = primVar[VELOCITY];
	lambda[PI] = (primVar[VELOCITY] + soundSpeed)/(1.0+primVar[VELOCITY]*soundSpeed);
	lambda[PHI] = (primVar[VELOCITY] - soundSpeed)/(1.0-primVar[VELOCITY]*soundSpeed);
	
}

/*
Impose a floor to keep the conserved variables from taking on unphysical values.
*/
int floorConserved(double* consVar){
	double floorValue = FLOOR_VALUE;
	int numVarFloored = 0;	
	if(consVar[CONS_DENSITY]<floorValue){
		consVar[CONS_DENSITY] = floorValue;
		numVarFloored++;
	} 
	if((consVar[CONS_DENSITY]+consVar[PI])<2.0*floorValue){
		consVar[PI] = 2.0*floorValue - consVar[CONS_DENSITY];
		numVarFloored++;
	}
	if((consVar[CONS_DENSITY]+consVar[PHI])<2.0*floorValue){
		consVar[PHI] = 2.0*floorValue - consVar[CONS_DENSITY];
		numVarFloored++;
	}
	/*
	if((consVar[PHI]+consVar[PI])<2.0*floorValue){
		double S = 0.5*(consVar[PI]-consVar[PHI]);
		consVar[PHI]=floorValue + S;
		consVar[PI]=floorValue - S;
		numVarFloored++;
	}
	*/	
	return numVarFloored;
}

int floorPrimitive(double* primVar){
	double floorValue = FLOOR_VALUE;
	int numVarFloored = 0;	
	if(primVar[REST_DENSITY]<floorValue){
		primVar[REST_DENSITY] = floorValue;
		numVarFloored++;
		if(primVar[SPECIFIC_ENERGY] > 100.0) primVar[SPECIFIC_ENERGY]=100.0;
	} 
	if(primVar[SPECIFIC_ENERGY]<floorValue){
		primVar[SPECIFIC_ENERGY] = floorValue;
		numVarFloored++;
	}
	if(fabs(primVar[VELOCITY])>=1.0){
		if(primVar[VELOCITY]>0.0){
			primVar[VELOCITY] = (1.0-floorValue);
		} else {
			primVar[VELOCITY] = (-1.0+floorValue);
		}
		numVarFloored++;
	}
	
	return numVarFloored;

}

/*
Compute the trace of the stress-energy tensor, T=3P-\rho
*/
void getStressEnergyTrace(double *T, double *primVar){
	*T=3.0*getPressure(primVar[REST_DENSITY],primVar[SPECIFIC_ENERGY])-(1.0+primVar[SPECIFIC_ENERGY])*primVar[REST_DENSITY];
}

