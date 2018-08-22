#define VERBOSE 1

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "flux.h"
#include "conservedPrimitive.h"
#include "matrixInverse.h"

double minmod(double a, double b){
	double product = a*b;
	if(product<0.0){
		return 0.0;
	}else {
		if(fabs(a)<fabs(b)){
			return a;
		} else{
			return b;
		}
	}
	
}

/*
Slope limiter interpolation
*/
void slopeLimiterInterp(double fminus, double f, double fplus, double fplusplus, double xminus, double x, double xplus, double xplusplus, double *left, double *right){

	double sminus = (f-fminus)/(x-xminus);
	double s = (fplus-f)/(xplus-x);
	double splus = (fplusplus-fplus)/(xplusplus-xplus);

	*left = f + minmod(sminus, s)*0.5*(xplus-x);
	*right = fplus + minmod(s, splus)*0.5*(x-xplus);

}

/*
Conservative values are first floored and primitive values are calcluated using conservative variables and guess of primitive values.
*/
void getPrimitiveArray(double *consVar, double *primVar, double *a, int length){
	
	int i;
	for(i=0; i<length*numVariables; i++){ if(isnan(consVar[i])){printf("c to p:cons nan at i=%d for length=%d\n",i,length); } }
	
	for(i=0; i<length; i++){
		int varFloored;
		varFloored = floorConserved((consVar+i*numVariables));
		if(varFloored>0){
			#if VERBOSE >= 2
			printf("%d Conservative variables floored at i=%d\n",varFloored,i);
			printf("cons %E %E %E prim %E %E %E\n",consVar[i*numVariables],consVar[i*numVariables+1],consVar[i*numVariables+2],primVar[i*numVariables],primVar[i*numVariables+1],primVar[i*numVariables+2]);
			#endif 
		}
		conservedToPrimitive((consVar+i*numVariables), a[i], (primVar+i*numVariables));

		varFloored = floorPrimitive((primVar+i*numVariables));

		if(varFloored>0){
			primitiveToConserved((primVar+i*numVariables), a[i], (consVar+i*numVariables));
//			printf("cons %E %E %E prim %E %E %E\n",consVar[i*numVariables],consVar[i*numVariables+1],consVar[i*numVariables+2],primVar[i*numVariables],primVar[i*numVariables+1],primVar[i*numVariables+2]);
			#if VERBOSE >= 2
			printf("%d Primitive variables floored at i=%d\n",varFloored,i);
			#endif 
		}
	}

	for(i=0; i<length*numVariables; i++){ 
		if(isnan(primVar[i])){
			int n = i/numVariables;
			printf("c to p: nan at i=%d\n\n\n",n); 
			printf("Conserved %E %E %E\n",consVar[n*numVariables+0],consVar[n*numVariables+1],consVar[n*numVariables+2]);
			printf("Primitive %E %E %E\n",primVar[n*numVariables+0],primVar[n*numVariables+1],primVar[n*numVariables+2]);
			abort(); 
		} 
	}
}

/*
Get numerical flux.  If fluxAlt is a valid pointer two flux terms will be calculated.  If it a Null pointer, only one will.
*/

void getNumericalFlux(double *consVar, double *primVar, double *a, double *r, int length, double *flux, double *fluxAlt){
	
	int i,j,k;
	
	//Use slope limiter method to reconstruct primitive variables at boundaries
	double *primVarLeft = malloc(sizeof(double)*(length-1)*numVariables);
	double *primVarRight = malloc(sizeof(double)*(length-1)*numVariables);
	i=0;{
		for(j=0; j<numVariables; j++){
			double s = (primVar[(i+1)*numVariables+j]-primVar[(i)*numVariables+j])/(r[i+1]-r[i]);
			double splus = (primVar[(i+2)*numVariables+j]-primVar[(i+1)*numVariables+j])/(r[i+2]-r[i+1]);
			primVarLeft[i*numVariables+j] = primVar[i*numVariables+j] + 0.5*s*(r[i+1]-r[i]);
			primVarRight[i*numVariables+j] = primVar[(i+1)*numVariables+j] + 0.5*minmod(s,splus)*(r[i]-r[i+1]);
		}
	}
	for(i=1; i<length-2; i++){
		for(j=0; j<numVariables; j++){
			slopeLimiterInterp(*(primVar + (i-1)*numVariables +j), *(primVar + (i)*numVariables +j), *(primVar + (i+1)*numVariables +j), *(primVar + (i+2)*numVariables +j), r[i-1], r[i], r[i+1], r[i+2], (primVarLeft + i*numVariables +j), (primVarRight + i*numVariables +j));
		}
	}
	i=length-2;{
		for(j=0; j<numVariables; j++){
			double sminus = (primVar[(i)*numVariables+j]-primVar[(i-1)*numVariables+j])/(r[i]-r[i-1]);
			double s = (primVar[(i+1)*numVariables+j]-primVar[(i)*numVariables+j])/(r[i+1]-r[i]);
			primVarLeft[i*numVariables+j] = primVar[i*numVariables+j] + 0.5*minmod(s,sminus)*(r[i+1]-r[i]);
			primVarRight[i*numVariables+j] = primVar[(i+1)*numVariables+j] + 0.5*s*(r[i]-r[i+1]);
		}
	}

	for(i=0; i<(length-1); i++){
		/*
		Calculate conservative left and right variables
		*/
		double consVarLeft[numVariables];
		double consVarRight[numVariables];
		double deltaConsVar[numVariables];
		double primVarLRAvg[numVariables];	
	
		primitiveToConserved((primVarLeft+i*numVariables), a[i], consVarLeft);
		primitiveToConserved((primVarRight+i*numVariables), a[i], consVarRight);
		
		for(j=0; j<numVariables; j++){
			deltaConsVar[j]=consVarRight[j]-consVarLeft[j];
			primVarLRAvg[j]=0.5*(primVarLeft[i*numVariables+j]+primVarRight[i*numVariables+j]);
		}

		/*
		Calculate eigenvalues and right eigenvectors of flux Jacobian
		*/
		double eigenVectors[numVariables][numVariables];
		double eigenValues[numVariables];
		getFluxJacobianEigen(primVarLRAvg, a[i], eigenVectors, eigenValues);
		
		/*
		Invert eigenvector matrix
		*/
		double invEigenVectors[numVariables][numVariables];
		matrix3Inverse(eigenVectors, invEigenVectors);

		/*
		Calculate jumps in the characteristic values
		*/
		double jumps[numVariables];
		for(j=0; j<numVariables; j++){
			jumps[j]=0.0;
			for(k=0; k<numVariables; k++){
				jumps[j]+=invEigenVectors[j][k]*deltaConsVar[k];
			}
		}
		
		/*
		Calculate flux at cell boundaries
		Notice that if fluxAlt!=Null then two different flux components will be calculated
		*/	
		double leftFlux[numVariables];
		double rightFlux[numVariables];
		double leftFluxAlt[numVariables];
		double rightFluxAlt[numVariables];		

		if(fluxAlt==NULL){
			getPhysicalFlux((primVarLeft+i*numVariables), consVarLeft, leftFlux);
			getPhysicalFlux((primVarRight+i*numVariables), consVarRight, rightFlux);
		} else {
			getPhysicalFluxAltOne((primVarLeft+i*numVariables), consVarLeft, leftFlux);
			getPhysicalFluxAltOne((primVarRight+i*numVariables), consVarRight, rightFlux);
			getPhysicalFluxAltTwo((primVarLeft+i*numVariables), consVarLeft, leftFluxAlt);
			getPhysicalFluxAltTwo((primVarRight+i*numVariables), consVarRight, rightFluxAlt);
		
		}
		
		for(j=0; j<numVariables; j++){
			flux[i*numVariables+j] = 0.5*(leftFlux[j]+rightFlux[j]);

			for(k=0; k<numVariables; k++){
				flux[i*numVariables+j]-= 0.5*fabs(eigenValues[k])*jumps[k]*eigenVectors[j][k];
			}
			
			if(fluxAlt!=NULL){
				fluxAlt[i*numVariables+j] = 0.5*(leftFluxAlt[j]+rightFluxAlt[j]);
			}
			
		}
	}//loop over i

	free(primVarLeft);
	free(primVarRight);		

}

/*
Compute source function for array of values.
*/
void getSourceArray(double *primVar, double *consVar, double* a, double* alpha, double* r, int length, double *source){

	int i;
	for(i=0; i<length; i++){
		getSource((primVar+i*numVariables), (consVar+i*numVariables), a[i], alpha[i], r[i], (source+i*numVariables));
	}
}

/*
Compute alternate source function for array of values.
*/
void getSourceArrayAlt(double *primVar, double *consVar, double* a, double* alpha, double* r, int length, double *source){

	int i;
	for(i=0; i<length; i++){
		getSourceAlt((primVar+i*numVariables), (consVar+i*numVariables), a[i], alpha[i], r[i], (source+i*numVariables));
	}
}

/*
Compute the trace of the stress-energy tensor, T=3P-\rho for an array
*/
void getStressEnergyTraceArray(double *T, double *primVar, int length){
	int i;
	for(i=0; i<length; i++){
		getStressEnergyTrace((T+i),(primVar+i*numVariables));	
	}

}

