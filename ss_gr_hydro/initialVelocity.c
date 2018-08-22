#define CONSTANT_PI 3.1415926535897932384626433

#include "initialVelocity.h"
#include "gravity.h"
#include "variableDef.h"
#include "equationOfState.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*
Given an initial coordinate velocity and primitive variables, solve for metric functions a and alpha using pointwise 2-d Newton iteration (and hence get velocity).
a(0) and alpha(0) should already be set
*/
void solveMetricCoordVelocity(double *coordVelocity, double *primVar, double *r, int length, double *a, double *alpha){

	
	int numNewtIter = NUM_NEWT_ITER; //Maximum number of allowed Newton iterations
	double relError = REL_ERROR; //Relative error tolerance in solution
	int i,j;

	//Solve Hamiltonian constraint and slicing condtion using pointwise Newtonian iteration
	//Recast equation using b=ln(a) and beta=ln(alpha)
	double currB = 0.0;
	double currBeta =0.0;
	double prevB = log(a[0]);
	double prevBeta = log(alpha[0]);
	double prevDeltaB = 0.0;
	double prevDeltaBeta = 0.0;

	
	for(i=1; i<length; i++){
		
		double rAvg = 0.5*(r[i]+r[i-1]);
		double dr = r[i]-r[i-1];
		
		double pressure = getPressure(primVar[(i-1)*numVariables+REST_DENSITY],primVar[(i-1)*numVariables+SPECIFIC_ENERGY]);
		double h = 1.0+primVar[(i-1)*numVariables+ SPECIFIC_ENERGY]+pressure/primVar[(i-1)*numVariables+REST_DENSITY];
		double rhoh = primVar[(i-1)*numVariables+REST_DENSITY]*h;
		double fourpiR = 4.0*CONSTANT_PI*rAvg;
		double term = 0.5/rAvg;
		double residual_b;
		double residual_beta;
		double deriv_b_b;
		double deriv_b_beta;
		double deriv_beta_b;
		double deriv_beta_beta;

		//Set initial guess
		currB = prevB + prevDeltaB;
		currBeta = prevBeta + prevDeltaBeta;

		//Use Newton iteration to improve answer
		for(j=0; j<numNewtIter; j++){
			double velTerm = coordVelocity[i-1]*coordVelocity[i-1]*exp(currB+prevB-currBeta-prevBeta);
			residual_b = (currB-prevB)/dr - exp(currB+prevB)*(fourpiR*(rhoh/(1.0-velTerm)-pressure)-term) - term;
			deriv_b_b = 1.0/dr - exp(currB+prevB)*(fourpiR*(rhoh/(1.0-velTerm)-pressure)-term)- exp(currB+prevB)*fourpiR*rhoh*velTerm/(1.0-velTerm)/(1.0-velTerm);
			deriv_b_beta = exp(currB+prevB)*fourpiR*rhoh*velTerm/(1.0-velTerm)/(1.0-velTerm);
			residual_beta =  (currBeta-prevBeta)/dr - exp(currB+prevB)*(fourpiR*(rhoh*velTerm/(1.0-velTerm)+pressure)+term) + term;
			deriv_beta_beta = 1.0/dr + exp(currB+prevB)*fourpiR*rhoh*velTerm/(1.0-velTerm)/(1.0-velTerm);
			deriv_beta_b = -1.0*exp(currB+prevB)*(fourpiR*(rhoh*velTerm/(1.0-velTerm)+pressure)+term)-exp(currB+prevB)*fourpiR*rhoh*velTerm/(1.0-velTerm)/(1.0-velTerm);
			double det = deriv_b_b*deriv_beta_beta-deriv_b_beta*deriv_beta_b;
			double deltaB = (deriv_beta_beta*residual_b-deriv_b_beta*residual_beta)/det;
			double deltaBeta = (deriv_b_b*residual_beta-deriv_beta_b*residual_b)/det;

			currB-=deltaB;
			currBeta-=deltaBeta;

			if((fabs(deltaB/currB)+fabs(deltaBeta/currBeta))<relError){
				break;
			}

		}

		if(j==numNewtIter){
			printf("Error! Coordinate velocity constraint solver did not converge at r=%f. Residual is %E %E for b=%E beta=%E\n", r[i], residual_b, residual_beta, currB, currBeta);
			abort();
		} 
		a[i]=exp(currB);
		alpha[i] = exp(currBeta);
		prevDeltaB = currB - prevB;
		prevB = currB;
		prevDeltaBeta = currBeta - prevBeta;
		prevBeta = currBeta;
	}
	
	//Set velocity
	for(i=0; i<length-1; i++){
		primVar[i*numVariables+VELOCITY] = 0.5*coordVelocity[i]*(a[i]/alpha[i]+a[i+1]/alpha[i+1]);
	}
}
