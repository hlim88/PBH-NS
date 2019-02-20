#define CONSTANT_PI 3.1415926535897932384626433

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "variableDef.h"
#include "equationOfState.h"
#include "gravity.h"

/*
Solve Hamiltonian constraint for a using pointwise Newtonian iteration.
It is assumed that a[0] is already set to the correct value.
*/
void solveHamiltonianConstraint(double *consVar, double *r, int length, double *a){
	int numNewtIter = NUM_NEWT_ITER; //Maximum number of allowed Newton iterations
	double relError = REL_ERROR; //Relative error tolerance in solution
	int i,j;
	

	
	//Solve Hamiltonian constraint using pointwise Newtonian iteration
	//Recast equation using A=ln(a)
	double currA = 0.0;
	double prevA = 0.0;
	double prevDeltaA = 0.0;
	for(i=1; i<length; i++){
		
		double E = 0.5*(consVar[(i-1)*numVariables+PI]+consVar[(i-1)*numVariables+PHI])+consVar[(i-1)*numVariables+REST_DENSITY];
		double rAvg = 0.5*(r[i]+r[i-1]);
		double dr = r[i]-r[i-1];
		double coef = 4.0*CONSTANT_PI*rAvg*E;
		double constant = 0.5/rAvg;
		double residual;
		double derivative;
		//Set initial guess
		currA = prevA + prevDeltaA;
		
		//Use Newton iteration to improve answer
		for(j=0; j<numNewtIter; j++){
			
			residual = (currA-prevA)/dr - coef*exp(currA+prevA) - constant*(1.0-exp(currA+prevA));
			derivative = 1.0/dr - (coef-constant)*exp(currA+prevA);
			double deltaA = residual/derivative;
				
			currA-=deltaA;
			if(fabs(deltaA/currA)<relError){
				break;
			}

		}
		if(j==numNewtIter){
			printf("Error! Hamiltonian constraint solver did not converge at r=%f. Residual is %E for A=%E, E=%E\n",r[i],residual,currA,E);
			printf("currA=%E prevA=%E dr=%E coef=%E constant=%E\n",currA,prevA,dr,coef,constant);
			abort();
		} 
		a[i]=exp(currA);
	
		prevDeltaA = currA - prevA;
		prevA = currA;
		
		//Black hole finder
		double Z = 1.0-1.0/a[i]/a[i];
		if(Z>0.99){
			printf("Black hole formation at r=%E\n",r[i]);
		}


	} //Calculated a[i]
}

/*
Calculate the residual of the Hamiltonian constraint
*/
void residualHamiltonianConstraint(double *consVar, double *r, int length, double *a, double *residual){
	
	double dr = r[1]-r[0];
	int i;
	for(i=1; i<length; i++){
		
		double E = 0.5*(consVar[(i-1)*numVariables+PI]+consVar[(i-1)*numVariables+PHI])+consVar[(i-1)*numVariables+REST_DENSITY];
		double rAvg = 0.5*(r[i]+r[i-1]);
		double aAvg = 0.5*(a[i]+a[i-1]);
		residual[i-1] = (a[i]-a[i-1])/dr-aAvg*aAvg*aAvg*(4.0*CONSTANT_PI*E*rAvg-0.5/rAvg)-0.5*aAvg/rAvg;
	}
}

/*
Solve polar-areal slicing condition for alpha.
It is assumed that the last entry of alpha is already set.
*/
void solveSlicingCondition(double *consVar, double *primVar, double *a, double *r, int length, double *alpha){

	int k;

	//Recast equations in terms of B=ln(alpha)
	double lnalpha = log(alpha[length-1]);
	
	for(k=length-1; k>0; k--){

		double pressure = getPressure(primVar[(k-1)*numVariables+REST_DENSITY],primVar[(k-1)*numVariables+SPECIFIC_ENERGY]);
		double coef = 0.5*(consVar[(k-1)*numVariables+PI]-consVar[(k-1)*numVariables+PHI])*primVar[(k-1)*numVariables+VELOCITY]+pressure;
		double rAvg = 0.5*(r[k]+r[k-1]);
		double aAvg = 0.5*(a[k]+a[k-1]);
		double dr = r[k]-r[k-1];
		double aSq = aAvg*aAvg;
		double psQuantity = aSq*4.0*CONSTANT_PI*rAvg*coef+0.5/rAvg*(aSq-1.0);
		lnalpha-=dr*psQuantity;
		alpha[k-1]=exp(lnalpha);
	}	

}

/*
Calculate the residual of the momentum equation (a evolution equation) using CN finite differencing.
*/
void calculateMomentumEqResidual(double *consVar_n, double *consVar_np1, double *a_n, double *a_np1, double *alpha_n, double *alpha_np1, double *r, double dt, int length, double *residual){
	int i;
	residual[0]=0.0;
	for(i=1; i<length; i++){
		double S_np1 = 0.25*(consVar_np1[(i-1)*numVariables+PI]-consVar_np1[(i-1)*numVariables+PHI])+0.25*(consVar_np1[i*numVariables+PI]-consVar_np1[i*numVariables+PHI]);
		double S_n = 0.25*(consVar_n[(i-1)*numVariables+PI]-consVar_n[(i-1)*numVariables+PHI])+0.25*(consVar_n[i*numVariables+PI]-consVar_n[i*numVariables+PHI]);
		residual[i]=(a_np1[i]-a_n[i])/dt+4.0*CONSTANT_PI*r[i]*0.5*(alpha_np1[i]*a_np1[i]*a_np1[i]*S_np1+alpha_n[i]*a_n[i]*a_n[i]*S_n);

	}

}

/*
Use momentum equation to evolve a
We let a_np1 = a_n + dt*F(a)
*/
void evolveMomentumEq(double *consVar, double *a,  double *a_n, double *a_np1, double *alpha, double *r, double dt, int length, int innerBdy, int outerBdy){
	int i;
	if(innerBdy) a_np1[0]=1.0;
	for(i=1; i<length-1; i++){
		double S = 0.25*(consVar[(i-1)*numVariables+PI]-consVar[(i-1)*numVariables+PHI])+0.25*(consVar[i*numVariables+PI]-consVar[i*numVariables+PHI]);
		a_np1[i] = a_n[i] - dt*4.0*CONSTANT_PI*r[i]*alpha[i]*a[i]*a[i]*S;
	}
	i=length-1;
	if(outerBdy) a_np1[i] = a_n[i] - dt*4.0*CONSTANT_PI*r[i]*alpha[i]*a[i]*a[i]*0.5*(consVar[(i-1)*numVariables+PI]-consVar[(i-1)*numVariables+PHI]);
}

#if 1
/*
  Slicing condition for alpha. 
  We choose 1+log slicing
*/
void solveAlphaSlicingCondition(double *consVar, double *primVar, double *r,
                             int length, double *alpha, double *trK) {

    int k;

	
	for(k=length-1; k>0; k--){

		double pressure = getPressure(primVar[(k-1)*numVariables+REST_DENSITY],
                                      primVar[(k-1)*numVariables+SPECIFIC_ENERGY]);
		double coef = 0.5*(consVar[(k-1)*numVariables+PI]-consVar[(k-1)*numVariables+PHI]) 
                         *primVar[(k-1)*numVariables+VELOCITY]+pressure;
		double rAvg = 0.5*(r[k]+r[k-1]);
		double trKAvg = 0.5*(trK[k]+trK[k-1]);
		double dr = r[k]-r[k-1];
		double trKSq = trKAvg*trKAvg;
		double psQuantity = trKSq*4.0*CONSTANT_PI*rAvg*coef+0.5/rAvg*(trKSq-1.0);
		alpha[k-1]=-dr*psQuantity;
	}	

}
/* 
 Slicing condition for beta
 We choose Gamma-freezing with auxilary variable
*/
void solveShiftSlicingCondition(double *consVar, double *primVar, double *r,
                              int length, double *betaR, double *Br, double *GamDelta){


    int i;
    for(i=length-1; i>0; i--) {
        Br[i] =3.0/4.0*GamDelta[i];
        betaR[i] = Br[i];
    }
}

#endif
