#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "variableDef.h"
#include "flux.h"
#include "gravity.h"
#include "boundaryConditions.h"
#include "fluxCorrection.h"
#include "fluidEvolve.h"


/*
2nd order Runge-Kutta evolution.  
If iter=1 we take the half step:
q_np1=q_n+0.5*dt*F(q_n)
If iter=2 we take the whole step:
q_np1=q_n+dt*F(q_np1)
We assume that a, alpha, rGravity have length entries.
We assume that consVar and primVar have numVariables*(length-1) entries.
if(boolInnerBdy) then the inner boundary is a physical boundary.
if(boolOuterBdy) then the outer boundary is a physical boundary.
*/
void timeStep(int iter, int boolInnerBdy, int boolOuterBdy, double dt, int length, double *rGravity, double* consVar_n, double* consVar_np1, double* primVar_n, double* primVar_np1, double *a_n, double *a_np1, double *phi_n, double *phi_np1, double *fluxCorrection, int *fluxCorrectionMask){

	

	int i,j;

	//DEBUG -- Check mass conservation
	static double initialMass = 0.0;
	if(boolInnerBdy && boolOuterBdy && length==4097){
		double totMass = 0.0;
		for(i=0; i<length-1; i++){
            totMass+= consVar_n[i*numVariables]*(rGravity[i+1]*rGravity[i+1]*rGravity[i+1]
                    - rGravity[i]*rGravity[i]*rGravity[i]);
        }
		if(initialMass==0.0){
            initialMass = totMass;
        }
		printf("totalMass on L=%d is %E loss fraction=%E\n",length,
                   totMass,(initialMass-totMass)/initialMass);
	}

	int lengthF = length-1;		
	double *consVar, *primVar, *phi;
	double a[length];
	double alpha[length];
	double flux[numVariables*(lengthF-1)];
	double rFluid[lengthF];
	double aFluid[lengthF];
	double alphaFluid[lengthF];
	double fluxAlt[numVariables*(lengthF-1)];
	double source[numVariables*(lengthF)];
	double deltaT;
	

	if(iter==1){
		consVar = consVar_n;
		primVar = primVar_n;
		for(i=0; i<length; i++) a[i] = a_n[i];
		phi = phi_n;
		deltaT=0.5*dt;

	} else if (iter==2){
		consVar = consVar_np1;
		primVar = primVar_np1;
		for(i=0; i<length; i++) a[i] = a_np1[i];
		phi = phi_np1;
		deltaT=dt;
	} else {
		//Make sure primitives arrays are in synch with conserved quantities
		for(i=0; i<lengthF; i++) aFluid[i]=0.5*(a_np1[i]+a_np1[i+1]);
		getPrimitiveArray(consVar_np1, primVar_np1, aFluid, lengthF);
		return;
	}
	

	/*
	Set alpha = exp(dphi/dr) 
	*/
	alpha[0] = exp((-1.5*phi[0]+2.0*phi[1]-0.5*phi[2])/(rGravity[1]-rGravity[0]));
	for(i=1; i<length-1; i++){
		alpha[i] = exp((phi[i+1]-phi[i-1])/(rGravity[i+1]-rGravity[i-1]));
	}
	alpha[length-1] = exp((1.5*phi[length-1]-2.0*phi[length-2]+0.5*phi[length-3])/(rGravity[length-1]-rGravity[length-2]));

	/*
	Evolve a using momentum equation
	*/
	evolveMomentumEq(consVar, a, a_n, a_np1, alpha, rGravity, deltaT, length, boolInnerBdy, boolOuterBdy);
	for(i=0; i<length; i++) if(a_np1[i]<1.0) a_np1[i]=1.0;
	
	/*
	These are the coordinates of the centers of the fluid cells
	*/
	for(i=0; i<(lengthF); i++){
		rFluid[i]=0.5*(rGravity[i+1]+rGravity[i]);
		aFluid[i]=0.5*(a[i+1]+a[i]);
		alphaFluid[i]=0.5*(alpha[i+1]+alpha[i]);
	}
	
	//Might want to recalculate primitives...
	getPrimitiveArray(consVar, primVar, aFluid, lengthF);
	
	//Get flux functions
	getNumericalFlux(consVar, primVar, aFluid, rFluid, lengthF, flux, fluxAlt);

	//Store flux corrections on the second iteration
	if(iter==2) storeFluxCorrection(flux, fluxAlt, rGravity, a, alpha, deltaT, lengthF, fluxCorrection, fluxCorrectionMask);
		
	//Get source functions
	getSourceArrayAlt(primVar, consVar, aFluid, alphaFluid, rFluid, lengthF, source);

	for(i=1; i<(lengthF-1); i++){

		double coeff = -3.0/(rGravity[i+1]*rGravity[i+1]*rGravity[i+1]-rGravity[i]*rGravity[i]*rGravity[i]);
		for(j=0; j<numVariables; j++){
			consVar_np1[i*numVariables+j] = consVar_n[i*numVariables+j] 
                                          + deltaT*(coeff*(rGravity[i+1]*rGravity[i+1]
                                          * alpha[i+1]/a[i+1]*flux[i*numVariables+j]
                                          - rGravity[i]*rGravity[i]*alpha[i]/a[i]*flux[(i-1)
                                          * numVariables+j])-(alpha[i+1]/a[i+1]
                                          * fluxAlt[i*numVariables+j]-alpha[i]/a[i]
                                          * fluxAlt[(i-1)*numVariables+j])/(rGravity[i+1]
                                          - rGravity[i]) + source[i*numVariables+j]);
			if(isnan(consVar_np1[i*numVariables+j])) {
				printf("timeStep: consVar_np1 nan at r=%E and" 
                        "consVar_n= %E a= %E alpha = %E flux =%E fluxAlt =%E\n",rFluid[i], 
                        consVar_n[i*numVariables+j], a[i], alpha[i], flux[i*numVariables+j], 
                        fluxAlt[i*numVariables+j]);
			}		
		}
	}
	
	//Apply boundary conditions if boundaries are physical
	if(boolInnerBdy){
		if(length<6){
			printf("Grid is not big enough for origin interpolation\n");
		} else{
			setOriginBoundaryCondition(consVar_np1);
		}
	}

	if(boolOuterBdy){

		if(length<4){
			printf("Grid is not big enough for outer boundary ghost cells\n");
		} else{
			setOuterBoundaryCondition(consVar_np1, lengthF);
		}
	}
	
	for(i=0; i<(lengthF); i++) aFluid[i]=0.5*(a_np1[i+1]+a_np1[i]);

	getPrimitiveArray(consVar_np1, primVar_np1, aFluid, lengthF);

	//residualHamiltonianConstraint(consVar_n, rGravity, length, a_n, alpha_np1);
	
}
	
