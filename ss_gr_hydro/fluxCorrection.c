#define PLUS_X 2
#define MINUS_X 4 
#include <stdio.h>
#include "variableDef.h"
#include "fluxCorrection.h"

/*
Store fluxes specified by flux mask
*/
void storeFluxCorrection(double *flux, double *fluxAlt, double *rGravity, double *a, double *alpha, double dt, int numCells, double *fluxCorrection, int *fluxMask){

	/*
	Need to also do fluxAlt corrections...
	
	*/
	int i,j;
	int sign;
	for(i=1; i<(numCells-1); i++){
		if(fluxMask[i]){
			if(fluxMask[i]%2) sign = 1;
			else sign = -1;
			if(fluxMask[i] & PLUS_X){
				for(j=0; j<numVariables; j++) fluxCorrection[i*numVariables+j] += -1.0*sign*dt*rGravity[i+1]*rGravity[i+1]*alpha[i+1]/a[i+1]*flux[i*numVariables+j];
			} else if(fluxMask[i] & MINUS_X){
				for(j=0; j<numVariables; j++) fluxCorrection[i*numVariables+j] += sign*dt*rGravity[i]*rGravity[i]*alpha[i]/a[i]*flux[(i-1)*numVariables+j];
			} 
		}
		
	}

}
/*
Apply flux corrections
*/
void applyFluxCorrection(double *consVar, double *fluxCorrection, double *rGravity, int numCells){
	
	int i,j;
	for(i=1; i<(numCells-1); i++){
		double coeff = 3.0/(rGravity[i+1]*rGravity[i+1]*rGravity[i+1]-rGravity[i]*rGravity[i]*rGravity[i]);
		for(j=0; j<numVariables; j++){
			consVar[i*numVariables+j] += coeff*fluxCorrection[i*numVariables+j];
		}
	}
}

/*
Clear flux correction variables of the type determined by typeFlag.
typeFlag is the parity of the flux mask value of the flux corrections to clear. 
*/
void clearFluxCorrection(int typeFlag, double *fluxCorrection, int *fluxMask, int numCells){
	int i,j;
	for(i=0; i<numCells; i++){
		if( fluxMask[i] && ((fluxMask[i]%2)==typeFlag) ){
			for(j=0; j<numVariables; j++) fluxCorrection[i*numVariables+j] = 0.0;
		}
	}
}
