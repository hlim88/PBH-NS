/* This code was generated by rnpl, a numerical programming language */
/* copyright (c) 1994-1998 by Robert L. Marsa and Matthew W. Choptuik */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <signal.h>
#include <sdf.h>
#include <librnpl.h>

#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

/*  This routine updates the following grid functions 
    q1 q2 q3 p1 p2 p3 flux1 flux2 flux3 fluxAlt1 fluxAlt2 fluxAlt3 source1 source2 source3 
*/
void rungeKuttaUpdater(int *rnpldone,double *q1_np1,double *q1_n,double *q2_np1,double *q2_n,
                       double *q3_np1,double *q3_n,double *p1_np1,double *p1_n,double *p2_np1,
                       double *p2_n,double *p3_np1,double *p3_n,double *source1_n,
                       double *source2_n,double *source3_n,double *flux1_n,double *flux2_n,
                       double *flux3_n,double *fluxAlt1_n,double *fluxAlt2_n,double *fluxAlt3_n, 
                       double *alpha_n,double *a_n,double *test_n,int g1_Nr,double *r,double dt)
{

/*
Second order Runge-Kutta.
*/
#include "flux.h"
#include "gravity.h"
#include "boundaryConditions.h"
int length = g1_Nr;
int nVar = 3;
double *primVar = malloc(nVar*length*sizeof(double));
double *consVar = malloc(nVar*length*sizeof(double));
double *flux  = malloc(nVar*length*sizeof(double));
double *fluxAlt  = malloc(nVar*length*sizeof(double));
double *source = malloc(nVar*length*sizeof(double));
double *rFluid = malloc(nVar*length*sizeof(double));
double *aFluid = malloc(nVar*length*sizeof(double));
double *alphaFluid = malloc(nVar*length*sizeof(double));

int i;
/*
These are the coordinates of the centers of the fluid cells
*/
for(i=0; i<(length-1); i++){
	rFluid[i]=0.5*(r[i+1]+r[i]);
}
rFluid[length-1]=2.0*r[length-1]-rFluid[length-2];

/*
This array shuffling needs to be fixed...
*/
for(i=0; i<length; i++){
	primVar[i*nVar] = p1_n[i];
	primVar[i*nVar+1] = p2_n[i];
	primVar[i*nVar+2] = p3_n[i];
}

	i = length/2;
	double rPlusHalf = r[i+1];
	double rMinusHalf = r[i];
	double aPlusHalf = a_n[i+1];
	double aMinusHalf = a_n[i];
	double alphaPlusHalf = alpha_n[i+1];
	double alphaMinusHalf = alpha_n[i];
	double coeff = -3.0/(rPlusHalf*rPlusHalf*rPlusHalf-rMinusHalf*rMinusHalf*rMinusHalf);


//Take half step 
for(i=1; i<(length-2); i++){
	double rPlusHalf = r[i+1];
	double rMinusHalf = r[i];
	double aPlusHalf = a_n[i+1];
	double aMinusHalf = a_n[i];
	double alphaPlusHalf = alpha_n[i+1];
	double alphaMinusHalf = alpha_n[i];
	double coeff = -3.0/(rPlusHalf*rPlusHalf*rPlusHalf-rMinusHalf*rMinusHalf*rMinusHalf);
		
	consVar[i*nVar] = q1_n[i] + 0.5*dt*(coeff*(rPlusHalf*rPlusHalf*alphaPlusHalf/aPlusHalf*flux1_n[i]-rMinusHalf*rMinusHalf*alphaMinusHalf/aMinusHalf*flux1_n[i-1])-(alphaPlusHalf/aPlusHalf*fluxAlt1_n[i]-alphaMinusHalf/aMinusHalf*fluxAlt1_n[i-1])/(rPlusHalf-rMinusHalf)  + source1_n[i]);
	consVar[i*nVar+1] = q2_n[i] + 0.5*dt*(coeff*(rPlusHalf*rPlusHalf*alphaPlusHalf/aPlusHalf*flux2_n[i]-rMinusHalf*rMinusHalf*alphaMinusHalf/aMinusHalf*flux2_n[i-1])-(alphaPlusHalf/aPlusHalf*fluxAlt2_n[i]-alphaMinusHalf/aMinusHalf*fluxAlt2_n[i-1])/(rPlusHalf-rMinusHalf) + source2_n[i]);
	consVar[i*nVar+2] = q3_n[i]+ 0.5*dt*(coeff*(rPlusHalf*rPlusHalf*alphaPlusHalf/aPlusHalf*flux3_n[i]-rMinusHalf*rMinusHalf*alphaMinusHalf/aMinusHalf*flux3_n[i-1])-(alphaPlusHalf/aPlusHalf*fluxAlt3_n[i]-alphaMinusHalf/aMinusHalf*fluxAlt3_n[i-1])/(rPlusHalf-rMinusHalf) + source3_n[i]);
		
}
i=0;{ 
	//Use regularity condition at the origin
	setOriginBoundaryCondition(consVar);
	
}
for(i=length-2; i<length; i++){
	//Ghost cells
	consVar[i*nVar] = consVar[(i-1)*nVar];
	consVar[i*nVar+1] = consVar[(i-1)*nVar+1];
	consVar[i*nVar+2] = consVar[(i-1)*nVar+2];
}



//These values are not used
flux[(length-1)*nVar] = 0.0;
flux[(length-1)*nVar+1] = 0.0;
flux[(length-1)*nVar+2] = 0.0;


solveHamiltonianConstraint(consVar, r, length, a_n);
//Interpolate metric to center of cells where fluid variables are defined
	for(i=0; i<length-1; i++){
		aFluid[i] = 0.5*(a_n[i]+a_n[i+1]);
	}
	aFluid[length-1] = aFluid[length-2];
getNumericalFlux(consVar, primVar, aFluid, rFluid, length, flux, fluxAlt);

solveSlicingCondition(consVar, primVar, a_n, r, length, alpha_n);
//Interpolate metric to center of cells where fluid variables are defined
	for(i=0; i<length-1; i++){
		alphaFluid[i] = 0.5*(alpha_n[i]+alpha_n[i+1]);
	}
alphaFluid[length-1] = alphaFluid[length-2];
getSourceArrayAlt(primVar, consVar, aFluid, alphaFluid, rFluid, length, source);

//Take full step
for(i=1; i<(length-2); i++){

	double rPlusHalf = r[i+1];
	double rMinusHalf = r[i];
	double aPlusHalf = a_n[i+1];
	double aMinusHalf = a_n[i];
	double alphaPlusHalf = alpha_n[i+1];
	double alphaMinusHalf = alpha_n[i];
	double coeff = -3.0/(rPlusHalf*rPlusHalf*rPlusHalf-rMinusHalf*rMinusHalf*rMinusHalf);	
		
	consVar[i*nVar] = q1_n[i] + dt*(coeff*(rPlusHalf*rPlusHalf*alphaPlusHalf/aPlusHalf*flux[i*nVar]-rMinusHalf*rMinusHalf*alphaMinusHalf/aMinusHalf*flux[(i-1)*nVar])-(alphaPlusHalf/aPlusHalf*fluxAlt[i*nVar]-alphaMinusHalf/aMinusHalf*fluxAlt[(i-1)*nVar])/(rPlusHalf-rMinusHalf) + source[i*nVar]);
	consVar[i*nVar+1] = q2_n[i] + dt*(coeff*(rPlusHalf*rPlusHalf*alphaPlusHalf/aPlusHalf*flux[i*nVar+1]-rMinusHalf*rMinusHalf*alphaMinusHalf/aMinusHalf*flux[(i-1)*nVar+1])-(alphaPlusHalf/aPlusHalf*fluxAlt[i*nVar+1]-alphaMinusHalf/aMinusHalf*fluxAlt[(i-1)*nVar+1])/(rPlusHalf-rMinusHalf) + source[i*nVar+1]);
	consVar[i*nVar+2] = q3_n[i] + dt*(coeff*(rPlusHalf*rPlusHalf*alphaPlusHalf/aPlusHalf*flux[i*nVar+2]-rMinusHalf*rMinusHalf*alphaMinusHalf/aMinusHalf*flux[(i-1)*nVar+2])-(alphaPlusHalf/aPlusHalf*fluxAlt[i*nVar+2]-alphaMinusHalf/aMinusHalf*fluxAlt[(i-1)*nVar+2])/(rPlusHalf-rMinusHalf) + source[i*nVar+2]);
		
}
i=0;{ 
	//Use regularity condition at the origin
	setOriginBoundaryCondition(consVar);
	
}
for(i=length-2; i<length; i++){
	//Ghost cells
	consVar[i*nVar] = consVar[(i-1)*nVar];
	consVar[i*nVar+1] = consVar[(i-1)*nVar+1];
	consVar[i*nVar+2] = consVar[(i-1)*nVar+2];
}

solveHamiltonianConstraint(consVar, r, length, a_n);
//Interpolate metric to center of cells where fluid variables are defined
	for(i=0; i<length-1; i++){
		aFluid[i] = 0.5*(a_n[i]+a_n[i+1]);
	}
	aFluid[length-1] = aFluid[length-2];

getNumericalFlux(consVar, primVar, aFluid, rFluid, length, flux, fluxAlt);

solveSlicingCondition(consVar, primVar, a_n, r, length, alpha_n);

//Interpolate metric to center of cells where fluid variables are defined
	for(i=0; i<length-1; i++){
		alphaFluid[i] = 0.5*(alpha_n[i]+alpha_n[i+1]);
	}
alphaFluid[length-1] = alphaFluid[length-2];

getSourceArrayAlt(primVar, consVar, aFluid, alphaFluid, rFluid, length, source);

for(i=0; i<length; i++){

	q1_np1[i] = consVar[i*nVar];
	q2_np1[i] = consVar[i*nVar+1];
	q3_np1[i] = consVar[i*nVar+2];	

	p1_np1[i] = primVar[i*nVar];
	p2_np1[i] = primVar[i*nVar+1];
	p3_np1[i] = primVar[i*nVar+2];

	flux1_n[i] = flux[i*nVar];
	flux2_n[i] = flux[i*nVar+1];
	flux3_n[i] = flux[i*nVar+2];
	
	fluxAlt1_n[i] = fluxAlt[i*nVar];
	fluxAlt2_n[i] = fluxAlt[i*nVar+1];
	fluxAlt3_n[i] = fluxAlt[i*nVar+2];

	source1_n[i] = source[i*nVar];
	source2_n[i] = source[i*nVar+1];
	source3_n[i] = source[i*nVar+2];
	
	
}
test_n[0]=q1_np1[0];

for(i=1; i<length; i++){
	test_n[i]=0.5*(q1_np1[i-1]+q1_np1[i]);
}

free(primVar);
free(consVar);
free(flux);
free(fluxAlt);
free(source);
free(rFluid);
free(aFluid);
free(alphaFluid);}

