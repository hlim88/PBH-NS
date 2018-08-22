/*
This array shuffling needs to be fixed...
*/
#include "initialData.h"
int length = g1_Nr;
int nVar = 3;
double *primVar = malloc(nVar*length*sizeof(double));
double *consVar = malloc(nVar*length*sizeof(double));
double *flux  = malloc(nVar*length*sizeof(double));
double *fluxAlt  = malloc(nVar*length*sizeof(double));
double *source = malloc(nVar*length*sizeof(double));
double *rFluid = malloc(nVar*length*sizeof(double));
	
int i;

/*
These are the coordinates of the centers of the fluid cells
*/
for(i=0; i<(length-1); i++){
	rFluid[i]=0.5*(r[i+1]+r[i]);
}
rFluid[length-1]=2.0*r[length-1]-rFluid[length-2];

//These values are not used
flux[(length-1)*nVar] = 0.0;
flux[(length-1)*nVar+1] = 0.0;
flux[(length-1)*nVar+2] = 0.0;


getInitialData(consVar, primVar, a_n, alpha_n, flux, fluxAlt, source, r, rFluid, length);

for(i=0; i<length; i++){
	
	q1_n[i] = consVar[i*nVar];
	q2_n[i] = consVar[i*nVar+1];
	q3_n[i] = consVar[i*nVar+2];
	
	p1_n[i] = primVar[i*nVar];
	p2_n[i] = primVar[i*nVar+1];
	p3_n[i] = primVar[i*nVar+2];

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

free(primVar);
free(consVar);
free(flux);
free(fluxAlt);
free(source);
free(rFluid);
