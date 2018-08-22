/*
This array shuffling needs to be fixed...
*/
#include "flux.h"
int length = g1_Nr;
int nVar = 3;
double *primVar = malloc(nVar*length*sizeof(double));
double *consVar = malloc(nVar*length*sizeof(double));
double *flux  = malloc(nVar*length*sizeof(double));
double *source = malloc(nVar*length*sizeof(double));
int i;
for(i=0; i<length; i++){
	consVar[i*nVar] = q1_np1[i];
	consVar[i*nVar+1] = q2_np1[i];
	consVar[i*nVar+2] = q3_np1[i];
		
	primVar[i*nVar] = p1_np1[i];
	primVar[i*nVar+1] = p2_np1[i];
	primVar[i*nVar+2] = p3_np1[i];
}
	//These values are not used
	flux[(length-1)*nVar] = 0.0;
	flux[(length-1)*nVar+1] = 0.0;
	flux[(length-1)*nVar+2] = 0.0;

getNumericalFlux(consVar, primVar, a_np1, r, length, flux);
getSourceArray(primVar, consVar, a_np1, alpha_np1, r, source);


for(i=0; i<length; i++){
	p1_np1[i] = primVar[i*nVar];
	p2_np1[i] = primVar[i*nVar+1];
	p3_np1[i] = primVar[i*nVar+2];

	flux1_np1[i] = flux[i*nVar];
	flux2_np1[i] = flux[i*nVar+1];
	flux3_np1[i] = flux[i*nVar+2];

	source1_np1[i] = source[i*nVar];
	source2_np1[i] = source[i*nVar+1];
	source3_np1[i] = source[i*nVar+2];
	
}

free(primVar);
free(consVar);
free(flux);
free(source);
