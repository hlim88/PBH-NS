/*
Conservative values are first floored and primitive values are calcluated using conservative variables and guess of primitive values.
*/
void getPrimitiveArray(double *consVar, double *primVar, double *a, int length);

/*
Get numerical flux.  If fluxAlt is a valid pointer two flux terms will be calculated.  If it a Null pointer, only one will.
*/

void getNumericalFlux(double *consVar, double *primVar, double *a, double *r, int length, double *flux, double *fluxAlt);
/*
Compute source function for array of values.
*/
void getSourceArray(double *primVar, double *consVar, double* a, double* alpha, double* r, int length, double *source);

/*
Compute alternate source function for array of values.
*/
void getSourceArrayAlt(double *primVar, double *consVar, double* a, double* alpha, double* r, int length, double *source);

/*
Compute the trace of the stress-energy tensor, T=3P-\rho of an array
*/
void getStressEnergyTraceArray(double *T, double *primVar, int length);
