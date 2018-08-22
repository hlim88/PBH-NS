/*
Store fluxes specified by flux mask
*/
void storeFluxCorrection(double *flux, double *fluxAlt, double *rGravity, double *a, double *alpha, double dt, int numCells, double *fluxCorrection, int *fluxMask);

/*
Apply flux corrections
*/
void applyFluxCorrection(double *consVar, double *fluxCorrection, double *rGravity, int numCells);

/*
Clear flux correction variables of the type determined by typeFlag.
typeFlag is the parity of the flux mask value of the flux corrections to clear. 
*/
void clearFluxCorrection(int typeFlag, double *fluxCorrection, int *fluxMask, int numCells);
