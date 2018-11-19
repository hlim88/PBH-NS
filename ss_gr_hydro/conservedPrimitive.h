#include "variableDef.h"
/*
Value below which nonnegative variables are not allowed to fall, used in floorConserved and floorPrimitive
*/
#define FLOOR_VALUE (1.0e-18)
/*
Velocity is only allowed to get this close to speed of light 
*/
#define EPS_VELOCITY (1.0e-16)

/*
Determine conserved variables from primitive variables
*/
void primitiveToConserved(double *primVar, double a, double *consVar);

/*
Determine primitive variables from conserved variables
*/
void conservedToPrimitive(double *consVar, double a, double *primVar);
//Old version of the above
void conservedToPrimitiveOld(double *consVar, double a, double *primVar);

/*
Compute physical flux from primitive and conservative variables
*/
void getPhysicalFlux(double *primVar, double *consVar, double *flux);

/*
Compute source function
*/
void getSource(double *primVar, double *consVar, double a, double alpha, double r, double *source);

/*
This is an alternate formulation of the above source and flux terms into two flux terms
*/

/*
Compute first physical flux for alternative form from primitive and conservative variables
*/
void getPhysicalFluxAltOne(double *primVar, double *consVar, double *flux);

/*
Compute second physical flux for alternative form from primitive and conservative variables
*/
void getPhysicalFluxAltTwo(double *primVar, double *consVar, double *flux);

/*
Compute alternate source function
*/
void getSourceAlt(double *primVar, double *consVar, double a, double alpha, double r, double *source);

/*
Compute right eigenvectors and eigenvalues of the Jacbian of the flux (w.r.t to conservative variables).
Return eigenvectors as matrix where the columns are the eigenvectors and eigenvalues as array.
*/
void getFluxJacobianEigen(double *primVar, double a, double A[3][3], double lambda[3]);

/*
Impose a floor to keep the conserved variables from taking on unphysical values.
*/
int floorConserved(double* consVar);

/*
Impose a floor to keep the primitive variables from taking on unphysical values.
*/
int floorPrimitive(double* primVar);

/*
Compute the trace of the stress-energy tensor, T=3P-\rho
*/
void getStressEnergyTrace(double *T, double *primVar);
