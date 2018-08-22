#include "variableDef.h"
/*
Value below which nonnegative variables are not allowed to fall, used in floorConserved and floorPrimative
*/
#define FLOOR_VALUE (1.0e-12)

/*
Determine conserved variables from primative variables
*/
void primativeToConserved(double *primVar, double a, double *consVar);

/*
Determine primative variables from conserved variables
*/
void conservedToPrimative(double *consVar, double a, double *primVar);

/*
Compute physical flux from primative and conservative variables
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
Compute first physical flux for alternative form from primative and conservative variables
*/
void getPhysicalFluxAltOne(double *primVar, double *consVar, double *flux);

/*
Compute second physical flux for alternative form from primative and conservative variables
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
Impose a floor to keep the primative variables from taking on unphysical values.
*/
int floorPrimative(double* primVar);

