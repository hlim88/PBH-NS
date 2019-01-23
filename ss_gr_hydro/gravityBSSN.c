#define CONSTANT_PI 3.1415926535897932384626433

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "variableDef.h"
#include "equationOfState.h"
#include "gravity.h"
#include "rungeKutta.h"

/*
Solve Einstein's equation in spherically symmetry with HPC that
couples to ideal fluid
*/

void bssnrhs(double *u){

  #if 0
  /*TODO : Need to clean up
    This is usual full 3D BSSN but we are not using
    all vars. I just keep it anyway.
    This also include 1+log slicing with Gamma driver condition
  */     

  const double *alpha = &uZipVars[U_ALPHA][offset];
  const double *chi = &uZipVars[U_CHI][offset];
  const double *trK = &uZipVars[U_TRK][offset];
  const double *a = &uZipVars[U_A][offset];
  const double *b = &uZipVars[U_B][offset];
  const double *Arr = &uZipVars[U_ARR][offset];
  const double *GamDelta = &uZipVars[U_GAMDELTA][offset];
  const double *betar = &uZipVars[U_BETAR][offset];
  const double *Br = &uZipVars[U_BR][offset];

  const double *alpha_rhs = &uZipVars[U_ALPHA][offset];
  const double *chi_rhs = &uZipVars[U_CHI][offset];
  const double *trK_rhs = &uZipVars[U_TRK][offset];
  const double *a_rhs = &uZipVars[U_A][offset];
  const double *b_rhs = &uZipVars[U_B][offset];
  const double *Arr_rhs = &uZipVars[U_ARR][offset];
  const double *GamDelta_rhs = &uZipVars[U_GAMDELTA][offset];
  const double *betar_rhs = &uZipVars[U_BETAR][offset];
  const double *Br_rhs = &uZipVars[U_BR][offset];
  #endif

  // Const var for BSSN gauge
  const unsigned int BSSN_LAMBDA[4];
  const unsigned int BSSN_LAMBDA_F[4];

  const unsigned int lambda[4] = {BSSN_LAMBDA[0], BSSN_LAMBDA[1],
                                    BSSN_LAMBDA[2], BSSN_LAMBDA[3]};
  const double lambda_f[2] = {BSSN_LAMBDA_F[0], BSSN_LAMBDA_F[1]};

  #include "gravityBSSNeqns.h"
}

//Solving this with RK routine

void solBSSN() {

  //Apply RK4 routine
  #if 0
  for(int i=1; i<length; i++) {

      fourthOrderRKStep(); //TODO : Reduce inputs?

  }
  #endif



}

#if 0
void computeConstraints(double *u){

  #include "gravityBSSNeqns.h"
  #include "consteqns.h"

}
#endif
