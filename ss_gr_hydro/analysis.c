/*

Brief : Some analysis tools for comparision to validate the code

*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "TOV.h"

void unpertTOV_comp(){

  //Put simple analytic form of TOV star


}

double computeConstraintL2(const double *constraintVec, const double *maskVec, 
                    unsigned int lbegin, unsigned int lend){

  for (unsigned int k = 3; k < nz - 3; k++) {
    double z = pmin[2] + k*hz;
    for (unsigned int j = 3; j < ny - 3; j++) {
      double y = pmin[1] + j*hy;
      for (unsigned int i = 3; i < nx - 3; i++) {
        double x = pmin[0] + i*hx;
        unsigned int pp = i + nx * (j + ny * k);

        double l2=0.0;
        double l2_g=0.0;
        const double MASK_THRESHOLD=1.0;
        for(unsigned int i=lbegin;i<lend;i++)
        {
            if(maskVec[i]<MASK_THRESHOLD)
                l2+=(constraintVec[i]*constraintVec[i]*maskVec[i]*maskVec[i]
                   *maskVec[i]*maskVec[i]);
            else
                l2+=(constraintVec[i]*constraintVec[i]);
        }

        return (sqrt(l2_g));

}

double printL2Constraint(const double **constraintVec, unsigned int timeSteps){

  const unsigned int numConstVars; //4?

  double constraintL2[numConstVars];

  for (unsigned int id = 0; id < numConstVars; id++) {
      constraintL2[id] = computeConstraintL2();
      printf("Constraint : %lf %lf\n",constraintL2[id]);
  }

  //TODO : Add routine to generate output using sprintf?
  #if 0
    ofstream fileConst;
    char fname[256];
    sprintf(fname,"%s_constraint.dat",val_name_
    if(timeSteps == 0)
    fileConst.open(fname)
    printf("Time steps : %d",timeSteps);
    printf("HAM : %lf",constraintL2[0]);
    printf("MOM1 : %lf",constraintL2[1]; //just 1d?
    fileConst.close();

    return 0;

  #endif


}

