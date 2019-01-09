/*

Brief : Some analysis tools for comparision to validate the code

*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void unpertTOV_comp(){


}

void constraintcheck(){

  for (unsigned int k = 3; k < nz - 3; k++) {
    double z = pmin[2] + k*hz;
    for (unsigned int j = 3; j < ny - 3; j++) {
      double y = pmin[1] + j*hy;
      for (unsigned int i = 3; i < nx - 3; i++) {
        double x = pmin[0] + i*hx;
        unsigned int pp = i + nx * (j + ny * k);

}

