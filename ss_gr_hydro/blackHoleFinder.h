
/*
 *  blackHoleFinder.h
 *  
 *
 *  Created by William East on 12/30/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */


/*
 Output char array to file.
 */
void outputString(char *output);

/*
 Finds black hole (defined as when 1-a^-2>BH_THRESHOLD) and mass.
 Returns 1 if black hole is found, 0 otherwise.
 */
int findBlackHole(double *a, double *rGravity, int length, double *rMax, double *aMax);

/*
 Returns the radius at which the maximum density occurs.
 */
double getMaxDensityRadius(double *consVar, double *r, int length);
