/*
 *  blackHoleFinder.c
 *  
 *
 *  Created by William East on 12/30/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
#define BH_THRESHOLD 0.995
#define FILE_NAME "critical_th"

#include <stdio.h>
#include <math.h>
#include "variableDef.h"
#include "blackHoleFinder.h"

/*
 Output char array to file.
 */
void outputString(char *output){
	char filename[] = FILE_NAME;
	FILE *filep = fopen(filename, "a");
	if(filep!=NULL){
		fprintf(filep,"%s", output);
		fclose(filep);
	} else{
		printf("Warning could not open %s!\n",filename);
	}
}

/*
 Finds black hole (defined as when 1-a^-2>BH_THRESHOLD) and mass.
 Returns 1 if black hole is found, 0 otherwise.
 The radius of the maximum of a is stored in rMax
 */
int findBlackHole(double *a, double *rGravity, int length, double *rMax, double *aMax){
	int i,j;
	double a_bh = sqrt(1.0/(1.0-BH_THRESHOLD));
	double a_max = a[0];
	int max_index = 0;
	for (i=1; i<length; i++) {
		if(a[i]>a_max){
			a_max=a[i];
			max_index = i;
		}
	}

	if(rMax!=NULL) *rMax=rGravity[max_index];	
	if(aMax!=NULL) *aMax=a_max;	

	if(a_max>a_bh){
		double mass = 0.5*rGravity[max_index]*(1.0-1.0/a_max/a_max);
		printf("Black hole detected at r=%E of mass %E\n",rGravity[max_index],mass);
		char output[80];
		snprintf(output, 79, "BH detected M=%E\n", mass);
		outputString(output);
		return 1;
	}
	
	return 0;
}

/*
 Returns the radius at which the maximum density occurs.
 */
double getMaxDensityRadius(double *consVar, double *r, int length){
	int i;
	int max_index=0;
	double max_density=consVar[CONS_DENSITY];
	for(i=1; i<length; i++){
		if(consVar[i*numVariables+CONS_DENSITY]>max_density){
			max_index = i;
			max_density = consVar[i*numVariables+CONS_DENSITY];
		}
	}
	return r[max_index];
}

