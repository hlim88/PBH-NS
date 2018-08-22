/*
 *  rungeKutta.c
 *  
 *
 *  Created by William East on 10/26/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include <stdlib.h>
#include <stdio.h>
/*
Take a fourth order Runge-Kutta step of size 'step' from startValue.  The dimension of the system is 'dimension'.
 deriv(x,y,dy/dx, auxFunc) stores the derivative at x, y to the pointer dy/dx and may use the aux function auxFunc.
 */
void fourthOrderRKStep(double* startValue, double x, double* endValue, double* startDeriv, int dimension, double step, void (*deriv) (double, double*, double*, double(*)(double)) , double (*auxFunc)(double)){
	
	double *initialDeriv;
	if(startDeriv==NULL){
		initialDeriv = malloc(dimension*sizeof(double));
		deriv(x, startValue, initialDeriv, auxFunc);
	} else {
		initialDeriv = startDeriv;
	}
	
	double a[dimension];
	double b[dimension];
	double c[dimension];
	double d[dimension];
	double argument[dimension];
	
	int i;
	for(i=0; i<dimension; i++){
		a[i] = step*initialDeriv[i];
		argument[i]=startValue[i]+0.5*a[i];
	}
	
	deriv((x+0.5*step), argument, b, auxFunc);
	
	for(i=0; i<dimension; i++){
		b[i] *= step;
		argument[i]=startValue[i]+0.5*b[i];
	}
	
	deriv((x+0.5*step), argument, c, auxFunc);
	
	for(i=0; i<dimension; i++){
		c[i] *= step;
		argument[i]=startValue[i]+c[i];
	}
	
	deriv((x+step), argument, d, auxFunc);
	
	for(i=0; i<dimension; i++){
		d[i] *= step;
		endValue[i] = startValue[i] + a[i]/6.0 + b[i]/3.0 + c[i]/3.0 + d[i]/6.0;
	}
	
	if(startDeriv==NULL){
		free(initialDeriv);
	}

}
