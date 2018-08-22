#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "brent.h"

/*
Bracket root by expanding geometrically.
Bracket boundaries must always be positive
If bracket is found return 0, otherwise return 1
Param is a pointer to some arbitrary fixed parameters that func() needs.
*/
int bracketRoot(double initialGuess, double* a, double* b, double (*func)(double, void *), void *param ){

	if(initialGuess<=0){
		printf("Error in bracketRoot().  Initial guess must be positive!\n");
		return 1;
	}
	
	int Nmax = MAX_ITER;
	double factor = BRACKET_FACTOR;
	int i;
	double fa, fb;	

	//First try method for roots that are close to initial guess
	double eps = log(0.5*ERROR_TOLERANCE);	
	for(i=1; i<=2; i++){
		*a=initialGuess/(1.0+exp(eps/i));
		*b=initialGuess*(1.0+exp(eps/i));
		fa = func(*a, param);
		fb = func(*b, param);
		if(fa*fb<0.0){
			return 0;
		}
	} 	

	//Now try method for roots that are further away
	*a = initialGuess/factor;
	*b = initialGuess*factor;
	fa = func(*a, param);
	fb = func(*b, param);
	
	for(i=0; i<Nmax; i++){
		if(fa*fb<0.0){
			return 0;
		}
		if(fabs(fa)<fabs(fb)){
			*a/=factor;
			fa = func(*a, param);
		} else {
			*b*=factor;
			fb = func(*b, param);
		}
	}

	//Failed to find bracket around root
	return 1;
}


/*
Use Brent's method to find the zero of function func which lies between k0 and k1
Param is a pointer to some arbitrary fixed parameters that func() needs.
*/
double findRoot(double k0, double k1, double (*func)(double, void *), void *param ){
  
  double errorTolerance = ERROR_TOLERANCE;
  int Nmax = MAX_ITER;
  int n=0;
  
  double a;
  double b;
  double c;
  double d=0;
  double f_a = func(k0, param);
  double f_b = func(k1, param);
  double f_c;
  double f_s;
  double s;
  int flag=1;
  
  if(f_a*f_b>=0.0){
  	printf("Error! No zero crossing.\n");
	printf("fa=%E fb=%E a=%E b=%E\n",f_a,f_b,k0,k1);
  	return (0.5*(k0+k1));
  }

  if(fabs(f_a)<fabs(f_b)){
	//Switch a and b
	a=k1;
	double buffer=f_a;
	f_a=f_b;
	b=k0;
	f_b=buffer;
  } else {
	a=k0;
	b=k1;
  }
  
  c=a;
  f_c=f_a;
  
  while(f_b!=0 && fabs(b-a)>(errorTolerance)*fabs(b) && n<Nmax){

  
    if(f_a!=f_c && f_b!=f_c){
	//Inverse quadratic interpolation
	s=a*f_b*f_c/(f_a-f_b)/(f_a-f_c)+b*f_a*f_c/(f_b-f_a)/(f_b-f_c)+c*f_a*f_b/(f_c-f_a)/(f_c-f_b);            
    } else{
	//Secant method
	s=b-f_b*(b-a)/(f_b-f_a);
    }
   
  
   if(((s<(3.0*a+b)/4.0) && (s<b)) || ((s>(3.0*a+b)/4.0) && (s>b)) || (flag && fabs(s-b)>=0.5*fabs(b-c)) || (!flag && fabs(s-b)>=0.5*fabs(c-d)) ){
	//Bisection method
	s=0.5*(a+b);
	flag=1;
   } else {
   	flag=0;
   }  
  
   d=c;
   c=b;
   f_c=f_b;
   
   f_s=func(s, param);
   
   if(f_a*f_s<0){
   	b=s;
   	f_b=f_s;
   }
   else {
   	a=s;
   	f_a=f_s;
   }
   
   if(fabs(f_a)<fabs(f_b)){
   	//Switch a and b
   	double buffer1 =a;
  	double buffer2=f_a;
  	a=b;
   	f_a=f_b;
   	b=buffer1;
   	f_b=buffer2;
   }
   
   n++;
   
  }//while
  if(n==Nmax){
  	printf("Desired precision not obtained! a=%e, b=%e, |a-b|=%e\n",a,b,fabs(b-a));
  }
  return b;
}//findRoot

