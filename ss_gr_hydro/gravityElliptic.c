#define CONSTANT_PI 3.1415926535897932384626433
#define CMASK_ON 1

#include <stdio.h>
#include <math.h>

#include "gravityElliptic.h"
#include "equationOfState.h"
#include "conservedPrimitive.h"

/*
Find L[phi] = phi''(r)-a^2(4(\pi)r(Sv-P)+m/r^2)
BC are phi(0)=0 and phi'(r_max)=-ln(a(r_max))
*/
void LPhi(double *phi, int *phyBdy, double *consVarVertex, double *primVarVertex, double *a, double *r, double *cmask, int length, double* LPhi){

	double dr = r[1]-r[0];

	int i;
	for(i=1; i<length-1; i++){
		if(((int)( cmask[i]+0.5))==CMASK_ON){
			double S = 0.5*(consVarVertex[(i)*numVariables+PI]-consVarVertex[(i)*numVariables+PHI]);
		//	conservedToPrimitive((consVarVertex+i*numVariables), a[i], (primVarVertex+i*numVariables));
			double pressure = getPressure(primVarVertex[i*numVariables+REST_DENSITY], primVarVertex[i*numVariables+SPECIFIC_ENERGY]);
		
			LPhi[i] = (phi[i+1]-2.0*phi[i]+phi[i-1])/dr/dr-a[i]*a[i]*4.0*CONSTANT_PI*r[i]*(S*primVarVertex[i*numVariables+VELOCITY] + pressure)-0.5*(a[i]*a[i]-1.0)/r[i]; 
		}
	}
	
	//Boundary Condtions
	if(phyBdy[0] && ((int)(cmask[0]+0.5))==CMASK_ON) LPhi[0] = phi[0];
	if(phyBdy[1] && ((int)(cmask[length-1]+0.5))==CMASK_ON) LPhi[length-1] = (1.5*phi[length-1]-2.0*phi[length-2]+0.5*phi[length-3])/dr +log(a[length-1]);
	
	
}

/*
Find residual = L[phi]-rhs
Also compute L2 norm of residual
*/
void residualPhi(double *phi, double *rhs, int *phyBdy, double *consVarVertex, double *primVarVertex, double *a, double *r, double *cmask, int length, double *norm, double* res){

	LPhi(phi, phyBdy, consVarVertex, primVarVertex, a, r, cmask, length, res);

	//Enforce physical boundaries if necessary
	int start, end;
	if(phyBdy[0]) start = 0;
	else start = 1; 	
	if(phyBdy[1]) end = length;
	else end = length-1; 	

	int N = 0;
	*norm = 0;


	int i;
	for(i=start; i<end; i++){
		if(((int) (cmask[i]+0.5))==CMASK_ON){
			res[i]-=rhs[i];
			*norm+=res[i]*res[i];
			N++;

		}		
	}

	if(N!=0) *norm = sqrt(*norm/N); else *norm = 0;

}

/*
Perform red and black Gauss-Seidel relaxation
*/
void relaxPhi(double *phi, double *rhs, int *phyBdy, double *consVarVertex, double *primVarVertex, double *a, double *r, double *cmask, int length, double *norm){

	double dr = r[1]-r[0];
	double jac = -2.0/dr/dr;

	double res; 
	*norm=0;
	int N=0;

	//DEBUG
	double max_norm=0.0;
	int max_index=0;
	
	//Boundary Condtions
	if(phyBdy[0] && ((int)(cmask[0]+0.5))==CMASK_ON){ 
		res = phi[0] - rhs[0];
		phi[0]-=res;
		*norm+=res*res;
		N++;
	}
	if(phyBdy[1] && ((int)(cmask[length-1]+0.5))==CMASK_ON){
		res = (1.5*phi[length-1]-2.0*phi[length-2]+0.5*phi[length-3])/dr +log(a[length-1]);
		phi[length-1]-=res/(1.5/dr);
		*norm+=res*res;
		N++; 
	}

	int parity,i;
	for(parity=0; parity<2; parity++){
		for(i=1; i<length-1; i++){
			if(((int)(cmask[i]+0.5))==CMASK_ON && (i+parity)%2==0){
				double S = 0.5*(consVarVertex[(i)*numVariables+PI]-consVarVertex[(i)*numVariables+PHI]);
			//	conservedToPrimitive((consVarVertex+i*numVariables), a[i], (primVarVertex+i*numVariables));
				double pressure = getPressure(primVarVertex[i*numVariables+REST_DENSITY], primVarVertex[i*numVariables+SPECIFIC_ENERGY]);

				res = (phi[i+1]-2.0*phi[i]+phi[i-1])/dr/dr-a[i]*a[i]*4.0*CONSTANT_PI*r[i]*(S*primVarVertex[i*numVariables+VELOCITY] + pressure)-0.5*(a[i]*a[i]-1.0)/r[i] -rhs[i]; 
				phi[i]-=res/jac;	
				*norm+=res*res;
				N++;
				if(res*res>max_norm){ max_norm=res*res; max_index=i;}
				
			}
		}
	}
	if(N!=0) *norm = sqrt(*norm/N); else *norm = 0;

	//DEBUG

	i = max_index;
	double density = primVarVertex[i*numVariables+REST_DENSITY];
	double specEnergy = primVarVertex[i*numVariables+SPECIFIC_ENERGY];
	double velocity = primVarVertex[i*numVariables+VELOCITY];
	double S = 0.5*(consVarVertex[(i)*numVariables+PI]-consVarVertex[(i)*numVariables+PHI]);
	double pressure = getPressure(primVarVertex[i*numVariables+REST_DENSITY], primVarVertex[i*numVariables+SPECIFIC_ENERGY]);
	//printf("Relax Mres =%E with norm %E i=%d out of %d r=%E with d=%E e=%E v=%E S %E P=%E phi = %E %E %E a= %E\n",max_norm,*norm,max_index,length,r[max_index],density,specEnergy,velocity,S,pressure,phi[i+1],phi[i],phi[i-1],a[i]);

}

