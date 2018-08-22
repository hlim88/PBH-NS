#define MAX_ITER (200)  //The maximum number of iterations the root bracketer or finder
#define BRACKET_FACTOR (1.4) //The factor determining the range around the given value a root bracket is searched for 
#define ERROR_TOLERANCE (1.0e-12)  //The allowed error when root finding

/*
Bracket root by expanding geometrically
*/
int bracketRoot(double initialGuess, double* a, double* b, double (*func)(double, void *), void *param );

/*
Use Brent's method to find the zero of function func which lies between k0 and k1
*/
double findRoot(double k0, double k1, double (*func)(double, void*), void  *param );
