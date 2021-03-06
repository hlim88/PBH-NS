/*
Take a fourth order Runge-Kutta step of size 'step' from startValue.  The dimension of the system is 'dimension'.
 deriv(x,y,dy/dx) stores the derivative at x, y to the pointer dy/dx.
 */
void fourthOrderRKStep(double* startValue, double x, double* endValue, double* startDeriv, 
                       int dimension, double step, 
                       void (*deriv) (double, double*, double*, double(*)(double)),
                       double (*auxFunc)(double));

/*
Partially Implicit Runge-Kutta Sovler 
This might be better to sovle G-BSSN in SS which is system we choose
*/
void PIRK2Step(double* startValue, double x, double* endValue, double* startDeriv, 
               int dimension, double step, 
               void (*deriv) (double, double*, double*, double(*)(double)),
               double (*auxFunc)(double));
