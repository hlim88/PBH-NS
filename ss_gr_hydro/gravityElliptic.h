
/*
Find L[phi] = phi''(r)-a^2(4(\pi)r(Sv-P)+m/r^2)
BC are phi(0)=0 and phi'(r_max)=ln(a(r_max))
*/
void LPhi(double *phi, int *phyBdy, double *consVarVertex, double *primVarVertex, double *a, double *r, double *cmask, int length, double* LPhi);

/*
Find residual = L[phi]-rhs
Also compute L2 norm of residual
*/
void residualPhi(double *phi, double *rhs, int *phyBdy, double *consVarVertex, double *primVarVertex, double *a, double *r, double *cmask, int length, double *norm, double* res);


/*
Perform red and black Gauss-Seidel relaxation
*/
void relaxPhi(double *phi, double *rhs, int *phyBdy, double *consVarVertex, double *primVarVertex, double *a, double *r, double *cmask, int length, double *norm);

