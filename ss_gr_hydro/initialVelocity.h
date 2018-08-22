

/*
Given an initial coordinate velocity and primitive variables, solve for metric functions a and alpha using pointwise 2-d Newton iteration (and hence get velocity).
a(0) and alpha(0) should already be set
*/
void solveMetricCoordVelocity(double *coordVelocity, double *primVar, double *r, int length, double *a, double *alpha);
