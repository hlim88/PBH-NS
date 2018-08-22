/*
Note by density we mean (1+specificEnergy)*restDensity.
The getDensity function takes pressure and returns density.
The value returned it the radius of the solution.
*/
double getTOVsolution(double centralPressure, double atmPressure, double (*getDensity)(double), double *pressure, double* a, double* alpha, double* r, int length);
