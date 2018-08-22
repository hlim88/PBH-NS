/*
2nd order Runge-Kutta evolution.  
If iter=1 we take the half step:
q_np1=q_n+0.5*dt*F(q_n)
If iter=2 we take the whole step:
q_np1=q_n+dt*F(q_np1)
We assume that a, alpha, rGravity have length entries.
We assume that consVar and primVar have numVariables*(length-1) entries.
if(boolInnerBdy) then the inner boundary is a physical boundary.
if(boolOuterBdy) then the outer boundary is a physical boundary.
*/
void timeStep(int iter, int boolInnerBdy, int boolOuterBdy, double dt, int length, double *rGravity, double* consVar_n, double* consVar_np1, double* primVar_n, double* primVar_np1, double *a_n, double *a_np1,  double *phi_n, double *phi_np1, double *fluxCorrection, int *fluxCorrectionMask);
