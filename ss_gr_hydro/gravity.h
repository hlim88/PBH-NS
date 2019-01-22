#define NUM_NEWT_ITER (400) //Maximum number of allowed Newton iterations
#define REL_ERROR (1.0e-9)  //Relative error tolerance in solution

/*
Solve Hamiltonian constraint for a using pointwise Newtonian iteration
*/
void solveHamiltonianConstraint(double *consVar, double *r, int length, double *a);

/*
Calculate the residual of the Hamiltonian constraint
*/
void residualHamiltonianConstraint(double *consVar, double *r, int length, double *a, double *residual);

/*
Solve polar-areal slicing condition for alpha
*/
void solveSlicingCondition(double *consVar, double *primVar, double *a, double *r, int length, double *alpha);

/*
Calculate the residual of the momentum equation (a evolution equation) using CN finite differencing.
*/
void calculateMomentumEqResidual(double *consVar_n, double *consVar_np1, double *a_n, double *a_np1, double *alpha_n, double *alpha_np1, double *r, double dt, int length, double *residual);

/*
Use momentum equation to evolve a
*/
void evolveMomentumEq(double *consVar, double *a, double *a_n, double *a_np1, double *alpha, double *r, double dt, int length, int innerBdy, int outerBdy);

// Aug.26.2018
// HL : Here we define different set of solving system 

/*
Solve 1+log slicing condition for alpha
*/
void solveAlphaSlicingCondition(double *consVar, double *primVar, double *r, 
                             int length, double *alpha, double *trK);

/*
Calculate the residual of the momentum equation (a evolution equation) using CN finite differencing.
*/
void calculateMomentumEqResidualHP(double *consVar_n, double *consVar_np1, double *a_n, 
                                 double *a_np1, double *alpha_n, double *alpha_np1, 
                                 double *r, double dt, int length, double *residual);

/*
Use momentum equation to evolve a
*/
void evolveMomentumEqHP(double *consVar, double *a, double *a_n, double *a_np1, 
                      double *alpha, double *r, double dt, int length, 
                      int innerBdy, int outerBdy);
