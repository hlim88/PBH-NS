#define NUM_NEWT_ITER (400) //Maximum number of allowed Newton iterations
#define REL_ERROR (1.0e-9)  //Relative error tolerance in solution

/*
   Define BSSN eqns with usual 1+log slicing and gamma driver
*/

void bssnrhs(double *u);

void solBSSN(double *u, double rkstage[]);

