#define NUM_NEWT_ITER (400) //Maximum number of allowed Newton iterations
#define REL_ERROR (1.0e-9)  //Relative error tolerance in solution

#define CONST_PI 3.1415926535897932384626433

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//Some numbers of vars
#define n_vars 11
#define n_gen 4
#define n_deriv 9
//#define n_const 2

//Define var names
#define U_A 0
#define U_B 1
#define U_TRK 2
#define U_ARR 3
#define U_CHI 4
#define U_GAMDELTA 5
#define U_BETAR 6
#define U_ALPHA 7
#define U_BR 8
#define U_HAMC 9
#define U_MOMC 10

#if 0
typedef struct coll_point_t {

    double *u[n_gen][n_vars];
    double *cst[n_gen][n_const];
    double *du[n_gen][n_deriv];
    double *ddu[n_gen][n_deriv];
    double *rhs[n_gen][n_vars];
    unsigned time_stamp; //TODO : Need it?

} coll_point_t;
#endif

typedef struct coll_point_t {

    double *u[n_vars];
    double *du[n_deriv];
    double *ddu[n_deriv];
    double *rhs[n_vars];
    const int lenght;
    unsigned time_stamp; //TODO : Need it?

} coll_point_t;

extern unsigned time_stamp;
extern const int length;

/*
   Define BSSN eqns with usual 1+log slicing and gamma driver
*/

void solBSSNeqns(coll_point_t *pfunc, const double t, const double dt, const int gen, int iter);

//2nd order central FD for first derivs
static double derivFirst(double *u, double dr, int n) {

    double du;
    return du = (u[n+1]-u[n-1])/(2.0*dr);

}

//2nd order central FD for sencond derivs
static double derivSecond(double *u, double dr, int n) {

    double ddu;
    return ddu = (u[n-1]-2.0*u[n]+u[n+1])/(dr*dr);

}

//RK routines
void rk4_helper1(coll_point_t *pfunc, const double dt);
void rk4_helper2(coll_point_t *pfunc, const double dt);
void rk4_helper3(coll_point_t *pfunc, const double dt);
void rk4_helper4(coll_point_t *pfunc, const double dt);
