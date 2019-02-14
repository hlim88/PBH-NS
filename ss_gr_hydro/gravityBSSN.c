
#include "gravityBSSN.h"

//Define BSSN eqns for our system and applying simper RK4 here

unsigned time_stamp;

void solBSSNeqns(coll_point_t *pfunc, const double t, const double dt, const int gen) {

    void(*rk4_helpers[4])(coll_point_t *, const double) = {rk4_helper1, rk4_helper2, rk4_helper3, rk4_helper4};

    //index
    int pp;
    int n_pp; //TODO : need to link it with previous defs for PAMR/AMRD


    for(pp = 0; pp < n_pp; pp++) {
        //Polar-areal radius in HPC
        double r[pp];

        //Define vars
        double a, b, trK, Arr, chi, GamDelta, alpha, betaR, Br;

        //RHS terms def
        double a_rhs, b_rhs, trK_rhs, Arr_rhs, chi_rhs, GamDelta_rhs, betaR_rhs, Br_rhs, alpha_rhs;

        //First Derivative with resect to radius term
        double dr_a, dr_b, dr_trK, dr_Arr, dr_GamDelta, dr_betaR, dr_alpha, dr_chi, dr_Br;

        //Source term from fluid : TODO Link it explicitly
        double rhoFlu, SrrFlu, SthetathetaFlu, SrFlu;    

        //Call arrays to start compute
        double *rhs = pfunc->rhs[pp];
        double *u = pfunc->u[pp];
        double *du = pfunc->du[pp];
        double *ddu = pfunc->ddu[pp];

        //Save it into array;
        a = u[U_A];
        b = u[U_B];
        trK = u[U_TRK];
        Arr = u[U_ARR];
        chi = u[U_CHI];
        GamDelta = u[U_GAMDELTA];
        betaR = u[U_BETAR];
        alpha = u[U_ALPHA];
        Br = u[U_BR];

        //Save first derivative
        dr_a = du[U_A];
        dr_b = du[U_B];
        dr_trK = du[U_TRK];
        dr_Arr = du[U_ARR];
        dr_chi = du[U_CHI];
        dr_GamDelta = du[U_GAMDELTA];
        dr_betaR = du[U_BETAR];
        dr_alpha = du[U_ALPHA];
        dr_Br = du[U_BR];
        
        //Grid spacing 
        //TODO : Do we need N?
        double dr;
        dr = r[pp]-r[pp-1];

        *du = derivFirst(u, dr, pp); //Correct?

        //Second derivative def
        double drdr_alpha, drdr_b, drdr_chi, drdr_a, drdr_betaR, drdr_trK, drdr_Arr, drdr_Br, 
               drdr_GamDelta;

        //Save second derivative
        drdr_a = ddu[U_A];
        drdr_b = ddu[U_B];
        drdr_trK = ddu[U_TRK];
        drdr_Arr = ddu[U_ARR];
        drdr_chi = ddu[U_CHI];
        drdr_GamDelta = ddu[U_GAMDELTA];
        drdr_betaR = ddu[U_BETAR];
        drdr_alpha = ddu[U_ALPHA];
        drdr_Br = ddu[U_BR];

        *ddu = derivSecond(u, dr, pp);

        //Save rhs expressions
        a_rhs = rhs[U_A];
        b_rhs = rhs[U_B];
        trK_rhs = rhs[U_TRK];
        Arr_rhs = rhs[U_ARR];
        chi_rhs = rhs[U_CHI];
        GamDelta_rhs = rhs[U_GAMDELTA];
        betaR_rhs = rhs[U_BETAR];
        alpha_rhs = rhs[U_ALPHA];
        Br_rhs = rhs[U_BR];

        //Fake zero rhs for testing code
        
        #if 0
        a_rhs = 0.0;
        b_rhs = 0.0;
        trK_rhs = 0.0;
        Arr_rhs = 0.0;
        chi_rhs = 0.0;
        GamDelta_rhs = 0.0;
        #endif

        //Non vanishing Ricci
        double Ricci_rr, Ricci_thetatheta, Ricci_scal;
        Ricci_rr = - 2.0/(a*r[pp]*b)*(b/a - r[pp]*b*dr_a/(a*a)
                   + r[pp]*dr_b/a);
        Ricci_thetatheta = (a - b*(b + r[pp]*dr_b)/a
                            - r[pp]*b*dr_a*(b+r[pp]*dr_b)/(a*a)
                            + r[pp]*dr_b*(b + r[pp]*dr_b)/a
                            + r[pp]*b*(2.0*dr_b + r[pp] * drdr_b)/a)
                            /(a*r[pp]*r[pp]*b*b);
        Ricci_scal = Ricci_rr + Ricci_thetatheta;
        
        
        //Define some common derivative terms for convenience
        double DivBetaR;
        DivBetaR = dr_betaR 
                     - betaR*((dr_a*b*b + 2.0*a*b*dr_b)
                                   /(2.0*a*b*b) + 2.0/r[pp]);
        double DelSqAlpha;
        DelSqAlpha = chi*chi/alpha
                       *(drdr_alpha - drdr_a/(2.0*a) + drdr_b/b
                         - drdr_chi/chi + dr_a*dr_a/(2.0*a*a)
                         - dr_b*dr_b/(b*b) + dr_chi*dr_chi/(chi*chi)
                         - 2.0/(r[pp]*r[pp]));
        double DelDelAlpha;
        DelDelAlpha = chi*chi/alpha 
                        *(drdr_alpha - drdr_a/(2.0*a) + drdr_chi/chi
                          + dr_a*dr_a/(2.0*a*a) 
                          - dr_chi*dr_chi/(chi*chi));

        //G-BSSN in SS with ideal fluid
        #if 1
        chi_rhs = betaR*dr_chi - chi/3.0*(alpha*trK - DivBetaR);
        a_rhs = betaR*dr_a + 2.0*a*dr_betaR 
                  - 2.0*a/3.0*DivBetaR - 2.0*alpha*a*Arr;
        b_rhs = betaR*dr_b + 2.0*b*betaR/r[pp] 
                  - 2.0*b/3.0*DivBetaR + alpha*b*Arr;
        trK_rhs = betaR*dr_trK - DelSqAlpha
                    + alpha * (3.0/2.0*Arr*Arr + trK*trK/3.0)
                    + 4*CONST_PI*alpha*(rhoFlu + SrrFlu + SthetathetaFlu);
        Arr_rhs = betaR*dr_Arr -(DelDelAlpha- DelSqAlpha/3.0) 
                    + alpha/3.0 * (2.0*Ricci_rr - Ricci_thetatheta)
                    + alpha*trK*Arr 
                    - 16.0*CONST_PI*alpha*(SrrFlu + SthetathetaFlu);
        GamDelta_rhs = betaR*dr_GamDelta - GamDelta*dr_betaR
                         + drdr_betaR/a + 2.0*(dr_betaR - betaR/r[pp])/(b*r[pp])
                         + ( (drdr_betaR 
                             + dr_betaR*( (dr_a*b*b + 2.0*a*b*dr_b)
                                              /(2.0*a*b*b) + 2.0/r[pp] 
                                            )
                             + betaR*( (b*drdr_b - dr_b*dr_b)/(b*b) 
                                           +(a*drdr_a - dr_a*dr_a)/(2.0*a*a)
                                           -2.0/(r[pp]*r[pp])
                                         )
                             )/a
                           + 2.0*GamDelta*DivBetaR/3.0)
                         - 2.0*(Arr*dr_alpha + alpha*dr_Arr)/a
                         + 2.0*alpha*(Arr*GamDelta + 3.0*Arr/(r[pp]*b))
                         + 2.0*alpha/a 
                           *( dr_Arr - 2.0/3.0*dr_trK - 3.0*Arr/chi*dr_chi
                             +3.0/2.0*Arr*(2.0/r[pp] + dr_b/b) - 8*CONST_PI*SrFlu
                            );
        
        #endif
        //Gauges. Here we use 1+log and Gamma driver

        alpha_rhs = - 2.0*alpha*trK;
        Br_rhs = 3.0/4.0*GamDelta_rhs;
        betaR_rhs = Br;
        
        // Applying solver here
        void(*rk4_helper)(coll_point_t *, const double) = rk4_helpers[pp];
        
        for (int i = 0; i<npts_in_array; i++) {
            //coll_point_t *pfunc = &coll_points->array[i];
            rk4_helper(pfunc, dt);
        }
        // TODO : Better way to put consteqns?

        // Constraint equations
        // Constraints are not evoling with respect to time. It receives the value for each
        // vars and evaluate the value
        double hamC_rhs, momC_rhs;

        //Fake zeros
        #if 0
        hamC_rhs = 0.0;
        momC_rhs = 0.0;
        #endif

        hamC_rhs = Ricci_scal - 3.0/2.0*Arr*Arr + 2.0/3.0*trK*trK 
                     - 16.0*CONST_PI*rhoFlu;
        momC_rhs = dr_Arr - 2.0/3.0*dr_trK - 3.0*Arr/(2.0*chi)*dr_chi
                     + 2.0/3.0*Arr*(2.0/r[pp] - dr_b/b) - 8.0*CONST_PI*rhoFlu;
    
    }   
    
}

void rk4_helper1(coll_point_t *pfunc, const double dt) {

    double *u_curr = pfunc->u[0];
    double *u_next = pfunc->u[1];
    double *rhs = pfunc->rhs[0];

    for (int i = 0; i < n_vars; i++)
        u_next[i] = u_curr[i] + 0.5*dt*rhs[i];

}


void rk4_helper2(coll_point_t *pfunc, const double dt) {

    double *u_curr = pfunc->u[0];
    double *u_next = pfunc->u[2];
    double *rhs = pfunc->rhs[1];

    for (int i = 0; i < n_vars; i++)
        u_next[i] = u_curr[i] + 0.5*dt*rhs[i];

}


void rk4_helper3(coll_point_t *pfunc, const double dt) {

    double *u_curr = pfunc->u[0];
    double *u_next = pfunc->u[3];
    double *rhs = pfunc->rhs[2];

    for (int i = 0; i < n_vars; i++)
        u_next[i] = u_curr[i] + 0.5*dt*rhs[i];

}


void rk4_helper4(coll_point_t *pfunc, const double dt) {

    double *u_curr = pfunc->u[0];
    double *k1 = pfunc->rhs[0];
    double *k2 = pfunc->rhs[1];
    double *k3 = pfunc->rhs[2];
    double *k4 = pfunc->rhs[3];

    for (int i = 0; i < n_vars; i++)
        u_curr[i] += dt*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i])/6.0;

}

