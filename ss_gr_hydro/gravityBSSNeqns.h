#define CONST_PI 3.1415926535897932384626433
//index
int pp;

//Polar-areal radius in HPC
double r[pp];

//Define vars
double a[pp], b[pp], trK[pp], Arr[pp], chi[pp], GamDelta[pp], betaR[pp], alpha[pp];

//RHS terms def
double a_rhs[pp], b_rhs[pp], trK_rhs[pp], Arr_rhs[pp], chi_rhs[pp],
       GamDelta_rhs[pp];

//First Derivative with respect to radius term
double dr_a[pp], dr_b[pp], dr_trK[pp], dr_Arr[pp], dr_GamDelta[pp], 
       dr_betaR[pp], dr_alpha[pp], dr_chi[pp];

//2nd order centeral FD for first dervis
// Meshing unit
double dr,h; 
dr = r[pp] - r[pp-1];
h = dr*dr;
dr_a[pp] = (a[pp+1]-a[pp-1])/(2.0*h);
dr_b[pp] = (b[pp+1]-b[pp-1])/(2.0*h);
dr_trK[pp] = (trK[pp+1]-trK[pp-1])/(2.0*h);
dr_Arr[pp] = (Arr[pp+1]-Arr[pp-1])/(2.0*h);
dr_GamDelta[pp] = (GamDelta[pp+1]-GamDelta[pp-1])/(2.0*h);
dr_betaR[pp] = (betaR[pp+1]-betaR[pp-1])/(2.0*h);
dr_alpha[pp] = (alpha[pp+1]-alpha[pp-1])/(2.0*h);
dr_chi[pp] = (chi[pp+1]-chi[pp-1])/(2.0*h);

//Second derivative def
double drdr_alpha[pp], drdr_b[pp], drdr_chi[pp], drdr_a[pp], drdr_betaR[pp];

//2nd order central FD for 2nd derivs
drdr_alpha[pp] = (alpha[pp-1] - 2.0*alpha[pp] + alpha[pp+1])/(h*h);
drdr_b[pp] = (b[pp-1] - 2.0*b[pp] + b[pp+1])/(h*h);
drdr_chi[pp] = (chi[pp-1] - 2.0*chi[pp] + chi[pp+1])/(h*h);
drdr_a[pp] = (a[pp-1] - 2.0*a[pp] + a[pp+1])/(h*h);
drdr_betaR[pp] = (betaR[pp-1] - 2.0*betaR[pp] + betaR[pp+1])/(h*h);

//Source term from fluid TODO : Link it explicitly
double rhoFlu[pp], SrrFlu[pp], SthetathetaFlu[pp], SrFlu[pp];

//Fake zero rhs for testing code
#if 0
a_rhs[pp] = 0.0;
b_rhs[pp] = 0.0;
trK_rhs[pp] = 0.0;
Arr_rhs[pp] = 0.0;
chi_rhs[pp] = 0.0;
GamDelta_rhs[pp] = 0.0;
#endif

//Non vanishing Ricci
double Ricci_rr, Ricci_thetatheta, Ricci_scal;
Ricci_rr = - 2.0/(a[pp]*r[pp]*b[pp]) 
         *(b[pp]/a[pp] - r[pp]*b[pp]*dr_a[pp]/(a[pp]*a[pp])
           + r[pp]*dr_b[pp]/a[pp]);
Ricci_thetatheta = (a[pp] - b[pp]*(b[pp] + r[pp]*dr_b[pp])/a[pp]
                    - r[pp]*b[pp]*dr_a[pp]*(b[pp]+r[pp]*dr_b[pp])/(a[pp]*a[pp])
                    + r[pp]*dr_b[pp]*(b[pp] + r[pp]*dr_b[pp])/a[pp]
                    + r[pp]*b[pp]*(2.0*dr_b[pp] + r[pp] * drdr_b[pp])/a[pp])
                    /(a[pp]*r[pp]*r[pp]*b[pp]*b[pp]);
Ricci_scal = Ricci_rr + Ricci_thetatheta;

//Define some common derivative terms for convenience
double DivBetaR[pp];
DivBetaR[pp] = dr_betaR[pp] 
             - betaR[pp]*((dr_a[pp]*b[pp]*b[pp] + 2.0*a[pp]*b[pp]*dr_b[pp])
                           /(2.0*a[pp]*b[pp]*b[pp]) + 2.0/r[pp]);
double DelSqAlpha[pp];
DelSqAlpha[pp] = chi[pp]*chi[pp]/alpha[pp]
               *(drdr_alpha[pp] - drdr_a[pp]/(2.0*a[pp]) + drdr_b[pp]/b[pp]
                 - drdr_chi[pp]/chi[pp] + dr_a[pp]*dr_a[pp]/(2.0*a[pp]*a[pp])
                 - dr_b[pp]*dr_b[pp]/(b[pp]*b[pp]) + dr_chi[pp]*dr_chi[pp]/(chi[pp]*chi[pp])
                 - 2.0/(r[pp]*r[pp]));
double DelDelAlpha[pp];
DelDelAlpha[pp] = chi[pp]*chi[pp]/alpha[pp] 
                *(drdr_alpha[pp] - drdr_a[pp]/(2.0*a[pp]) + drdr_chi[pp]/chi[pp]
                  + dr_a[pp]*dr_a[pp]/(2.0*a[pp]*a[pp]) 
                  - dr_chi[pp]*dr_chi[pp]/(chi[pp]*chi[pp]));

//G-BSSN in SS with ideal fluid
#if 1
chi_rhs[pp] = betaR[pp]*dr_chi[pp] - chi[pp]/3.0*(alpha[pp]*trK[pp] - DivBetaR[pp]);
a_rhs[pp] = betaR[pp]*dr_a[pp] + 2.0*a[pp]*dr_betaR[pp] 
          - 2.0*a[pp]/3.0*DivBetaR[pp] - 2.0*alpha[pp]*a[pp]*Arr[pp];
b_rhs[pp] = betaR[pp]*dr_b[pp] + 2.0*b[pp]*betaR[pp]/r[pp] 
          - 2.0*b[pp]/3.0*DivBetaR[pp] + alpha[pp]*b[pp]*Arr[pp];
trK_rhs[pp] = betaR[pp]*dr_trK[pp] - DelSqAlpha[pp] 
            + alpha[pp] * (3.0/2.0*Arr[pp]*Arr[pp] + trK[pp]*trK[pp]/3.0)
            + 4*CONST_PI*alpha[pp]*(rhoFlu[pp] + SrrFlu[pp] + SthetathetaFlu[pp]);
Arr_rhs[pp] = betaR[pp]*dr_Arr[pp] -(DelDelAlpha[pp]- DelSqAlpha[pp]/3.0) 
            + alpha[pp]/3.0 * (2.0*Ricci_rr - Ricci_thetatheta)
            + alpha[pp]*trK[pp]*Arr[pp] 
            - 16.0*CONST_PI*alpha[pp]*(SrrFlu[pp] + SthetathetaFlu[pp]);
GamDelta_rhs[pp] = betaR[pp]*dr_GamDelta[pp] - GamDelta[pp]*dr_betaR[pp]
                 + drdr_betaR[pp]/a[pp] + 2.0*(dr_betaR[pp]- betaR[pp]/r[pp])/(b[pp]*r[pp])
                 + ( (drdr_betaR[pp] 
                     + dr_betaR[pp]*( (dr_a[pp]*b[pp]*b[pp] + 2.0*a[pp]*b[pp]*dr_b[pp])
                                      /(2.0*a[pp]*b[pp]*b[pp]) + 2.0/r[pp] 
                                    )
                     + betaR[pp]*( (b[pp]*drdr_b[pp]-dr_b[pp]*dr_b[pp])/(b[pp]*b[pp]) 
                                   +(a[pp]*drdr_a[pp]-dr_a[pp]*dr_a[pp])/(2.0*a[pp]*a[pp])
                                   -2.0/(r[pp]*r[pp])
                                 )
                     )/a[pp]
                   + 2.0*GamDelta[pp]*DivBetaR[pp])/3.0
                 - 2.0*(Arr[pp]*dr_alpha[pp] + alpha[pp]*dr_Arr[pp])/a[pp]
                 + 2.0*alpha[pp]*(Arr[pp]*GamDelta[pp] + 3.0*Arr[pp]/(r[pp]*b[pp]))
                 + 2.0*alpha[pp]/a[pp] 
                   *( dr_Arr[pp] - 2.0/3.0*dr_trK[pp] - 3.0*Arr[pp]/chi[pp]*dr_chi[pp]
                     +3.0/2.0*Arr[pp]*(2.0/r[pp] + dr_b[pp]/b[pp]) - 8*CONST_PI*SrFlu[pp]
                    );
#endif

//Gauges. Here we use 1+log and Gamma driver
double alpha_rhs[pp], Br_rhs[pp], betaR_rhs[pp];

alpha_rhs[pp] = - 2.0*alpha[pp]*trK[pp];
Br_rhs[pp] = 3.0/4.0*GamDelta_rhs[pp];


// TODO : Better way to put consteqns?

// Constraint equations
// Constraints are not evoling with respect to time. It receives the value for each
// vars and evaluate the value
double hamC_rhs[pp], momC_rhs[pp];

//Fake zeros
#if 0
hamC_rhs[pp] = 0.0;
momC_rhs[pp] = 0.0;
#endif

hamC_rhs[pp] = Ricci_scal - 3.0/2.0*Arr[pp]*Arr[pp] + 2.0/3.0*trK[pp]*trK[pp] 
             - 16.0*CONST_PI*rhoFlu[pp];
momC_rhs[pp] = dr_Arr[pp] - 2.0/3.0*dr_trK[pp] - 3.0*Arr[pp]/(2.0*chi[pp])*dr_chi[pp]
             + 2.0/3.0*Arr[pp]*(2.0/r[pp] - dr_b[pp]/b[pp]) - 8.0*CONST_PI*rhoFlu[pp];
