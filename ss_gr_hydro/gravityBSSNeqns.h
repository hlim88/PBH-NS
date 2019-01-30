//index
int pp;

double a_rhs[pp], b_rhs[pp], trK_rhs[pp], Arr_rhs[pp], chi_rhs[pp],
       GamDelta_rhs[pp];

//Fake zero rhs for testing code
#if 1
a_rhs[pp] = 0.0;
b_rhs[pp] = 0.0;
trK_rhs[pp] = 0.0;
Arr_rhs[pp] = 0.0;
chi_rhs[pp] = 0.0;
GamDelta_rhs[pp] = 0.0;
#endif

//G-BSSN in SS with ideal fluid
#if 0

#endif

// TODO : Better way to put consteqns?

#include "consteqns.h"
