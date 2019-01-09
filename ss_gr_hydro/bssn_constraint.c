/*
    Enforce physical constraints on BSSN vars
    det(gt) = 1, tr(At) = 0, alpha > 0, and chi > 0
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void force_bssn_constraints(){

  const double one_third = 1.0 / 3.0;
  //TODO : Change var name associated with code var name
  double gtd[3][3], Atd[3][3];

  double det_gtd  =  gtd[0][0]*( gtd[1][1]*gtd[2][2] - gtd[1][2]*gtd[1][2])
                  -  gtd[0][1]*gtd[0][1]*gtd[2][2]
                  +  2.0*gtd[0][1]*gtd[0][2]*gtd[1][2]
                  -  gtd[0][2]*gtd[0][2]*gtd[1][1];

  if (det_gtd < 0.0) {
    det_gtd = 1.0;
  }
  
  double det_gtd_to_neg_third = 1.0 / pow(det_gtd, one_third);

  for (unsigned int j = 0; j < 3; j++)
  {
    for (unsigned int i = 0; i < 3; i++)
    {
      gtd[i][j] *= det_gtd_to_neg_third;
    }
  }

  det_gtd =   gtd[0][0]*( gtd[1][1]*gtd[2][2] - gtd[1][2]*gtd[1][2])
            - gtd[0][1]*gtd[0][1]*gtd[2][2]
            + 2.0*gtd[0][1]*gtd[0][2]*gtd[1][2]
            - gtd[0][2]*gtd[0][2]*gtd[1][1];

  double detgt_m1 = det_gtd - 1.0;

  // Error checking
  if (fabs(detgt_m1) > 1.0e-6) {
    printf("enforce_bssn_constraint: det(gtd) != 1. det=%lf\n",det_gtd);
    printf("gtd_11 = %lf\n",gtd[0][0]);
    printf("gtd_12 = %lf\n",gtd[0][1]);
    printf("gtd_13 = %lf\n",gtd[0][2]);
    printf("gtd_11 = %lf\n",gtd[1][1]);
    printf("gtd_12 = %lf\n",gtd[1][2]);
    printf("gtd_22 = %lf\n",gtd[2][2]);
  }

  double gtu[3][3];
  double idet_gtd = 1.0/det_gtd;
  gtu[0][0] = idet_gtd*(gtd[1][1]*gtd[2][2]-gtd[1][2]*gtd[1][2]);
  gtu[0][1] = idet_gtd*(-gtd[0][1]*gtd[2][2]+gtd[0][2]*gtd[1][2]);
  gtu[0][2] = idet_gtd*(gtd[0][1]*gtd[1][2]-gtd[0][2]*gtd[1][1]);
  gtu[1][0] = gtu[0][1];
  gtu[1][1] = idet_gtd*(gtd[0][0]*gtd[2][2]-gtd[0][2]*gtd[0][2]);
  gtu[1][2] = idet_gtd*(-gtd[0][0]*gtd[1][2]+gtd[0][1]*gtd[0][2]);
  gtu[2][0] = gtu[0][2];
  gtu[2][1] = gtu[1][2];
  gtu[2][2] = idet_gtd*(gtd[0][0]*gtd[1][1]-gtd[0][1]*gtd[0][1]);   

        /* Require Atd to be traceless. */
  double one_third_trace_Atd =   one_third * (
                        Atd[0][0]*gtu[0][0]
                      + Atd[1][1]*gtu[1][1]
                      + Atd[2][2]*gtu[2][2]
                      + 2.0 * (   Atd[0][1]*gtu[0][1]
                                + Atd[0][2]*gtu[0][2]
                                + Atd[1][2]*gtu[1][2]  )
                      );

  Atd[0][0] -= one_third_trace_Atd * gtd[0][0];
  Atd[0][1] -= one_third_trace_Atd * gtd[0][1];
  Atd[0][2] -= one_third_trace_Atd * gtd[0][2];
  Atd[1][1] -= one_third_trace_Atd * gtd[1][1];
  Atd[1][2] -= one_third_trace_Atd * gtd[1][2];
  Atd[2][2] -= one_third_trace_Atd * gtd[2][2];

  double tr_A =    Atd[0][0]*gtu[0][0]
                 + Atd[1][1]*gtu[1][1]
                 + Atd[2][2]*gtu[2][2]
                 + 2.0 * (   Atd[0][1]*gtu[0][1]
                           + Atd[0][2]*gtu[0][2]
                           + Atd[1][2]*gtu[1][2]  );
  
  /* apply a floor to chi and alpha */
  //TODO : Same, change vars

  double alpha, chi
  double floor_val = 1e-6;

  if ( chi <= 0.0) {
    chi = floor_val;
  }

  if ( alpha <= 0.0) {
    alpha = floor_val;
  }

}
