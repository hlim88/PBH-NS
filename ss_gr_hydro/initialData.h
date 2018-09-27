//Initial Data Parameters
#define CENTRAL_DENSITY (22.77e-4) 
#define ATM_PRESSURE (1e-14) //Pressure of the 'atmosphere'
#define VELOCITY_AMP (0.4) //Parameter determining magnitude of velcity profile
#define X_ATM (1.2) //Parameter determining the radius (in units of the star radius) above which the velocity is zero

/*
Set inital data.  All the input arrays are assumed to be empty with the exception of r.
velocityAmp is parameter determining magnitude of velcity profile
*/
double getInitialData(double centralPressure, double velocityAmp, double *consVar, double *primVar, double *a, double *alpha, double *rGravity, double *rFluid, int length);

double getID2(double centralPressure, double velocityAmp,
              double *consVar, double *primVar, double *a,
              double *alpha, double *rGravity, double *rFluid,
              int length){

/*
 For partial domains, this function extends the domain from r=0 to r=rMaxExt
 Then it uses getInitialData() to find ID over this extended domain and then
 copies in the required part.
 */
double getInitialDataPartialDomain(double centralPressure, double velocityAmp, double rMaxExt, double *consVar, double *primVar, double *a, double *alpha, double *rGravity, double *rFluid, int length);
