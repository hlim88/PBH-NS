
enum eosTypes {
	GAMMA_LAW,
    ADIABATIC,
	SHEN_EOS
};

/*
Call first to set the equation of state that will be used.
*/
void setEquationOfState(int eos_type, double param);

//Use equation of state to return pressure
double getPressure(double restDensity, double specificEnergy);

//Return sound speed as function of specific energy and rest density
double getSoundSpeed(double restDensity, double specificEnergy);

/*
This equation is used in addition to the EOS for example in finding initial TOV solutions
*/
double getPolytropicRestDensity(double pressure);

//Use equation of state to return specific energy
double getSpecificEnergy(double restDensity, double pressure);
