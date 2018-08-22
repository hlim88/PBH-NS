# define DEFAULT_GAMMA (2.0)

#include "shenInterface.h"

# include "equationOfState.h"
# include <math.h>
# include <stdio.h>

static int eosType = 0;
static double GAMMA_INDEX = DEFAULT_GAMMA;

/*
Call first to set the equation of state that will be used.
*/
void setEquationOfState(int eos_type, double param){
	eosType = eos_type;
	if(eosType==GAMMA_LAW){
		GAMMA_INDEX = param;	
	} else if(eosType==SHEN_EOS){
		//Load Shen EOS with Ye=param
		loadShenTable(param);
	}
}


//Use equation of state to return pressure
double getPressure(double restDensity, double specificEnergy){
	double pressure;
	
	if(eosType == GAMMA_LAW){
		pressure = (GAMMA_INDEX-1.0)*restDensity*specificEnergy;
	} else if(eosType == SHEN_EOS){
		pressure = getShenPressure(restDensity, specificEnergy);
	}

	return pressure; 
}

//Use equation of state to return specific energy
double getSpecificEnergy(double restDensity, double pressure){
	double specificEnergy;
	
	if(eosType == GAMMA_LAW){
		specificEnergy = pressure/restDensity/(GAMMA_INDEX-1.0);
	} else if(eosType == SHEN_EOS){
		specificEnergy = getShenSpecificEnergy(restDensity, pressure);
	}
	
	return specificEnergy;
}

//Return sound speed as function of specific energy and rest density
double getSoundSpeed(double restDensity, double specificEnergy){
	double soundSpeed;
	
	if(eosType == GAMMA_LAW){
		soundSpeed = sqrt(GAMMA_INDEX*(GAMMA_INDEX-1.0)*specificEnergy/(1.0+GAMMA_INDEX*specificEnergy));
	} else if(eosType == SHEN_EOS){ 
		soundSpeed = getShenSoundSpeed(restDensity, specificEnergy);
	}	

	return soundSpeed;
}


//Return temperature in MeV, for use with tabulated Shen EOS
double getTemperature(double restDensity, double specificEnergy){
	double temperature;
	if(eosType == GAMMA_LAW){
		temperature = 0.0;
	} else if(eosType == SHEN_EOS){
		temperature = getShenTemperature(restDensity, specificEnergy);	
	}
	return temperature;
}

/*
This equation is used in addition to the EOS for example in finding initial TOV solutions
*/
double getPolytropicRestDensity(double pressure){

		double restDensity = pow(pressure, (1.0/GAMMA_INDEX));
		return restDensity;
}

