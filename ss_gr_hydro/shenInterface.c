# define DEFAULT_ELECTRON_FRAC (0.2)    //Electron fraction Ye
# define DENSITY_UNIT (1.0e13) // Units of rest density in gm/cm^3
# define VERBOSE (0) 
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include "loadTable.h"
# include "shenInterface.h"



static double *pressureTable;
static double *energyTable;
static double *soundSpeedTable;
static double *temperatureTable;


/*
Do bisection search.  Assumes data contains monotonically increasing values.
*/
int bisectionSearch(double* data, int length, double value){
	int lower_index = 0;
	int upper_index = length-1;
	int index;
	
	if(data[lower_index]>=value){
		return -1;	
	} else if(data[upper_index]<value){
		return length-1;
	} 
	while((upper_index-lower_index)>1){
		index = lower_index+(upper_index-lower_index)/2;
		if(data[index]>value){
			upper_index = index;
		} else {
			lower_index = index;
		}
	}
	return lower_index;
}
/*
 Do bilinear interpolation to (x,y) given 
 f00=f(x0,y0)
 f01=f(x0,y1)
 f10=f(x1,y0)
 f11=f(x1,y1)
 */
double bilinearInterpReg(double x, double y, double x0, double x1, double y0, double y1, double f00, double f01, double f10, double f11){
	double deltaX = (x1-x0);
	double deltaY = (y1-y0);
	double f = (f11+f00-f01-f10)*(x-x0)*(y-y0)/deltaX/deltaY+(f01-f00)*(x-x0)/deltaX+(f10-f00)*(y-y0)/deltaY+f00 ;
	return f;
}

/*
 Do bilinear interpolation to get f(x,y) given 
 f00=f(x0,y00)
 f01=f(x0,y01)
 f10=f(x1,y10)
 f11=f(x1,y11)
 */
double bilinearInterpHalfReg(double x, double y, double x0, double x1, double y00, double y01, double y10, double y11, double f00, double f01, double f10, double f11){
	double a = ((f11-f10)/(y11-y10)-(f01-f00)/(y01-y00))/(x1-x0);
	double b = (f11-f00-(y11-y00)*(f11-f10)/(y11-y10))/(x1-x0);
	double f = a*(x-x0)*(y-y00)+b*(x-x0)+(f01-f00)/(y01-y00)*(y-y00)+f00;

	return f;
}

/*
Load the tabulated EOS into memory.
We only load the tables for one Ye value.
*/
void loadShenTable(double Ye){
	static int init = 0;
	if(!init) {
		
		# include "shenTableParam.h"
		
		int Ye_index = ((Ye-Ye_min)*(num_Ye-1))/(Ye_max-Ye_min)+0.5;
		if(Ye_index<0){
			Ye_index=0; 
		} else if(Ye_index >= num_Ye){
			Ye_index=num_Ye-1;
		}
		
		loadTable( &energyTable, &pressureTable, &soundSpeedTable, Ye_index);
		
		temperatureTable = malloc(sizeof(double)*num_energy);
		int i;
		for(i=0; i<num_energy; i++){
			temperatureTable[i] = log_temp_min + ((double) i)/(num_energy-1)*(log_temp_max-log_temp_min);
		}

		init=1;
	}
}

/*
Use rest density and input to find output.
First we search in inputTable for the closest points to inputVal.  
Then we linearly interpolate in outputTable to get the output value.
Note that outputTableDim is the dimension of the output table.  
It should either be 1 if the output table is a function of the input variable index only,
Or 2 if the output table is a function of rest density index and input variable index.
*/
double eosLookup(double restDensity, double inputVal, double **pInputTable, double **pOutputTable, int outputTableDim){
	
	# include "shenTableParam.h"

	/*
	Load tables if loadShenTable() has not already been called.
	For non-default Ye value, loadShenTable() should have already been called.
	*/
	loadShenTable(DEFAULT_ELECTRON_FRAC);

	double *outputTable = *pOutputTable;
	double *inputTable = *pInputTable;


	int offTable = 0;
	double log_density = log10(restDensity)+log10(DENSITY_UNIT);
	double log_density_low, log_density_high;
	int index_rho = ((log_density-log_rho_min)/(log_rho_max-log_rho_min))*(num_rho-1);
	if(index_rho>num_rho-2){
		offTable = 1;
		index_rho = num_rho-2;
	} else if(index_rho<0){
		offTable = 1;
		index_rho = 0;
	}
	log_density_low = log_rho_min + ((double) index_rho)/(num_rho-1)*(log_rho_max-log_rho_min);
	log_density_high = log_rho_min + ((double) (index_rho+1))/(num_rho-1)*(log_rho_max-log_rho_min);

	int index_input_m = bisectionSearch((inputTable+num_energy*index_rho), num_energy, inputVal);
	int index_input_p = bisectionSearch((inputTable+num_energy*(index_rho+1)), num_energy, inputVal);

	if(index_input_m>num_energy-2){
		offTable = 1;
		index_input_m = num_energy-2;
	} else if(index_input_m<0){
		offTable = 1;
		index_input_m = 0;
	}
	if(index_input_p>num_energy-2){
		offTable = 1;
		index_input_p = num_energy-2;
	} else if(index_input_p<0){
		offTable = 1;
		index_input_p = 0;
	}

	//interpolate
	double p00, p01, p10, p11;
	if(outputTableDim==1){
		p00 = outputTable[index_input_m];
		p01 = outputTable[index_input_m+1];
		p10 = outputTable[index_input_p];
		p11 = outputTable[index_input_p+1];
	} else{
		if(outputTableDim!=2) printf("Error in output table dim\n");
		p00 = outputTable[index_rho*num_energy+index_input_m];
		p01 = outputTable[index_rho*num_energy+index_input_m+1];
		p10 = outputTable[(index_rho+1)*num_energy+index_input_p];
		p11 = outputTable[(index_rho+1)*num_energy+index_input_p+1];
	} 

	double e00 = inputTable[index_rho*num_energy+index_input_m];
	double e01 = inputTable[index_rho*num_energy+index_input_m+1];
	double e10 = inputTable[(index_rho+1)*num_energy+index_input_p];
	double e11 = inputTable[(index_rho+1)*num_energy+index_input_p+1];
	

	double output = bilinearInterpHalfReg(log_density, inputVal, log_density_low, log_density_high, e00, e01, e10, e11, p00, p01, p10, p11);
		
	//if(offTable) printf("Off table at rho=%E input=%E output=%e\n",restDensity, inputVal, output);

	return output;
}


/*
Get the pressure according to the Shen EOS as a function of rest density and specific energy. 
*/
double getShenPressure(double restDensity, double specificEnergy){
		
	double log_energy = log10(specificEnergy);
	double pressure = eosLookup(restDensity, log_energy, &energyTable, &pressureTable, 2); 
	pressure = pressure/DENSITY_UNIT;
	
	if(pressure<0.0){
		pressure = 0.0;
	#if VERBOSE >= 2
		printf("Pressure floored!\n");
	#endif
	}

	return pressure;
}

/*
Get the sound speed according to the Shen EOS as a function of rest density and pressure 
*/
double getShenSoundSpeed(double restDensity, double specificEnergy){
	
	double log_energy = log10(specificEnergy);
	double soundSpeed = eosLookup(restDensity, log_energy, &energyTable, &soundSpeedTable, 2); 
	
	if(soundSpeed<0.0){
		soundSpeed = 0.0;
	#if VERBOSE >=2	
		printf("Sound speed floored!\n");
	#endif
	} else if (soundSpeed>1.0){
		soundSpeed = 1.0;
	#if VERBOSE >=2
		printf("Sound speed capped!\n");
	#endif
	}

	return soundSpeed;
}

/*
Get the specific energy according to the Shen EOS as a function of rest density and pressure 
*/
double getShenSpecificEnergy(double restDensity, double pressure){
	pressure=pressure*DENSITY_UNIT;
	double log_energy = eosLookup(restDensity, pressure, &pressureTable, &energyTable, 2); 
	double specificEnergy = exp(log(10.0)*log_energy);

	return specificEnergy;
}

/*
Get the temperature (in MeV) according to the Shen EOS as a function of rest density and specific energy
*/
double getShenTemperature(double restDensity, double specificEnergy){
	
	double log_energy = log10(specificEnergy);
	double log_temp = eosLookup(restDensity, log_energy, &energyTable, &temperatureTable, 1); 
	double temp = exp(log(10.0)*log_temp);
	
	return temp;
}



