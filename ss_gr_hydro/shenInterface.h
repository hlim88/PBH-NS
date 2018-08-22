/*
Load the tabulated EOS into memory.
We only load the tables for one Ye value.
*/
void loadShenTable(double Ye);

/*
Get the pressure according to the Shen EOS as a function of rest density and specific energy. 
*/
double getShenPressure(double restDensity, double specificEnergy);

/*
Get the sound speed according to the Shen EOS as a function of 
*/
double getShenSoundSpeed(double restDensity, double specificEnergy);

/*
Get the specific energy according to the Shen EOS as a function of rest density and pressure 
*/
double getShenSpecificEnergy(double restDensity, double pressure);

/*
Get the temperature (in MeV) according to the Shen EOS as a function of rest density and specific energy
*/
double getShenTemperature(double restDensity, double specificEnergy);

