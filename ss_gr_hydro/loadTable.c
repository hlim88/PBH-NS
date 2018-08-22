#include <stdio.h>
#include <stdlib.h> 
#include <string.h>
#include <math.h>
#include "loadTable.h"

//Read in a binary (of byte size entrySize) array to data from a binary file fileName located in filePath
void* readBinaryArray(char *filePath, char*fileName, int entrySize, int *dataLength){
        void *data;
        char *inputFileName = malloc((strlen(filePath)+strlen(fileName)+strlen("/")+1)*sizeof(char));
        strcpy(inputFileName, filePath);
        strcat(inputFileName,"/");
        strcat(inputFileName, fileName);
        FILE *file = fopen(inputFileName, "rb");
        if(file!=NULL){
                fseek(file, 0, SEEK_END);
                int size = ftell(file);
                fseek(file, 0, SEEK_SET);
                if(size%entrySize!=0){
                        printf("Error %s is not an array of %d byte entries!\n",inputFileName, entrySize);
                        return 0;
                }
                int length = size/entrySize;
                data = malloc(size);
                int n = fread(data, entrySize, length, file);
                if(n>0){
                        printf("%d samples read in from %s\n", n, inputFileName);
                } else {
                        printf("Error reading from %s\n", inputFileName);
                }
                fclose(file);
                *dataLength = length;
                return data;
        } else {
                printf("Error could not open %s\n", inputFileName);
                *dataLength = 0;
        }
	return NULL;
}

//Read in a double array to data from a binary file fileName located in filePath
int readBinaryDoubleArray(char *filePath, char*fileName, double *(data[])){
        int dataLength;
        *data = readBinaryArray(filePath, fileName, sizeof(double), &dataLength);
        return dataLength;
}

/*
Load Shen EOS table from data file.
Only load table for one specified Ye value.
*/
void loadTable(double **energyTable, double **pressureTable, double **soundSpeedTable, int Ye_index){
	#include "shenTableParam.h"

	double *table;

	int numEntries = readBinaryDoubleArray(filepath, filename, &table);	
	int sizeVar = num_rho*num_energy*num_Ye;

	if(numEntries!=sizeVar*num_elements) printf("Error table parameters %d and number of file entries %d do not match!\n",sizeVar,numEntries);

	int invalid_table_flag = 0; //Set to zero if there are nans in the table
	int size = num_rho*num_energy;
	*energyTable = malloc(sizeof(double)*size);
	*pressureTable = malloc(sizeof(double)*size);
	*soundSpeedTable = malloc(sizeof(double)*size);

	double *gmx = malloc(sizeof(double)*size);
	
	int var_index = 0;
	int i;
	//Load specific energy
	for(i=0; i<size; i++) {
		//log of specific energy in units of c^2
		(*energyTable)[i]=log10((table[i+Ye_index*size+var_index*sizeVar]+9.3)*(1.06431302e-3));
		if(isnan((*energyTable)[i])) {
			invalid_table_flag = 1;
		}
	}

	if(invalid_table_flag) printf("Error in energy table\n");
	
	//Load pressure
	var_index = 1;
	for(i=0; i<size; i++) {
		//pressure in units of c^2*gm/cm^3
		(*pressureTable)[i]=(table[i+Ye_index*size+var_index*sizeVar])*(1.7826627e12);
		if(isnan((*pressureTable)[i])) {
			invalid_table_flag = 1;
		} else if((*pressureTable)[i]<0.0) {
			invalid_table_flag = 1;
		}
	}
	
	if(invalid_table_flag) printf("Error in pressure table\n");

	//Load gmx
	int k;
	var_index = 11;
	for(i=0; i<size; i++) {
		gmx[i]=table[i+Ye_index*size+var_index*sizeVar];
		//Interpolate over nan's
		if(isnan(gmx[i])){
			int j = 1;
			//Find closest non-nan
			while((i+j)<size && isnan(table[i+j+Ye_index*size+var_index*sizeVar])){
				j++;
			}
			if((i%num_energy)>((i+j)%num_energy) || j>=num_energy){
				if((i%num_energy)+j==num_energy){
					for(k=0; k<j; k++){
						gmx[i+k] = gmx[i-2] + (gmx[i-1]-gmx[i-2])*((double) (k+1));
					}
					i = i+j-1; 
				} else {
					printf("Error in interpolating over gmx nan's\n");
				}
			} else if((i+j)<size && i!=0){
				for(k=0; k<j; k++){
					gmx[i+k] = gmx[i-1] + (table[i+j+Ye_index*size+var_index*sizeVar]-gmx[i-1])*((double) (k+1))/(j+1);
				}
				i = i+j-1;
			} else {
				printf("Error, nan's in gmx at edge of table\n");
			}
		}
	}

	//Ensure gamma is between 1 and 3
	for(i=0; i<size; i++){
		if(gmx[i]<1.0){ 
			gmx[i]=1.0;
		} else if(gmx[i]>3.0){
			gmx[i]=3.0;
		}
	}	

	//Load speed of sound
	for(i=0; i<size; i++) {
		int index_rho = (i%(num_energy*num_rho))/num_energy;
		//if(isnan(gmx[i]) || gmx[i]<=0.0) printf("Error in gmx=%e at i=%d rho_index=%d\n",gmx[i],i,index_rho);
		double log_density = log_rho_min + ((double) index_rho)/(num_rho-1)*(log_rho_max-log_rho_min);
		//Speed of sound in units of c
		(*soundSpeedTable)[i] = sqrt((*pressureTable)[i]/(exp(log(10.0)*log_density)*(1.0+exp(log(10.0)*(*energyTable)[i]))+(*pressureTable)[i])*fabs(gmx[i]));

		if(isnan((*soundSpeedTable)[i])) {
			invalid_table_flag = 1;
		}

		if((*soundSpeedTable)[i]>1.0) {
			(*soundSpeedTable)[i]=1.0;
		} else if((*soundSpeedTable)[i]<0.0){
			(*soundSpeedTable)[i]=0.0;	
		}  
	}

	if(invalid_table_flag) {	
		printf("Warning, EOS table is invalid!\n");
		abort();
	}

	printf("done\n");
	free(gmx);
	free(table);
}
