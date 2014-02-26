#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>



/* This function return the Cumulative Average of a number samples of files of the form: filename_num.dat, with 
num from 1 to samples */
/* The first row contains information about the geometry and the origin of the sample. It has the form: 
   nx ny nz (global index of first cell) (global index of last cell) timestep */
/* Subsequent rows correspond each to one particular cell in the simulation. This can be a particular x_slice in the simulation box,
   a set of cells or the whole simulation box. */
/* The array average is the output */
/* The factor depends on the particular quantity being CAM-averaged: e.g. 1/3 for temperature, 1.0 for velocity, 1/3a^3 for pressure (possibly, double check this) */
/* WARNING: This function does not check that all the files correspond to the same simulation, or whether they are "averageable" or not */
/* This function HAS BEEN TESTED, compared to velocity_average.c everything OK */
inline void CAM_average(char * filename, int samples, double factor )
{
	/*Read header of first file to obtain geometry data */
	FILE * fp;
	char filename1[50]="";
	strcpy(filename1,filename);
	strcat(filename1,"_1.dat");	
	
	fp=fopen(filename1,"r"); 
	if( fp == NULL ){ perror("CAM_average: Error while opening the file.\n"); exit(EXIT_FAILURE);}
 
	int nx=0, ny=0, nz=0;
	long int step=0;
	int cell_idx_start=0 ,cell_idx_end= 0;
	
	fscanf( fp, "%d \t %d \t %d \t %d \t %d \t %ld \n", &nx, &ny, &nz, &cell_idx_start, &cell_idx_end, &step);
	
	/* Identify volume to study (set of collision cells) */
	int num_cells = cell_idx_end - cell_idx_start + 1;
	
	
	/* Allocate memory */
	double * average;
	int * part_num;
	
	average = malloc( num_cells * sizeof(double) );
	if (average==NULL) {printf("Error allocating average in CAM_average call \n"); exit(EXIT_FAILURE);}
	part_num = malloc( num_cells * sizeof(int) );
	if (part_num==NULL) {printf("Error allocating part_num in slice_average_CAM call \n"); exit(EXIT_FAILURE);}
	
	/* Initialization */
	int i,j,k;
	
	for(i=0; i<num_cells; i++ )
	{
		average[i] = 0.0;
		part_num[i] = 0;
	}
	
	/* Read data from all files and start the average */
	double row_length; // this is the occupation number for a given cell. It is stored as a double because the data file will be a 2D array of doubles, only the first element will be "special"
	double aux;
	char intaschar[20];

	
	for(i=1; i<=samples; i++) //file 1 is open the first time that we pass through here	
	{
		for(j=0;j<num_cells;j++)
		{
			fscanf(fp, "%lf", &row_length);
			
			part_num[j] += (int)row_length;
		
			for(k=1; k<=(int)row_length; k++)
			{
				fscanf(fp, "%lf", &aux);
				average[j]+=aux;
			}
			
		}
		/* close this file...*/
		fclose(fp);
		/* ... prepare the name of the next file...*/
		if(i==samples)
		{
			break; 
		}
		strcpy(filename1,filename);		
		sprintf(intaschar, "_%d.dat", i+1);
		strcat(filename1,intaschar);
		
		/*...and open the next file */		
		fp=fopen(filename1,"r"); 
		if( fp == NULL ){ perror("CAM_average: Error while opening the file.\n"); exit(EXIT_FAILURE);}

		/* read fist line to set the pointer in the right position. Checking of parameters for robustness is possible here */
		fscanf( fp, "%d \t %d \t %d \t %d \t %d \t %ld \n", &nx, &ny, &nz, &cell_idx_start, &cell_idx_end, &step);
		
	}
	
	
	/* Complete the average and export to file */
	strcpy(filename1,filename);
	strcat(filename1,"_CAM_averaged.dat");
	
	fopen(filename1,"w");
	
	fprintf( fp, "%d \t %d \t %d \t %d \t %d \n", nx, ny, nz, cell_idx_start, cell_idx_end);
	for(i=0;i<num_cells;i++)
	{
		if( part_num[i]!=0)
		{
			average[i] = factor * (average[i]/(double)part_num[i]) ;
		}else{
			average[i] = 0.0; // no particles in this collision cell, in any of the samples
		}
		fprintf(fp,"%lf\n", average[i]);
	}
	
	//fprintf(fp,"\n");
	
	/* Clean up */
	free(average);
	free(part_num);
	fclose(fp);
	
	return; /* back to main */
}




/* This function takes a different type of files compared to CAM_average. The first line is a header (same as CAM_average), but 
   the line afterwards contains the values of the quantity to be averaged at every cell, in succession. I.e.
   header
   value_for_cell_start \t value_for_cell_start+1 \t ...
*/
/* This different data file can be generated directly by the main code, or it can be created by using the other data file formats. An
   auxiliary function would be necessary in that case. 	
*/
/* This function HAS BEEN TESTED against temperature_average. It works OK. */
inline void SAM_average(char * filename, int samples, double factor )
{
	/*Read header of first file to obtain size of slice */
	FILE * fp;
	char filename1[50]="";
	strcpy(filename1,filename);
	strcat(filename1,"_1.dat");	
	
	fp=fopen(filename1,"r"); 
	if( fp == NULL ){ perror("SAM_average: Error while opening the file.\n"); exit(EXIT_FAILURE);}
 
	int nx=0, ny=0, nz=0;
	long int step=0;
	int cell_idx_start=0 ,cell_idx_end= 0;
	
	fscanf( fp, "%d \t %d \t %d \t %d \t %d \t %ld \n", &nx, &ny, &nz, &cell_idx_start, &cell_idx_end, &step);
	
	/* Identify volume to study (set of collision cells) */
	int num_cells = cell_idx_end - cell_idx_start + 1;
	
	/* Allocate memory */
	double * average;
	
	average = malloc( num_cells * sizeof(double) );
	if (average==NULL) {printf("Error allocating average in average_SAM call \n"); exit(EXIT_FAILURE);}
	
	/* Initialization */
	int i,j;
	
	for(i=0; i<num_cells; i++ )
	{
		average[i] = 0.0;
	}
	
	/* Read data from all files and start the average */
	double aux;
	char intaschar[20];
		
	for(i=1; i<=samples; i++) //file 1 is open the first time that we pass through here	
	{
		for(j=0;j<num_cells;j++)
		{
			fscanf(fp, "%lf", &aux);
			if(aux >=0.0 )
			{
				average[j] += aux;
			}
			
		}
		/* close this file...*/
		fclose(fp);
		/* ... prepare the name of the next file...*/
		if(i==samples)
		{
			break; 
		}
		strcpy(filename1,filename);		
		sprintf(intaschar, "_%d.dat", i+1);
		strcat(filename1,intaschar);
		
		/*...and open the next file */		
		fp=fopen(filename1,"r"); 
		if( fp == NULL ){ perror("SAM_average: Error while opening the file.\n"); exit(EXIT_FAILURE);}

		/* read fist line to set the pointer in the right position. Checking of parameters for robustness is possible here */
		fscanf( fp, "%d \t %d \t %d \t %d \t %d \t %ld \n", &nx, &ny, &nz, &cell_idx_start, &cell_idx_end, &step);
		
	}
	
	
	/* Complete the average and export to file */
	strcpy(filename1,filename);
	strcat(filename1,"_SAM_averaged.dat");
	fopen(filename1,"w");
	
	fprintf( fp, "%d \t %d \t %d \t %d \t %d \n", nx, ny, nz, cell_idx_start, cell_idx_end);
		
	double fact = factor * (1.0/(double)samples);
	
	for(i=0; i<num_cells; i++)
	{
		fprintf(fp, "%lf\t", fact*average[i] );
	}
		
	/* Clean up */
	free(average);
	fclose(fp);
	
	return; /* back to main */
}




/*-------------------------------*/
/*		MAIN		 */
/*-------------------------------*/

/* Applicable to CAM averages only, for now */
/* call with arguments filename numberOFsamples number, with number being 1 for velocities averages, 3.0 for temperature averages and 3*a^3 for pressure averages */
int main(int argc, char **argv) {

	/*
	argv[0] is the program name
	argv[1] is the type of average: CAM or SAM
	argv[2] is name of files
	argv[3] is the number of samples 
	argv[4] is a factor (see note above)
	*/
	
	if(argc!=5)
	{
		printf("Wrong number of arguments. Call as:\n");
		printf("./average <av.type> <filename> <n.samples> <factor>\n<av.type>\tSAM or CAM\n<filename>\tdatafile basic name\n<n.samples>\tnumber of datafiles\n<factor>\tcorrection factor\n");
		printf("Aborting...\n");
		exit(EXIT_FAILURE);
	}
	
	
	char filename[50]="./DATA/";
	strcat(filename,argv[2]);
	
	double number;
	number = strtod(argv[4], NULL);
	
	
	if(strcmp(argv[1], "CAM")==0)
	{
		//printf("you chose CAM\n");
		CAM_average(filename, atoi(argv[3]), 1.0/number );
	}else if(strcmp(argv[1], "SAM")==0)
	{
		//printf("you chose SAM\n");
		SAM_average(filename, atoi(argv[3]), 1.0/number);
	}else{
		printf("Type of average not recognized. Aborting...\n");
		exit(EXIT_FAILURE);
	}

		
	return 0;
}
