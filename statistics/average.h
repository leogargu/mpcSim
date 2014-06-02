#ifndef AVERAGE_H
#define AVERAGE_H


#include <stdlib.h>
#include <stdio.h>
//#include <math.h>
//#include <time.h>
#include <string.h>




/* To do: */
/* NOTE: CAM_average and SAM_average contain very similar code: it should be possible to merge in one function, with the option of calculating CAM or SAM */
/* These two functons do not work for vector values */



////////////////////////////////////////////////////////////////////////////////////////
/// This function returns the Cumulative Average of a SCALAR over a number of files of the form: filename_num.dat, with num from first_file to last_file.
/// The first row of the output datafile contains information about the geometry and the origin of the sample. It has the form:  
/// nx ny nz (global index of first cell) (global index of last cell) timestep
/// Subsequent rows correspond each to one particular cell in the simulation. This can be a particular x_slice in the simulation box, a set of cells or the whole simulation box.
/// This is determined by the values of the cell_idx_start and cell_idx_end
/// The output data is stored in the array average.
/// The factor depends on the particular quantity being CAM-averaged: 
/// e.g. 1/3 for temperature, 1.0 for velocity, 1/3a^3 for pressure (possibly, double check this) 
/// WARNING: This function does not check that all the files correspond to the same simulation, or whether they are "averageable" or not
/// verbose: 0 for running silently, 1 for printing on screen the file being processed
////////////////////////////////////////////////////////////////////////////////////////
inline void CAM_average(char * filename, int first_file, int last_file, double factor, int verbose )
{
	int flag;
	/*Read header of first file to obtain geometry data */
	FILE * fp;
	char intaschar[20];
	char filename1[50]="";
	strcpy(filename1,filename);
	sprintf(intaschar, "_%d.dat", first_file);
	strcat(filename1,intaschar);	
	
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

	
	for(i=first_file; i<=last_file; i++) //file 1 is open the first time that we pass through here	
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
		flag = fclose(fp);
		if( flag !=0 )
		{
			fprintf(stderr,"export_vtk_plasma: error closing file\n Aborting...\n");
			exit(EXIT_FAILURE);
		}
		/* ... prepare the name of the next file...*/
		if(i==last_file)
		{
			break; 
		}
		strcpy(filename1,filename);		
		sprintf(intaschar, "_%d.dat", i+1);
		strcat(filename1,intaschar);
		
		if( verbose ){	printf("Processing file %s\n",filename1); }
		
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
		
	/* Clean up */
	free(average);
	free(part_num);
	flag = fclose(fp);
	if( flag !=0 )
	{
		fprintf(stderr,"CAM_average: error closing file\n Aborting...\n");
		exit(EXIT_FAILURE);
	}
	
	return; /* back to main */
}




////////////////////////////////////////////////////////////////////////////////////////
/// This function returns the Sample Average of a SCALAR over a number of files of the form: filename_num.dat, with num from fist_file to last_file.
///This function takes a different type of files compared to CAM_average. The first line is a header (same as CAM_average), but the line afterwards contains the values of the quantity to be averaged at every cell, in succession. I.e.
///   header
///   value_for_cell_start \t value_for_cell_start+1 \t ...
/// This different data file can be generated with export_SAM_data_scalar(), or it can be created by using the CAM data file formats. 
/// An auxiliary function would be necessary in that case (implemented, not tested). 
/// The factor depends on the particular quantity being SAM-averaged: 
/// e.g. 1/3 for temperature, 1.0 for velocity, 1/3a^3 for pressure (possibly, double check this) 
/// WARNING: This function does not check that all the files correspond to the same simulation, or whether they are "averageable" or not
/// verbose: 0 for running silently, 1 for printing on screen the file being processed
////////////////////////////////////////////////////////////////////////////////////////
inline void SAM_average(char * filename, int first_file, int last_file, double factor, int verbose )
{
	int flag;
	/*Read header of first file to obtain size of slice */
	FILE * fp;	
	char intaschar[20];
	char filename1[50]="";
	strcpy(filename1,filename);
	sprintf(intaschar, "_%d.dat", first_file);
	strcat(filename1,intaschar);
	
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
	double aux=0.0;
	int read_flag=0;
		
	if(verbose){ printf("Processing file %s\n",filename1); }
	for(i=first_file; i<=last_file; i++) //file 1 is open the first time that we pass through here	
	{
		for(j=0;j<num_cells;j++)
		{
			read_flag=fscanf(fp, "%lf", &aux);
			if(read_flag!=1)
			{
				fprintf(stderr,"SAM_average: Error reading number from file. Aborting...\n");
				exit(EXIT_FAILURE);
			}
				
			average[j] += aux;
			
			
		}
		/* close this file...*/
		flag = fclose(fp);
		if( flag !=0 )
		{
			fprintf(stderr,"export_vtk_plasma: error closing file\n Aborting...\n");
			exit(EXIT_FAILURE);
		}
	
		/* ... prepare the name of the next file...*/
		if(i==last_file)
		{
			break; 
		}
		strcpy(filename1,filename);		
		sprintf(intaschar, "_%d.dat", i+1);
		strcat(filename1,intaschar);
		
		if(verbose){ printf("Processing file %s\n",filename1); }
		
		/*...and open the next file */		
		fp=fopen(filename1,"r"); 
		if( fp == NULL ){ perror("SAM_average: Error while opening the file.\n"); exit(EXIT_FAILURE);}

		/* read fist line to set the pointer in the right position. Checking of parameters for robustness is possible here */
		read_flag=fscanf( fp, "%d \t %d \t %d \t %d \t %d \t %ld \n", &nx, &ny, &nz, &cell_idx_start, &cell_idx_end, &step);
		/* code repetition!*/
		if(read_flag!=6)
		{
			fprintf(stderr,"SAM_average: Error reading numbers from file. Aborting...\n");
			exit(EXIT_FAILURE);
		}
		
	}
	
	
	/* Complete the average and export to file */
	strcpy(filename1,filename);
	strcat(filename1,"_SAM_averaged.dat");
	fopen(filename1,"w");
	
	fprintf( fp, "%d \t %d \t %d \t %d \t %d \n", nx, ny, nz, cell_idx_start, cell_idx_end);
		
	int samples = last_file - first_file + 1;
	
	double fact = factor * (1.0/(double)samples);
	
	for(i=0; i<num_cells; i++)
	{
		fprintf(fp, "%lf\n", fact*average[i] );
	}
		
	/* Clean up */
	free(average);
	flag = fclose(fp);
	if( flag !=0 )
	{
		fprintf(stderr,"export_vtk_plasma: error closing file\n Aborting...\n");
		exit(EXIT_FAILURE);
	}
		
	return; /* back to main */
}




/////////////////////////////////////////////////////////////////////////////
/// Transforms a CAM file filename into a SAM file
/////////////////////////////////////////////////////////////////////////////
/* not tested */
inline void CAM_to_SAM(char * filename)
{	
	int flag, i,j;

	FILE * fp;
	FILE * fp_out;	
	fp = fopen(filename,"r");
	if( fp==NULL ){printf("CAM_to_SAM: Error opening file. Aborting...\n");exit(EXIT_FAILURE);}
	
	
	char filename1[50]="";
	strcpy(filename1,filename);
	strcat(filename1,"_SAM");
	

	fp_out = fopen(filename1,"w");
	if( fp_out == NULL ){printf("CAM_to_SAM: error opening file. Aborting...\n"); exit(EXIT_FAILURE);}
	
	/* Read header from input file and parse to output file */
	int nx=0, ny=0, nz=0;
	long int step=0;
	int cell_idx_start=0 ,cell_idx_end= 0;
	
	fscanf( fp, "%d \t %d \t %d \t %d \t %d \t %ld \n", &nx, &ny, &nz, &cell_idx_start, &cell_idx_end, &step);
	fprintf( fp_out, "%d \t %d \t %d \t %d \t %d \t %ld \n", nx, ny, nz, cell_idx_start, cell_idx_end, step);
	
	
	/* Identify volume to study (set of collision cells) */
	int num_cells = cell_idx_end - cell_idx_start + 1;
	
	/* Parsing from CAM to SAM */
	double average=0;
	double aux;
	double row_length;
	
	for(i=0; i<num_cells; i++)
	{
		/* get local density */
		fscanf(fp, "%lf", &row_length);
		fprintf(fp_out,"%lf", row_length);
		for(j=1; j<=row_length; j++)
		{
			fscanf(fp,"%lf",&aux);
			average += aux;
		}
		if( row_length !=0 )
		{
			fprintf(fp_out,"\t %lf",average/((double)row_length) );
		}
		average = 0.0;
		
	}
	
	/* Closing files */
	flag = fclose(fp);
	if(flag!=0)
	{
		printf("CAM_to_SAM: Error closing file. Aborting...\n");
		exit(EXIT_FAILURE);
	}
	flag = fclose(fp_out);
	if(flag!=0)
	{
		printf("CAM_to_SAM: Error closing file. Aborting...\n");
		exit(EXIT_FAILURE);
	}
	
	
	return;/* back to main */
}






#endif /* header guard*/
