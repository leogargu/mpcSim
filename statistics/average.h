#ifndef AVERAGE_H
#define AVERAGE_H


#include <stdlib.h>
#include <stdio.h>
#include <math.h> //isnan lives here
//#include <time.h>
#include <string.h>
#include <assert.h>




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
/// The output file has two columns: the first column contains the number of particles, over all the samples used, used in the average.
/// The second column has the actual value of the average, and it is empty is the corresponding cell was empty across the samples used to calculate the average.
////////////////////////////////////////////////////////////////////////////////////////
inline void CAM_average(char * filename, int first_file, int last_file, int stride, double factor, int verbose )
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
	int actual_last_file;
	
	//this loop should be a do-while instead
	if(verbose){ printf("Processing file %s\n",filename1); }
	for(i=first_file; i<=last_file; i=i+stride) //file 1 is open the first time that we pass through here	
	{
		actual_last_file = i;
		
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
		if(i+stride>last_file)
		{
			break; 
		}
		strcpy(filename1,filename);		
		sprintf(intaschar, "_%d.dat", i+stride);
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
	
	/* Export expanded header */
	fprintf( fp, "%d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \n", nx, ny, nz, cell_idx_start, cell_idx_end, first_file, stride, actual_last_file);
	for(i=0;i<num_cells;i++)
	{
		fprintf(fp, "%d \t ",part_num[i]);
		if( part_num[i]!=0 )
		{
			average[i] = factor * (average[i]/(double)part_num[i]) ;
			fprintf(fp,"%lf\n", average[i]);
		}else{
			fprintf(fp,"\n");
		}
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
/// This function returns the Sample Average of a **SCALAR** over a number of files of the form: filename_num.dat, with num from fist_file to last_file.
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
inline void SAM_average(char * filename, int first_file, int last_file, int stride, double factor, int verbose )
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
	if (average==NULL) {printf("SAM_average: Error allocating average array. \n"); exit(EXIT_FAILURE);}
	
	
	/* Need to discard empty cells in their average */
	int * samples;
	samples=malloc(num_cells*sizeof(int));
	if (samples==NULL) {printf("SAM_average: error allocating samples array. \n"); exit(EXIT_FAILURE);}
	
	//int samples_baseline = last_file - first_file +	 1;//maximum number of samples that will be used in any cell
	int samples_baseline = 1+(int)((last_file - first_file)/stride );
	
	/* Initialization */
	int i,j;
	
	for(i=0; i<num_cells; i++ )
	{
		average[i] = 0.0;
		samples[i] = samples_baseline;
	}
		
	/* Read data from all files and start the average */
	double aux=0.0;
	int read_flag=0;
	int actual_last_file;
		
	if(verbose){ printf("Processing file %s\n",filename1); }
	for(i=first_file; i<=last_file; i=i+stride) //file 1 is open the first time that we pass through here	
	{
		actual_last_file = i;
		
		for(j=0;j<num_cells;j++)
		{
			read_flag=fscanf(fp, "%lf", &aux);
			if(read_flag!=1){ fprintf(stderr,"SAM_average: Error reading number from file. Aborting...\n"); exit(EXIT_FAILURE); }
			
			if(aux > 0.0)//then we have some particle(s) in that cell, so we accept this sample
			{
				read_flag=fscanf(fp, "%lf", &aux);
				if(read_flag!=1){ fprintf(stderr,"SAM_average: Error reading number from file. Aborting...\n"); exit(EXIT_FAILURE); }
				
				average[j]+=aux;
				
			}else{//this cell (j) is empty in this sample. We discard this data *for this cell only*
				samples[j] -= 1 ;
			}		
			
		}
		/* close this file...*/
		flag = fclose(fp);
		if( flag !=0 ){ fprintf(stderr,"export_vtk_plasma: error closing file\n Aborting...\n"); exit(EXIT_FAILURE); }
	
		/* ... prepare the name of the next file...*/
		if(i+stride>last_file)
		{
			break; 
		}
		strcpy(filename1,filename);		
		sprintf(intaschar, "_%d.dat", i+stride);
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
	
	/* Export extended header */
	fprintf( fp, "%d \t %d \t %d \t %d \t %d \t %d \t %d \t %d\n", nx, ny, nz, cell_idx_start, cell_idx_end, first_file, stride, actual_last_file);
		

	for(i=0; i<num_cells; i++)
	{
		assert( samples[i]>=0 );
		if( samples[i]!=0 ){
			fprintf(fp, "%d \t %lf\n", samples[i], average[i]/(double)samples[i] );
		}else{ //cells that are empty (e.g. external cells) should ideally get an empty space or NAN. genfromtxt needs the second column to have some value, so 0 would be needed here if not using pandas for visualization
			fprintf(fp, "0 \t \n" );
		}
	}
		
	/* Clean up */
	free(average);
	free(samples);
	flag = fclose(fp);
	if( flag !=0 ){ fprintf(stderr,"export_vtk_plasma: error closing file\n Aborting...\n"); exit(EXIT_FAILURE); }
		
	return; /* back to main */
}




/////////////////////////////////////////////////////////////////////////////
/// Transforms a CAM file filename into a SAM file
/// filename - name of the CAM-formated file, without the .dat extension, to be converted into SAM format
/// path - relative path where the file is located
/////////////////////////////////////////////////////////////////////////////
inline void CAM_to_SAM(char * path, char * filename)
{	
	int flag, i,j;

	FILE * fp;
	FILE * fp_out;

	char input_file[50]="";
	strcpy(input_file,path);
	strcat(input_file,filename);
	strcat(input_file,".dat");
	
	char output_file[50]="";
	strcpy(output_file,path);
	
	char tag[50]="SAMconverted_";
	strcat(tag,filename);
	strcat(tag,".dat");
	strcat(output_file,tag);
	
	
	fp = fopen(input_file,"r");
	if( fp==NULL ){printf("CAM_to_SAM: Error opening file. Aborting...\n");exit(EXIT_FAILURE);}
	
	fp_out = fopen(output_file,"w");
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
	double average=0.0;
	double aux;
	double row_length;
	
	for(i=0; i<num_cells; i++)
	{
		/* get local density */
		fscanf(fp, "%lf", &row_length);
		/* export local density */
		fprintf(fp_out,"%.0lf\t", row_length);
		for(j=1; j<=row_length; j++)
		{
			fscanf(fp,"%lf",&aux);
			average += aux;
		}
		if( row_length !=0 )
		{
			fprintf(fp_out,"%lf\n",average/((double)row_length) );
		}else{   //else, leave the entry in this second column empty
			fprintf(fp_out,"\n");
		}
		
		/* reset for the next average line */
		average = 0.0;
		
	}
	
	/* Closing files */
	flag = fclose(fp);
	if(flag!=0)
	{
		printf("CAM_to_SAM: Error closing input file. Aborting...\n");
		exit(EXIT_FAILURE);
	}
	flag = fclose(fp_out);
	if(flag!=0)
	{
		printf("CAM_to_SAM: Error closing output file. Aborting...\n");
		exit(EXIT_FAILURE);
	}
	
	
	return;/* back to main */
}





/////////////////////////////////////////////////////////////////////////////
/// Finds cells that are symmetrically equivalent to the given cell in the cylinder slice. 
/// All indices are relative to the slice.
/// Input:
/// ny - Number of cells in the y direction in the slice
/// nz - Number of cells in the z direction in the slice
/// cell_idx - reative index of the cell of interest. It must be in the first quadrant (upper right)
/// Output:
/// sym_cells - array (int, length 4) containng the relative indices of the cells in the slice that are 
///             symmetrically equivalent to cell_idx, including cell_idx. The order is:
///	  	lower right quadrant (0), upper right quadrant (1), lower left (2), upper left (3)
/// 
/// Rationale behind the algorithm:
/// q1 = given cell index
/// q0 = given cell index - twice the relative position respect to the start of the first quadrant, which is lambda=cell_idx%(nz/2), -1 to correct for indices starting at 0
/// q2 = given cell index + twice the distance to midpoint (mu), -1 to correct for indices starting at 0
/// q3 = q2 + the quantity that was substracted from q1 to get to q0
/////////////////////////////////////////////////////////////////////////////
inline void find_sym_cells( int ny, int nz, int cell_idx, int * sym_cells )
{
	/* Defensive programming */ //find a way to do this only once, not at every call?
	if( ny%2!=0 || nz%2!=0 )
	{
		fprintf(stderr,"find_sym_cells: ny,nz need to be even numbers. Aborting..\n");
		exit(EXIT_FAILURE);
	}
	
	int nz_half = nz/2;
	int alpha = 2*(cell_idx % nz_half) + 1;
	int mu = ny*nz_half;

	
	sym_cells[0] = cell_idx - alpha;
	sym_cells[1] = cell_idx;
	sym_cells[2] = 2 * mu - cell_idx - 1; // q2=cell_idx+2*(mu-cell_idx)-1
	sym_cells[3] = sym_cells[2] + alpha;	
	
	return;
}








#endif /* header guard*/
