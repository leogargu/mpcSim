#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <string.h>

#include "mpc.h"    /*  Struct definitions, function prototypes and global variables that must be made available to other files  */


#include "macros.h"




//////////////////////////////////////////////////////////////////////////
/// This function exports the data contained in the array data in a file format suitable for calculating a SAM average with average.c 
/// Input: 
/// data - array of length l, containing the value of a scalar property. l = cell_end-cell_start+1
/// This function will export the section of the data array specified by cell_start and cell_end (both inclusive)
/// header: array of ints containing the following information: file_number nx ny nz step
/// Output:
/// The file has the format: header \n data with:
/// header= nx \t ny \t nz \t cell_start \t cell_end \t current time step
/// data= each number corresponds to a cell. There is only one line, numbers are separated by \t 
//////////////////////////////////////////////////////////////////////////
//inline void export_SAM_data_scalar(double * data, char * filename, Geometry geometry, int cell_start, int cell_end, int file_number, int step)
inline void export_SAM_data_scalar(double * data, char * filename, int * header, int cell_start, int cell_end)
{
	int i;
	FILE * file_fp;
	char file_name[50];
	char filename1[50]="";
	
	/* Assemple the file name */
	strcpy(filename1,filename);
	strcat(filename1,"_%d.dat");
	snprintf(file_name,sizeof(char)*30,filename1,header[0]);
	
	/* Prepare file pointer */
	file_fp=fopen(file_name,"w");
	if( file_fp ==NULL )
	{
		printf("export_SAM_data_scalar: Error opening file. Aborting...\n");
		exit(EXIT_FAILURE);
	}
	
	/* Print header */
	fprintf(file_fp,"%d \t %d \t %d \t %d \t %d \t %d \n", header[1], header[2], header[3], cell_start, cell_end, header[4]);
	for(i=cell_start; i<=cell_end; i++)
	{
		fprintf(file_fp,"%lf\t",data[i]);
	}
	fprintf(file_fp,"\n");

	/* Clean up*/
	fclose(file_fp);		
		
	return; /* back to main*/
}


//////////////////////////////////////////////////////////////////////////
/// This function exports the data contained in the array data in a file format suitable for calculating a SAM average with average.c 
/// Input: 
/// data - array of length l, containing the value of a vector property. l = cell_end-cell_start+1
/// This function will export the section of the data array specified by cell_start and cell_end (both inclusive).
/// header: array of ints containing the following information: file_number nx ny nz step
/// Output:
/// The file has the format: header \n data with:
/// header= nx \t ny \t nz \t cell_start \t cell_end \t current time step
/// data= each row corresponds to a cell. For each cell (row), the x, y, z coordinates of the vector quantity are separated by tabs:
/// v0x \t v0y \t v0z \n
/// v1x \t v1y \t v1z \n
/// ...
//////////////////////////////////////////////////////////////////////////
inline void export_SAM_data_vector(double ** data, char * filename, int * header, int cell_start, int cell_end)
{// Shame there's no overloading in C...
	int i;
	FILE * file_fp;
	
	/* Prepare file name*/
	char file_name[50];
	char filename1[50]="";
	strcpy(filename1,filename);
	strcat(filename1,"_%d.dat");
	snprintf(file_name,sizeof(char)*30,filename1,header[0]);
	
	/* Open file */
	file_fp=fopen(file_name,"w");
	if( file_fp ==NULL )
	{
		printf("export_SAM_data_vector: Error opening file. Aborting...\n");
		exit(EXIT_FAILURE);
	}
	
	/* Print header*/
	fprintf(file_fp,"%d \t %d \t %d \t %d \t %d \t %d \n",header[1],header[2],header[3], cell_start, cell_end, header[4]);
	
	/* Export data */
	for(i=cell_start; i<=cell_end; i++)
	{
		fprintf(file_fp,"%lf\t%lf\t%lf\n",data[i][0],data[i][1],data[i][2]);
	}

	/* Clean up*/
	fclose(file_fp);		
		
	return; /* back to main*/
}

//////////////////////////////////////////////////////////////////////////
/// This function exports the data contained in the array data in a file format suitable for calculating a CAM average with average.c.
/// The array data is a 2D array with each row containing the scalar/vector values of interest for the particles in the corresponding collision cell. 
/// Cells are stored in consecutive order. The occupation numbers (local instantaneous densities), ie length of the rows, are stored in the array of 
/// ints local_density.
/// Input: 
/// data: (cells)x(variable+1) array containing the data to be exported to file. Each row corresponds to a cell. data[ci][0] is the local density in cell ci. alpha is 
/// the number of particles in each cell (data[ci][0]=alpha). If data is scalar, an example could be: 3 e1 e2 e3, for row ci. If the quantity is a vector, for example the velocities, 
/// the the row for cell ci would look like: 3 v1x v1y v1z v2x v2y v2z v3x v3y v3z.
/// cell_start: index (within the data) of the first cell that will be exported. 
/// cell_end: index (within the data) of the last cell that will be exported.
/// file_number: number to include in the file name
/// step: timestep at which this data was gathered
//////////////////////////////////////////////////////////////////////////
inline void export_CAM_data(int is_scalar, double ** data, char * filename, int * header, int cell_start, int cell_end)
{
	FILE * fp;
	
	/* Prepare file name*/
	char file_name[50];
	char filename1[50]="";
	strcpy(filename1,filename);
	strcat(filename1,"_%d.dat");
	snprintf(file_name,sizeof(char)*30,filename1,header[0]);
	
	int i,j,local_density;
	
	/* Open file name */
	fp=fopen(file_name,"w");
	if(fp==NULL)
	{
		printf("export_CAM_data: Error opening file. Aborting...\n");
		exit(EXIT_FAILURE);
	}
	
	/* Print header */
	fprintf(fp,"%d \t %d \t %d \t %d \t %d \t %d \n", header[1], header[2], header[3], cell_start, cell_end, header[4]);
	
	/* Export to file */
	for(i=cell_start; i<=cell_end; i++) // i is the cell in the slice of interest
	{	
		local_density = (int)data[i][0];
		fprintf(fp, "%d\t", local_density);
		if( local_density > 0 )
		{
			for(j=1; j<=local_density; j++ ) // j is the particle in cell i
			{
				if( is_scalar ==1 )
				{
					fprintf(fp,"%lf\t",data[i][j]);
				}else if(is_scalar == 0)
				{
					fprintf(fp,"%lf\t%lf\t%lf\t",data[i][3*(j-1)+1],data[i][3*(j-1)+2],data[i][3*(j-1)+3]);
					
				}else{
					fprintf(stderr, "export_CAM_data: Option not recognised is_scalar should be 1 or 0. Aborting...\n");
					exit(EXIT_FAILURE);
				}
			}
			fprintf(fp,"\n");
		}else{
			fprintf(fp,"\n");
	  	}
	}
	
	/* Clean up and release memory*/
	fclose(fp);	
		
	return; /* back to main*/
}



#define TEST_CAM 1
#define TEST_SAM 1


int main(int argc, char **argv) {

	srand(time(NULL));
	
	
	int n_cells=8;
	int maxocc=3;
	double number_range = 100.0;

	n_cells=atoi(argv[1]);
	maxocc =atoi(argv[2]);
	number_range=atof(argv[3]);
	
	//printf("N_cells=%d,maxocc=%d,number_range=%lf\n",n_cells,maxocc,number_range);

	int cell_start = n_cells*(1.0*rand()/RAND_MAX);
	int cell_end = cell_start + (int)(n_cells-cell_start)*(1.0*rand()/RAND_MAX);
	
	printf("start: %d, end: %d\n",cell_start,cell_end);
	
	double a= 1.0;
	double L=8.0;
	double Lx=12.0;
	double radius=0.5*(L-a);
	double verlet_epsilon = (1e-5)*radius; 
	double verlet_min_sq = radius - verlet_epsilon;
	verlet_min_sq *=verlet_min_sq;
	double verlet_max_sq = radius + verlet_epsilon;
	verlet_max_sq *= verlet_max_sq;
	Geometry cylinder={a,Lx,L,0.5*L,0.5*( L - a ), verlet_epsilon, verlet_min_sq, verlet_max_sq, n_cells, { (int)(Lx/a), (int)(L/a), (int)(L/a)} }; 
	


	int header[5] = { 0, cylinder.n_cells_dim[0],cylinder.n_cells_dim[1],cylinder.n_cells_dim[2], 1};

	double * SAM_scalar;
	double ** SAM_vector;
	double ** array_scalar;
	double ** array_vector;
	int i,j;
	int density = 0;
	array_scalar = malloc(n_cells*sizeof(double*));
	if( array_scalar == NULL){ printf("malloc error\n");exit(EXIT_FAILURE);}
	array_vector = malloc(n_cells*sizeof(double*));
	if( array_vector == NULL){ printf("malloc error\n");exit(EXIT_FAILURE);}
	SAM_scalar = malloc(n_cells*sizeof(double));
	if( SAM_scalar == NULL){ printf("malloc error\n");exit(EXIT_FAILURE);}
	SAM_vector = malloc(n_cells*sizeof(double*));
	if( SAM_vector == NULL){ printf("malloc error\n");exit(EXIT_FAILURE);}

	for(i=0; i<n_cells; i++)
	{
		SAM_vector[i] = malloc(3*sizeof(double));
		if( SAM_vector[i] ==NULL){printf("malloc error\n");exit(EXIT_FAILURE);}
		array_scalar[i] = malloc((maxocc+1)*sizeof(double));
		if( array_scalar[i] ==NULL){printf("malloc error\n");exit(EXIT_FAILURE);}
		array_vector[i] = malloc((1+maxocc*3)*sizeof(double));
		if( array_vector[i] ==NULL){printf("malloc error\n");exit(EXIT_FAILURE);}
	}
	
	/* Initialise */
	for(i=0; i<n_cells; i++)
	{
		density=(int)((maxocc+1)*(1.0*rand()/RAND_MAX)); //number between 0 and maxocc, inclusive
		array_scalar[i][0] = density;
		array_vector[i][0] = density;
		
		//counter=1;
		for(j=1; j<=density; j++)
		{
			array_scalar[i][j] = number_range*rand()/RAND_MAX;
			array_vector[i][3*(j-1)+1] = number_range*rand()/RAND_MAX;
			array_vector[i][3*(j-1)+2] = number_range*rand()/RAND_MAX;
			array_vector[i][3*(j-1)+3] = number_range*rand()/RAND_MAX;
		}
		
		SAM_scalar[i] = 10.0*rand()/RAND_MAX;
		
		for(j=0; j<3; j++)
		{
			SAM_vector[i][j] = number_range*rand()/RAND_MAX;
		}
	}
	
	/* Print to screen */
#if TEST_CAM
	for(i=cell_start; i<cell_end; i++)
	{
		printf("%.0lf",array_scalar[i][0]);
			
		for(j=1; j<=array_scalar[i][0]; j++)
		{
			printf("\t%lf",array_scalar[i][j]);
		}
		printf("\n");
	}	
	

	for(i=cell_start; i<cell_end; i++)
	{
		printf("%.0lf",array_vector[i][0]);
			
		for(j=1; j<=array_vector[i][0]; j++)
		{
			printf("\t%lf\t%lf\t%lf",array_vector[i][3*(j-1)+1],array_vector[i][3*(j-1)+2],array_vector[i][3*(j-1)+3]);
		}
		printf("\n");
	}
	
	/* Use export routines */
	export_CAM_data(1, array_scalar, "CAM_scalar", header, cell_start, cell_end);
	export_CAM_data(0, array_vector, "CAM_vector", header, cell_start, cell_end); //defensive programming in case it is called like this?
#endif
#if TEST_SAM
	for(i=cell_start; i<=cell_end; i++)
	{
		printf("%lf\t", SAM_scalar[i]);
	}
	printf("\n");
	for(i=cell_start; i<=cell_end; i++)
	{
		printf("%lf\t%lf\t%lf\n",SAM_vector[i][0],SAM_vector[i][1],SAM_vector[i][2]);
	}
	
	/* Use export routines */
	export_SAM_data_scalar(SAM_scalar, "SAM_scalar", header, cell_start, cell_end);
	export_SAM_data_vector(SAM_vector, "SAM_vector", header, cell_start, cell_end);
#endif
	
	
	
	/* Clean up */
	for(i=0; i<n_cells; i++)
	{
		free(array_scalar[i]);
	}
	free(array_scalar);
	

	for(i=0; i<n_cells; i++)
	{
		free(array_vector[i]);
	}
	free(array_vector);
	

	free(SAM_scalar);
	
	for (i=0; i<n_cells; i++)
	{
		free(SAM_vector[i]);
	}
	free(SAM_vector);
	
	return 0;
}

		
