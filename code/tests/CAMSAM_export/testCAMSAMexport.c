#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <string.h>

#include "./../../mpc.h"    /*  Struct definitions, function prototypes and global variables that must be made available to other files  */


#include "./../../macros.h"

#include "./../../export_routines.h"



// Compile as:
// gcc testCAMSAMexport.c -o test
// Call as:
// ./test 10, 5, 20 a > screen.out
// where a=0,1,2,3 for SAM scalar, SAM vector, CAM scalar or CAM vector respectively
// Compare output with exported files by running:
// bash test_export.sh
// Modify variables in the .sh script to do stress testing
int main(int argc, char **argv) {

	srand(time(NULL));
	
	
	int n_cells=8;
	int maxocc=3;
	double number_range = 100.0;

	n_cells=atoi(argv[1]);
	maxocc =atoi(argv[2]);
	number_range=atof(argv[3]);
	
	int test_case=atoi(argv[4]);
	assert(test_case==0 || test_case==1 || test_case==2 || test_case==3);
	

	//int cell_start = n_cells*(1.0*rand()/RAND_MAX);
	//int cell_end = cell_start + (int)(n_cells-cell_start)*(1.0*rand()/RAND_MAX);
	
	int cell_start = 0; 
	int cell_end = n_cells-1;
	
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
	
	/*----------------------------------------------------------------------------------------*/
	/* Print to screen */
	if(test_case==2)
	{
		for(i=cell_start; i<=cell_end; i++)
		{
			printf("%.0lf",array_scalar[i][0]);
				
			for(j=1; j<=array_scalar[i][0]; j++)
			{
				printf("\t%lf",array_scalar[i][j]);
			}
			printf("\n");
		}
		/* Use export routines */
		export_CAM_data(1, array_scalar, "CAM_scalar", header, cell_start, cell_end);
	}else if(test_case==3)
	{

		for(i=cell_start; i<=cell_end; i++)
		{
			printf("%.0lf",array_vector[i][0]);
			
			for(j=1; j<=array_vector[i][0]; j++)
			{
				printf("\t%lf\t%lf\t%lf",array_vector[i][3*(j-1)+1],array_vector[i][3*(j-1)+2],array_vector[i][3*(j-1)+3]);
			}
			printf("\n");
		}
		/* Use export routines */
		export_CAM_data(0, array_vector, "CAM_vector", header, cell_start, cell_end); //defensive programming in case it is called like this?
	}
	
	if(test_case==0)
	{
		for(i=cell_start; i<=cell_end; i++)
		{
			printf("%lf\t", SAM_scalar[i]);
		}
		printf("\n");
		/* Use export routines */
		export_SAM_data_scalar(SAM_scalar, "SAM_scalar", header, cell_start, cell_end);
	
	}else if(test_case==1)
	{
		for(i=cell_start; i<=cell_end; i++)
		{
			printf("%lf\t%lf\t%lf\n",SAM_vector[i][0],SAM_vector[i][1],SAM_vector[i][2]);
		}
		/* Use export routines */
		export_SAM_data_vector(SAM_vector, "SAM_vector", header, cell_start, cell_end);
	}
	
	/*----------------------------------------------------------------------------------------*/	
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

		
