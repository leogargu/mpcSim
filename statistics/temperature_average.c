#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>


/* need to modify the function that exports the temperature to include a header with nx, ny, nz and the timestep */
inline void average_SAM(char * filename, int samples )
{
	/*Read header of first file to obtain size of slice */
	FILE * fp;
	char filename1[50]="";
	strcpy(filename1,filename);
	strcat(filename1,"_1.dat");	
	
	fp=fopen(filename1,"r"); 
	if( fp == NULL ){ perror("Error while opening the file.\n"); exit(EXIT_FAILURE);}
 
	int nx=0;
	int ny=0;
	int nz=0;
	long int step=0;
	
	fscanf( fp, "%d \t %d \t %d \t %ld \n", &nx, &ny, &nz, &step );
	
	int n_cells= nx*ny*nz;	
	
	/* Allocate memory */
	double * average;
	
	average = malloc( n_cells * sizeof(double) );
	if (average==NULL) {printf("Error allocating average in average_SAM call \n"); exit(EXIT_FAILURE);}
	
	/* Initialization */
	int i,j,k;
	
	for(i=0; i<n_cells; i++ )
	{
		average[i] = 0.0;
	}
	
	/* Read data from all files and start the average */
	double aux;
	char intaschar[20];
	
	//printf("ncells=%d\n",n_cells);
	
	for(i=1; i<=samples; i++) //file 1 is open the first time that we pass through here	
	{
		for(j=0;j<n_cells;j++)
		{
			//printf("hello!\n");
			fscanf(fp, "%lf", &aux);
			//printf("I read this %lf \n",aux);
			if(aux >=0.0 )
			{
				average[j] += aux;
				//printf("adding %lf to average\n",aux);
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
		if( fp == NULL ){ perror("Error while opening the file.\n"); exit(EXIT_FAILURE);}

		/* read fist line to set the pointer in the right position. Checking of parameters for robustness is possible here */
		fscanf( fp, "%d \t %d \t %d \t %ld \n", &nx, &ny, &nz, &step );
		
	}
	
	
	/* Complete the average and export to file */
		
	fopen("./DATA/temperature_SAM_averaged.dat","w");
	
	fprintf( fp, "%d \t %d \t %d \n", nx, ny, nz );
	
	double factor=1.0/(double)samples;
	for(i=0; i<n_cells; i++)
	{
		fprintf(fp, "%lf\t", factor*average[i] );
	}
		
	/* Clean up */
	free(average);
	fclose(fp);
	
	return; /* back to main */
}



/*-------------------------------*/
/*		MAIN		 */
/*-------------------------------*/

int main(int argc, char **argv) {

	/*
	argv[0] is the program name
	argv[1] is name of files
	argv[2] is the number of samples 
	*/
	char filename[50]="./DATA/";
	strcat(filename,argv[1]);
	
	average_SAM(filename, atoi(argv[2]) );

		
	return 0;
}
