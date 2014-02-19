#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>


/*nynz is the product of the number of cells in the y direction and in the z direction*/
/* a is the length of a side of each cubic collision cell*/
/* That slice includes the cells of indices from slice_start_index return to that value-1+nynz, inclusive*/
/*inline int slice_start_index(double x, int nynz, double a)
{	
	return nynz*(int)floor(x/a);
}*/



/* This function return the Cumulative Average of a number samples of files of the form: filename_num.dat, with 
num from 1 to samples */
/* This function does not check that all the files correspond to the same simulation, or whether they are "averageable" or not */
inline void slice_average_CAM(char * filename, int samples )
{
	/*Read header of first file to obtain size of slice */
	FILE * fp;
	char filename1[50]="";
	strcpy(filename1,filename);
	strcat(filename1,"_1.dat");	
	
	fp=fopen(filename1,"r"); 
	if( fp == NULL ){ perror("Error while opening the file.\n"); exit(EXIT_FAILURE);}
 
	double x_slice=0.0;
	int nz=0;
	int nynz=0;
	int step=0;
	
	fscanf( fp, "%lf \t %d \t %d \t %d \n", &x_slice, &nynz, &nz, &step );
	
	
	/* Allocate memory */
	double * vel_ave;
	int * part_num;
	
	vel_ave = malloc( nynz * sizeof(double) );
	if (vel_ave==NULL) {printf("Error allocating vel_ave in slice_average_CAM call \n"); exit(EXIT_FAILURE);}
	part_num = malloc( nynz * sizeof(int) );
	if (part_num==NULL) {printf("Error allocating part_num in slice_average_CAM call \n"); exit(EXIT_FAILURE);}
	
	/* Initialization */
	int i,j,k;
	
	for(i=0; i<nynz; i++ )
	{
		vel_ave[i] = 0.0;
		part_num[i] = 0;
	}
	
	/* Read data from all files and start the average */
	float row_length;
	double vel_aux;
	char intaschar[20];

	
	for(i=1; i<=samples; i++) //file 1 is open the first time that we pass through here	
	{
		for(j=0;j<nynz;j++)
		{
			fscanf(fp, "%f", &row_length);
			
			part_num[j] += (int)row_length;
		
			for(k=1; k<=(int)row_length; k++)
			{
				fscanf(fp, "%lf", &vel_aux);
				vel_ave[j]+=vel_aux;
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
		fscanf( fp, "%lf \t %d \t %d \t %d \n", &x_slice, &nynz, &nz, &step );
		
	}
	
	
	/* Complete the average and export to file */
		
	fopen("./DATA/vel_CAM_averaged.dat","w");
	
	fprintf( fp, "%lf \t %d \t %d \t %d \n", x_slice, nynz, nz, step );
	for(i=0;i<nynz;i++)
	{
		if( part_num[i]!=0)
		{
			vel_ave[i] = vel_ave[i]/(double)part_num[i];
		}else{
			vel_ave[i] = 0.0;
		}
		fprintf(fp,"%lf\n", vel_ave[i]);
	}
	
	//fprintf(fp,"\n");
	
	/* Clean up */
	free(vel_ave);
	free(part_num);
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
	
	slice_average_CAM(filename, atoi(argv[2]) );

		
	return 0;
}
