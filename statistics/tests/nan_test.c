#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//#include <time.h>
#include <assert.h>
#include <string.h>

#include "../../code/export_routines.h"


int main(int arg, char ** argsv){
	
	double * p;
	int size=4;
	p=malloc(size*sizeof(double));

	int i;
	for(i=0;i<size;i++)
	{
		p[i]=i*0.2;
	}
	
	/* Storing nan's*/
	p[1]=NAN;
	
	for(i=0;i<size;i++)
	{
		printf("%lf ",p[i]);
	}
	printf("\n");

	for(i=0;i<size;i++)
	{
		printf("%d ",isnan(p[i]));
	}
	printf("\n");	
	
	/* exporting nan's*/
	int header[6]={0, 2, 2, 2 ,0, 0};
	export_SAM_data_scalar(p, "nan_test", header, 0, size-1);

	/* Reading nan's*/
	FILE * fp;
	fp=fopen("nan_test_0.dat","r");
	
	double * q;
	q=malloc(size*sizeof(double));

	int nx,ny,nz,cell_idx_start,cell_idx_end,step;	
	//read header:
	fscanf( fp, "%d \t %d \t %d \t %d \t %d \t %d \n", &nx, &ny, &nz, &cell_idx_start, &cell_idx_end, &step);

	
	for(i=0;i<size;i++)
	{
		fscanf(fp, "%lf",&q[i]);
	}
	
	for(i=0;i<size;i++)
	{
		printf("%lf ",q[i]);
	}
	printf("\n");
	
	/*Cleaning up*/
	fclose(fp);
	free(p);
	free(q);
	return 0;
}
