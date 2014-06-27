#ifndef EXPORT_ROUTINES_H
#define EXPORT_ROUTINES_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//#include <time.h>
#include <assert.h>
#include <string.h>




//////////////////////////////////////////////////////////////////////////
/// This function exports the data contained in the array data in a file format suitable for calculating a SAM average with average.c 
/// Input: 
/// data - array of length l, containing the value of a scalar property. l = cell_end-cell_start+1
/// This function will export the section of the data array specified by cell_start and cell_end (both inclusive)
/// header: array of ints containing the following information: file_number nx ny nz first_cell_index step
/// where the first_cell_index is the global index of the first cell recorded in data (inclusive), or -1 if the data array constains information on the whole simulation box (i.e. all slices)
/// Output:
/// The file has the format: header \n data with:
/// header= nx \t ny \t nz \t first_cell_index (global) \t current time step
/// data= each number corresponds to a cell. There is only one line, numbers are separated by \t 
/// IMPORTANT NOTE: There is a flaw in this routine. If a cell is empty, what should the data array contain in the entry corresponding to that empty cell? 0 is not a good option.
/// A better option would be a nan. Routines handling SAM-formatted data files should expect nan's and handle them appropriately.
/// This method is "blind" to this issue. If nan's are present, it will print them to file.
//////////////////////////////////////////////////////////////////////////
// MODIFIED BUT NOT TESTED (with the new SAM format, it was tested with teh old one)!!!! TO DO
inline void export_SAM_data_scalar(double * data, char * filename, int * header, int cell_start, int cell_end)
{
	int i, flag;
	FILE * fp;
	char file_name[100];
	char filename1[100]="";
	
	/* Assemble the file name */
	strcpy(filename1,filename);
	strcat(filename1,"_%d.dat");
	snprintf(file_name,sizeof(char)*120,filename1,header[0]);
	
	/* Prepare file pointer */
	fp=fopen(file_name,"w");
	if( fp ==NULL )
	{
		printf("export_SAM_data_scalar: Error opening file. Aborting...\n");
		exit(EXIT_FAILURE);
	}
	
	int cell_idx_start, cell_idx_end;
	if( header[4] == -1 )
	{
		cell_idx_start = 0;
		cell_idx_end = cell_end - cell_start + 1;
	}else{
		cell_idx_start = cell_start;
		cell_idx_end = cell_idx_start + cell_end - cell_start;
	}
	
	/* Print header */
	fprintf(fp,"%d \t %d \t %d \t %d \t %d \t %d \n", header[1], header[2], header[3], cell_idx_start, cell_idx_end ,header[5]);

	for(i=cell_start; i<=cell_end; i++)
	{
		fprintf(fp,"%lf\t",data[i]);
	}
	fprintf(fp,"\n");

	/* Clean up*/
	flag = fclose(fp);
	if( flag !=0 )
	{
		fprintf(stderr,"export_SAM_data_scalar: error closing file\n Aborting...\n");
		exit(EXIT_FAILURE);
	}	
	
	return; /* back to main*/
}


//////////////////////////////////////////////////////////////////////////
/// This function exports the data contained in the array data in a file format suitable for calculating a SAM average with average.c 
/// Input: 
/// data - array of length l, containing the value of a vector property. l = cell_end-cell_start+1
/// This function will export the section of the data array specified by cell_start and cell_end (both inclusive).
/// header: array of ints containing the following information: file_number nx ny nz first_cell_index step
/// where the first_cell_index is the global index of the first cell recorded in data (inclusive), or -1 if the data array constains information on the whole simulation box (i.e. all slices)
/// Output:
/// The file has the format: header \n data with:
/// header= nx \t ny \t nz \t first_cell_index (global) \t current time step
/// data= each row corresponds to a cell. For each cell (row), the x, y, z coordinates of the vector quantity are separated by tabs:
/// v0x \t v0y \t v0z \n
/// v1x \t v1y \t v1z \n
/// ...
//////////////////////////////////////////////////////////////////////////
// NOT TESTED WITH NEW FORMAT
inline void export_SAM_data_vector(double ** data, char * filename, int * header, int cell_start, int cell_end)
{// Shame there's no overloading in C...
	int i, flag;
	FILE * fp;
	
	/* Prepare file name*/
	char file_name[100];
	char filename1[100]="";
	strcpy(filename1,filename);
	strcat(filename1,"_%d.dat");
	
	snprintf( file_name, sizeof(char)*120, filename1, header[0] );
	
	/* Open file */
	fp=fopen(file_name,"w");
	if( fp ==NULL )
	{
		printf("export_SAM_data_vector: Error opening file. Aborting...\n");
		exit(EXIT_FAILURE);
	}
	
	int cell_idx_start, cell_idx_end;
	if( header[4] == -1 )
	{
		cell_idx_start = 0;
		cell_idx_end = cell_end - cell_start + 1;
	}else{
		cell_idx_start = cell_start;
		cell_idx_end = cell_idx_start + cell_end - cell_start;
	}
	
	/* Print header */
	fprintf(fp,"%d \t %d \t %d \t %d \t %d \t %d \n", header[1], header[2], header[3], cell_idx_start, cell_idx_end ,header[5]);

	/* Export data */
	for(i=cell_start; i<=cell_end; i++)
	{
		fprintf(fp,"%lf\t%lf\t%lf\n",data[i][0],data[i][1],data[i][2]);
	}

	/* Clean up*/
	flag = fclose(fp);
	if( flag !=0 )
	{
		fprintf(stderr,"export_SAM_data_vector: error closing file\n Aborting...\n");
		exit(EXIT_FAILURE);
	}	
		
	return; /* back to main*/
}

//////////////////////////////////////////////////////////////////////////
/// This function exports the data contained in the array data in a file format suitable for calculating a CAM average with average.c.
/// The array data is a 2D array with each row containing the scalar/vector values of interest for the particles in the corresponding collision cell. 
/// Cells are stored in consecutive order. The occupation numbers (local instantaneous densities), ie length of the rows, are stored in the array of 
/// ints local_density.
/// The header of the output file will contain the fillowing information:
/// nx \t ny \t nz \t (global index of the first cell) \t (global index of the last cell) \t step
/// Input: 
/// data: (cells)x(variable+1) array containing the data to be exported to file. Each row corresponds to a cell. data[ci][0] is the local density in cell ci. alpha is 
/// the number of particles in each cell (data[ci][0]=alpha). If data is scalar, an example could be: 3 e1 e2 e3, for row ci. If the quantity is a vector, for example the velocities, 
/// the the row for cell ci would look like: 3 v1x v1y v1z v2x v2y v2z v3x v3y v3z.
/// cell_start: index (within the data) of the first cell that will be exported. 
/// cell_end: index (within the data) of the last cell that will be exported.
/// header: array of ints containing the following information: file_number nx ny nz first_cell_index step
/// where the first_cell_index is the global index of the first cell recorded in data (inclusive), or -1 if the data array contains information on the whole simulation box (i.e. all slices)
//////////////////////////////////////////////////////////////////////////
inline void export_CAM_data(int is_scalar, double ** data, char * filename, int * header, int cell_start, int cell_end)
{
	FILE * fp;
	
	/* Prepare file name*/
	char file_name[100];
	char filename1[100]="";
	strcpy(filename1,filename);
	strcat(filename1,"_%d.dat");
	snprintf(file_name,sizeof(char)*120,filename1,header[0]);
	
	int i,j,local_density, flag;
	
	/* Open file name */
	fp=fopen(file_name,"w");
	if(fp==NULL)
	{
		printf("export_CAM_data: Error opening file. Aborting...\n");
		exit(EXIT_FAILURE);
	}
	
	int cell_idx_start, cell_idx_end;
	if( header[4] == -1 )
	{
		cell_idx_start = 0;
		cell_idx_end = cell_end - cell_start + 1;
	}else{
		cell_idx_start = cell_start;
		cell_idx_end = cell_idx_start + cell_end - cell_start;
	}
	
	/* Print header */
	//fprintf(fp,"%d \t %d \t %d \t %d \t %d \t %d \n", header[1], header[2], header[3], cell_start, cell_end, header[4]);
	fprintf(fp,"%d \t %d \t %d \t %d \t %d \t %d \n", header[1], header[2], header[3], cell_idx_start, cell_idx_end ,header[5]);

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
	flag = fclose(fp);	
	if( flag !=0 )
	{
		fprintf(stderr,"export_CAM_data: error closing file\n Aborting...\n");
		exit(EXIT_FAILURE);
	}	
	return; /* back to main*/
}


#endif /* header guard */
