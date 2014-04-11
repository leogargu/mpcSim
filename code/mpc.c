/** MPC simulations.
 *
 * Author: Leonor Garcia-Gutierrez, 2013
 * Version: mpcSim-1.0
 *
 *
 * METHODS:
 * populate
 * export_data
 * export_vtk_plasma
 * export_vtk_gid
 * export_vessel_geometry
 *
 * mpc_paraview.py IS NOT CREATED BY THIS CODE.
 *
 * To lauch paraview from the python script: paraview --script=mpc_paraview.py &
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include<time.h>
#include<assert.h>
#include <string.h>

#include "mpc.h"    /*  Struct definitions, function prototypes and global variables that must be made available to other files  */

/* GNU GSL*/
#include <gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
/* See available algorithms for GSL in: */
/*  http://www.gnu.org/software/gsl/manual/html_node/Random-number-generator-algorithms.html   */


#include "stream.h"
#include "collide.h"
#include "macros.h"

/*-------------------------------*/
/* 		METHODS   	 */
/*-------------------------------*/

// 0,1,2 equivalent to x, y, z

//////////////////////////////////////////////////////////////////////////
/// Returns the norm squared of a 3D vector 
//////////////////////////////////////////////////////////////////////////
inline double norm_sq(double * v)
{
	return v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
}

//////////////////////////////////////////////////////////////////////////
/// Given a value of x, this function returns the global index of the first cell (bottom right) in the slice corresponding to that x.
/// nynz is the product of the number of cells in the y direction and in the z direction. a is the length of a side of each cubic collision cell. 
/// That slice includes the cells of indices from slice_start_index return to that value-1+nynz, inclusive.
//////////////////////////////////////////////////////////////////////////
inline int slice_start_index(double x, int nynz, double a)
{	
	return nynz*(int)floor(x/a);
}

//////////////////////////////////////////////////////////////////////////
/// Generates a random position vector uniformly distributed inside of the cylinder given by the Geometry struct
/// Output: pos
/// Adapted from: http://stackoverflow.com/questions/5837572/generate-a-random-point-within-a-circle-uniformly
//////////////////////////////////////////////////////////////////////////
inline void populate(gsl_rng * r, Geometry cylinder, double * pos) //could be modified using the fact that the sum of two uniformly distributed numbers is symmetrically triangularly distributed
{
	double t,u,rvar;
	
	t = 2*M_PI*gsl_rng_uniform_pos(r); 
	u = cylinder.radius*( gsl_rng_uniform_pos(r) + gsl_rng_uniform_pos(r) ); // triangular distribution?
	
	if( u > cylinder.radius )
	{
		/* fold the triangle up */
		rvar = 2*cylinder.radius - u;
	}else{
		/* just accept */
		rvar = u;
	}		
		
	pos[0] = cylinder.Lx * gsl_rng_uniform_pos(r);
	pos[1] = cylinder.L_half + rvar * cos(t);
	pos[2] = cylinder.L_half + rvar * sin(t);
	
	return; /* back to main */
}

/* Notice this can be called without undoing the shift after collide */
/* NOT TESTED */
/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Control routine that checks for compressibility effects: it stopsthe simulation if, at the time of being invoked, the density in any 
/// collision cell is greater than density+X or smaller than density-X with X=DENSITY_TOL*density.
/////////////////////////////////////////////////////////////////////////////////////////////////////
inline void check_compressibility(double * cell_mass, int n_cells, double m_inv, double density)
{
	int ci;
	double tolerance = DENSITY_TOL*density;
	for(ci=0; ci< n_cells; ci++)
	{
		//if(m_inv*cell_mass[ci] > DENSITY_TOL*density)
		if(m_inv*cell_mass[ci] > (density+tolerance) || m_inv*cell_mass[ci] < (density-tolerance) )
		{
			fprintf(stderr,"Compressibility effects exceeded tolerance: *shifted* cell %d has instantaneous density of %.2lf\n",ci,cell_mass[ci]*m_inv);
			exit(EXIT_FAILURE);
		}
	}
	
	return; /* Back to main*/
}


/* TODO: at the moment the export is done in ASCII only, binary option has not been implemented yet*/
/* possible to do by keeping the lines in ASCII and writing the data (only) in binary to the same file stream*/
//////////////////////////////////////////////////////////////////////////
/// It writes a VTK file with positions and velocities of all particles 
//////////////////////////////////////////////////////////////////////////
inline void export_vtk_plasma(int id, double ** pos, double ** vel,double ** acc, int n_part)
{
	int i;
	
	FILE * fp;
	char pos_name[30];
	//snprintf(pos_name,sizeof(char)*30,"plasma%04d.vtk",id); //leading zeroes: not pvpython-friendly
	snprintf(pos_name,sizeof(char)*30,"./DATA/plasma%d.vtk",id);

	fp=fopen(pos_name,"w");
	
	fprintf(fp,"# vtk DataFile Version 3.0\nPlasma snapshot\n");
	#if BINARY_EXPORT		//AQUI
		fprintf(stderr,"Feature not yet implemented\n");
	        exit(EXIT_FAILURE);
		fprintf(fp,"BINARY\n");
	#else    
	
	fprintf(fp,"ASCII\n");	
	fprintf(fp,"DATASET UNSTRUCTURED_GRID\n\n");
	fprintf(fp,"POINTS %d double\n",n_part);
	
	for(i=0; i<n_part; i++)
	{
		fprintf(fp,"%.6lf %.6lf %.6lf\n",pos[i][0],pos[i][1],pos[i][2]);
	}
	
	fprintf(fp,"\nCELLS 1 %d\n%d",n_part+1,n_part);
	for(i=0; i<n_part; i++)
	{
		fprintf(fp," %d",i);
	}
	fprintf(fp,"\n\n");
	
	fprintf(fp,"CELL_TYPES 1\n2\n\nPOINT_DATA %d\nVECTORS velocity double\n",n_part);
	
	for(i=0; i<n_part; i++)
	{
		fprintf(fp,"%.6lf %.6lf %.6lf\n",vel[i][0],vel[i][1],vel[i][2]);
	}
	
	fprintf(fp,"\n\n");
	
	fprintf(fp,"VECTORS acceleration double\n");
	
	for(i=0; i<n_part; i++)
	{
		fprintf(fp,"%.6lf %.6lf %.6lf\n",acc[i][0],acc[i][1],acc[i][2]);
	}
	
	#endif
	
	fclose(fp);
	
	return; /* back to main */
}


/* 
/* No need to do in in binary, unless the geometry becomes much more complicated (i.e. some triangulated surface)*/
/* NOTE that if the directory DATA is not present, the program crashes with a Segmentation Fault */
//////////////////////////////////////////////////////////////////////////
/// This exports the VTK file defining the grid used (not the vessel) 
/// This is called only once at the start of the simulation. It exports in ASCII.
//////////////////////////////////////////////////////////////////////////
inline void export_vtk_gid(Geometry cylinder)
{
	FILE * fp;
	fp = fopen("./DATA/collision_grid.vtk","w");
	if( fp == NULL )
	{
		printf("export_vtk_grid: Error opening ./DATA/collision_grid.vtk. Check directory exists. Aborting...\n");
		exit(EXIT_FAILURE);
	}
	fprintf(fp,"# vtk DataFile Version 2.0\nGrid representation\nASCII \nDATASET STRUCTURED_POINTS\nDIMENSIONS %d %d %d\nORIGIN 0 0 0\nSPACING %.1lf %.1lf %.1lf",cylinder.n_cells_dim[0]+1,cylinder.n_cells_dim[1]+1,cylinder.n_cells_dim[2]+1,cylinder.a,cylinder.a,cylinder.a);
	fclose(fp);	
}

//////////////////////////////////////////////////////////////////////////
/// This function exports the data necessary for the python script to generate a cylinder of the right dimensions on paraview.
//////////////////////////////////////////////////////////////////////////					
inline void export_vessel_geometry(Geometry cylinder, int num_steps)
{
	FILE * fp;
	fp = fopen("./DATA/vessel_geometry.py","w");
	if( fp == NULL )
	{
		printf("export_vessel_geometry: Error opening ./DATA/vessel_geometry.vtk. Check directory exists. Aborting...\n");
		exit(EXIT_FAILURE);
	}
	fprintf(fp,"# Cylinder data\n");
	fprintf(fp,"Lx=%lf\n",cylinder.Lx);
	fprintf(fp,"L=%lf\n",cylinder.L);
	//fprintf(fp,"a=%lf\n",cylinder.a);
	fprintf(fp,"radius=%lf\n",cylinder.radius);
	fprintf(fp,"L_half=%lf\n",cylinder.L_half);
	fprintf(fp,"num_frames=%d\n",num_steps);
	fclose(fp);
}

/* export is done in ASCII or BINARY depending on the value of the constant
   BINARY_EXPORT defined in macros.h					*/
/* This function is a helpful auxiliary tool, even if not used in a particular version of mpcSIM: DO NOT DELETE*/
//////////////////////////////////////////////////////////////////////////
/// Exports the lengthx3 array data into a file given by the pointer fp.
/// Can be used to export the positions/velocities/accelerations of all particles at a given timestep. 
//////////////////////////////////////////////////////////////////////////
inline void export_data(FILE * fp, double ** data ,int length )//length=n_part or n_cells
{
	int i;
	
	for(i=0; i<length; i++)
	{
		#if BINARY_EXPORT
			fwrite(&(data[i][0]), sizeof(double), 3, fp);
		#else
			fprintf(fp,"%.6lf \t %.6lf \t %.6lf \t",data[i][0],data[i][1],data[i][2]);
		#endif
	}
	#if BINARY_EXPORT
		char str[]="\n";
		fwrite(str,sizeof(char),sizeof(str),fp);
	#else
		fprintf(fp,"\n");
	#endif
	
	return; /* back to main */
}



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



//AQUI

/* Calculates the total linear momentum of the system per particle. Returns a 3D vector */
/*momentum_output contains the 3 components of the total momentum, followed by the three variances */
inline void total_momentum(int n_part, double m, double ** vel, double * momentum_output)
{
	/* Define variables */
	int i,j;
	double factor = 1.0/(double)n_part;
	
	/* Initialise output vector */
	for(j=0; j<6; j++)
	{
		momentum_output[j]=0.0;
	}
	
	/* Add momentums */
	for( i=0; i<n_part; i++) //loop over particles
	{
		for(j=0; j<3; j++) // loop over cartesian components of the velocity
		{
			momentum_output[j]+=m*vel[i][j];
		}
	}
	
	/* Calculate momentum per particle */
	for(j=0; j<3; j++)
	{
		momentum_output[j]*=factor; //momentum per particle
	}
	
	for(i=0; i<n_part; i++)
	{
		for(j=3; j<6; j++)
		{
			momentum_output[j] += (m*vel[i][j-3] - momentum_output[j-3])*(m*vel[i][j-3] - momentum_output[j-3]);
		}
	}
	
	for(j=3; j<6; j++)
	{
		momentum_output[j]*=factor; //momentum per particle
	}
	
	return; /* Back to main */
}



/* Calculates the total energy per particle, over the whole system, at the time of the call*/
/* output vector contains the total kinetic energy of the whole system (output[0]), and the 
variance (output[1])	*/
inline void total_kinetic_energy(int n_part, double m, double ** vel, double * output)
{
	double energy = 0.0;
	int i;
	double sigma = 0.0;
	
	for(i=0; i<n_part; i++)
	{
		//energy += m * ( vel[i][0]*vel[i][0] + vel[i][1]*vel[i][1] + vel[i][2]*vel[i][2] );
		energy += m * norm_sq(vel[i]);
	}
	energy = 0.5*energy/(double)n_part;
	
	for(i=0; i<n_part; i++)
	{
		sigma += (0.5*m*norm_sq(vel[i])-energy)*(0.5*m*norm_sq(vel[i])-energy);
	}
	
	/*debugging*/
	/*double aux=0.0;
	for(i=0; i<n_part; i++)
	{
		aux += (0.5*m*norm_sq(vel[i]))*(0.5*m*norm_sq(vel[i]));
	}
	
	aux=aux/(double)(n_part);*/
	
	
	
	sigma = sigma /(double)n_part;
	
	/*assert( is_approx_zero((aux-energy*energy)-sigma,1e-5) );*/
	
	output[0] = energy; //energy per particle
	output[1] = sigma; //variance (sigma^2)
	
	return;
}




////////////////////////////////////////////////////////////////////////////////////////////
///
////////////////////////////////////////////////////////////////////////////////////////////		
/* Maps the particle positions to the cells they ar at, after the grid has been shifted by a vector shift.
   It fills in the vector c_p, cell_mass and cell_occupation (CAM-like) 		*/
/* MODYDIFY THIS */
//inline void encage(double ** pos, double * shift, int n_part, Geometry cylinder, int * c_p, double * cell_mass, double m, int calculate_cell_occupation , int ** cell_occupation)
inline void encage(double ** pos, double * shift, int n_part, Geometry cylinder, int * c_p, int calculate_cell_occupation , int ** cell_occupation)
{
	int i;
	int cell_idx;
	
	/* Clean arrays*/
	if(calculate_cell_occupation)
	{
		for(i=0; i<cylinder.n_cells; i++)
		{
			cell_occupation[i][0] =0;
		}
	}

	
	for(i=0; i<n_part; i++)
	{
		cell_idx = pos2cell_idx(cylinder, pos[i], shift );// canonize is included here, and check of being outside of the cylinder too
		c_p[i] = cell_idx;
		if(calculate_cell_occupation)
		{
			cell_occupation[cell_idx][0]++;
			cell_occupation[cell_idx][cell_occupation[cell_idx][0]] = i;
		}
	}

	return; /* back to main */
}




/* cell_temp is output array of the local (collision cell) temperature */
/* To do SAM average, call with factor 1.0 */
/* This function calculates teh temperature in each cell over the whole simulation box */
/* Use for SAM averaging */
/* cell_temp is the output array. cell_occupation is input. Call encage first to generate it. */
inline void check_temperature(int ** cell_occupation, double ** vel ,double m, int n_cells, double * cell_temp)
{
	int ci,j,k;
	double cell_vel[3];
	
	for(ci=0; ci<n_cells; ci++)
	{
		cell_temp[ci] = 0.0;
	}
	
	for(ci=0; ci<n_cells; ci++)
	{
		if( cell_occupation[ci][0] == 0)
		{
			cell_temp[ci] = 0.0;  //-1.0  right choice?
			continue;
		}
		cell_vel[0] = 0.0;
		cell_vel[1] = 0.0;
		cell_vel[2] = 0.0;
		
		for(j=1; j<=cell_occupation[ci][0]; j++)
		{
			for(k=0; k<3; k++)
			{
				cell_vel[k] += vel[ cell_occupation[ci][j] ][k];
			}
		}

		for(k=0; k<3; k++)
		{
			cell_vel[k]/=(double)cell_occupation[ci][0];
		}
		
		for(j=1; j<=cell_occupation[ci][0]; j++)
		{
			for(k=0; k<3; k++)
			{
				cell_temp[ci] += ((vel[ cell_occupation[ci][j] ][k]-cell_vel[k]) * (vel[ cell_occupation[ci][j] ][k]-cell_vel[k]));
				//cell_temp[ci] += ((vel[ cell_occupation[ci][j] ][k] ) * (vel[ cell_occupation[ci][j] ][k]));

			}
		}
		cell_temp[ci]*= (m/((double)3*cell_occupation[ci][0]));
	}
	
	//cell_temp[ci] contains the cvalue of k_BT_cell for each cell ci
	
	
	return; /* return to main */
}









/* It exports the data to do the averages to calculate the 3Dvelocity profile at a given */
/* point x, at a particular time step */
/* The first row in the exported file contains, in this order: the x value at which the slice is taken, */
/* the number of cells in the slice (ny*nz), and the number of cells in each column in that slice (nz)  */
/* IDEA TO INCREASE PERFORMANCE: DO NOT ALLOCATE AND DEALOCATE MEMORY IN EVERY CALL, EITHER USE STATIC OR pass the arrays as arguments, allocate once in main*/
/* PING - check this after changing export_CAM_data */
inline void export_vel_profile(int n_part, double density, double ** vel, double ** pos, const Geometry cylinder, double x_slice, int file_number ,int step)
{
	/* Declare variables */
	char file_name[30]="./DATA/velprofile";
	int i,j;
	double ** slice_vel;
	int * local_density;
	int max_part_density = (int)density * DENSITY_TOL; // This is the maximum particle density in any given cell that is tolerated (WARNING: only controlled if the function check_incompressibility is called in the main simulation loop!)
	int idx;
	int nynz = cylinder.n_cells_dim[1] * cylinder.n_cells_dim[2];
	double nullshift[3]={0.0,0.0,0.0};
	
	/* Allocate memory and initialize */
	/* http://c-faq.com/aryptr/dynmuldimary.html */
	slice_vel = malloc(nynz * sizeof(double *)); //allocate array of pointers to rows
	if (slice_vel==NULL) {printf("Error allocating slice_vel in mpc.c\n"); exit(EXIT_FAILURE);}	
	slice_vel[0] = malloc(nynz * max_part_density * sizeof(double)); // allocate the whole memory of the 2D array to store the data, and initialize to 0
	if (slice_vel[0]==NULL) {printf("Error allocating slice_vel[0] in mpc.c\n"); exit(EXIT_FAILURE);}
	
	local_density = malloc(nynz * sizeof(int)); 
	if (local_density==NULL) {printf("Error allocating local_density in mpc.c\n"); exit(EXIT_FAILURE);}
	
	for(i = 0; i<nynz; i++)
	{
		slice_vel[i] = slice_vel[0] + i * max_part_density; // assign the pointers to the correct start of their rows
	}
	
	for(i=0; i<nynz; i++) 
	{
		local_density[i] = 0;
		
		for(j=0; j< max_part_density; j++) //Optimising: the rightmost index in the inner loop.
		{
			slice_vel[i][j] = 0.0;
		}
	}
	
	/* Find slice of interest  */
	int first_cell = slice_start_index(x_slice, nynz, cylinder.a);
	

	/* Find out where each particle is, and gather the data */
	/* This information could be stored in c_p, if supplied to this method. I do not consider it necessary, in any case a separate call would do*/
	/* For example, in case another property that needs c_p is going to be calculated after the velocity profile */
	int cell_idx = 0;
	for(i=0; i<n_part; i++) //loop over particles
	{
		/* Which cell is this particle in?*/
		cell_idx = pos2cell_idx(cylinder, pos[i], nullshift);
		if( cell_idx >= first_cell && cell_idx < first_cell+nynz)//if particle i is in a cell in the slice of interest...
		{
			idx = cell_idx - first_cell; //relative cell index
			local_density[idx]+=1; // update the counter
			if(local_density[idx]>max_part_density) //this prevents array overflow/out of bounds
			{
			  	  printf("export_vel_profile in mpc.c: Local density exceeded max tolerance (%d). Aborting...\n",max_part_density);
			  	  exit(EXIT_FAILURE);
			}
			slice_vel[ idx ][ local_density[idx] ] = vel[i][0];  //export x-velocity only
		}
	}
	

	/* now, export the data */
	export_CAM_data(slice_vel, local_density, file_name, cylinder, first_cell, first_cell + nynz - 1, file_number, step);
	
	
	/* Clean up and release memory*/
	free(slice_vel[0]);
	free(slice_vel);
	free(local_density);

}



////////////////////////////////////////////////////
/// Calculates the mean velocity of each cell at the time of being called (instantaneous mean velocity).
/// Input: cell_occupation, 2D array cointaining in row j the total number of particles found in that cell, followed by 
///        the indices of those particles. This array can be produced calling encage()
/// Output: cell_vel_output, n_cells X 3 array
////////////////////////////////////////////////////
inline void calculate_cell_velocity(int n_cells, int n_part, double ** vel, int ** cell_occupation, double ** cell_vel_output )
{
	int i,j,k;
	int local_density;
	
	
	/* Initialization */
	for(i=0; i<n_cells; i++)
	{
		for(k=0; k<3; k++)
		{
			cell_vel_output[i][k] = 0.0;
		}
	}
	
	
	for(i=0; i<n_cells; i++ )
	{
		local_density = cell_occupation[i][0];
		for(k=1; k<=local_density; k++)
		{
			for(j=0; j<3; j++)
			{
				cell_vel_output[i][j] += vel[ cell_occupation[i][k] ][ j ];
			}
		}
		for(j=0;j<3;j++)
		{
			if( local_density != 0 )
			{
				cell_vel_output[i][j] /= (double)local_density;
			}else{
				cell_vel_output[i][j]=0;
			}
		}
	}
	
	return; /* back to main */
}




/* Calculates the distribution of particles in the whole domain. The output vector densities has length=n_cells,
 and contains the number of particles found in each cell, for the given timestep.		*/
inline void density_distribution(int n_part, double m, double ** pos, Geometry geometry, double * shift, double * densities )
{
	int i;
	int cell_idx;
	double factor=1.0/geometry.a;
	
	/* Clean array*/
	for(i=0; i<geometry.n_cells; i++)
	{
		densities[i] = 0.0;
	}

	for(i=0; i<n_part; i++)
	{
		cell_idx = pos2cell_idx(geometry, pos[i], shift );// canonize is included here, and check of being outside of the cylinder too
		assert( cell_idx>= 0 && cell_idx <= geometry.n_cells);
		densities[cell_idx] += factor * m;
	}
	
	return; /* Back to main*/
}




/*-------------------------------*/
/*		MAIN		 */
/*-------------------------------*/

int main(int argc, char **argv) {
	
	int i=0,j=0,k=0; 		 	/* Generic loop counters */
	
	FILE * input_fp;
	input_fp = fopen(argv[1],"r");
	//input_fp = fopen("input.dat","r");
	if(input_fp == NULL){ fprintf(stderr,"Cannot open input file.\n Aborting...\n"); exit(EXIT_FAILURE);}
	
	
	char variable_name[15];
	double value=0.0;
	double variable_array[10];
	while (fscanf(input_fp, "%s\t%lf", variable_name, &value) == 2) {
		
		//printf("reading variable %s = %.3lf\n",variable_name, value);
		variable_array[i] = value;
		i++;
  	}

  	const double a = variable_array[0];		/// Size of a collision cell (cubic cell)
	const double Lx = variable_array[1];		/// Length of the simulation box in the direction of the flow //10
	const double L = variable_array[2];		/// Width and height of the simulation box //4
	const int n_part = variable_array[3];   	/// Number of particles in the simulation box//density=10
	const double T = variable_array[4];  		/// Temperature, in kT units
	const double m = variable_array[5]; 		/// Mass of a particle (plasma)
	const double dt = variable_array[6];
	const double g = variable_array[7];		/// "Gravity" (acceleration that drags the flow in the x-direction)
	const int steps = variable_array[8];	/// Number of steps in the simulation
	int equilibration_time = variable_array[9]; 
	
	
	int equilibration_export = 10;
	fclose(input_fp);
	/*------------------------------------------------------------------*/
	/*	SETUP							    */
	/* 	Definition of constants, initialization of variables,       */
	/*      memory allocation, I/O setup			            */
	/*								    */
	/*	 See description of units in Watari2007 		    */
	/*------------------------------------------------------------------*/
	/* Parameters, constants and variables */
	/*const double a = 1.0;		/// Size of a collision cell (cubic cell)
	const double Lx = 30.0;		/// Length of the simulation box in the direction of the flow //10
	const double L = 10.0;		/// Width and height of the simulation box //4
	const int n_part = 30000;   	/// Number of particles in the simulation box//density=10
	const double T = 1.0;  		/// Temperature, in kT units
	const double m = 1.0;*/ 		/// Mass of a particle (plasma)
	const double m_inv = 1.0/m;	/// Inverse mass of a particle (plasma)
	const double radius = 0.5*( L - a );	   	/// Radius of the cylindrical vessel. It is a/2 smaller than half the box side to allow for grid shifting 
	const int n_cells = (int)((Lx*L*L)/(a*a*a));  	/// Number of cells in the simulation box	
	//const double dt = 0.01;		/// Timestep
	const double density = (1.0*n_part)/(L*L*Lx);  	/// Number of particles per collision cell //const int
	double null_shift[3] = {0.0, 0.0, 0.0};
	
	/* TODO: check Whitmer's implementation of Verlet step and why he does not have a lambda */
	//const double lambda = 0.65;  /// Verlet's algorithm lambda parameter. 0.65 is a "magic number" for DPD, see page 5 of \cite{verlet}," Novel Methods in Soft Matter Simulations", ed. Karttunen et. al.
	/* consider setting lambda inside stream()? maybe as a define? */
	
	//const double g = 1e-4;	/// "Gravity" (acceleration that drags the flow in the x-direction)
	double t = 0.0;			/// Time
	/*const int steps = 10000;	/// Number of steps in the simulation
	int equilibration_time = 300; */
	
	#if CHECK_EQUILIBRATION
		equilibration_time = 0;
	#endif
	
	double verlet_epsilon = (1e-5)*radius; // i.e., tolerance zone thickness 2*epsilon = 2% of the radius.
	double verlet_min_sq = radius - verlet_epsilon;
	verlet_min_sq *=verlet_min_sq;
	double verlet_max_sq = radius + verlet_epsilon;
	verlet_max_sq *= verlet_max_sq;
	
	/* Struct initialization */
	const Geometry cylinder={a,Lx,L,0.5*L,0.5*( L - a ), verlet_epsilon, verlet_min_sq, verlet_max_sq, n_cells, { (int)(Lx/a), (int)(L/a), (int)(L/a)} }; 
	
	/* Particle quantities lists */
	double * pos_rmo;
	double * vel_rmo;
	double * acc_rmo;
	double * cell_vel_rmo;
	double * cell_rnd_vel_rmo;
	double ** pos;
	double ** vel;
	double ** acc;
	double ** cell_vel;
	double ** cell_rnd_vel;
	double * densities;
		
	/* Check parameter values are compatible */
	if( fmod(Lx,a) != (double)0.0 ){
		printf(" Lx must be a multiple of a. Exiting...\n");
		exit(EXIT_FAILURE);
	}else if( fmod(L,a) != (double)0.0 ){
		printf("L must be a multiple of a. Exiting...\n");
		exit(EXIT_FAILURE);
	}

	printf("--------------------------------\n");
	printf("Parameter \t Value\n");
	printf("--------------------------------\n");
	printf("density \t %.2lf\n",density);
	printf("n_part \t\t %d\n",n_part);
	printf("a \t\t %.3lf\n",a);
	printf("n_cells \t %d\n",n_cells);
	printf("Lx \t\t %.3lf\n",L);
	printf("L \t\t %.3lf\n",Lx);
	printf("radius \t\t %.3lf\n",radius);
	printf("m \t\t %.3lf\n",m);
	printf("T \t\t %.3lf\n",T);
	printf("grid \t\t %.0lfx%.0lfx%.0lf \n",L/a,L/a,Lx/a);
	printf("simulation box \t %.2lfx%.2lfx%.2lf\n",Lx,L,L);
	printf("Mean free path \t %.2lf\n",dt * sqrt(T*m_inv));
	printf("Kinetic viscosity /eta/ \t %.2lf\n", density * T * dt * ((density/(density-1+exp(-density)))-0.5)/(a*a*a) );
	printf("Collision viscosity /eta/\t %.2lf\n", m * (density-1+exp(-density))/(12*a*dt));
	printf("Total viscosity /eta/ \t\t %.2lf\n", density * T * dt * ((density/(density-1+exp(-density)))-0.5)/(a*a*a) + m * (density-1+exp(-density))/(12*a*dt));
	printf("Knudsen number Kn=%.3lf\n",dt * sqrt(T*m_inv)/a); // Kn>>1 is free molecular regime, Kn<<1 is hydrodynamic regime, Kn around 1 is gas rarefaction
	// export hydrodynamic numbers?
	printf("--------------------------------\n");
 

	/* Export collision grid data for visualization, and the vessel geometry data */
	#if EXPORT_STATES
		export_vtk_gid(cylinder);
		export_vessel_geometry(cylinder,steps);
	#endif
	
	/* Set up RNG */
	const gsl_rng_type * rngtype; 	 /* Type of rng generator   */
	gsl_rng * r;			 /* Instance of a generator */
    
	//gsl_rng_env_setup();		/* the default seed is used (0) */ //-------------------------------------->TODO: REMEMBER TO CHANGE THE SEED WHEN WORKING CODE IS RUNNING
	
	rngtype = gsl_rng_default; 	/* this is gsl_rng_mt19937 */
	r = gsl_rng_alloc(rngtype);	/* create pointer to instance of a RNG of type rngtype */
	
	long int seed = (long int)time(NULL);
	//seed = 1379012447 ;
	//seed = 1379520223 ;
	//seed = 1380619236; //Report 1/10/13
	//#if VERBOSE
	printf("RNG seed \t %ld\n",seed);
	//#endif
	gsl_rng_set(r , seed); //sets the seed for the RNG r

	/* Allocate memory for particles and set up initial values */
	pos_rmo = malloc( 3 * n_part * sizeof(double) );
	if (pos_rmo==NULL) {printf("Error allocating pos_rmo in mpc.c\n"); exit(EXIT_FAILURE);}
	vel_rmo = malloc( 3 * n_part * sizeof(double) );
	if (vel_rmo==NULL) {printf("Error allocating vel_rmo in mpc.c\n"); exit(EXIT_FAILURE);}
	acc_rmo = malloc( 3 * n_part * sizeof(double) );
	if (acc_rmo==NULL) {printf("Error allocating acc_rmo in mpc.c\n"); exit(EXIT_FAILURE);}

	cell_vel_rmo = malloc( 3 * n_cells * sizeof(double) );
	if (cell_vel_rmo==NULL) {printf("Error allocating cell_rnd_vel_rmo in mpc.c\n"); exit(EXIT_FAILURE);}
	cell_rnd_vel_rmo = malloc( 3 * n_cells * sizeof(double) );
	if (cell_rnd_vel_rmo==NULL) {printf("Error allocating cell_rnd_vel_rmo in mpc.c\n"); exit(EXIT_FAILURE);}
	
	pos = malloc( n_part * sizeof(double*) );
	if (pos==NULL) {printf("Error allocating pos in mpc.c\n"); exit(EXIT_FAILURE);}
	vel = malloc( n_part * sizeof(double*) );
	if (vel==NULL) {printf("Error allocating vel in mpc.c\n"); exit(EXIT_FAILURE);}
	acc = malloc( n_part * sizeof(double*) );
	if (acc==NULL) {printf("Error allocating acc in mpc.c\n"); exit(EXIT_FAILURE);}
	
	cell_vel = malloc( n_cells * sizeof(double*) );
	if (cell_vel==NULL) {printf("Error allocating cell_vel in mpc.c\n"); exit(EXIT_FAILURE);}
	cell_rnd_vel = malloc( n_cells * sizeof(double*) );
	if (cell_rnd_vel==NULL) {printf("Error allocating cell_rnd_vel in mpc.c\n"); exit(EXIT_FAILURE);}
	
	densities = malloc(n_cells * sizeof(double));
	if (densities==NULL){ printf("Error allocating densities in mpc.c\n"); exit(EXIT_FAILURE); }
	
	#if CHECK_MOMENTUM_CONSERVATION
		// auxiliary arrays for checking collide 
		double ** cell_vel_beforecollide;
		double * cell_vel_beforecollide_rmo;
		double ** cell_vel_aftercollide;
		double * cell_vel_aftercollide_rmo;
	
		cell_vel_beforecollide_rmo=malloc( 3 * n_cells * sizeof(double) );
		if (cell_vel_beforecollide_rmo==NULL) {printf("Error allocating cell_vel_beforecollide_rmo in mpc.c\n"); exit(EXIT_FAILURE);}
		cell_vel_beforecollide = malloc( n_cells * sizeof(double*) );
		if (cell_vel_beforecollide==NULL) {printf("Error allocating cell_vel_beforecollide in mpc.c\n"); exit(EXIT_FAILURE);}
		
		cell_vel_aftercollide_rmo=malloc( 3 * n_cells*sizeof(double) );
		if (cell_vel_aftercollide_rmo==NULL) {printf("Error allocating cell_vel_aftercollide_rmo in mpc.c\n"); exit(EXIT_FAILURE);}
		cell_vel_aftercollide = malloc( n_cells * sizeof(double*) );
		if (cell_vel_aftercollide==NULL) {printf("Error allocating cell_vel_aftercollide in mpc.c\n"); exit(EXIT_FAILURE);}
	#endif
	// If not checking the temperature this is not needed: create a dummy pointer for encage() and carry on
	int ** cell_occupation;
	int * cell_occupation_rmo;
	int max_oc = (int)(density * DENSITY_TOL);
	cell_occupation_rmo=malloc(max_oc*n_cells*sizeof(int));
	if (cell_occupation_rmo==NULL) {printf("Error allocating cell_occupation_rmo in mpc.c\n"); exit(EXIT_FAILURE);}
	cell_occupation = malloc(n_cells*sizeof(int*));
	if (cell_occupation==NULL) {printf("Error allocating cell_occupation in mpc.c\n"); exit(EXIT_FAILURE);}
	// -------------------------------------------------------

	for(i=0; i<n_cells; i++)
	{
		cell_vel[i] = &cell_vel_rmo[3*i];
		cell_rnd_vel[i] = &cell_rnd_vel_rmo[3*i];
		/* Initialization is done at the beginning of collide() */
		cell_occupation[i] = &cell_occupation_rmo[max_oc*i];
		
		#if CHECK_MOMENTUM_CONSERVATION
			cell_vel_beforecollide[i]=&cell_vel_beforecollide_rmo[3*i];
			cell_vel_aftercollide[i]=&cell_vel_aftercollide_rmo[3*i];
		#endif
	}
		
	int * c_p;				/// c_p[i] is the cell index of particle i 
	c_p = malloc( n_part*sizeof(int) );     /// cell idx from particle idx
	if (c_p==NULL) {printf("Error allocating c_p in mpc.c\n"); exit(EXIT_FAILURE);}

    
	/* cell_mass[ci] is the total mass of the cell ci */
	double * cell_mass;
	cell_mass = malloc( n_cells*sizeof(double) ); // mass of particles in cell
	if (cell_mass==NULL) {printf("Error allocating cell_mass in mpc.c\n"); exit(EXIT_FAILURE);}

	
	double * shift;			/// random grid shift vector (Galilean invariance)
	shift = malloc(3*sizeof(double));
	if(shift == NULL) {printf("Error allocating shift in mpc.c\n"); exit(EXIT_FAILURE);}
	

	/*----------------------------------------------------*/
		
	/* Initialise */
	/* Position is uniformly distributed in space within vessel: geometry dependent */
	/* Velocity is normally distributed (Maxwell-Boltzmann) 			*/
	/* Note: Components of the velocity are chosen uncorrelated   			*/
	for(i=0; i<n_part; i++)
	{
		/* Complete memory allocation */
		pos[i] 	 = &pos_rmo[ 3*i ];
		
		/* Initialise */
		/* For a cylindrical geometry, random points are chosen uniformly within the circle. */
		populate( r, cylinder, &(pos[i][0]) );
		
		/* Complete memory allocation */
		vel[i] = &vel_rmo[ 3*i ];
		acc[i] = &acc_rmo[ 3*i ];
		
		/* Initialise */
		for(j=0; j<3; j++)
		{
			vel[i][j] = gsl_ran_gaussian_ziggurat(r, sqrt(T * m_inv) ); 	/* mean zero */ //std deviation= sqrt(k_B T /m)
			acc[i][j] = 0.0; 
		}
		acc[i][0] = g; 
		
		/*-----------------------------------------------*/
		// compute_acc(g,&(acc[i][0])); //This line is included in the previous two lines, so it is redundant.
		//it would, however, be necessary if there were forces acting on the plasma particles
		 
		 /* These lines were used in the testing of the fractioning of the timestep in streaming step */
		//acc[i][0] = g;// -((double)rand()/(RAND_MAX));
		//acc[i][1] = 40.0*(1-2*((double)rand()/(RAND_MAX)));
		//acc[i][2]=  40.0*(1-2*((double)rand()/(RAND_MAX)));
		/*-----------------------------------------------*/
		
	}	
	
	for( i=0; i<n_cells; i++)
	{
		/* Complete memory allocation */
		cell_vel[i] = &cell_vel_rmo[ 3*i ];
		
		/* Initialise */
		cell_vel[i][0] =0.0;
		cell_vel[i][1] =0.0;
		cell_vel[i][2] =0.0;
	}

	
	
	/*------------------------------------------------------------------*/
	/*	SIMULATION						    */
	/*------------------------------------------------------------------*/
	
	/* Export initial data */
	#if EXPORT_STATES
		export_vtk_plasma(0, pos, vel,acc, n_part);
	#endif
	
	/* Export velocity profile setup*/
	int file_counter = 0;
	
	/* Choose x-position at which the velocity profile will be exported (see implementation in collide.c)*/
	double x_slice = 15; // check that it is less than Lx (periodic conditions probably compensate anyway...)
	
	if( x_slice >= Lx || x_slice <=0 )
	{
		fprintf(stderr,"Slice of interest is outside cylinder. Exiting...\n");
		exit(EXIT_FAILURE);
	}
	
	int counter = 0;
	
	#if CHECK_EQUILIBRATION
		/* Variables to monitor equilibration*/
		double momentumPpart[6] = {0.0,0.0,0.0,0.0,0.0,0.0}; //linear momentum per particle, at any given timestep, followed by the variances
		double energyPpart = 0.0;
		double systemTemp = 0.0;
		int equilibration_counter=0;
		FILE * equilibration_fp;
 		//fprintf(equilibration_fp,"Timestep \t p_x \t p_y \t p_z \t e_k \n");
 		equilibration_fp=fopen("equilibration.dat","w");
 		double total_energy[2]={0.0,0.0};
 		
 		fprintf(equilibration_fp,"%s %6s %8s %8s %8s %8.s %8s %8s %8s %8s\n","t","px","sigma(px)","py","sigma(py)","pz","sigma(pz)","ek","sigma(ek)","T");
 	#endif
	
		
	
	
	#if DEBUGGING_STREAMCOLLIDE
		//FILE * debug_fp;
		debug_fp=fopen("debug_data.dat","w");
	#endif
	#if CHECK_TEMPERATURE
		int file_temp_counter = 0;
	#endif
	
	
	#if GALILEAN_SHIFT == 0
		shift[0] = 0.0;
		shift[1] = 0.0;
		shift[2] = 0.0;
	#endif
	
	/* Loop for a fixed number of times */
	for(i=1; i<= steps; i++)
	{
		#if DEBUGGING_STREAMCOLLIDE
			export_data(debug_fp, pos , n_part );
			export_data(debug_fp, vel , n_part );
			export_data(debug_fp, acc , n_part );	
		#endif
		/* Streaming step */
		stream( dt, n_part, g, cylinder, pos, vel, acc);
		 
		/* debugging */
		if(TEST_all_particles_in_lumen(n_part, cylinder, pos)!=1)
		{
			fprintf(stderr,"mpc.c: particles escaped the lumen after stream step in timestep %d\n",i);
			exit(EXIT_FAILURE);
		}

		/* Update timestep */
		t += dt; 

		#if GALILEAN_SHIFT
			shift[0] = a * ( -0.5 + gsl_rng_uniform_pos(r) ) ;
			shift[1] = a * ( -0.5 + gsl_rng_uniform_pos(r) );
			shift[2] = a * ( -0.5 + gsl_rng_uniform_pos(r) );
		#endif
		
		/* Collision step */
		
		
		#if CHECK_MOMENTUM_CONSERVATION
			encage(pos, shift, n_part, cylinder, c_p, 1 , cell_occupation);
			calculate_cell_velocity(n_cells, n_part, vel, cell_occupation, cell_vel_beforecollide );
		#endif	
		collide(n_part, T, m, m_inv, cylinder, r, c_p, cell_mass, cell_vel, cell_rnd_vel, shift, pos, vel); // BEWARE! the velocity of cells after this is not to be trusted, due to the Galilean shift: it is necessary to undo it!
		// c_p contains the occupation information, but with the Galilean shift: if we wanted to use it, we would need to undo it.		
		
		#if CHECK_MOMENTUM_CONSERVATION
			calculate_cell_velocity(n_cells, n_part, vel, cell_occupation, cell_vel_aftercollide );

			for(j=0;j<n_cells;j++)
			{
				for(k=0;k<3;k++)
				{
					assert( is_approx_zero(cell_vel_aftercollide[j][k] - cell_vel_beforecollide[j][k],1e-6 ) );
					/*if( is_approx_zero(cell_vel_aftercollide[j][k] - cell_vel_beforecollide[j][k],1e-5 )  )
					{ //do nothing 
					}else{
						printf("Not zero: %lf\n Aborting...\n",cell_vel_aftercollide[j][k] - cell_vel_beforecollide[j][k]);
						exit(EXIT_FAILURE);
					}*/
				}
			}
		#endif
		
		
		
		#if DEBUGGING_STREAMCOLLIDE
			export_data(debug_fp, pos , n_part );
			export_data(debug_fp, vel , n_part );
			export_data(debug_fp, acc , n_part );	
		#endif
		
		#if CHECK_COMPRESSIBILITY
			check_compressibility(cell_mass, n_cells, m_inv, density);
		#endif

		/* Export particle data for a 3D paraview visualisation */
		#if EXPORT_STATES
			export_vtk_plasma(i, pos, vel,acc, n_part);		
		#endif
		
	
		/* Export velocity data every 50 timesteps, after the first equilibration_time timesteps*/
		#if EXPORT_VEL_PROFILE
			counter++;
			if( counter==equilibration_time + 50 )
			{
				//fprintf(stderr,"Exporting velocity profile data at step number %d\n", i);
				file_counter++;
				//DEPRECATED: export_vel_profile(n_part, density, vel, pos, cylinder, x_slice, file_counter, i); 
				export_vel_profile(n_part, density, vel, pos, cylinder, x_slice, file_counter, i);
				counter = equilibration_time;
			}
		#endif
		
		

		#if 0
			equilibration_counter++;
			if(equilibration_counter>20) // export data on total momentum per particle every 10 timesteps
			{
				equilibration_counter = 0;
				
				#if CHECK_TEMPERATURE
					file_temp_counter++;
					encage(pos, null_shift, n_part, cylinder, c_p, 1, cell_occupation);
					check_temperature(cell_occupation, vel, m, n_cells, cell_mass);//recycling vector cell_mass as cell_temp. It is reset to 0 at the start of collide()
					//export_temp_profile(cell_mass, cylinder, file_temp_counter, i);
					/* the following line exports exactly the same as above - TESTED OK */
					export_SAM_data(cell_mass, "./DATA/temperatureSAM", cylinder, 0, cylinder.n_cells-1, file_temp_counter, i);
				#endif
				
			}
		#endif
		
		#if CHECK_EQUILIBRATION
			equilibration_counter++;
			if( equilibration_counter >= equilibration_export ){
			//total_momentum(n_part, m, vel, momentumPpart);
			//total_kinetic_energy(n_part, m, vel, total_energy); 
			//systemTemp = (2.0*energyPpart - m_inv*norm_sq(momentumPpart) )/3.0;
			//systemTemp = (2*total_energy[0]-
			//	(momentumPpart[0]*momentumPpart[0]+momentumPpart[1]*momentumPpart[1]+momentumPpart[2]*momentumPpart[2])*m_inv)/3.0;
			/*fprintf(equilibration_fp,"%d %8.4lf  %8.4lf  %8.4lf  %8.4lf  %8.4lf  %8.4lf  %8.4lf  %8.4lf  %8.4lf\n", 
				i, momentumPpart[0], sqrt(momentumPpart[3]), momentumPpart[1], sqrt(momentumPpart[4]), momentumPpart[2], 
				sqrt(momentumPpart[5]), total_energy[0], sqrt(total_energy[1]), systemTemp);*/
			
			/* Density distribution */	
			//density_distribution(n_part, m, pos, cylinder, null_shift, densities );
			//export_SAM_data(densities, "./DATA/density", cylinder, 0, cylinder.n_cells-1, i, i);
			/* Velocity profile */
			/*encage(pos, shift, n_part, cylinder, c_p, 1 , cell_occupation);
			calculate_cell_velocity(n_cells, n_part, vel, cell_occupation, cell_vel );
			export_data(cell_vel,"./DATA/vel_equi",cylinder,0,cylinder.n_cells-1,i,i);*/
			equilibration_counter=0;
			}
			
		#endif
		
		

	}/* end for loop (stream-collide steps) */


	

	
	/*------------------------------------------------------------------*/
	/*	EXITING							    */
	/*------------------------------------------------------------------*/
	/* Clean up */
	
	#if DEBUGGING_STREAMCOLLIDE
		fclose(debug_fp);//global file pointer
	#endif
	#if CHECK_EQUILIBRATION		
 		fclose(equilibration_fp);
 	#endif
	
	free(pos_rmo);
	free(vel_rmo);
	free(acc_rmo);
	free(cell_vel_rmo);
	free(cell_rnd_vel_rmo);
	free(pos);
	free(vel);
	free(acc);
	free(cell_vel);
	free(cell_rnd_vel);
	free(c_p);
	free(cell_mass);
	free(shift);
	free(cell_occupation_rmo);
	free(cell_occupation);
	free(densities);
	#if CHECK_MOMENTUM_CONSERVATION
		free(cell_vel_beforecollide_rmo);
		free(cell_vel_beforecollide);
		free(cell_vel_aftercollide_rmo);
		free(cell_vel_aftercollide);
	#endif

	gsl_rng_free(r); /* free memory associated with the rng r */

	return(EXIT_SUCCESS);
}





