/** MPC simulations.
 *
 * Author: Leonor Garcia-Gutierrez, 2014
 *
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
#include <time.h>
#include <assert.h>
#include <string.h>

#include "mpc.h"    /*  Struct definitions, function prototypes and global variables that must be made available to other files  */
#include "export_routines.h"

/* GNU GSL*/
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
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

/* DO NOT DELETE !!*/
/* Notice this can be called without undoing the shift after collide */
/* NOT TESTED : THIS ALSO NEEDS TO INCLUDE THE VOLUME INFORMATION, IT NEEDS TO COMPARE DENSITIES, NOT PARTICLE OCCUPATIONS */
/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Control routine that checks for compressibility effects: it stops the simulation if, at the time of being invoked, the density in any 
/// collision cell is greater than density+X or smaller than density-X with X=DENSITY_TOL*density.
/////////////////////////////////////////////////////////////////////////////////////////////////////
/*inline void check_compressibility(double * cell_mass, int n_cells, double m_inv, double density)
{
	int ci;
	double tolerance = DENSITY_TOL*density;
	for(ci=0; ci< n_cells; ci++)
	{
		if(m_inv*cell_mass[ci] > (density+tolerance) || m_inv*cell_mass[ci] < (density-tolerance) )
		{
			fprintf(stderr,"Compressibility effects exceeded tolerance: *shifted* cell %d has instantaneous density of %.2lf\n",ci,cell_mass[ci]*m_inv);
			exit(EXIT_FAILURE);
		}
	}
	
	return; /* Back to main*
}
*/

/* TODO: at the moment the export is done in ASCII only, binary option has not been implemented yet*/
/* possible to do by keeping the lines in ASCII and writing the data (only) in binary to the same file stream*/
//////////////////////////////////////////////////////////////////////////
/// It writes a VTK file with positions and velocities of all particles 
//////////////////////////////////////////////////////////////////////////
inline void export_vtk_plasma(int id, double ** pos, double ** vel,double ** acc, int n_part)
{
	int i, flag;
	
	FILE * fp;
	char pos_name[30];
	//snprintf(pos_name,sizeof(char)*30,"plasma%04d.vtk",id); //leading zeroes: not pvpython-friendly
	snprintf(pos_name,sizeof(char)*30,"./../experiments/plasma%d.vtk",id);

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
	
	flag = fclose(fp);
	if( flag !=0 )
	{
		fprintf(stderr,"export_vtk_plasma: error closing file\n Aborting...\n");
		exit(EXIT_FAILURE);
	}
	
	return; /* back to main */
}



/* No need to do in in binary, unless the geometry becomes much more complicated (i.e. some triangulated surface)*/
/* NOTE that if the directory DATA is not present, the program crashes with a Segmentation Fault */
//////////////////////////////////////////////////////////////////////////
/// This exports the VTK file defining the grid used (not the vessel) 
/// This is called only once at the start of the simulation. It exports in ASCII.
//////////////////////////////////////////////////////////////////////////
inline void export_vtk_grid(Geometry cylinder)
{
	int flag;
	FILE * fp;
	fp = fopen("./../experiments/collision_grid.vtk","w");
	if( fp == NULL )
	{
		printf("export_vtk_grid: Error opening ./../experiments/collision_grid.vtk. Check directory exists. Aborting...\n");
		exit(EXIT_FAILURE);
	}
	fprintf(fp,"# vtk DataFile Version 2.0\nGrid representation\nASCII \nDATASET STRUCTURED_POINTS\nDIMENSIONS %d %d %d\nORIGIN 0 0 0\nSPACING %.1lf %.1lf %.1lf",cylinder.n_cells_dim[0]+1,cylinder.n_cells_dim[1]+1,cylinder.n_cells_dim[2]+1,cylinder.a,cylinder.a,cylinder.a);
	flag = fclose(fp);	
	if( flag !=0 )
	{
		fprintf(stderr,"export_vtk_grid: error closing file\n Aborting...\n");
		exit(EXIT_FAILURE);
	}
}

//////////////////////////////////////////////////////////////////////////
/// This function exports the data necessary for the python script to generate a cylinder of the right dimensions on paraview.
//////////////////////////////////////////////////////////////////////////					
inline void export_vessel_geometry(Geometry cylinder, int num_steps)
{
	int flag;
	FILE * fp;
	fp = fopen("./../experiments/vessel_geometry.py","w");
	if( fp == NULL )
	{
		printf("export_vessel_geometry: Error opening ./../experiments/vessel_geometry.vtk. Check directory exists. Aborting...\n");
		exit(EXIT_FAILURE);
	}
	fprintf(fp,"# Cylinder data\n");
	fprintf(fp,"Lx=%lf\n",cylinder.Lx);
	fprintf(fp,"L=%lf\n",cylinder.L);
	//fprintf(fp,"a=%lf\n",cylinder.a);
	fprintf(fp,"radius=%lf\n",cylinder.radius);
	fprintf(fp,"L_half=%lf\n",cylinder.L_half);
	fprintf(fp,"num_frames=%d\n",num_steps);
	flag = fclose(fp);
	if( flag !=0 )
	{
		fprintf(stderr,"export_vessel_geometry: error closing file\n Aborting...\n");
		exit(EXIT_FAILURE);
	}
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






////////////////////////////////////////////////////////////////////////////////////////////
/// Maps the particle positions to the cells they are at, after the grid has been shifted by a vector shift
/// Input:
/// pos - Position vector containing the position vectors of all the (plasma) particles in teh simulation box
/// shift - 3D shift vector. 
/// n_part - Total number of (plasma) particles in the simulation box
/// cylinder - Geometry object (struct)
/// calculate_cell_occupation -  Boolean to switch between calculating the cell_occupation array (1) or not (0)
/// Output:
/// c_p - Vector of length n_part containing, for each particle (row), the global index of the cell the particle is at
/// cell_occupation - 2D array in which each row represents a cell. Within each row, the first element is the local denisty. Teh following numbers are the indices of
/// the particles at that cell. For example:
/// 4 \t 120 \t 146 \t 135 \t 27
/// 2 \t 189 \t 100
/// 0
/// 3 \t 110 \t 199 \t 55
/// ....
////////////////////////////////////////////////////////////////////////////////////////////		
inline void encage(double ** pos, double * shift, int n_part, Geometry cylinder, int * c_p, int calculate_cell_occupation , int ** cell_occupation, int max_oc)
{
	int i;
	int cell_idx;
	
	/* Clean arrays */
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
			
			if( cell_occupation[cell_idx][0] > max_oc )
			{
				fprintf(stderr,"Cell occupation for cell %d > density*DENSITY_MAX.",cell_idx);
				fprintf(stderr,"Increase DENSITY_MAX (now = %d), or check compressibility conditions. Aborting...\n",DENSITY_MAX);
				exit(EXIT_FAILURE);
			}
			
			
			cell_occupation[cell_idx][cell_occupation[cell_idx][0]] = i;
		}
	}

	return; /* back to main */
}


//////////////////////////////////////////////////////////////////////////
/// Calculates the total linear momentum of the system (only plasma particles).
/// Input:
/// n_part - Total number of particles in the simulation box
/// m - Mass of plasma particles
/// vel - Velocity vector of all the particles (n_part x 3)
/// Output:
/// momentum_output - 3D vector with total momentum of the simulation box (sum of the momentum of all particles)
//////////////////////////////////////////////////////////////////////////
inline void total_momentum(int n_part, double m, double ** vel, double * momentum_output)
{
	/* Define variables */
	int i,j;
	
	/* Initialise output vector */
	for(j=0; j<3; j++)
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
	
	
	return; /* Back to main */
}


//////////////////////////////////////////////////////////////////////////
/// Calculates the total linear momentum of the cell (only plasma particles).
/// Input:
/// m - Mass of plasma particles
/// vel - Velocity vector of all the particles (n_part x 3)
/// cell_particles - 1D array containing the indices of the particles currently at this cell, preceded by the local particle density (length of cell_particles)
/// For example: cell_particles = {4, 102, 150, 164, 98}. This can be cell_occupation[ci], with cell_occupation being the array generated by encage()
/// Output:
/// momentum_output - 3D vector with total momentum of the cell (sum of the momentum of the particles present in the cell)
//////////////////////////////////////////////////////////////////////////
inline void cell_momentum(int * cell_particles , double m, double ** vel, double * momentum_output)
{
	/* Define variables */
	int i,j;
	
	/* Initialise output vector */
	for(j=0; j<3; j++)
	{
		momentum_output[j]=0.0;
	}
	
	/* Add momentums */
	for( i=1; i<=cell_particles[0]; i++) //loop over particles
	{
		for(j=0; j<3; j++) // loop over cartesian components of the velocity
		{			
			momentum_output[j] += m*vel[ cell_particles[i] ][j];
		}
	}
	
	return; /* Back to main */
}


//////////////////////////////////////////////////////////////////////////
/// Calculates the total kinetic energy of the simulation box (plasma particles only)
/// Input:
/// n_part - Total number of (plasma) particles in the simulation box
/// m - Mass of the particles
/// vel - Velocity array for the system
///Output:
/// Total kinetic energy of the system of n_part particles
//////////////////////////////////////////////////////////////////////////
inline double total_kinetic_energy(int n_part, double m, double ** vel)
{
	double energy = 0.0;
	int i;
	
	for(i=0; i<n_part; i++)
	{
		energy += m * norm_sq(vel[i]);
	}

	
	return 0.5*energy;
}


//////////////////////////////////////////////////////////////////////////
/// Calculates the total kinetic energy of a given cell (plasma particles only)
/// Input:
/// m - Mass of the particles
/// vel - Velocity array for the system
/// cell_particles - 1D array containing the indices of the particles currently at this cell, preceded by the local particle density (length of cell_particles)
/// For example: cell_particles = {4, 102, 150, 164, 98}. This can be cell_occupation[ci], with cell_occupation being the array generated by encage()
///Output:
/// Total kinetic energy of the system of n_part particles
//////////////////////////////////////////////////////////////////////////
inline double cell_kinetic_energy(int * cell_particles, double m, double ** vel)
{
	double energy = 0.0;
	int i;
	
	for(i=1; i<=cell_particles[0]; i++)
	{
		energy += m * norm_sq(vel[i]);
	}

	
	return 0.5*energy;
}


// DO NOT DELETE!
/*  This function assumes the cell_volume array in the struct Geometry has been defined.
It is also necessary to implement a function to calculate these volumes (algorithm already done)
//////////////////////////////////////////////////////////////////////////
/// Calculates the instantaneous density distribution in the simulation box
/// Input:
/// cell_occupation - 2D array generated by encage() 
/// n_cells - Total number of collision cells in teh simulation box
/// Output:
/// densities - 1D array containing the 
//////////////////////////////////////////////////////////////////////////
inline void density_distribution(int ** cell_occupation, double m, Geometry cylinder, int * densities)
{
	int i;
	
	for(i=0; i<cylinder.n_cells; i++)
	{
		densities[i] = m*cell_occupation[i][0]/cylinder.cell_volume[i];
	}
	
	return;/* Back to main*
}
*/



//////////////////////////////////////////////////////////////////////////
/// This function builds the CAM array of velocities for a particular slice (at position float x), at a particular time
/// Input:
/// x - x position of the slice that is being considered.
/// cylinder - Geometry struct containing the boundary data
/// cell_occupation - 2D array generated by encage()
/// vel - 2D array containing the velocity vectors of all the particles in the system
/// only_x - Export only the x-component of the velocities if set to 1, export the three components if 0. 
/// IMPORTANT: This option affects the output vector to be called. If set to 1, pass slice_CAM_scalar, if set to 0 pass slice_CAM_vector
/// Output:
/// slice_CAM_output - Output array of size (ny*nz) x (3*max_density) if only_x=0, or (ny*nz)x(max_density), in the CAM format, 
/// with the indices of the cells being referred to the local slice. There is no header. For example: 2 v1x v1y v1z v2x v2y v2z
//////////////////////////////////////////////////////////////////////////
inline void slice_CAM_velocities(double x, Geometry cylinder, int ** cell_occupation, double ** vel, int only_x, double ** slice_CAM_output)
{
	int i,j,k;
	int nynz = cylinder.n_cells_dim[1] * cylinder.n_cells_dim[2];//I know this in advance because I am passing the output vector
	int first_cell = slice_start_index(x, nynz, cylinder.a);
	int last_cell =  first_cell + nynz - 1;
	
	assert( only_x == 0 || only_x == 1 );
	
	int counter = 0;
	int particle_idx;
	int local_density;
	
	for(i=first_cell; i<=last_cell; i++)
	{
		local_density =  cell_occupation[i][0];
		slice_CAM_output[counter][0] = 1.0*local_density;
		
		for(j=1; j<=local_density; j++)
		{
			particle_idx = cell_occupation[i][j];
			
			if( only_x == 0 )
			{
				for(k=0; k<3; k++)
				{
					slice_CAM_output[counter][3*j-2+k] = vel[ particle_idx ][k];//3*(j-1)+k+1
				}
			}else{
				slice_CAM_output[counter][j] = vel[ particle_idx ][0];
			}
		}
		
		counter++;
		
	}
	assert( counter == nynz );
	
	return; /* Back to main */
}





/*
/* It exports the data to do the averages to calculate the 3Dvelocity profile at a given *
/* point x, at a particular time step */
/* The first row in the exported file contains, in this order: the x value at which the slice is taken, *
/* the number of cells in the slice (ny*nz), and the number of cells in each column in that slice (nz)  *
/* IDEA TO INCREASE PERFORMANCE: DO NOT ALLOCATE AND DEALOCATE MEMORY IN EVERY CALL, EITHER USE STATIC OR pass the arrays as arguments, allocate once in main*/
/* PING - check this after changing export_CAM_data *
inline void export_vel_profile(int n_part, double density, double ** vel, double ** pos, const Geometry cylinder, double x_slice, int file_number ,int step)
{
	/* Declare variables *
	char file_name[30]="./DATA/velprofile";
	int i,j;
	double ** slice_vel;
	int * local_density;
	int max_part_density = (int)density * DENSITY_TOL; // This is the maximum particle density in any given cell that is tolerated (WARNING: only controlled if the function check_incompressibility is called in the main simulation loop!)
	int idx;
	int nynz = cylinder.n_cells_dim[1] * cylinder.n_cells_dim[2];
	double nullshift[3]={0.0,0.0,0.0};
	
	/* Allocate memory and initialize *
	/* http://c-faq.com/aryptr/dynmuldimary.html *
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
	
	/* Find slice of interest  *
	int first_cell = slice_start_index(x_slice, nynz, cylinder.a);
	

	/* Find out where each particle is, and gather the data */
	/* This information could be stored in c_p, if supplied to this method. I do not consider it necessary, in any case a separate call would do*/
	/* For example, in case another property that needs c_p is going to be calculated after the velocity profile *
	int cell_idx = 0;
	for(i=0; i<n_part; i++) //loop over particles
	{
		/* Which cell is this particle in?*
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
	

	/* now, export the data *
	export_CAM_data(slice_vel, local_density, file_name, cylinder, first_cell, first_cell + nynz - 1, file_number, step);
	
	
	/* Clean up and release memory*
	free(slice_vel[0]);
	free(slice_vel);
	free(local_density);

}
*/


//AQUI


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





/*-------------------------------*/
/*		MAIN		 */
/*-------------------------------*/

int main(int argc, char **argv) {
	
	int i=0,j=0,k=0; 		 	/* Generic counters */
	int fclose_flag;			/* Flag to check file close */
	/*------------------------------------------------------------------*/
	/*	SETUP: Import						    */
	/* 	Reading variables from input data file.			    */
	/*								    */
	/*	 See description of units in Watari2007 		    */
	/*------------------------------------------------------------------*/	
	
	// Include macros in input data?
	
	FILE * input_fp;
	input_fp = fopen(argv[1],"r");
	if(input_fp == NULL){ fprintf(stderr,"Cannot open input file.\n Aborting...\n"); exit(EXIT_FAILURE);}
	
	
	char variable_name[15];
	double value=0.0;
	double variable_array[10];
	while (fscanf(input_fp, "%s\t%lf", variable_name, &value) == 2) {
		variable_array[i] = value;
		i++;
  	}

  	const double a = variable_array[0];		/// Size of a collision cell (cubic cell)
	const double Lx = variable_array[1];		/// Length of the simulation box in the direction of the flow //10
	const double L = variable_array[2];		/// Width and height of the simulation box //4
	const int n_part = variable_array[3];   	/// Number of particles in the simulation box//density=10
	const double T = variable_array[4];  		/// Temperature, in kT units
	const double m = variable_array[5]; 		/// Mass of a particle (plasma)
	const double dt = variable_array[6];		/// Timestep
	const double g = variable_array[7];		/// "Gravity" (acceleration that drags the flow in the x-direction)
	const int steps = variable_array[8];		/// Number of steps in the simulation
	int equilibration_time = variable_array[9]; 	/// Number of timesteps to run before the production phase starts
	
	fclose_flag = fclose(input_fp);
	if( fclose_flag !=0 )
	{
		fprintf(stderr,"Error closing input file\n Aborting...\n");
		exit(EXIT_FAILURE);
	}
	
	
	/*------------------------------------------------------------------*/
	/*	SETUP: Complete						    */
	/* 	Definition of variables derived from the input data file,   */
	/*      and initialization of other independent variables.  	    */
	/*------------------------------------------------------------------*/
	const double m_inv = 1.0/m;			/// Inverse mass of a particle (plasma)
	const double radius = 0.5*( L - a );	   	/// Radius of the cylindrical vessel. It is a/2 smaller than half the box side to allow for grid shifting 
	const int n_cells = (int)((Lx*L*L)/(a*a*a));  	/// Number of cells in the simulation box	
	const double density = (1.0*n_part)/(L*L*Lx);  	/// Average number of particles per collision cell: note bulk cells have a greater density, and boundary cells have a smaller density than this
	double null_shift[3] = {0.0, 0.0, 0.0};		/// Null 3D vector
	double shift[3] = {0.0,0.0,0.0};			/// random grid shift vector (Galilean invariance)
	
	/* TODO: check Whitmer's implementation of Verlet step and why he does not have a lambda */
	//const double lambda = 0.65;  /// Verlet's algorithm lambda parameter. 0.65 is a "magic number" for DPD, see page 5 of \cite{verlet}," Novel Methods in Soft Matter Simulations", ed. Karttunen et. al.
	/* consider setting lambda inside stream()? maybe as a define? */
	
	/* Epsilon at the boundary */
	double verlet_epsilon = (1e-5)*radius; // i.e., tolerance zone thickness 2*epsilon = 2% of the radius.
	double verlet_min_sq = radius - verlet_epsilon;
	verlet_min_sq *=verlet_min_sq;
	double verlet_max_sq = radius + verlet_epsilon;
	verlet_max_sq *= verlet_max_sq;
	
	/* Struct initialization */
	const Geometry cylinder={a,Lx,L,0.5*L,0.5*( L - a ), verlet_epsilon, verlet_min_sq, verlet_max_sq, n_cells, { (int)(Lx/a), (int)(L/a), (int)(L/a)} }; 
	
	double t = 0.0;					/// Time
	int equilibration_export = 10;			/// Number of timesteps in between data exports during the equilibration phase
	/* data_header contains this info: file_number nx ny nz timestep x_slice_1stCell 
	*  Being the x_slice1stCell the global index of the first cell in that slice (equivalent to specifying x_slice -double- )
	*/
	int data_header_size = 7;
	int * data_header;
	data_header=malloc(data_header_size * sizeof(int));
	if(data_header==NULL){printf("mpc.c: Error allocating memory for data_header\nAborting...\n");exit(EXIT_FAILURE);}
	data_header[0] = 0;				// file number
	data_header[1] = cylinder.n_cells_dim[0]; 	// nx
	data_header[2] = cylinder.n_cells_dim[1];	// ny
	data_header[3] = cylinder.n_cells_dim[2];	// nz
	data_header[4] = -1;				// GLOBAL index first cell included in datafile (first cell in slice, if data is a slice)
	data_header[5] = -1;				// GLOBAL index last cell included in datafile
	data_header[6] = -1;				// timestep
	
	
	double aux3Dvector[3] = {0.0, 0.0, 0.0};

	
	/*------------------------------------------------------------------*/
	/*	SETUP: Check						    */
	/*------------------------------------------------------------------*/
	/* Check parameter values are compatible */
	if( fmod(Lx,a) != (double)0.0 ){
		printf(" Lx must be a multiple of a. Exiting...\n");
		exit(EXIT_FAILURE);
	}else if( fmod(L,a) != (double)0.0 ){
		printf("L must be a multiple of a. Exiting...\n");
		exit(EXIT_FAILURE);
	}

	/* Timestamp log file */
	time_t now;
	time(&now);
	printf("Start: %s\n",ctime(&now));
	
	/* Echo parameter values to log file */
	printf("========================================\n");
	printf("%-15s \t %-15s\n","Parameter","Value");
	printf("----------------------------------------\n");
	printf("%-15s \t %.4lf\n","density",density);
	printf("%-15s \t %d\n","n_part",n_part);
	printf("%-15s \t %.4lf\n","a",a);
	printf("%-15s \t %d\n","n_cells",n_cells);
	printf("%-15s \t %.4lf\n","Ly,Lz",L);
	printf("%-15s \t %.4lf\n","Lx",Lx);
	printf("%-15s \t %.4lf\n","radius",radius);
	printf("%-15s \t %.4lf\n","m",m);
	printf("%-15s \t %.4lf\n","T",T);
	printf("%-15s \t %.4lf\n","dt",dt);
	printf("%-15s \t %.4lf\n","g",g);
	printf("%-15s \t %d\n","simulation steps",steps);
	printf("%-15s \t %d\n","equilibration steps",equilibration_time);
	printf("%-15s \t %.0lfx%.0lfx%.0lf \n","grid",Lx/a,L/a,L/a);
	printf("%-15s \t %.2lfx%.2lfx%.2lf\n","simulation box",Lx,L,L);
	// export hydrodynamic numbers?
		
	/* Echo relevant macro values to log file*/
	printf("========================================\n");
	printf("%-15s \t %-15s\n","Macro","Value");
	printf("----------------------------------------\n");
	printf("%-15s \t %lf\n","EQN_EPS (not used)",EQN_EPS);
	printf("%-15s \t %e\n","TIME_EPS",TIME_EPS);
	printf("%-15s \t %d \n","CHECK_COMPRESSIBILITY",CHECK_COMPRESSIBILITY); //not fully funtional yet. TO DO
	printf("%-15s \t %lf \n","DENSITY_TOL",DENSITY_TOL);
	printf("%-15s \t %d \n","DENSITY_MAX",DENSITY_MAX);
	printf("%-15s \t %d \n","GALILEAN_SHIFT",GALILEAN_SHIFT);
	printf("----------------------------------------\n");
	printf("%-15s \t %d\n","MONITOR_EQUILIBRATION",MONITOR_EQUILIBRATION);
	printf("%-15s \t %d\n","EXPORT_STATES",EXPORT_STATES);
	printf("%-15s \t %d\n","EXPORT_VEL_PROFILE",EXPORT_VEL_PROFILE);
	printf("%-15s  %d\n","EXPORT_VEL_PROFILE_SKIP",EXPORT_VEL_PROFILE_SKIP); 
	printf("%-15s \t %d\n","EXPORT_TEMPERATURE",EXPORT_TEMPERATURE);

	
	/* Echo additional info to log file */
	printf("========================================\n");
	printf("Other information \n");
	printf("----------------------------------------\n");
	printf("%-30s \t %.4lf\n","Mean free path",dt * sqrt(T*m_inv));
	printf("%-30s \t %.4lf\n", "Kinetic viscosity /eta/",density * T * dt * ((density/(density-1+exp(-density)))-0.5)/(a*a*a) );
	printf("%-30s \t %.4lf\n","Collision viscosity /eta/", m * (density-1+exp(-density))/(12*a*dt));
	printf("%-30s \t %.4lf\n","Total viscosity /eta/ ", density * T * dt * ((density/(density-1+exp(-density)))-0.5)/(a*a*a) + m * (density-1+exp(-density))/(12*a*dt));
	printf("%-30s \t %.4lf\n","Knudsen number (Kn)",dt * sqrt(T*m_inv)/a); // Kn>>1 is free molecular regime, Kn<<1 is hydrodynamic regime, Kn around 1 is gas rarefaction
	printf("----------------------------------------\n");
	
	/*------------------------------------------------------------------*/
	/*	RNG INITIALIZATION					    */
	/*------------------------------------------------------------------*/
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
	gsl_rng_set(r , seed); //sets the seed for the RNG r
	
	printf("%-30s \t %ld\n","RNG seed", seed);

	
	/*------------------------------------------------------------------*/
	/*	MEMORY ALLOCATION AND INITIALIZATION OF ARRAYS		    */
	/*------------------------------------------------------------------*/
	/* Particle quantities lists */
	double * pos_rmo;
	double * vel_rmo;
	double * acc_rmo;
	double * cell_vel_rmo;
	double * cell_rnd_vel_rmo;
	double * cell_mass;	/// cell_mass[ci] is the total mass of the cell ci 
	double ** pos;
	double ** vel;
	double ** acc;
	double ** cell_vel;
	double ** cell_rnd_vel;
	
	/* auxiliary lists */
	//double * densities;  /// array of length n_cells which ... volumes?
	/*densities = malloc(n_cells * sizeof(double));
	if (densities==NULL){ printf("Error allocating densities in mpc.c\n"); exit(EXIT_FAILURE); }
	*/
	int ** cell_occupation;
	int * cell_occupation_rmo;
	int * c_p;				/// c_p[i] is the cell index of particle i 
	
	
	int max_oc = (int)(DENSITY_MAX * density);
	//(int)density + (int)(density * DENSITY_TOL); //rounding down
	int slice_size = cylinder.n_cells_dim[1] * cylinder.n_cells_dim[2];

	printf("%-30s \t %d\n","cell_occupancy storage limit", max_oc);

	//printf("Minimum cell occupancy: %d\n", (int)density-(int)(DENSITY_TOL*density));
	printf("========================================\n");
	
	cell_occupation_rmo=malloc((max_oc+1)*n_cells*sizeof(int));
	if (cell_occupation_rmo==NULL) {printf("Error allocating cell_occupation_rmo in mpc.c\n"); exit(EXIT_FAILURE);}
	cell_occupation = malloc(n_cells*sizeof(int*));
	if (cell_occupation==NULL) {printf("Error allocating cell_occupation in mpc.c\n"); exit(EXIT_FAILURE);}
	
	c_p = malloc( n_part*sizeof(int) );     /// cell idx from particle idx
	if (c_p==NULL) {printf("Error allocating c_p in mpc.c\n"); exit(EXIT_FAILURE);}


	
	/* Data export lists */
	double * slice_SAM_scalar; //density or particle occupation (casted to double), kinetic energy, mass, temperature
	double * slice_SAM_vector_rmo;
	double * slice_CAM_scalar_rmo;
	double * slice_CAM_vector_rmo;
	double ** slice_SAM_vector; // velocity
	double ** slice_CAM_scalar; //density or particle occupation (casted to double), kinetic energy, mass, temperature
	double ** slice_CAM_vector; // velocity	
	
	
 
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
	
	cell_mass = malloc( n_cells*sizeof(double) ); 
	if (cell_mass==NULL) {printf("Error allocating cell_mass in mpc.c\n"); exit(EXIT_FAILURE);}
	
	slice_SAM_vector_rmo = malloc( 3 * slice_size * sizeof(double) ); //any x-slice has nynz cells
	if (slice_SAM_vector_rmo==NULL) {printf("Error allocating slice_SAM_vector_rmo in mpc.c\n"); exit(EXIT_FAILURE);}
	slice_CAM_scalar_rmo = malloc( (max_oc + 1) *  slice_size * sizeof(double) );//CAM data includes local densities, hence +1
	if (slice_CAM_scalar_rmo==NULL) {printf("Error allocating slice_CAM_scalar_rmo in mpc.c\n"); exit(EXIT_FAILURE);}
	slice_CAM_vector_rmo = malloc( (3 * max_oc + 1) *  slice_size * sizeof(double) );
	if (slice_CAM_vector_rmo==NULL) {printf("Error allocating slice_CAM_vector_rmo in mpc.c\n"); exit(EXIT_FAILURE);}
	
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
	
	
	for(i=0; i<n_cells; i++)
	{
		cell_vel[i] = &cell_vel_rmo[3*i];
		cell_rnd_vel[i] = &cell_rnd_vel_rmo[3*i];
		/* Initialization is done at the beginning of collide() */
		cell_occupation[i] = &cell_occupation_rmo[(max_oc+1)*i];
	}
	
	slice_SAM_vector = malloc(  slice_size * sizeof(double*) );
	if (slice_SAM_vector==NULL) {printf("Error allocating slice_SAM_vector in mpc.c\n"); exit(EXIT_FAILURE);}
	slice_CAM_scalar = malloc(  slice_size * sizeof(double*) );
	if (slice_CAM_scalar==NULL) {printf("Error allocating slice_CAM_scalar in mpc.c\n"); exit(EXIT_FAILURE);}
	slice_CAM_vector = malloc(  slice_size * sizeof(double*) );
	if (slice_CAM_vector==NULL) {printf("Error allocating slice_CAM_vector in mpc.c\n"); exit(EXIT_FAILURE);}
	
	for(i=0; i<slice_size; i++)
	{
		slice_SAM_vector[i] = &slice_SAM_vector_rmo[ 3*i ];
		slice_CAM_scalar[i] = &slice_CAM_scalar_rmo[ (max_oc+1)*i ];
		slice_CAM_vector[i] = &slice_CAM_vector_rmo[ (3*max_oc+1)*i ];
	}
	
	slice_SAM_scalar = malloc(  slice_size * sizeof(double) );
	if (slice_SAM_scalar==NULL) {printf("Error allocating slice_SAM_scalar in mpc.c\n"); exit(EXIT_FAILURE);}
	
		
	/* Initialize and complete memory allocation */
	/* Position is uniformly distributed in space within vessel: geometry dependent */
	/* Velocity is normally distributed (Maxwell-Boltzmann) 			*/
	/* Note: Components of the velocity are chosen uncorrelated   			*/
	for(i=0; i<n_part; i++)
	{
		/* Complete memory allocation */
		pos[i] 	 = &pos_rmo[ 3*i ];
		
		/* Initialize */
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
		
		/*---------------------------------------------------------------------------------------------------*/
		// compute_acc(g,&(acc[i][0])); //This line is included in the previous two lines, so it is redundant.
		//it would, however, be necessary if there were forces acting on the plasma particles
		 
		 /* These lines were used in the testing of the fractioning of the timestep in streaming step */
		//acc[i][0] = g;// -((double)rand()/(RAND_MAX));
		//acc[i][1] = 40.0*(1-2*((double)rand()/(RAND_MAX)));
		//acc[i][2]=  40.0*(1-2*((double)rand()/(RAND_MAX)));
		/*---------------------------------------------------------------------------------------------------*/
		
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
	/*	EXPORT pre-simulation data				    */
	/*------------------------------------------------------------------*/

	/* Export collision grid data for visualization, and the vessel geometry data */
	#if EXPORT_STATES
		export_vtk_grid(cylinder);
		export_vessel_geometry(cylinder,steps);
		export_vtk_plasma(0, pos, vel,acc, n_part);
	#endif
	
	/*------------------------------------------------------------------*/
	/*	PREPARE optional pre-simulation variables	            */
	/*------------------------------------------------------------------*/
	#if MONITOR_EQUILIBRATION
		FILE * eq_fp;
		eq_fp = fopen("./../experiments/equilibration.dat","w");
		if(eq_fp == NULL){printf("mpc.c: Error opening eq_fp file. Aborting...\n "); exit(EXIT_FAILURE);}
	#endif
	
	
	
	/*------------------------------------------------------------------*/
	/*     PREPARE  x-positions where the velocity profile is exported    */
	/*------------------------------------------------------------------*/
	/* Choose x-position at which the velocity profile will be exported (see implementation in collide.c)<--mysterious comment */
	/* Memory allocation */ 
	#if EXPORT_VEL_PROFILE
	
		/* Check we know at what x-positions the velocity profile will be exported */
		if( argc<3 )
		{
			printf("mpc.c: Missing .dat file containing positions of x_slices for velocity profile study\nCall as:\n./mpcSim input.dat x_slices.dat\nAborting...\n");
			exit(EXIT_FAILURE);
		}
	
		/* Generous allocation of memory */
		double * x_slices;
		x_slices = malloc( cylinder.n_cells_dim[0] * sizeof(double) );
		if( x_slices==NULL ){ printf("mpc.c: Allocating memory for x_slices failed.\nAborting..."); exit(EXIT_FAILURE);}
		int * x_slices_idxs;
		x_slices_idxs = malloc( cylinder.n_cells_dim[0] * sizeof(int) );
		if( x_slices_idxs == NULL){ printf("mpc.c: Allocating memory for x_slices_idxs failes\nAborting...|n"); exit(EXIT_FAILURE);}
	
		/* Prepare to read in the values */
		int num_slices=0;
	
		/* Open the file with the x-positions of the slices */
		FILE * x_slices_fp;
		x_slices_fp = fopen(argv[2],"r");
		if(x_slices_fp == NULL){ fprintf(stderr,"Cannot open x_slices file.\n Aborting...\n"); exit(EXIT_FAILURE);}
		
		/* Read the x-positions, count how many we have */
		while ( fscanf(x_slices_fp, "%lf", &value) == 1) { //value was defined above as an auxiliary variable
			x_slices[num_slices] = value;
			num_slices++;
		}
		
		
		/*Health check*/
		assert( num_slices <= cylinder.n_cells_dim[0] && num_slices >=0 ); 
		for(i=0; i<num_slices;i++){
			if( x_slices[i] >= Lx || x_slices[i] <0 ) {
				fprintf(stderr,"Slice of interest num %d is outside cylinder. Exiting...\n",i);
				exit(EXIT_FAILURE);
			}
		}
	
		/* Calculate indices of fist cell in each of the slices */
		for(i=0; i<num_slices; i++)
		{
			x_slices_idxs[i] = slice_start_index(x_slices[i], cylinder.n_cells_dim[1]*cylinder.n_cells_dim[2], cylinder.a);
			/* Health check */
			assert( x_slices_idxs[i] < cylinder.n_cells );
			assert( x_slices_idxs[i] >=0 );
		}
		
	#endif
	

	
	/*------------------------------------------------------------------*/
	/*	SIMULATION						    */
	/*------------------------------------------------------------------*/
	// Import state?
	// Equilibration run?
	// production run?
	
	
	/* Export velocity profile setup*/
	int file_counter = 0;
	int counter = 0;
	

	#if GALILEAN_SHIFT == 0
		shift[0] = 0.0;
		shift[1] = 0.0;
		shift[2] = 0.0;
	#endif
	
		

	/*----------------------------------*/
	/* STREAM - COLLIDE ALGORITHM 	    */
	/* Loop for a fixed number of times */
	/*----------------------------------*/
	for(i=1; i <= steps; i++)
	{
	
		/*--------------------*/
		/*  MPC-AT algorithm  */
		/*--------------------*/
		
		/* Streaming step */
		stream( dt, n_part, g, cylinder, pos, vel, acc);
		 
		assert( TEST_all_particles_in_lumen(n_part, cylinder, pos) == 1 );
		

		/* Update timestep */
		t += dt; 
		
		
		/* Draw a random 3D vector for the Galilean shift, if required */
		#if GALILEAN_SHIFT
			shift[0] = a * ( -0.5 + gsl_rng_uniform_pos(r) ) ;
			shift[1] = a * ( -0.5 + gsl_rng_uniform_pos(r) );
			shift[2] = a * ( -0.5 + gsl_rng_uniform_pos(r) );
		#endif
		
		/* Collision step */
		collide(n_part, T, m, m_inv, cylinder, r, c_p, cell_mass, cell_vel, cell_rnd_vel, shift, pos, vel); 
		// BEWARE! the velocity of cells after this is not to be trusted, due to the Galilean shift: it is necessary to undo it!
		// c_p contains the occupation information, but with the Galilean shift: if we wanted to use it, we would need to undo it.		
		
		#if CHECK_COMPRESSIBILITY
			// substitute or combine with an assert?
			//check_compressibility(cell_mass, n_cells, m_inv, density);
		#endif
		
		/*--------------------*/
		/*  DATA EXPORT	      */
		/*--------------------*/

		/* Export particle data for a 3D paraview visualisation */
		#if EXPORT_STATES
			export_vtk_plasma(i, pos, vel,acc, n_part);		
		#endif
		
	

		
		/* Export velocity data every 10 timesteps, after the first equilibration_time timesteps*/
		#if EXPORT_VEL_PROFILE
			const char data_filename[60]="./../experiments/velprof_slice%d";
			char slice_filename[100]="";
		
			counter++;
			if( counter==equilibration_time + EXPORT_VEL_PROFILE_SKIP )
			{
				encage(pos, null_shift, n_part, cylinder, c_p, 1, cell_occupation, max_oc);
				/* Update data header (1/2) */
				data_header[0] = file_counter;
				data_header[6] = i;
				
				for(j=0; j<num_slices; j++)
				{
					/* Gather data for the slice*/
					slice_CAM_velocities(x_slices[j], cylinder, cell_occupation, vel, 1, slice_CAM_scalar); // only x-component of the velocities
					/* Update data header (2/2) */
					data_header[4] = x_slices_idxs[j];
					data_header[5] = x_slices_idxs[j] + cylinder.n_cells_dim[1]*cylinder.n_cells_dim[2] - 1;
					assert( data_header[1] == cylinder.n_cells_dim[0] );
					assert( data_header[2] == cylinder.n_cells_dim[1] );
					assert( data_header[3] == cylinder.n_cells_dim[2] );
					/* build export filename */
					snprintf(slice_filename,100,data_filename,j);
					/* export the data */
					export_CAM_slice(1, slice_CAM_scalar, slice_filename, data_header_size, data_header);
				}
				file_counter++;
				counter = equilibration_time;
			}
		#endif
		
		#if MONITOR_EQUILIBRATION
			total_momentum(n_part, m, vel, aux3Dvector);
			fprintf(eq_fp,"%lf \t %lf \t %lf \t %lf\n", aux3Dvector[0], aux3Dvector[1], aux3Dvector[2], total_kinetic_energy(n_part, m, vel) );
		#endif
		
	
		

	}/* end for loop (stream-collide steps) */
	

	/*------------------------------------------------------------------*/
	/*	FINISH BUSINESS						    */
	/*------------------------------------------------------------------*/
	
	#if MONITOR_EQUILIBRATION
		fclose_flag = fclose(eq_fp);
		if(fclose_flag != 0){printf("mpc.c: Error closing eq_fp file. Aborting...\n "); exit(EXIT_FAILURE);}
	#endif
	
	time(&now);
	printf("\nFinish: %s\n",ctime(&now));
	
	/*------------------------------------------------------------------*/
	/*	EXITING							    */
	/*	Free memory.						    */
	/*------------------------------------------------------------------*/
	/* Clean up */
	
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
	free(cell_occupation_rmo);
	free(cell_occupation);
	//free(densities);
	free(slice_SAM_scalar);
	free(slice_SAM_vector_rmo);
	free(slice_SAM_vector);
	free(slice_CAM_scalar_rmo);
	free(slice_CAM_scalar);
	free(slice_CAM_vector_rmo);
	free(slice_CAM_vector);
	free(x_slices);
	
	gsl_rng_free(r); /* free memory associated with the rng r */

	return(EXIT_SUCCESS);
}





