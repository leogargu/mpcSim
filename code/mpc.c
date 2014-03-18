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

/* Adapted from: */
// http://stackoverflow.com/questions/5837572/generate-a-random-point-within-a-circle-uniformly
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


/* Exports the positions/velocities/accelerations of all particles at a given timestep */
/* export is done in ASCII or BINARY depending on the value of the constant
   BINARY_EXPORT defined in macros.h					*/
/* This function is a helpful auxiliary tool, even if not used in a particular version of mpcSIM: DO NOT DELETE*/
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





/* It writes a VTK file with positions and velocities of all particles */
/* TODO: at the moment the export is done in ASCII only, binary option has not been implemented yet*/
/* possible to do by keeping the lines in ASCII and writing the data (only) in binary to the same file stream*/
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



/* This exports the VTK file defining the grid used (not the vessel) */
/* This is called only once at the start of the simulation. It exports in ASCII*/
/* No need to do in in binary, unless the geometry becomes much more complicated (i.e. some triangulated surface)*/
/* NOTE that if the directory DATA is not present, the program crashes with a Segmentation Fault */
inline void export_vtk_gid(Geometry cylinder)
{
	FILE * fp;
	fp = fopen("./DATA/collision_grid.vtk","w");
	fprintf(fp,"# vtk DataFile Version 2.0\nGrid representation\nASCII \nDATASET STRUCTURED_POINTS\nDIMENSIONS %d %d %d\nORIGIN 0 0 0\nSPACING %.1lf %.1lf %.1lf",cylinder.n_cells_dim[0]+1,cylinder.n_cells_dim[1]+1,cylinder.n_cells_dim[2]+1,cylinder.a,cylinder.a,cylinder.a);
	fclose(fp);	
}


/* This function exports the data necessary for the python script to generate a
 cylinder of the right dimensions on paraview					*/
inline void export_vessel_geometry(Geometry cylinder, int num_steps)
{
	FILE * fp;
	fp = fopen("./DATA/vessel_geometry.py","w");
	fprintf(fp,"# Cylinder data\n");
	fprintf(fp,"Lx=%lf\n",cylinder.Lx);
	fprintf(fp,"L=%lf\n",cylinder.L);
	//fprintf(fp,"a=%lf\n",cylinder.a);
	fprintf(fp,"radius=%lf\n",cylinder.radius);
	fprintf(fp,"L_half=%lf\n",cylinder.L_half);
	fprintf(fp,"num_frames=%d\n",num_steps);
	fclose(fp);
}

/*nynz is the product of the number of cells in the y direction and in the z direction*/
/* a is the length of a side of each cubic collision cell*/
/* That slice includes the cells of indices from slice_start_index return to that value-1+nynz, inclusive*/
inline int slice_start_index(double x, int nynz, double a)
{	
	return nynz*(int)floor(x/a);
}



/* Calculates the total linear momentum of the system per particle. Returns a 3D vector */
inline void total_momentum_ppart(int n_part, double m, double ** vel, double * momentum)
{
	/* Define variables */
	int i,j;
	double factor = 1.0/n_part;
	
	/* Initialise output vector */
	for(j=0; j<3; j++)
	{
		momentum[j]=0.0;
	}
	
	/* Add momentums */
	for( i=0; i<n_part; i++) //loop over particles
	{
		for(j=0; j<3; j++) // loop over cartesian components of the velocity
		{
			momentum[j]+=m*vel[i][j];
		}
	}
	
	/* Calculate momentum per particle */
	for(j=0; j<3; j++)
	{
		momentum[j]*=factor;
	}
	
	return; /* Back to main */
}


/* returnd the square norm of a 3D vector v */
inline double norm_sq(double * v)
{
	return v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
}

/* Calculates the total energy per particle, over the whole system, at the time of the call*/
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
	energy = 0.5*energy;
	
	for(i=0; i<n_part; i++)
	{
		sigma += (0.5*m*norm_sq(vel[i])-energy)*(0.5*m*norm_sq(vel[i])-energy);
	}
	sigma = sigma /(double)n_part;
	
	output[0] = energy;
	output[1] = sigma;
	
	return;
}


/* Control routine that checks for compressibility effects: it stopsthe simulation if, at the time of being invoked, the density in any 
   collision cell is greater than DENSITY_TOL x initial density */
/* Notice this can be called without undoing the shift after collide */
/* NOT TESTED */
inline void check_compressibility(double * cell_mass, int n_cells, double m_inv, double density)
{
	int ci;
	for(ci=0; ci< n_cells; ci++)
	{
		if(m_inv*cell_mass[ci] > DENSITY_TOL*density)
		{
			fprintf(stderr,"Compressibility effects exceeded tolerance: *shifted* cell %d has instantaneous density of %.2lf\n",ci,cell_mass[ci]*m_inv);
			exit(EXIT_FAILURE);
		}
	}
	
	return; /* Back to main*/
}

		
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
	
	/*for(i=0; i<cylinder.n_cells; i++)
	{
		cell_mass[i] = 0.0;
	}*/
	
	for(i=0; i<n_part; i++)
	{
		cell_idx = pos2cell_idx(cylinder, pos[i], shift );// canonize is included here, and check of being outside of the cylinder too
		c_p[i] = cell_idx;
		//cell_mass[cell_idx]+=m; //this line...
		if(calculate_cell_occupation)
		{
			cell_occupation[cell_idx][0]++;//...and this line are equivalent
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



/* This function exports the data contained in the array data in a file format suitable for calculating a SAM average with average.c */
/* The data in data should contain a list of the average value of the property under study in each cell considered, in consecutive order. */
inline void export_SAM_data(double * data, char * filename, Geometry geometry, int cell_start, int cell_end, int file_number, int step)
{
	FILE * file_fp;
	char file_name[50];
	char filename1[50]="";
	strcpy(filename1,filename);
	strcat(filename1,"_%d.dat");
	snprintf(file_name,sizeof(char)*30,filename1,file_number);
	
	
	file_fp=fopen(file_name,"w");
	int i;		
	int num_cells = cell_end - cell_start +1;
	fprintf(file_fp,"%d \t %d \t %d \t %d \t %d \t %d \n",geometry.n_cells_dim[0],geometry.n_cells_dim[1],geometry.n_cells_dim[2], cell_start, cell_end, step);
	for(i=0; i<num_cells; i++)
	{
		fprintf(file_fp,"%lf\t",data[i]);
	}

	fclose(file_fp);		
		
	return; /* back to main*/
}

/* This function exports the data contained in the array data in a file format suitable for calculating a CAM average with average.c */
/* The data in data should contain a 2D array with each row containing teh values of interest for the particles in the corresponding collision cell. Cells are stored
   in consecutive order. The occupation numbers, ie length of the rows, ie the local instantaneous densities are stored in the array of ints local_density.*/
/* not tested yet*/
/* Any 2D array ontaining CAM data can be exported to a file calling this function */
inline void export_CAM_data(double ** data, int * local_density, char * filename, Geometry geometry, int cell_start, int cell_end, int file_number, int step)
{
	FILE * fp;
	char file_name[50];
	char filename1[50]="";
	strcpy(filename1,filename);
	strcat(filename1,"_%d.dat");
	snprintf(file_name,sizeof(char)*30,filename1,file_number);
	
	int num_cells = cell_end - cell_start + 1;
	int i,j;
	
	fp=fopen(file_name,"w");
	fprintf(fp,"%d \t %d \t %d \t %d \t %d \t %d \n", geometry.n_cells_dim[0], geometry.n_cells_dim[1], geometry.n_cells_dim[2], cell_start, cell_end, step);
	
	/* Export to file */
	for(i=0; i<num_cells; i++) // i is the cell in the slice of interest
	{	
		fprintf(fp, "%d\t", local_density[i]);
		if( local_density[i]>0 )
		{
			for(j=1; j<=local_density[i]; j++ ) // j is the particle in cell i
			{
				fprintf(fp,"%lf\t",data[i][j]);
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









/*-------------------------------*/
/*		MAIN		 */
/*-------------------------------*/

int main(int argc, char **argv) {
	
	int i=0,j=0; 		 	/* Generic loop counters */
	
	FILE * input_fp;
	input_fp = fopen(argv[1],"r");
	//input_fp = fopen("input.dat","r");
	if(input_fp == NULL){ fprintf(stderr,"Cannot open input file.\n Aborting...\n"); exit(EXIT_FAILURE);}
	
	
	char variable_name[15];
	double value=0.0;
	double variable_array[10];
	while (fscanf(input_fp, "%s\t%lf", variable_name, &value) == 2) {
		
		printf("reading variable %s = %.3lf\n",variable_name, value);
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
		
	// If not checking the temperature this is not needed: create a dummy pointer for encage() and carry on
	int ** cell_occupation;
	int * cell_occupation_rmo;
	int max_oc = (int)(density * DENSITY_TOL);
	cell_occupation_rmo=malloc(max_oc*cylinder.n_cells*sizeof(int));
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
		double momentumPpart[3] = {0.0,0.0,0.0}; //linear momentum per particle, at any given timestep
		double energyPpart = 0.0;
		double systemTemp = 0.0;
		int equilibration_counter=0;
		FILE * equilibration_fp;
 		//fprintf(equilibration_fp,"Timestep \t p_x \t p_y \t p_z \t e_k \n");
 		equilibration_fp=fopen("equilibration.dat","w");
 		double output[2]={0.0,0.0};
 	#endif
	
		
	
	
	#if DEBUGGING_STREAMCOLLIDE
		//FILE * debug_fp;
		debug_fp=fopen("debug_data.dat","w");
	#endif
	#if CHECK_TEMPERATURE
		int file_temp_counter = 0;
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
	
		/* Set up grid shift */ 
		//shift[0] = a * ( -0.5 + gsl_rng_uniform_pos(r) ) ;
		//shift[1] = a * ( -0.5 + gsl_rng_uniform_pos(r) );
		//shift[2] = a * ( -0.5 + gsl_rng_uniform_pos(r) );
		
		shift[0]=0.0;
		shift[1]=0.0;
		shift[2]=0.0;

		/* Collision step */
		collide(n_part, T, m, m_inv, cylinder, r, c_p, cell_mass, cell_vel, cell_rnd_vel, shift, pos, vel); // BEWARE! the velocity of cells after this is not to be trusted, due to the Galilean shift: it is necessary to undo it!
		// c_p contains the occupation information, but with the Galilean shift: if we wanted to use it, we would need to undo it.		
		
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
				//total_momentum_ppart(n_part, m, vel, momentumPpart);
				//energyPpart = total_kinetic_energy_ppart(n_part, m, vel); 
				//printf("%d \t %.4lf \t %.4lf \t %.4lf \t %.4lf\n", i, momentumPpart[0],momentumPpart[1],momentumPpart[2], energyPpart);
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
			total_momentum_ppart(n_part, m, vel, momentumPpart);
			total_kinetic_energy(n_part, m, vel, output); 
			energyPpart = output[0] /(double)n_part;
			//systemTemp = (2.0*energyPpart - m_inv*norm_sq(momentumPpart) )/3.0;
			systemTemp = (2*energyPpart-(momentumPpart[0]*momentumPpart[0]+momentumPpart[1]*momentumPpart[1]+momentumPpart[2]*momentumPpart[2])*m_inv)/3.0;
			fprintf(equilibration_fp,"%d \t %.4lf \t %.4lf \t %.4lf \t %.4lf \t %.4lf \t %.4lf\n", i, momentumPpart[0],momentumPpart[1],momentumPpart[2], energyPpart,sqrt(output[1])/(double)n_part, systemTemp);
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

	gsl_rng_free(r); /* free memory associated with the rng r */

	return(EXIT_SUCCESS);
}





