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


/* Exports the positions/velocities of all particles at a given timestep */
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
			fprintf(fp,"%.6lf \t %.6lf \t %.6lf \n",data[i][0],data[i][1],data[i][2]);
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


inline int pos2cell_idx(const Geometry geometry, const double * pos, const double * shift );






/* It exports the data to do the averages to calculate the 3Dvelocity profile at a given */
/* point x, at a particular time step */
/* The first row in the exported file contains, in this order: the x value at which the slice is taken, */
/* the number of cells in the slice (ny*nz), and the number of cells in each column in that slice (nz)  */
inline void export_vel_profile(int n_part, double density, double ** vel, double ** pos, const Geometry cylinder, double x_slice, int file_number ,int step)
{
	/* Declare variables */
	FILE * fp;
	char file_name[30];
	int i,j;
	double ** slice_vel;
	int max_part_density = (int)10*density; // If the particle density in any given cell exceeds this value, the simulation stops (high compressibility effects)
	int index;
	int nynz = cylinder.n_cells_dim[1] * cylinder.n_cells_dim[2];
	double nullshift[3]={0.0,0.0,0.0};
	
	/* Set up file */
	snprintf(file_name,sizeof(char)*30,"./DATA/velprofile_%d.dat",file_number);

	fp=fopen(file_name,"w");
	fprintf(fp,"%lf \t %d \t %d \t %d  \n", x_slice, nynz, cylinder.n_cells_dim[2], step);//nynz in the file is the number of rows, after the first one, that represent real data (it might happen that the last cells in the slice of interest do not have particles)
	
	/* Allocate memory and initialize */
	/* http://c-faq.com/aryptr/dynmuldimary.html */
	slice_vel = malloc(nynz * sizeof(double *)); //allocate array of pointers to rows
	if (slice_vel==NULL) {printf("Error allocating slice_vel in mpc.c\n"); exit(EXIT_FAILURE);}	
	slice_vel[0] = malloc(nynz * max_part_density * sizeof(double)); // allocate the whole memory of the 2D array to store the data, and initialize to 0
	if (slice_vel[0]==NULL) {printf("Error allocating slice_vel[0] in mpc.c\n"); exit(EXIT_FAILURE);}
	
	for(i = 1; i<nynz; i++)
	{
		slice_vel[i] = slice_vel[0] + i * max_part_density; // assign the pointers to the correct start of their rows
	}
	
	for(i=0; i<nynz; i++) 
	{
		for(j=0; j< max_part_density; j++) //Optimising: the rightmost index in the inner loop.
		{
			slice_vel[i][j] = 0.0;
		}
	}
	
	/* Find slice of interest  */
	int first_cell = nynz*(int)floor(x_slice/cylinder.a);
	
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
			index = cell_idx - first_cell;
			slice_vel[index][0]+=1; // update the counter
			if(slice_vel[index][0]>max_part_density)
			  {
			  	  printf("export_vel_profile in mpc.c: Local density exceeded max tolerance (%d). Aborting...\n",max_part_density);
			  	  exit(EXIT_FAILURE);
			    
			  }
			slice_vel[ index ][ (int)slice_vel[index][0] ] = vel[i][0];  //export x-velocity only
		}
	}
	
	/* DO NOT DELETE THIS: it will be useful if the above comment is implemented
	//---------------------------
	// Gather data 
	for(i=0; i<n_part; i++) //loop over particles
	{
		if( c_p[i] >= first_cell && c_p[i] < first_cell+nynz)//if particle i is in a cell in the slice of interest...
		{
			index = c_p[i] - first_cell;
			slice_vel[index][0]+=1; // update the counter
			if(slice_vel[index][0]>max_part_density)
			  {
			  	  printf("export_vel_profile in mpc.c: Local density exceeded max tolerance (%d). Aborting...\n",max_part_density);
			    exit(EXIT_FAILURE);
			    
			  }
			slice_vel[ index ][ (int)slice_vel[index][0] ] = vel[i][0];  //export x-velocity only
		}
	}
	*/
	
	int local_density = 0;
	/* Export to file */
	for(i=0; i<nynz; i++) // i is the cell in the slice of interest
	{
		local_density = (int)slice_vel[i][0];
		fprintf(fp, "%d\t", local_density);
		if( local_density>0 )
		{
			for(j=1; j<=local_density; j++ ) // j is the particle in cell i
			{
				fprintf(fp,"%lf\t",slice_vel[i][j]);
			}
			fprintf(fp,"\n");
		}else{
			fprintf(fp,"\n");
	  	}
	}
	
	/* Clean up and release memory*/
	fclose(fp);
	free(slice_vel[0]);
	free(slice_vel);

}



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

inline double total_kinetic_energy_ppart(int n_part, double m, double ** vel)
{
	double energy = 0.0;
	int i;
	
	for(i=0; i<n_part; i++)
	{
		energy += m * ( vel[i][0]*vel[i][0] + vel[i][1]*vel[i][1] + vel[i][2]*vel[i][2] );
	}
	
	return 0.5*energy/(double)n_part;
}

/*-------------------------------*/
/*		MAIN		 */
/*-------------------------------*/

int main(int argc, char **argv) {
	/*------------------------------------------------------------------*/
	/*	SETUP							    */
	/* 	Definition of constants, initialization of variables,       */
	/*      memory allocation, I/O setup			            */
	/*								    */
	/*	 See description of units in Watari2007 		    */
	/*------------------------------------------------------------------*/
	/* Parameters, constants and variables */
	const double a = 1.0;		/// Size of a collision cell (cubic cell)
	const double Lx = 30.0;		/// Length of the simulation box in the direction of the flow //10
	const double L = 10.0;		/// Width and height of the simulation box //4
	const int n_part = 30000;   	/// Number of particles in the simulation box//density=10
	const double T = 1.0;  		/// Temperature, in kT units
	const double m = 1.0; 		/// Mass of a particle (plasma)
	const double m_inv = 1.0/m;	/// Inverse mass of a particle (plasma)
	const double radius = 0.5*( L - a );	   	/// Radius of the cylindrical vessel. It is a/2 smaller than half the box side to allow for grid shifting 
	const int n_cells = (int)((Lx*L*L)/(a*a*a));  	/// Number of cells in the simulation box	
	const double dt = 1.0;		/// Timestep
	const double density = (1.0*n_part)/(L*L*Lx);  	/// Number of particles per collision cell //const int
	
	/* TODO: check Whitmer's implementation of Verlet step and why he does not have a lambda */
	//const double lambda = 0.65;  /// Verlet's algorithm lambda parameter. 0.65 is a "magic number" for DPD, see page 5 of \cite{verlet}," Novel Methods in Soft Matter Simulations", ed. Karttunen et. al.
	/* consider setting lambda inside stream()? maybe as a define? */
	
	const double g = 0.5;		/// "Gravity" (acceleration that drags the flow in the x-direction)
	double t = 0.0;			/// Time
	int i,j; 		 	/* Generic loop counters */
	const long int steps = 2000;	/// Number of steps in the simulation
	long int equilibration_time = 300; 
	
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
	int * cell_idx;
		
	/* Check parameter values are compatible */
	if( fmod(Lx,a) != (double)0.0 ){
		printf(" Lx must be a multiple of a. Exiting...\n");
		exit(EXIT_FAILURE);
	}else if( fmod(L,a) != (double)0.0 ){
		printf("L must be a multiple of a. Exiting...\n");
		exit(EXIT_FAILURE);
	}

	fprintf(stderr,"--------------------------------\n");
	fprintf(stderr,"Parameter \t Value\n");
	fprintf(stderr,"--------------------------------\n");
	fprintf(stderr,"density \t %.2lf\n",density);
	fprintf(stderr,"n_part \t\t %d\n",n_part);
	fprintf(stderr,"a \t\t %.3lf\n",a);
	fprintf(stderr,"n_cells \t %d\n",n_cells);
	fprintf(stderr,"Lx \t\t %.3lf\n",L);
	fprintf(stderr,"L \t\t %.3lf\n",Lx);
	fprintf(stderr,"radius \t\t %.3lf\n",radius);
	fprintf(stderr,"m \t\t %.3lf\n",m);
	fprintf(stderr,"T \t\t %.3lf\n",T);
	fprintf(stderr,"grid \t\t %.0lfx%.0lfx%.0lf \n",L/a,L/a,Lx/a);
	fprintf(stderr,"simulation box \t %.2lfx%.2lfx%.2lf\n",Lx,L,L);
	fprintf(stderr,"Mean free path \t %.2lf\n",dt * sqrt(T*m_inv));
	fprintf(stderr,"--------------------------------\n");
 
	

	/* Export collision grid data for visualization, and the vessel geometry data */
	export_vtk_gid(cylinder);
	export_vessel_geometry(cylinder,steps);
	
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
	#if VERBOSE
	printf("RNG seed: %ld\n",seed);
	#endif
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
	
	for(i=0; i<n_cells; i++)
	{
		cell_vel[i] = &cell_vel_rmo[3*i];
		cell_rnd_vel[i] = &cell_rnd_vel_rmo[3*i];
		/* Initialization is done at the beginning of collide() */
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
	//int file_idx = 0;
	
	/* Variables to monitor equilibration*/
	double momentumPpart[3] = {0.0,0.0,0.0}; //linear momentum per particle, at any given timestep
	int equilibration_counter=0;
	printf("Timestep \t p_x \t\t p_y \t\t p_z \t\t e_k \n");
	double energyPpart = 0.0;
	
	/* Loop for a fixed number of times */
	for(i=1; i<= steps; i++)
	{
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
		shift[0] = a * ( -0.5 + gsl_rng_uniform_pos(r) ) ;
		shift[1]= a * ( -0.5 + gsl_rng_uniform_pos(r) );
		shift[2]= a * ( -0.5 + gsl_rng_uniform_pos(r) );
		
		/* Collision step */
		collide(n_part, T, m, m_inv, cylinder, shift, r, c_p, cell_mass, cell_vel, cell_rnd_vel, pos, vel); // BEWARE! the velocity of cells after this is not to be trusted, due to the Galilean shift: it is necessary to undo it!
		// c_p contains the occupation information, but with the Galilean shift: if we wanted to use it, we would need to undo it.		
		
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
				export_vel_profile(n_part, density, vel, pos, cylinder, x_slice, file_counter, i); 
				counter = equilibration_time;
			}
		#endif


		equilibration_counter++;
		if(equilibration_counter>20) // export data on total momentum per particle every 10 timesteps
		{
			total_momentum_ppart(n_part, m, vel, momentumPpart);
			energyPpart = total_kinetic_energy_ppart(n_part, m, vel); 
			printf("%d \t %.4lf \t %.4lf \t %.4lf \t %.4lf\n", i, momentumPpart[0],momentumPpart[1],momentumPpart[2], energyPpart);
			equilibration_counter = 0;
		}

	}/* end for loop (stream-collide steps) */

	
	/*------------------------------------------------------------------*/
	/*	EXITING							    */
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
	free(cell_idx);
	free(c_p);
	free(cell_mass);
	free(shift);

	gsl_rng_free(r); /* free memory associated with the rng r */
	

	return(EXIT_SUCCESS);
}





