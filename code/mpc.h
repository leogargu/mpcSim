/** MPC simulations. Header file
 *
 * Author: Leonor Garcia-Gutierrez, 2013
 * Version: mpcSim-1.0
 *
 */
 
#ifndef MPC_H	/* header guard */
#define MPC_H

/*------------------------------*/
/* 	GLOBAL VARIABLES 	*/
/*------------------------------*/


/*-------------------------------*/
/*		STRUCTS		 */
/*-------------------------------*/
struct geometry
{
	double a;		/// Size of a collision cell (cubic cell)
	double Lx;		/// Length of the simulation box in the direction of the flow
	double L;		/// Width and height of the simulation box
	double L_half;		/// L/2
	double radius;	   	/// Radius of the cylindrical vessel. It is a/2 smaller than half the box side to allow for grid shifting 
	double epsilon;		/// Half-width of the tolerance zone (i.e., 'numerical thickness' of the boundary). Set to a percentage of the radius.
	double min_sq;		/// Set to (redius - epsilon)^2
	double max_sq;		/// Set to (radius + epsilon)^2
	int n_cells;		/// Number of cells in the simulation box
	int n_cells_dim[3];	/// Number of cells in each axis direction (array of 3 elements x,y,z)
	
};

typedef struct geometry Geometry;


/*//////////////////////////////////////////*/
/*-----------------*/
/* TEST FUNCTIONS  */
/*-----------------*/



/*//////////////////////////////////////////*/
/* This tests whether the given position of a particle is within the lumen */
/*inline int TEST_particle_in_lumen(Geometry cylinder, double * pos)
{
	if(pos[0] > cylinder.Lx || pos[0] < 0. )
	{
		return 0;
	}

	double distance_sq = (pos[1]-cylinder.L_half)*(pos[1]-cylinder.L_half)+(pos[2]-cylinder.L_half)*(pos[2]-cylinder.L_half);
	double radius_sq = cylinder.radius * cylinder.radius;
	
	if( distance_sq >= radius_sq )
	{
		//printf("distance: %lf > radius %lf\n",distance_sq, radius_sq);
		//printf("outside circle!\n");
		return 0;
	}
	
	return 1;
}


/* Second function of a class of auxiliary tester functions. Not part of the program. Used for debugging *
inline int TEST_all_particles_in_lumen(int n_part, Geometry cylinder, double ** pos)
{
	int i;
	for(i=0;i<n_part;i++)
	{
		if( TEST_particle_in_lumen(cylinder, &(pos[i][0])) )
		{
			//do nothing
		}else{
			printf("Particle %d outside lumen\n",i);
			return 0;
		}
	}
	return 1;
}
*/


/*------------------------------*/
/* 	METHODS PROTOTYPES 	*/
/*------------------------------*/

#endif	/* header guard */



