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




/*------------------------------*/
/* 	METHODS PROTOTYPES 	*/
/*------------------------------*/

#endif	/* header guard */



