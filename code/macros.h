#ifndef MACROS_H /* header guard */
#define MACROS_H

/* Uncomment to deactivate all assert statements */
//#define NDEBUG

#define VERBOSE 1
#define DEBUGGING 1
/* Expots positions, velocities, accelerations and rnd velocities of all particles, to debig the stream-collide alternating algorithm*/
#define BINARY_EXPORT 0
#define DEBUGGING_QUARTIC_SOLVER 0   // put this in quartic_solver.h, or change back to DEBUGGING

/* used in quartic_solver.h, only in the function solve_quartic, which is not used by stream (it relies of gsl entirely after testing) */
#define EQN_EPS 1e-5
/*  A NOTE ON MACHINE EPSILON VALUES:

	Running find_machine_epsilon in Godzilla and Minerva yields:

	FLT_EPSILON: 1.19209e-07
	FLT_MIN: 1.17549e-38
	DBL_EPSILON: 2.22045e-16
	DBL_MIN: 2.22507e-308
	LDBL_EPSILON: 1.0842e-19
	LDBL_MIN: 3.3621e-4932
*/



#define TIME_EPS 1e-10

/* Exports vtk files for every timestep. These vtk files are snapshots of the simulation, containing positions, velocities 
and accelerations of all particles. It also exports auxiliary vtk files: collision_grid.vtk and vessel_geometry.vtk */
#define EXPORT_STATES 0

/* Exports files to do CAM average over a x_slice, ready to plot the parabolic profile of POiseuille flow */
#define EXPORT_VEL_PROFILE 1

/* Defines how many timesteps to skip between exports. Set to 1 to export the vel profile at every timestep (with EXPORT_VEL_PROFILE 1).*/
#define EXPORT_VEL_PROFILE_SKIP 1

/* Checks, at every timestep, that the instantaneous density at any collision cell is less than DENSITY_TOLxInitial density */
#define CHECK_COMPRESSIBILITY 0

/* At setup, there is an avarege density per cell (bulk cells have a higher density than this, boundary cells have a lower density)*/
/* DENSITY_TOL defines the % of variation allowed from the average density. For example, 50% variation is DENSITY_TOL = 0.5 */
#define DENSITY_TOL 0.5
/* This defines the maximum occupancy allowed in any given cell as density*DENSITY_MAX. This is not related to complressibility chekcs. */
/* Thsi is solely motivated by the definition of the cell_occupation data structure */
#define DENSITY_MAX 10

/* Exports a file, at regular timesteps, containing the temperature in each collision cell. */
#define EXPORT_TEMPERATURE 0

/* Set to 1 when running equilibration tests. It exports the momentum (3 components) and kinetic energy per particle (average over all particles at that timestep),
at every timestep. Redo calculations using the script.pbs.oXXX file data, which includes all geometrical data and the RNG seed. */
#define MONITOR_EQUILIBRATION 1

/* Turns the shift of the grid prior to each collision step on/off*/
#define GALILEAN_SHIFT 1


/* This function checks whether a double is zero. This function is not necessary and its code can be changed
within the calling functions*/
inline static int is_zero(double x) 
{
	return  x == (double)0.0;
}




/* This function checks whether a double is zero, within the precision specified by eps */
inline static int is_approx_zero(double x, double eps)
{
	return x > -eps && x < eps;

}


#endif /* header guard */





/* Note: max values for ints */
/*
INT_MAX: 2147483647		(~2e09)
UINT_MAX: 4294967295		(~4e09)
LONG_MAX: 9223372036854775807   (~9e18)
ULONG_MAX: 18446744073709551615 ~(18e18)
*/
