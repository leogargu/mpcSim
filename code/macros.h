#ifndef MACROS_H /* header guard */
#define MACROS_H


#define VERBOSE 1
#define DEBUGGING 1
#define BINARY_EXPORT 0

#define DEBUGGING_QUARTIC_SOLVER 0   // put this in quartic_solver.h, or change back to DEBUGGING

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

#define EXPORT_STATES 0

#define EXPORT_VEL_PROFILE 1



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

#endif /* header guard */





/* Note: max values for ints */
/*
INT_MAX: 2147483647		(~2e09)
UINT_MAX: 4294967295		(~4e09)
LONG_MAX: 9223372036854775807   (~9e18)
ULONG_MAX: 18446744073709551615 ~(18e18)
*/
