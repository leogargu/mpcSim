
#include "collide.h"
#include<gsl/gsl_randist.h>
#include<gsl/gsl_poly.h>
#include "macros.h"
#include <math.h>


///////////////////////////////////////////////////////////////////////////////
/// Gives the index of the cell from the cell coordinates in the grid.
/// call with 3d vectors!
///
///////////////////////////////////////////////////////////////////////////////
inline int cell_coord2cell_idx(int * cell_coord, const int * n_cells_dim)
{
	return cell_coord[2] + n_cells_dim[2] * ( cell_coord[0] * n_cells_dim[1] + cell_coord[1] );
}

/*converts from the cell index to the cell coordinates*/
/* Coded for completeness, but it is better not to use: involves two divisions */
/*     NOT TESTED !!!!!!!!!     */
inline void cell_idx2cell_coord(int cell_idx, const int * n_cells_dim, int * cell_coords)
{
	int aux;
	cell_coords[0] =  cell_idx % (n_cells_dim[1]*n_cells_dim[2]);
	aux = cell_idx - cell_coords[0]*n_cells_dim[1]*n_cells_dim[2];
	cell_coords[1] = aux % n_cells_dim[2];
	cell_coords[2] = aux - n_cells_dim[2]*cell_coords[1];
	
}



///////////////////////////////////////////////////////////////////////////////
/// \brief Gives the index of the cell a particle in position pos is at, with the grid
/// shifted by vector shift
///
///////////////////////////////////////////////////////////////////////////////
inline int pos2cell_idx(const Geometry geometry, const double * pos, const double * shift )
{
	int i;
	int cell_coord[3];
	int ci;
	
	/* The compiler should unroll this */
	for(i=0; i<3; i++)
	{
		ci = (int)( (pos[i] - shift[i]) / geometry.a ); 
		cell_coord[i] = ci;
		
		/* Checking particles remain in simulation box, and applying periodic boundary conditions */
		if( ci >= geometry.n_cells_dim[i] )
		{
			if( i==0 ) // if in the x-direction (flow direction)
			{
				cell_coord[0] = ci - geometry.n_cells_dim[0];
			}else{
				//this should never happen if the correct geometry and reasonable shift vector are used 
				fprintf(stderr,"A particle escaped the simulation box during shift in pos2cell_idx\n"); 
				fprintf(stderr,"It was at (%.3lf,%.3lf,%.3lf), cell %d in direction %d when max_cell in that direction is %d\n",
					pos[0],pos[1],pos[2],ci,i,geometry.n_cells_dim[i]);
				fprintf(stderr,"Check geometry chosen and magnitude of the shift vector\n");
				exit(EXIT_FAILURE);
			}

			
		}
	
		
		if( ci < 0 )
		{
			if( i==0 ) //if in the direction of the flow (x-direction)
			{
				cell_coord[0] = ci + geometry.n_cells_dim[0];
			}else{
				fprintf(stderr,"A particle escaped the simulation box during shift in pos2cell_idx\n"); 
				fprintf(stderr,"It was at (%.3lf,%.3lf,%.3lf), cell %d in direction %d \n",
					pos[0],pos[1],pos[2],ci,i);
				exit(EXIT_FAILURE);
			}
		}
		
		
	}
	
	
	return cell_coord2cell_idx( cell_coord, geometry.n_cells_dim );
}



/* The collide step does not change the positions of the particles. */
void collide(const int n_part, const double T, const double m, const double m_inv, Geometry geometry, gsl_rng * r, 
	int * c_p, double * cell_mass ,double ** cell_vel, double ** cell_rnd_vel, double * shift, double ** pos, double ** vel)
{
	/*-------------------------------------------*/
	/*	STEP 1: variables setup 	     */
	/*-------------------------------------------*/
	int n_cells = geometry.n_cells;
	int i,j;	
	int ci;
	double cell_inv_mass = 0.0;
    	
    	/* reset cell arrays */
    	for (ci = 0; ci < n_cells; ci++) {
    		cell_mass[ci] = 0.0;
    		for( j=0; j<3; j++)
    		{
    			cell_vel[ci][j]=0.0;
    			cell_rnd_vel[ci][j]=0.0;
    		}
	}
	/* There is no need to reset c_p because it will be changed in the loop below */
	
	/*-------------------------------------------*/
	/*	STEP 2: scan through the particles   */
	/*-------------------------------------------*/
	/* From this point on, the grid should be considered shifted by the vector shift */

	/* Loop over the particles */
	for(i=0; i<n_part; i++)
	{

		/* Generate a new random relative velocity vector from Maxwell-Boltzmann distribution */
		/* Vector components are uncorrelated */
		double rnd_vel[3] ={ gsl_ran_gaussian_ziggurat(r, sqrt(T * m_inv) ) ,
				     gsl_ran_gaussian_ziggurat(r, sqrt(T * m_inv) ) ,
			   	     gsl_ran_gaussian_ziggurat(r, sqrt(T * m_inv) ) };    /* zero mean */
   	   			
		#if DEBUGGING_STREAMCOLLIDE
			fprintf(debug_fp,"%lf \t %lf \t %lf \t",rnd_vel[0],rnd_vel[1],rnd_vel[2]);	
		#endif	   	     
		/* Get the index of the cell particle i is at */
		c_p[i] = pos2cell_idx( geometry, pos[i], shift ); // canonize is included here, and check of being outside of the cylinder too
		//c_p[i] = ci; /* store the index of the cell that particle i belongs to */
				
		/* Add mass and momentum to that cell */
		for(j=0; j<3; j++)
		{
			cell_vel[c_p[i]][j] += m * vel[i][j]; // update momentum of the cell
			cell_rnd_vel[c_p[i]][j] += m * rnd_vel[j]; // update random momentum of the cell
		}
		cell_mass[c_p[i]] += m;
				
		/* assign new velocity to particle: step 1/2 */
		for(j=0; j<3; j++)
		{
			vel[i][j] = rnd_vel[j];
		}
	}
	#if DEBUGGING_STREAMCOLLIDE
		fprintf(debug_fp,"\n");	
	#endif	
		#if DEBUGGING_STREAMCOLLIDE
			int p;
			for(p=0;p<n_part;p++)
			{
				fprintf(debug_fp,"%d \t",c_p[p]);
			}
			fprintf(debug_fp,"\n");
			for(p=0;p<n_cells;p++)
			{
				fprintf(debug_fp,"%d \t",(int)(cell_mass[p]*m_inv));
			}
			fprintf(debug_fp,"\n");
		#endif
	
	/*-------------------------------------------*/
	/*	STEP 3: loop over the cells 	     */
	/*-------------------------------------------*/
	for(ci=0; ci<n_cells; ci++)
	{
		if(cell_mass[ci] != 0.0)
		{
			cell_inv_mass=1.0/cell_mass[ci];
			//double cell_inv_mass = 1.0/cell_mass[ci];
			for(j=0; j<3; j++)
			{
				cell_vel[ci][j] = cell_inv_mass*( cell_vel[ci][j]-cell_rnd_vel[ci][j] );	
			}

		}else{
			//do nothing, cell_vel[ci] is already set to {0,0,0}, as is cell_rnd_vel[ci]
		}
		
		
	}
	
	
	/*-------------------------------------------*/
	/*	STEP 4: loop over the particles	     */
	/*-------------------------------------------*/
	/* assign new velocities */
	
	for(i=0; i<n_part; i++)
	{
		for(j=0; j<3; j++)
		{
			vel[i][j] += cell_vel[ c_p[i] ][j]; /* assign new velocity: step 2/2 */
		}
	
	}
	
	/* cell_vell still contains the average momentum (velocity if all masses of plasma particles are equal) when collide() returns to main */
	/* Likewise, c_p is up to date and available outside collide. */
}



