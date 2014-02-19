
#include "stream.h"
#include "quartic_solver.h"
#include <stdio.h>
#include "macros.h"



#define GATHER_OFFENDERS_STATISTICS 0
#define TRACKING 0

#if TRACKING
	FILE * fp;
	FILE * fp2;
#endif


/*=======================*/
/* Global variable, needed for counting offenders statistics */
int stuck;
/* another global variable for debugging/tracking of particles  */
//int display = 0;
/*======================*/

			
inline int TEST_particle_in_lumen(Geometry cylinder, double * pos)
{
	if(pos[0] > cylinder.Lx || pos[0] < 0. )
	{
		fprintf(stderr,"TEST_particle_in_lumen: Particle out of bounds in x direction\n");
		return 0;
	}

	double distance_sq = (pos[1]-cylinder.L_half)*(pos[1]-cylinder.L_half)+(pos[2]-cylinder.L_half)*(pos[2]-cylinder.L_half);
	
	if( distance_sq >= cylinder.max_sq)
	{
		fprintf(stderr,"distance^2: %lf > (radius+epsilon)^2 %lf\n",distance_sq, cylinder.max_sq);
		fprintf(stderr,"outside tolerance zone!\n");
		fprintf(stderr,"TEST_particle_in_lumen: Increase value of epsilon\n");
		return 0;
	}
	
	return 1;
}


/* Second function of a class of auxiliary tester functions. Not part of the program. Used for debugging */
inline int TEST_all_particles_in_lumen(int n_part, Geometry cylinder, double ** pos)
{
	int i;
	for(i=0;i<n_part;i++)
	{
		if( TEST_particle_in_lumen(cylinder, &(pos[i][0])) )
		{
			//do nothing
		}else{
			fprintf(stderr,"TEST_all_particles_in_lumen: Particle %d outside lumen\n",i);
			return 0;
		}
	}
	return 1;
}

///////////////////////////////////////////////////////////////////////////////
/// \brief Imposes periodic boundary condition in the x-direction.
/// Other directions are unnaffected.
///
///////////////////////////////////////////////////////////////////////////////
inline void canonize(double Lx, double * pos)
{
	int testvar=(int)( pos[0] / Lx );
	
	if( pos[0] > Lx )
	{
		if( testvar > 1 )
		{
			fprintf(stderr,"canonize: Particle covered %d times the length of the capillary in 1 timestep\n",testvar);
			fprintf(stderr,"pos[0]: %lf, Lx: %lf\n",pos[0],Lx);
			fprintf(stderr,"canonize: Aborting...\n");
			exit(EXIT_FAILURE);
			//pos[0] -= ((int)( pos[0] / Lx ))*Lx; <--This would be the way to implement it, if we tolerated this behaviour!!
		}
		pos[0] -= Lx;
		
	}
	if( pos[0] < 0.0 )
	{
		
		if( testvar < -1 )
		{
			fprintf(stderr,"canonize: Particle covered %d times the length of the capillary in 1 timestep\n",testvar);
			fprintf(stderr,"pos[0]: %lf, Lx: %lf\n",pos[0],Lx);
			fprintf(stderr,"canonize: Aborting...\n");
			exit(EXIT_FAILURE);
			//pos[0] += +((int) ( pos[0] / Lx))*Lx; <--This would be the way to implement it, if we tolerated this behaviour!!
		}
		pos[0] += Lx;
	}
}


/* TODO: This needs adaptation (when RBC enter in the picture) */
/* Basically, this is compute_force */
/* acc is the output variable */
inline void compute_acc(double g, double * acc) 
{
	acc[0] = g;
	acc[1] = 0.0;//10.0*(1-2*((double)rand()/(RAND_MAX)));//5000;10
	acc[2] = 0.0;//10.0*(1-2*((double)rand()/(RAND_MAX)));//5000; 20
	return; /* back to main */
}

/* Bounce-back the velocity */
inline void bb_velocity(double * vel)
{
	int j;
	
	for( j=0; j<3; j++ )
	{
		vel[j] = -vel[j];
	}
	return; /* back to main */
}


/* Calculate the coefficients of the quartic defining the impact times of the particle trajectory
   with the wall*/
inline void calculate_coeffs( double * pos, double * vel, double * acc, const Geometry cylinder, double * coeff )
{
	double acc_z, acc_y, aux_y, aux_z;
	
	acc_z = acc[2]  ;
	acc_y = acc[1]  ;
	aux_y = pos[1] - cylinder.L_half;
	aux_z = pos[2] - cylinder.L_half;
	
	coeff[0] = aux_y*aux_y + aux_z*aux_z - cylinder.radius*cylinder.radius;
	coeff[1] = 2 *( aux_y * vel[1] + aux_z * vel[2] ) ;
	coeff[2] = vel[1] * vel[1] + vel[2]*vel[2] + aux_y*acc_y + aux_z*acc_z ;
	coeff[3] = vel[1]*acc_y + vel[2]*acc_z ;
	coeff[4] = 0.25*( acc_y*acc_y + acc_z*acc_z) ;
	
	return; /* back to main */
}



/* This function performs the Verlet update in pos and vel IF the particle is not streamed outside of the capillary.
 "Outside" is defined as the area beyond the tolerance zone, ie a circle of radius (R + epsilon)
 If the particle escapes, this function does not modify the position nor the velocity of the particle, and returns 1, otherwise returns 0.
 Epsilon (encoded in the struct cylinder, defined in mpc.h) is the half-width of the tolerance zone. The tolerance zone extends a distance 
 epsilon outside the real boundary and a distance -epsilon inside the lumen.							*/

/* Note that this function does not account for the case in which a particle escapes and comes back to the capillary within dt.
   A way to avoid this would be to solve the quartic for every particle, or every particle "sufficiently close" to the wall: if no solutions are 
   found, that is to be understood as the particle never intersecting the wall. <----------TO DO. */
inline int vel_verlet(double g, double * pos, double * vel, double * acc, const Geometry cylinder, double dt)
{
	int i;
	double acc_temp[3];
	/* Note we work with the projection of the vectors onto the cross-section plane. */
	
	/* -------------------------------------------- */
	/* VERLET STEP 1:  tentatively update positions */
	/* -------------------------------------------- */
		
	double aux_pos[3];
	double aux_y, aux_z;

	for(i=0; i<3; i++)    	
	{
		aux_pos[i] = pos[i] + dt * vel[i] + 0.5*dt*dt*acc[i];
	}
	
	aux_y = aux_pos[1] - cylinder.L_half ;
	aux_z = aux_pos[2] - cylinder.L_half ;
	double aux_circunf_sq = aux_y*aux_y + aux_z*aux_z;

	/* If the particle is out, do not stream at all and return 1 */
	if( aux_circunf_sq >= cylinder.max_sq  ) 
	{
		return 1;
	}
	
	/* ...Else, we accept the streaming step: actually update positions  */
	for(i=0; i<3; i++)
	{
		pos[i] = aux_pos[i];
	}
	
	
	/* correct the positions above: impose periodic boundary conditions */
	canonize(cylinder.Lx, pos);
	
	
	/* ---------------------------------------------------- */
	/* VERLET STEP 2: update velocity (and acceleration) 	*/	
	/* ---------------------------------------------------- */	
	
	compute_acc( g, &(acc_temp[0]) );
	for(i=0; i<3; i++)
	{
		//vel[i] += acc[i]*dt; // For forward Euler scheme
		vel[i] += 0.5 * (acc[i] + acc_temp[i] ) * dt ;	/* update velocities*/
		acc[i] = acc_temp[i];	 			/* update accelerations (forces)*/
	}
	
	
	return 0; /* back to main */
}




/* This function returns the amount of time that the particle still needs to be streamed for */
/* Time tau must be positive. Time dt is the maximum time */
/* NOTE ON INLINING ------------
 * This function is not inlined. The keyword inline has been removed to eliminate warwinings at compile time.
 * A way to get the same results bypassing this is to define a new variable tau2 in boundary_fractional_stream, and
 * set its value to either +tau or -tau. Then, copy the code of advance_rajectory into the body of boundary_fractional_stream.
 * However, this might result in boundary_fractional_stream not being inlined anymore. Whether that is a more efficient approach
 * or not, can only be decided by the use of a profiler. <---- check later
 * --------------------------------------
 */
double advance_trajectory(double g, double * pos, double * vel, double * acc, const Geometry cylinder, double tau, double dt)
{
	double remaining_time;
	double has_escaped;
	
	if( tau < dt ) /*then the particle impacts on boundary before its streaming is done */
	{
		has_escaped = vel_verlet(g, pos, vel, acc, cylinder, tau );
		if( has_escaped > 0 ) /* catch errors */
		{
			fprintf(stderr, "advance_trajectory: Something went wrong with dt=%.3lf>tau=%.3lf branch. Aborting...\n",dt,tau);
			fprintf(stderr,"advance_trajectory: Something went wrong. Testing for verlet_epsilon...\n");
			/*----------------------------------------------------------*/
			double aux_pos[3];
			double aux_y, aux_z;
			int j;
			for(j=0; j<3; j++)    	
			{
				aux_pos[j] = pos[j] + tau * vel[j] + 0.5*tau*tau*acc[j];
			}
			fprintf(stderr,"Code attempted to move particle to (%.15lf, %.15lf, %.15lf)\n",aux_pos[0],aux_pos[1],aux_pos[2]);
			aux_y = aux_pos[1] - cylinder.L_half ;
			aux_z = aux_pos[2] - cylinder.L_half ;
			double aux_circunf_sq = aux_y*aux_y + aux_z*aux_z;
			fprintf(stderr,"This position is at a radius %.15lf, tolerance zone is between %.15lf and %.15lf\n",sqrt(aux_circunf_sq),sqrt(cylinder.min_sq),sqrt(cylinder.max_sq));			
			fprintf(stderr,"Aborting now...\n");							
			exit(EXIT_FAILURE);
		}
		
		#if TRACKING
			fprintf(fp2,"%.20lf \t",tau);		
		#endif
		
		bb_velocity( vel );
		remaining_time = dt - tau;
		
	}else{ /* then the particle never actually impacts against the wall in this advance along its trajectory */
		has_escaped = vel_verlet(g, pos, vel, acc, cylinder, dt );
		if( has_escaped > 0 ) /* catch errors */
		{
			fprintf(stderr,"tau=%.15lf, dt=%.15lf\n",tau,dt);
			fprintf(stderr,"advance_trajectory: Something went wrong with dt<tau branch. Aborting...\n");
			fprintf(stderr,"advance_trajectory: Something went wrong. Testing for verlet_epsilon...\n");
			/*----------------------------------------------------------*/
			double aux_pos[3];
			double aux_y, aux_z;
			int j;
			for(j=0; j<3; j++)    	
			{
				aux_pos[j] = pos[j] + dt * vel[j] + 0.5*dt*dt*acc[j];
			}
			fprintf(stderr,"Code attempted to move particle to (%.15lf, %.15lf, %.15lf)\n",aux_pos[0],aux_pos[1],aux_pos[2]);
			aux_y = aux_pos[1] - cylinder.L_half ;
			aux_z = aux_pos[2] - cylinder.L_half ;
			double aux_circunf_sq = aux_y*aux_y + aux_z*aux_z;
			fprintf(stderr,"This position is at a radius %.15lf, tolerance zone is between %.15lf and %.15lf\n",sqrt(aux_circunf_sq),sqrt(cylinder.min_sq),sqrt(cylinder.max_sq));
			fprintf(stderr,"Aborting now...\n");		
			exit(EXIT_FAILURE);
		}
		
		#if TRACKING
			fprintf(fp2,"%.20lf \t",dt);	
		#endif
		
		//bb_velocity( vel ); //why was this here?!
		remaining_time = 0.0;
	}
	

	return remaining_time;
}


/* This function attempts to stream a particle initially at the boundary for a time dt, and if it cannot do it it will
continue fractioning the timestep and streaming the particle from a point at the boundary to another point at the 
boundary until the streaming time is finished. Then the particle will end up either in the lumen or at the boundary. */
/* It returns the number of times the timestep was fractioned per function call (i.e. for a given particle) */
/* This function is to be called ONLY for particles AT THE BOUNDARY */
inline int boundary_fractional_stream(double g,  double * pos, double * vel, double * acc, const Geometry cylinder, double * coeff, double dt)
{
	int continue_streaming = 1;
	double tau;
	int j;	
	int count = 0 ;
		
	
	while( continue_streaming )
	{
		
		/* Calculate the coefficients */
		calculate_coeffs( pos, vel, acc, cylinder, coeff);
		tau = find_trajectory_boundary_gsl( coeff );
				
		/*if( display ) /* display is a global variable, initially set to 0 and changed by the function stream()*
		{
			//printf("Coeffs for particle: %.14lf + %.14lf t + %.14lf t^2+ %.14lf t^3+ %.14lf t^4\n",coeff[0],coeff[1],coeff[2],coeff[3],coeff[4]);
			printf("boundary_fractional_stream: tau=%.14lf\n",tau);		
		}*/
				
		if( is_zero(tau) )
		{
			/* then particle is stuck at the boundary */
			for(j=0;j<3;j++)
			{
				vel[j] = 0.0;
			}
			continue_streaming = 0;
			stuck++;
			fprintf(stderr,"boundary_fractional_stream: A particle is stuck at the boundary\n");
		}else if( tau > 0.0)
		{
			#if TRACKING
				fprintf(fp,"%.20lf \t %.20lf \t %.20lf \t",pos[0],pos[1],pos[2]);
				fprintf(fp,"%.20lf \t %.20lf \t %.20lf \t",vel[0],vel[1],vel[2]);
				fprintf(fp,"%.20lf \t %.20lf \t %.20lf \t",acc[0],acc[1],acc[2]);
				
			#endif
			
			dt = advance_trajectory( g, pos, vel, acc, cylinder, tau, dt);
	
			if( dt > 0.0 )
			{
				count++;
			}
			
		}else{
			bb_velocity( vel );
			#if TRACKING
				fprintf(fp,"%.20lf \t %.20lf \t %.20lf \t",pos[0],pos[1],pos[2]);
				fprintf(fp,"%.20lf \t %.20lf \t %.20lf \t",vel[0],vel[1],vel[2]);
				fprintf(fp,"%.20lf \t %.20lf \t %.20lf \t",acc[0],acc[1],acc[2]);
			#endif
			
			dt = advance_trajectory( g, pos, vel, acc, cylinder, -1.0*tau, dt);
			if( dt > 0.0 )
			{
				count++;
			}
		}
		
		if( dt ==(double)0.0 )
		{
			continue_streaming = 0;
		}
	}
	
	#if TRACKING
		fprintf(fp,"%.20lf \t %.20lf \t %.20lf \t",pos[0],pos[1],pos[2]);
		fprintf(fp,"%.20lf \t %.20lf \t %.20lf \t",vel[0],vel[1],vel[2]);
		fprintf(fp,"%.20lf \t %.20lf \t %.20lf \t",acc[0],acc[1],acc[2]);		
	#endif
	
	
	return count; /* back to main */
}
		




///////////////////////////////////////////////////////////////////////////////
/// \brief Streaming step 
/// COMPLETE FOR MPC PARTICLES ONLY (only force is gravity in flow direction)
///////////////////////////////////////////////////////////////////////////////
/* Returns number of offenders   */
void stream(double dt, const int n_part, const double g, const Geometry cylinder, double ** pos, double ** vel, double ** acc)
{	
	#if TRACKING
 		fp=fopen("particle_track.dat","wr");
 		fp2=fopen("particle_times.dat","wr");
	#endif

	int i;
	int count;
	double has_escaped = 1.0; 
	double tau=0.0;
	double aux_y,aux_z;
	double circunf_sq;
	double * coeff;
	coeff = malloc( 5 * sizeof(double) );
	if (coeff==NULL) {fprintf(stderr,"Error allocating coeff in stream.c\n"); exit(EXIT_FAILURE);}

	#if GATHER_OFFENDERS_STATISTICS
		int * dt_fractions;
		dt_fractions = malloc( n_part * sizeof(int) );
		if (dt_fractions==NULL) {fprintf(stderr,"Error allocating dt_fractions in stream.c\n"); exit(EXIT_FAILURE);}
		
	#endif
	
	
	stuck = 0;
	/* Stream all particles a timestep dt */
	for(i=0; i<n_part; i++) /* For each plasma particle... */
	{
	
		
		/*if( i==0 )
		{
			display = 1; /* display is a global variable *
		}else{
			display = 0;
		}*/
		
		
		/* reset counter */
		count = 0;
		
		/*...figure out where the particle is. */
		aux_y = pos[i][1]-cylinder.L_half;
		aux_z = pos[i][2]-cylinder.L_half;
		circunf_sq = aux_y*aux_y + aux_z*aux_z;
	
		
		
		if( TEST_all_particles_in_lumen(n_part, cylinder, pos)!=1 )
		{
			fprintf(stderr,"stream: Aborting...\n");
			exit(EXIT_FAILURE);
		}
		
		/* this bit of code is duplicated, also appears in TEST_particle_in_lumen*/ 
		/* If particle is outside the wall beyond the tolerance zone, something went wrong: abort the simulation */
		//if( circunf_sq >= cylinder.max_sq )// this check probably can be eliminated because vel_verlet also checks <--- TO DO
		//{
		//	fprintf(stderr, "A particle escaped beyond the tolerance zone: aborting. ");
		//	exit(EXIT_FAILURE);	
			
		/* FOR PARTICLES IN THE LUMEN */
		//}else if( circunf_sq < cylinder.min_sq)
		
		if (circunf_sq < cylinder.min_sq)
		{
			#if TRACKING
				fprintf(fp,"%.20lf \t %.20lf \t %.20lf \t",pos[i][0],pos[i][1],pos[i][2]);
				fprintf(fp,"%.20lf \t %.20lf \t %.20lf \t",vel[i][0],vel[i][1],vel[i][2]);
				fprintf(fp,"%.20lf \t %.20lf \t %.20lf \t",acc[i][0],acc[i][1],acc[i][2]);
			#endif
			
			/* attempt streaming */
			has_escaped = vel_verlet(g, pos[i], vel[i], acc[i], cylinder, dt );
			if( has_escaped > 0 )
			{
				/* Fraction the timestep: find impact time */	
				calculate_coeffs( pos[i], vel[i], acc[i], cylinder, coeff );
				tau = find_impact_time_gsl(coeff,dt);
				
				/*if( display ) /* display is a global variable *
				{
					//calculate_coeffs( pos[i], vel[i], acc[i], cylinder, coeff );
					//double aux1=sqrt( (vel[i][1])*(vel[i][1])+(vel[i][2])*(vel[i][2]) );
					//double aux2= sqrt( (pos[i][1]-cylinder.L_half)*(pos[i][1]-cylinder.L_half)+(pos[i][2]-cylinder.L_half)*(pos[i][2]-cylinder.L_half) );
					//printf("aux1=%lf, aux2=%lf\n",aux1,aux2);
					//printf("Angle of particle %d is: %.14lf\n",i,(1.0/M_PI)*180*acos(coeff[1]/(aux1*aux2)   ));
					//printf("Coeffs for particle %d: %.14lf + %.14lf t + %.14lf t^2+ %.14lf t^3+ %.14lf t^4\n",i,coeff[0],coeff[1],coeff[2],coeff[3],coeff[4]);
					//printf("Coeffs accounting for tau: %.14lf + %.14lf + %.14lf+ %.14lf + %.14lf\n",coeff[0],coeff[1]*tau,coeff[2]*tau*tau,coeff[3]*tau*tau*tau,coeff[4]*tau*tau*tau*tau);
					//printf("Solution is tau=%.14lf\n",tau);
				}*/
				

				if( tau < 0.0 )
				{
					fprintf(stderr,"stream.c: something went wrong, aborting...\n");
					exit(EXIT_FAILURE);
				}
		    		
				/* stream for a time tau: this places the particle at the boundary */
		    		has_escaped = vel_verlet(g, pos[i], vel[i], acc[i], cylinder, tau );
		    		
		    		
		    		if( has_escaped > 0 ) /* catch errors */
		    		{
		    			fprintf(stderr,"stream.c: Something went wrong when taking a lumen particle to the boundary. Testing for verlet_epsilon...\n");
			    		/*----------------------------------------------------------*/
			    		double aux_pos[3];
			    		double aux_var_y, aux_var_z;
					int j;
					for(j=0; j<3; j++)    	
					{
						aux_pos[j] = pos[i][j] + dt * vel[i][j] + 0.5*dt*dt*acc[i][j];
					}
					fprintf(stderr,"Code attempted to move particle %d to (%.15lf, %.15lf, %.15lf)\n",i,aux_pos[0],aux_pos[1],aux_pos[2]);
					aux_var_y = aux_pos[1] - cylinder.L_half ;
					aux_var_z = aux_pos[2] - cylinder.L_half ;
					double aux_circunf_sq = aux_var_y*aux_var_y + aux_var_z*aux_var_z;
					fprintf(stderr,"This position is at a radius %.15lf, tolerance zone is between %.15lf and %.15lf\n",sqrt(aux_circunf_sq),sqrt(cylinder.min_sq),sqrt(cylinder.max_sq));
					fprintf(stderr,"Aborting now...\n");
					exit(EXIT_FAILURE);
		    		}
				
		    		/* No-stick boundary conditions */ 
		    		bb_velocity( vel[i] );
		    		
		    		#if TRACKING
					fprintf(fp2,"%.20lf \t",tau);
				#endif
		    		
				//count should be 1 at this point...
				
				/* Now try to stream from the boundary. */
		    		count = boundary_fractional_stream(g, pos[i], vel[i], acc[i], cylinder, coeff, dt-tau);
		    		count++;//...so we correct count here
			}else{
				/* If particle did not escape, it has already been updated by the call to vel_verlet: go to next particle*/
				 //do nothing
				 #if TRACKING
				 	fprintf(fp,"%.20lf \t %.20lf \t %.20lf \t",pos[i][0],pos[i][1],pos[i][2]);
					fprintf(fp,"%.20lf \t %.20lf \t %.20lf \t",vel[i][0],vel[i][1],vel[i][2]);
					fprintf(fp,"%.20lf \t %.20lf \t %.20lf \t",acc[i][0],acc[i][1],acc[i][2]);
					fprintf(fp2,"%.20lf \t",dt);
					 
				 #endif
			}
		/* FOR PARTICLES AT THE BOUNDARY (i.e. WITHIN TOLERANCE ZONE)*/	
		}else{	
			/*if( display ) /*display is a global variable *
			{
				calculate_coeffs( pos[i], vel[i], acc[i], cylinder, coeff );
				printf("Coeffs for particle %d: %.14lf + %.14lf t + %.14lf t^2+ %.14lf t^3+ %.14lf t^4\n",i,coeff[0],coeff[1],coeff[2],coeff[3],coeff[4]);
				printf("Coeffs accounting for tau: %.14lf + %.14lf + %.14lf+ %.14lf + %.14lf\n",coeff[0],coeff[1]*tau,coeff[2]*tau*tau,coeff[3]*tau*tau*tau,coeff[4]*tau*tau*tau*tau);

			}*/
			count = boundary_fractional_stream(g, pos[i], vel[i], acc[i], cylinder, coeff, dt);
		}
		
		#if GATHER_OFFENDERS_STATISTICS
			dt_fractions[i] = count;
		#endif
		
		#if TRACKING
			fprintf(fp ,"\n");
			fprintf(fp2,"\n");				
		#endif
		
	}/* end for -loop over plasma particles- */
		
		
	// Now post-process the bounce-back statistics
	#if GATHER_OFFENDERS_STATISTICS
		int total_bounces = 0;
		int bounces[10]={0,0,0,0,0,0,0,0,0,0};
		int max_bounces=0;
		int num_multiple_bouncers=0;
		
		for(i=0; i<n_part; i++)
		{	
			if( dt_fractions[i] > max_bounces )
			{
				max_bounces = dt_fractions[i];
			}
			total_bounces += dt_fractions[i];
			switch( dt_fractions[i] )
			{
				case 0 : bounces[0]++; break;
				case 1 : bounces[1]++; /* printf("Particle %d bounced once\n",i);*/    break;
				case 2 : bounces[2]++; /* printf("Particle %d bounced twice\n",i);*/ break;
				case 3 : bounces[3]++; /* printf("Particle %d bounced 3 times\n",i);*/ break;
				case 4 : bounces[4]++; printf("Particle %d bounced 4 times\n",i);break;
				case 5 : bounces[5]++; printf("Particle %d bounced 5 times\n",i); break;
				case 6 : bounces[6]++; printf("Particle %d bounced 6 times\n",i); break;
				case 7 : bounces[7]++; printf("Particle %d bounced 7 times\n",i); break;
				case 8 : bounces[8]++; printf("Particle %d bounced 8 times\n",i); break;
				default : bounces[9]++; printf("Particle %d bounced %d times\n",i, dt_fractions[i]);
			}
					
		}
		#if DEBUGGING
			printf("Part/Bounced \t x0 \t x1 \t x2 \t x3 \t x4 \t x5 \t x6  \t x7 \t x8 \t >8 \t Stuck \t Total Bounces\n");
			printf("---------------------------------------------------------------------------------------------------------------------------\n");
			printf("%d \t\t %d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t %d\n",n_part, bounces[0], bounces[1], bounces[2], bounces[3], bounces[4],bounces[5], bounces[6], bounces[7],bounces[8],bounces[9],stuck,total_bounces);
		#endif
		
		#if GATHER_OFFENDERS_STATISTICS
			free(dt_fractions);
		#endif
	#endif
	
	
	free(coeff);
	
	#if TRACKING
 		fclose(fp);
 		fclose(fp2);
	#endif
	
	
	return; /* back to main*/

} /* end of stream */		
 		
		
	









			
		
		
		

