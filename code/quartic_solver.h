#ifndef QUARTIC_SOLVER_H
#define QUARTIC_SOLVER_H


#include"macros.h"
#include<stdio.h>
#include<math.h>
#include<gsl/gsl_poly.h>


//--------------------------------------------------------------------------------------
// Based on source code found at:
// http://tog.acm.org/resources/GraphicsGems/gemsiv/vec_mat/ray/solver.c
// Jochen Schwarze. "Cubic and Quartic Roots". In Graphics Gems, Academic Press, 1990, pp. 404?407.



/* this function returns the number of real solutions, which are stored in the array sols */
/* equation is coeff[3] x^3 + coeff[2] x^2 + coeff[1] x  + coeff[0] == 0 */
int solve_cubic(double c[4], double s[3])
{
	int	i, num;
	double	sub, A, B, C, sq_A, p, q, cb_p, D;
	double aux = 1.0 / c[3];
	
	static const double inv_3 = 1.0 / 3.0 ;
	static const  double inv_3_sq = 1.0 / 9.0; //inv_3 * inv_3;
	
	
// normalize the equation:x ^ 3 + Ax ^ 2 + Bx  + C = 0
A = c[2] * aux;
B = c[1] * aux;
C = c[0] * aux;

// substitute x = y - A / 3 to eliminate the quadric term: x^3 + px + q = 0

sq_A = A * A;
p = inv_3 * ( B - inv_3 * sq_A );
sub = inv_3 * A;
q = 0.5 * ( sub * (2 * inv_3_sq *sq_A - B) + C );

// use Cardano's formula

cb_p = p * p * p;
D = q * q + cb_p;

if ( is_zero(D) )
    {
    if ( is_zero(q) )
	{
	// one triple solution
	s[0] = 0.0;
	num = 1;
	}
    else
	{
	// one single and one double solution
	double u = cbrt(-q);
	s[0] = 2.0 * u;
	s[1] = - u;
	num = 2;
	}
    }
else
    if (D < 0.0)
	{
	// casus irreductibilis: three real solutions
	double phi = inv_3 * acos(-q / sqrt(-cb_p));
	double t = 2.0 * sqrt(-p);
	s[0] = t * cos(phi);
	s[1] = -t * cos(phi + M_PI * inv_3);
	s[2] = -t * cos(phi - M_PI * inv_3);
	num = 3;
	}
    else
	{
	// one real solution
	double sqrt_D = sqrt(D);
	double u = cbrt(sqrt_D + fabs(q));
	aux = p / u;
	if (q > 0.0)
	    s[0] = - u + aux ;
	else
	    s[0] = u - aux ;
	num = 1;
	}

// resubstitute
//sub = inv_3 * A;   
for (i = 0; i < num; i++)
    s[i] -= sub;
return num;
}



/* This function determines the real roots of a linear equation.*/						
/* It outputs the number of roots found.			*/
int solve_linear(double c[2], double s[1])
{
	if (is_zero(c[1]))
	{
		return 0;
	}
	s[0] = - c[0] / c[1];
	return 1;
}




/* This function determines the real roots of a quadric equation.*/						
/* It outputs the number of roots found.			*/
int solve_quadric(double c[3], double s[2])
{
	double p, q, D;
	
	// make sure we have a d2 equation
	if (is_zero(c[2]))
	{
		return solve_linear(c, s);
	}
	
	double aux = 1.0/c[2];
	
	// normal for: x^2 + px + q
	p = 0.5 * c[1] *aux ;
	q = c[0] * aux;
	D = p * p - q;

	if (is_zero(D))
	{
		// one double root
		s[0] = s[1] = -p;
		return 1;
	}

	if (D < 0.0)
		// no real root
	return 0;

	else
	{
		// two real roots
		double sqrt_D = sqrt(D);
		s[0] = sqrt_D - p;
		s[1] = -sqrt_D - p;
		return 2;
  	}
}



/* My custom function to treat degenerated cases, based NOT on gsl functions */
inline int solve_quartic_degenerate(double * coeff, double * sols)
{
	int num = -1;
		
		if( is_zero(coeff[3]) )
		{
			num = solve_quadric(coeff, sols);
		}else{
			num = solve_cubic( coeff, sols);
		}
	
	
	return num;
}


inline int solve_quartic_degenerate_gsl(double * coeff, double * sols)
{
	int num;
		
		if( is_zero(coeff[3]) )
		{
			num = gsl_poly_solve_quadratic(coeff[2], coeff[1], coeff[0],&sols[0],&sols[1]);
		}else{
			double aux = 1.0 / coeff[3] ;
			num = gsl_poly_solve_cubic( coeff[2] * aux, coeff[1] * aux, coeff[0] * aux, &sols[0], &sols[1], &sols[2]);
		}	
	
	return num;
}




// Returns the number of real solutions of the quartic. Stores such solutions in sols, which
// must be a pointer to 4 doubles.
// This function only works if the quartic is non-degenerate!
inline int solve_quartic(double * coeff, double * sols )
{
	int num;

	/* STEP 1: check degenerated cases */
	/*----------------------------------*/
	// make sure we have a d4 equation
	if( is_zero(coeff[4]) )
	{	
		/*solve_quartic_degenerate_gsl is also possible, but it is slower: */
		/* 345.450 VS 468.730 secs for cubics, (584627 equations solved) */
		/* 55.630  VS 75.300 secs for cuadrics (602482 equations solved) */
		num = solve_quartic_degenerate(coeff, sols);

	}else{ 	/* solve the not degenerated quartic */
		
		double z, u, v, sub, A, B, C, D, A_2, p, q, r;
		int i;	
		double xi[3]; 	/* array to store solutions of the resolvent cubic*/
		double aux;

		/* STEP 2: handle general case.     */
		/*----------------------------------*/
	
		aux = 1.0/coeff[4]; // this change in the computation of A, B, C, D made the computation of 2166 solutions go down from 2.08 secs to 1.95 secs

		/* Normalize the equation to obtain:   x ^ 4 + Ax ^ 3 + Bx ^ 2 + Cx + D = 0   */
		
		A = coeff[3] * aux;
		B = coeff[2] * aux;
		C = coeff[1] * aux;
		D = coeff[0] * aux;
	
		/* Substitute x = y - A / 4 to depress the quartic to:  x^4 + px^2 + qx + r = 0 (eliminate cubic term) */

		A_2 = A * A;
		aux = 0.0625 * A_2;   // recycling of variables.

		p = -0.375 * A_2 + B;
		q = 0.5 * A * ( 0.25* A_2 - B ) + C; 			
		r = aux *(  B - 3 * aux ) - 0.25 * A* C + D; 	
		
		if ( is_zero(r) )  
		{
			num = gsl_poly_solve_cubic( (double)0.0, p,q, &xi[0], &xi[1], &xi[2]);
			/* The number of real roots (either one or three) is stored in num, and their 
			* locations are stored in xi. If one real root is found then 
			* only xi[0] is modified. When three real roots are found they are stored in 
			* xi in ascending order. As in the quadratic case, finite precision 
			* may cause equal or closely-spaced real roots to move off the real axis into 
			* the complex plane, leading to a discrete change in the number of real roots.*/
    
			if( num > 0 )
			{
				sols[0] = xi[0]; 
 	   	    	}
    
 	   	}else{
 	   		/* solve the resolvent cubic... */
 	   		 	   		
 	   		num = gsl_poly_solve_cubic( -0.5* p, -r, 0.5*r*p-0.125*q*q, &xi[0], &xi[1], &xi[2]);
 	
 	   		/* and take the smallest real solution */ 
 	   		z = xi[0];

 	   		// ...to build two quadratic equations
 	   		u = z * z - r;
 	   		v = 2.0 * z - p;
 	   		
 	   	// try modifying this as if u or v are negative but very small, set to 0, otherwise respect?
 	   		if( u <= (double)0.0 && u > -EQN_EPS ) //is_approx_zero(u)
 	   		{
 	   			u = (double)0.0 ; 
 	   		}else if( u > 0.0 )
 	   		{
 	   			u = sqrt(u);
    	    		}else
    	    		{
    	    			#if DEBUGGING_QUARTIC_SOLVER
    	    				printf("Custom: u is negative! No real sols found\n");
    	    				printf("u=%.20lf\n",u);
    	    			#endif
    	    			return 0;
    	    		}
    
    
    	    		if( v <=(double)0.0 && v >-EQN_EPS  ) //is_approx_zero(v)
    	    		{
    	    			v = (double)0.0;
    	    		}else if( v > 0.0 )
    	    		{
    	    			v=sqrt(v);
    	    		}else
    	    		{
    	    			#if DEBUGGING_QUARTIC_SOLVER
    	    				printf("Custom: v is negative! no real sols found (?)\n");
    	    				printf("v=%.20lf\n",v);
    	    				printf("These are the coeffs I am trying to solve:\n");
    	    				int k;
    	    				for(k=0;k<5;k++)
    	    				{
    	    					printf(" %.30lf\n",coeff[k]);
    	    				}
    	    			#endif
    	    			return 0;	    
		    	}


		    	aux = q < 0.0 ? -v : v; // more recycling of variables
		    	
		 
		    	/* This only looks for real roots*/
		    	num = gsl_poly_solve_quadratic(1.0, aux, z-u , &sols[0], &sols[1]);
		   
		    	/* This only looks for real roots */
		    	num += gsl_poly_solve_quadratic(1.0, -aux, z+u, &sols[2], &sols[3]);
		}

		// resubstitute
		sub = 0.25 * A;

		for(i=0; i<num; i++)
		{
			sols[i] -= sub;
		}

	}
	
	return num;

}


/* sol has to be a pointer to an array of 4 doubles. This function returns the number of real solutions found */
/* This function does not return the complex solutions */
inline int solve_quartic_gsl( double * coeff, double * sol  )
{
	int i;
	int num = 0;

       /* handle degenerate cases */
       if( is_zero(coeff[4]) )
       {
       	       if( is_zero(coeff[3]) )
       	       {
       	       	       sol[0] = 0.0;
       	       	       sol[1] = 0.0;
		       num = gsl_poly_solve_quadratic(coeff[2], coeff[1], coeff[0], &sol[0], &sol[1]);
		       #if DEBUGGING_QUARTIC_SOLVER
		       		printf("GSL: Real solutions are: \n");
		       		for (i = 0; i < num; i++)
		       		{
		       			printf("z%d = %+.18f \n", i, sol[i]);
                   		}
                       #endif
                       return num;
		
		       	       
	       }else{
			double aux = 1.0 / coeff[3] ;
			sol[0] = 0.0;
       	       	        sol[1] = 0.0;
       	       	        sol[2] = 0.0;
			num = gsl_poly_solve_cubic( coeff[0] * aux, coeff[1] * aux, coeff[0] * aux, &sol[0], &sol[1], &sol[2]);
			#if DEBUGGING_QUARTIC_SOLVER
		       		printf("GSL: Real solutions are: \n");
		       		for (i = 0; i < num; i++)
		       		{
		       			printf("z%d = %+.18f \n", i, sol[i]);
                   		}
                       #endif
                       
                       return num;
 
	       }
       }else{
       	
       	       double sol_complex[8];
       	       /* Solve the quartic */
       	       gsl_poly_complex_workspace * w = gsl_poly_complex_workspace_alloc(5);
       	       gsl_poly_complex_solve(coeff, 5, w, sol_complex);
       	       gsl_poly_complex_workspace_free(w);
     
       	      
       	       #if DEBUGGING_QUARTIC_SOLVER
       	       		printf("GSL: Solutions are: \n");
       	       		for (i = 0; i < 4; i++)
       	       		{
       	       			printf("z%d = %+.25lf %+.25lf i \n", 
       	       				i, sol_complex[2*i], sol_complex[2*i+1]);
       	       		}
       	       	#endif	
       	       	
       	       	num = 0; //peace of mind
       	       	for(i=0; i<4; i++)
       	       	{
       	       		if( is_zero(sol_complex[2*i+1]) )
       	       		{
       	       			sol[num] = sol_complex[2*i];
       	       			num++;
       	       		}
       	       	}
       	       	

	
       }
	
	return num;
}



inline double find_impact_time(double coeff[5], const double dt)
{
     int i, num;
     int flag = 0;
     double tentative = dt;
     double sols[4];
            
      num = solve_quartic(coeff, &sols[0]);
      
      if( num == 0 )
      {
      	      #if DEBUGGING_QUARTIC_SOLVER
      	     		 printf("find_impact_time: The quartic doesn't have any real solutions: parameter values are not realistic\n");
      	      #endif
      	      return -1.0;
      }
      #if DEBUGGING_QUARTIC_SOLVER
      		for (i = 0; i < num; i++)
      		{	        
      		     	 printf("find_impact_time: Real root number %d is: %+.18lf\n",i,sols[i]);      
      		}
      #endif
      
     
      /* Select right tau */
      for(i=0; i<num; i++)
      {
      	      if( sols[i]>0 && sols[i]<tentative)
      	      {
      	      	      tentative = sols[i];
      	      	      flag=1;
      	      }
      	      
      }
      
      if(flag){
      	      #if DEBUGGING_QUARTIC_SOLVER
      	     	 	printf("find_impact_time: Solution is: %+.18lf<%.8lf\n",tentative,dt);
              #endif
              
              return tentative;
      }else{
      	      #if DEBUGGING_QUARTIC_SOLVER
      	      		printf("find_impact_time: None of the solutions is acceptable.\n");
	      #endif
	      
	      return -1.0;
      }
      		       
}



inline double find_impact_time_gsl(double * coeff, const double dt)
{
       double sol[4];
       int flag = 0;
       double tentative = dt;
       int i;
       int num;
       
       
       num = solve_quartic_gsl( coeff,  &(sol[0])  );

       
       for( i=0; i<num; i++)			
       {
       	       if( sol[i]>0.0 && sol[i]<tentative )
       	       {
       	       	       tentative = sol[i];
		       flag = 1;
	       }
	}
	
	if(flag)
	{
		return tentative;
	}else{
		#if DEBUGGING_QUARTIC_SOLVER
			printf("No solution found for impact time!\n");
		#endif
       	        return -1.0;
	}
	
	
}



/* takes a positive real number. If it is 0, it returns 0. */
inline int order_of_magnitude(double a)
{	
	int order_a = 0;
	
	if( a < 1 )
	{
		while( a < 1 )
		{
			a *= 10;
			order_a--;
		}
	}else{
		while ( a >= 1 )
		{
			a *= 0.1;
			order_a++;
		}
		order_a--;
	}
	
	return order_a;
}



/* This function is to be called for a particle at the boundary */
/* coeff is the array containing the coefficient of the full quartic.*/
inline double find_trajectory_boundary_gsl( double * coeff )
{
     int i, sol_num;
     int size_non_0_sols=0;
     double non_0_sols[4]={ 0.0, 0.0, 0.0, 0.0 };
     double tentative;
     int num_positive = 0;
     
     /* We are here because the particle is at the boundary and we want to know what the trajectory looks like. */
     /* This means that coeff[0] is basically 0. We'll use it to compare other numbers, when necessary */ 
     //double time_eps = fabs(coeff[0]); //this was not necessary
     

     /* Since particle is at the boundary, we know one solution is tau =0.0 (approx). So we don't need to solve the quartic */
     /* We'll simplify it first and solve the corresponding cubic: */
     /*  coeff[1] + coeff[2] tau + coeff[3] tau^2 + coeff[4] tau^3==0 */
     
     double a,b,c; 
     gsl_complex * z;
     z = malloc(3*sizeof(gsl_complex));
     
     if( is_zero(coeff[4]) )
     {
     	     //this is a cuadric. Degenerated cases are automatically handled by function below
    	     sol_num = gsl_poly_complex_solve_quadratic(coeff[3],coeff[2],coeff[1],&(z[0]),&(z[1]));
     	     
     }else{
     	     double aux = 1.0/coeff[4];
     	     a = coeff[3] * aux;
     	     b=  coeff[2] * aux;
     	     c = coeff[1] * aux;
     	     
     	     sol_num= gsl_poly_complex_solve_cubic( a,  b,  c, &(z[0]),  &(z[1]),  &(z[2]) );	     
     }
     
     
     /* Trajectory depends on the number of non-zero real solutions */
     /* We first find the number of real solutions...*/
     for(i=0;i<sol_num;i++)
     {
     	     if( is_zero( GSL_IMAG(z[i]) ) && !is_zero( GSL_REAL( z[i] ) )  ) 
     	     {
     	     	     //this is a real, non-zero solution
     	     	     non_0_sols[size_non_0_sols] = GSL_REAL(z[i]);
     	     	     size_non_0_sols++;
     	     }     	     
     }

     //printf("This particle has %d non 0 real solutions\n",size_non_0_sols);

     /* Now we handle each case */
     double tau;
     
     switch( size_non_0_sols )
     {
     	case 0:
     		//sticking the particle:
     		tau = 0.0;
     		break;
     		
     	case 1: tau = non_0_sols[0];
     		break;
     		
     	case 2: 
     		/* If roots have same sign, particle is stuck at boundary */
     		if( non_0_sols[0]*non_0_sols[1] > 0)
     		{
     			tau = 0.0;
     		}else if( non_0_sols[0] > non_0_sols[1] ) /* Otherwise, return the positive element (the greater of the two)*/
     		{
     			tau = non_0_sols[0];
     		}else{
     			tau = non_0_sols[1];
     		}
     		break;
     		
     	case 3:
     		/* In this case we need to find the number of positive solutions and then select the maximum or minimum element	  */
     		for( i=0; i<size_non_0_sols; i++ )
     		{
     			if( non_0_sols[i] > 0 )
     			{
     				num_positive++;
     			}    
     		}
     	     
     		tentative = non_0_sols[0];
     	     
     		if( num_positive >=2 )
     		{
     			/* return the minimum */
     			for( i=1; i<size_non_0_sols; i++ )
     			{
     				if( non_0_sols[i] < tentative )
     				{
     					tentative = non_0_sols[i];
     				}//else do nothing
     			}
		    
     		}else{
     			/* return the maximum */
     			for( i=1; i<size_non_0_sols; i++ )
     			{
     				if( non_0_sols[i] > tentative )
     				{
     					tentative = non_0_sols[i];
     				}//else do nothing
     			}
     		}
     	     
     		tau = tentative;
     		break;
     		
     	default:
     		printf("Unexpected error in find_trajectory gsl. Aborting...\n");
     		exit(EXIT_FAILURE);	 
     		
     }/* end switch */
     

     free(z);
     
     return tau;

}
	


#endif

