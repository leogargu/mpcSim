#include<gsl/gsl_randist.h>
#include<gsl/gsl_poly.h>
#include <math.h>
#include <time.h>
#include <assert.h>



int num_elements(int * bins, int length)
{
	int i;
	int total=0;
	for(i=0; i<length; i++)
	{
		total += bins[i];
	}
	
	return total;
}

/* The collide step does not change the positions of the particles. */
int main(int argc, char **argv) {
	/* Set up RNG */
	const gsl_rng_type * rngtype; 	 /* Type of rng generator   */
	gsl_rng * r;			 /* Instance of a generator */
    
	//gsl_rng_env_setup();		/* the default seed is used (0) */ //-------------------------------------->TODO: REMEMBER TO CHANGE THE SEED WHEN WORKING CODE IS RUNNING
	
	rngtype = gsl_rng_default; 	/* this is gsl_rng_mt19937 */
	r = gsl_rng_alloc(rngtype);	/* create pointer to instance of a RNG of type rngtype */
	
	long int seed = (long int)time(NULL);
	//seed = 1379012447 ;
	//seed = 1379520223 ;
	printf("RNG seed \t %ld\n",seed);
	gsl_rng_set(r , seed); //sets the seed for the RNG r

	double m_inv = 1.0;
	double T = 1.0;
	
	int test;
	
	if( atoi(argv[1]) == 1 )
	{
	int bin_num = 200; //this must be even
	
	assert( bin_num%2 == 0.0 );
	
	int bins[bin_num];
	int min = -5;
	int max = 5;
	double bin_size = (max-min)/(double)bin_num;
	int i,idx;
	double rnd;
	int middle = bin_num/2;
	
	int samples=1e7;
	int counter=0;
	
	for(i=0; i<bin_num; i++)
	{
		bins[i]=0;
	}
	
	for(i=0; i<samples; i++)
	{
		rnd = gsl_ran_gaussian(r, sqrt(T * m_inv) );
		if(rnd <= min || rnd >= max)
		{
			counter++;
			continue;
		}
		idx = (int)(rnd/bin_size) + middle;
		
		if( rnd<0 )
		{
			idx-=1;
		}
		
		if( rnd==0.0 )
		{
			printf("Zero value obtained\n");
			if( gsl_ran_gaussian_ziggurat(r, sqrt(T * m_inv) ) >0 )
			{
				idx = middle;
			}else{
				idx = middle;//-1
			}
		}
		

		assert( idx>=0 && idx<bin_num );
		
		bins[idx]++;
		
	}
	
	
	assert(num_elements(bins,bin_num)==samples-counter);
	
	FILE * fp;
	fp = fopen("maxwell.dat","w");
	
	for(i=0; i<bin_num; i++)
	{
		fprintf(fp,"%lf\t%lf\n",min+i*bin_size,(bin_num/(max-min))*(bins[i]/(double)(samples-counter)));
	}
		
	fclose(fp);
	}else if(atoi(argv[1]) == 2 )
	{
		double sigma_sq = 0.0;
		int i;
		double aux;
		int max=10000;
		
		for(i=0; i<max; i++)
		{
			aux = gsl_ran_gaussian_ziggurat(r, sqrt(T * m_inv) );
			sigma_sq += aux*aux;
		}
		
		printf("Sigma set to %lf, sigma estimated as %lf\n",T*m_inv,sigma_sq/(double)max);
		
		
		
	}else{
		printf("Argument was %d, call with 1 or 2\n",atoi(argv[1]));
	}
	
	
	
	return 0;
				     
}
