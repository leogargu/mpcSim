#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "average.h"


/*-------------------------------*/
/*		MAIN		 */
/*-------------------------------*/

/* call with arguments: <average type> filename (first file index) (last file index) number */
/* with number being 1 for velocities averages, 3.0 for temperature averages and 3*a^3 for pressure averages (perhaps, needs checking) */
int main(int argc, char **argv) {

	/*
	argv[0] is the program name
	argv[1] is the type of average: CAM or SAM
	argv[2] is name of files
	argv[3] is the index of the first file to be processed
	argv[4] is the index of the last file to be processed
	argv[5] is the skip parameter/stride (e.g. if it is 1=take all files, 2=take every other file, etc.)
	argv[6] is a factor (see note above)
	argv[7] is 1 or 0 (verbose mode)
	*/
	
	printf("WARNING: average.c only calculates CAM or SAM averages of scalar quantities.\n");
	
	
	/*---------------------------*/
	/* Main  		     */
	/*---------------------------*/
	if(argc!=8)
	{
		printf("Wrong number of arguments. Call as:\n");
		printf("./average <av.type> <filename> <first file index> <last file index> <factor> <verbose>\n<av.type>\tSAM, CAM or CAM2SAM\n<filename>\tdatafile basic name in averages folder (do not write path)\n<first file index>\t Index of first file (inclusive)\n<last file index>\t Index of last file (inclusive)\n<stride>\t number of files to skip after each sample\n<factor>\tcorrection factor\n<verbose>\tverbose mode 1/0\n");
		printf("Aborting...\n");
		exit(EXIT_FAILURE);
	}
	
	char directory_name[70]="./../experiments/";
	char filename[70]="";
	strcpy(filename,directory_name);
	strcat(filename,argv[2]);
	
	double number;
	number = strtod(argv[6], NULL);
	
	int stride = atoi(argv[5]);
	
	if(strcmp(argv[1], "CAM")==0)
	{
		CAM_average(filename, atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), 1.0/number, atoi(argv[7]) ); //strtol + error checking desirable here

	}else if(strcmp(argv[1], "SAM")==0)
	{
		SAM_average(filename, atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), 1.0/number, atoi(argv[7]));

	}else if(strcmp(argv[1], "CAM2SAM")==0)
	{
		printf("Converting CAM data files to SAM data format: no average is done. Factor value is not used.\n");
		char filename_numbered[50]="";
		char intaschar[20];
		int i, first, last, verbose;
		first = atoi(argv[3]);
		last = atoi(argv[4]);
		verbose = atoi(argv[7]);
		
		for(i=first; i<=last; i=i+stride)
		{
			strcpy(filename_numbered,argv[2]);		
			sprintf(intaschar, "_%d", i);
			strcat(filename_numbered,intaschar);
						
			CAM_to_SAM(directory_name, filename_numbered);
			if(verbose)
			{
				printf("File %s.dat processed\n",filename_numbered);
			}
		}
	}else{
		printf("Type of average not recognized. Aborting...\n");
		exit(EXIT_FAILURE);
	}
	
			
	return 0;
}
