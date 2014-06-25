#include <stdlib.h>
#include <stdio.h>
//#include <math.h>
//#include <time.h>
#include <string.h>

#include "./../average.h"


/*-------------------------------*/
/*		MAIN		 */
/*-------------------------------*/

/* call with arguments: filename (first file index) (last file index) number */
/* with number being 1 for velocities averages, 3.0 for temperature averages and 3*a^3 for pressure averages (perhaps, needs checking) */
int main(int argc, char **argv) {

	/*
	argv[0] is the program name
	argv[1] is the type of average: CAM or SAM
	argv[2] is name of files
	argv[3] is the index of the first file to be processed
	argv[4] is the index of the last file to be processed
	argv[5] is a factor (see note above)
	*/

	SAM_average("./SAM1/testSAM1", 0, 2, 1.0, 0);
	SAM_average("./SAM2/testSAM2", 0, 2, 1.0, 0);
	SAM_average("./SAM3/testSAM3", 0, 2, 1.0, 0);
	SAM_average("./SAM4/testSAM4", 0, 2, 1.0, 0);
	
	CAM_average("./CAM1/testCAM1", 0, 2, 1.0, 0);
	CAM_average("./CAM2/testCAM2", 0, 2, 1.0, 0);

	CAM_to_SAM("./CAMtoSAM/", "testCAM1");
	CAM_to_SAM("./CAMtoSAM/", "testCAM2");


	return 0;
}