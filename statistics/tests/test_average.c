#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "./../average.h"


/*-------------------------------*/
/*		MAIN		 */
/*-------------------------------*/

int main(int argc, char **argv) {

	
	int stride = 1;
	/*----------------------------*/
	assert(stride >0 );

	printf("Test SAM 1: stride %d, all cells occupied \n",stride);
	SAM_average("./SAM1/testSAM1", 0, 2, stride, 1.0, 0);
	printf("Test SAM 2: stride %d, one cell empty in one sample \n",stride);
	SAM_average("./SAM2/testSAM2", 0, 2, stride, 1.0, 0);
	printf("Test SAM 3: stride %d, multiple cells empty in multiple samples \n",stride);
	SAM_average("./SAM3/testSAM3", 0, 2, stride, 1.0, 0);
	printf("Test SAM 4: stride %d, multiple cells empty in multiple sample, including first and last cells in all files \n",stride);
	SAM_average("./SAM4/testSAM4", 0, 2, stride, 1.0, 0);

	printf("Test CAM1: stride %d, some empty cells \n",stride);	
	CAM_average("./CAM1/testCAM1", 0, 2, stride, 1.0, 0);
	printf("Test CAM2: stride %d, various empty cells \n",stride);
	CAM_average("./CAM2/testCAM2", 0, 2, stride, 1.0, 0);
	
	
	
	stride = 2;
	/*----------------------------*/
	printf("Test SAM5: stride %d, skipped files are empty \n",stride);
	SAM_average("./SAM5/testSAM5", 0, 4, stride, 1.0, 1);

	printf("Test CAM3: stride %d, skipped files are empty \n",stride);
	CAM_average("./CAM3/testCAM3", 0, 4, stride, 1.0, 1);
	
	
	
	stride = 10;
	/*----------------------------*/
	printf("Test CAM4: stride %d, skipped files are empty \n",stride);
	CAM_average("./CAM4/testCAM4", 0, 25, stride, 1.0, 1);
	
	
	
	/* The way CAMtoSAM is implemented, it only transforms 1 file at a time. It is average.c main() function what transforms multiple files at once*/
	/* Therefore, to test the stride functionality it is necessary to refactor teh code or to duplicate it in here.*/
	CAM_to_SAM("./CAMtoSAM/", "testCAM1");
	CAM_to_SAM("./CAMtoSAM/", "testCAM2");


	printf("----------------------------\n");
	return 0;
}