#include <stdlib.h>
#include <stdio.h>
//#include <math.h>
//#include <time.h>
#include <string.h>

#include "./../../average.h"


//Compared two int arrays of length 4. Length is not checked.
int are_equal(int * a, int * b)
{
	int equal = 1;
	int i;
	
	for(i=0;i<4;i++)
	{
		//printf("%d\t%d\n",a[i],b[i]);
		if(a[i]!=b[i])
		{
			equal=0;
			break;
		}
	}
	
	return equal;
}


/*-------------------------------*/
/*		MAIN		 */
/*-------------------------------*/

int main(int argc, char **argv) {

	/* General set up */
	int i;
	int sym_cells[4] = {0,0,0,0};
	
	/* Set up test cases */
	int input[7][3] = {
		{6,6,5},
		{6,6,9},
		{6,6,15},
		{4,4,2},
		{4,4,7},
		{6,6,10},
		{6,6,17}
	};
	
	
	/* Set up solutions to tests cases */
	int solution[7][4] = {
		{0,5,30,35},
		{8,9,26,27},
		{14,15,20,21},
		{1,2,13,14},
		{4,7,8,11},
		{7,10,25,28},
		{12,17,18,23}
	};
	
	/* Remember to update this when adding new test cases */
	int num_tests = 7;
	
	
	/* Run tests */
	for(i=0;i<num_tests;i++)
	{
		printf("Test %d ... ",i);
		find_sym_cells( input[i][0], input[i][1], input[i][2], sym_cells );
		if(are_equal(sym_cells,solution[i]))
		{
			printf("PASS\n");
		}else{
			printf("*FAIL*\n");
		}
	}
	
	
	
	return 0;
}