/*
 *  normalize.c
 *  
 *
 *  Created by Joydeep Singha on 01/04/15.
 *  Copyright 2015 __MyCompanyName__. All rights reserved.
 *
 Given a matrix U each whose columns are considered to be vectors with 2n number of elements this function normalizes each of those column vectors 
 */

 #include<stdio.h>
 #include<stdlib.h>
 #include<math.h>
 #include<gsl/gsl_math.h>
 #include<gsl/gsl_matrix.h>
double normalize(gsl_matrix * U, int n);					/*function prototype*/
double norm(double *vector, int n);							/*function prototype normalizing a vector*/
double normalize(gsl_matrix * U, int n)						/*function definition*/
{
	int i,j;												/*interger variables to run loops*/
	double *row;											/*array to store each column from the matrix U*/
	row = (double *) malloc((2 * n) * sizeof(double));		/*dynamic allocation of the array "row"*/
	if(row == NULL)											
	{
		printf("cannot allocate memory for row\n");
		exit(1);
	}
	
	for(i = 0;i < (2 * n);i++)								/*reading columns from the matrix U*/
	{
		for(j = 0;j < (2 * n);j++)
		{ 
			row[j] = gsl_matrix_get (U, j, i);				/*copying the column i from the matrix U to the array "row"*/
		}
		norm(row, n);										/*calling the function to calculate the unit vector*/
		for(j = 0;j < (2 * n);j++)
		{
			gsl_matrix_set (U, j, i, row[j]);				/*replacing U with the normalized column vectors*/
		}
	}
	free(row);												/*freeing the arrays*/
}
