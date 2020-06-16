/*
 *  ellipse.c
 *  
 *
 *  Created by Joydeep Singha on 31/03/15.
 *  Copyright 2015 __MyCompanyName__. All rights reserved.
 *
 basically this program multiplies two matrices which are real. Given the matrices namely "data" and "U", both of which are of order 2n, the code multiplies them a and stores them in the matrix U
 */
 #include<stdio.h>
 #include<stdlib.h>
 #include<math.h>
 #include<gsl/gsl_math.h>
 #include<gsl/gsl_matrix.h>

double ellipse( gsl_matrix * data, gsl_matrix * U, int n);	/*function prototype*/
double ellipse( gsl_matrix * data, gsl_matrix * U, int n)	/*function definition*/
{
	int i, j, k;											/*variable to run loops*/
	double row;											/*variables to calculate intermediate values*/
	double column;
	double value;
	gsl_matrix * ellipse;									/*a temporary gsl matrix to store the result of the matrix multiplication*/
	
	ellipse = gsl_matrix_alloc ((2 * n), (2 * n));			/*allocating the gsl matrix*/
	
	for(i = 0;i < (2 * n);i++)								/*calculating the (i, j)th element of the matrix ellipse*/
	{
		for(j = 0;j < (2 * n);j++)
		{
			value = 0.0;									/*intializing*/
			for(k = 0;k < (2 * n);k++)
			{
				row = gsl_matrix_get (data, i, k);			/*accessing the values of the ith row and jth column*/
				column = gsl_matrix_get (U, k, j);
				value = value + (row * column);			/*multiplying the (i, k)th element and the (k, j)th element and summing*/
			}
			gsl_matrix_set (ellipse, i, j, value);			/*filling up the (i, j)th element in the matrix "ellipse"*/
		}
	}
	gsl_matrix_memcpy (U, ellipse);							/*copying the matrix elements from "ellipse" to 'U'*/
	gsl_matrix_free (ellipse);								/*freeing the variables*/
//	free(value);
//	free(row);
//	free(column);
}
